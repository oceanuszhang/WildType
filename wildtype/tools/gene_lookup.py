"""Deterministic gene lookup — no Claude involved — backed by an
assembly's own official GFF3 gene annotation file, fetched once per
assembly and cached in memory. Factored out of
wildtype/agent/iterative_agent.py's search_genes_near_position tool so
the same real annotation also drives genome_pipeline.py's triage stage
(broad, cheap scan first, full LLM-driven investigation only on
survivors) and wildtype/agent/report_validation.py's gene-name lookups.

History, all from 2026-08-15, because this file went through two real,
live-caught correctness bugs in one day and the reasoning matters:

  1. First version queried NCBI's Gene database (db="gene", esearch +
     esummary) per position, filtering out any record NCBI marked
     `currentid` (superseded) on the theory those were just stale
     duplicates of the live record.

  2. That broke the moment the input data (Ollie's Embark TPED) turned
     out to be in CanFam3.1 coordinates, not UU_Cfam_GSD_1.0 (NCBI's
     current "reference genome" pick for dog — see reference_genome.py's
     resolve_assembly_by_accession() docstring for that whole story).
     NCBI's Gene database only reports the CURRENT assembly's coordinates
     on the live record; older coordinates, when they exist at all, only
     show up on records marked superseded. Fixed by matching a record's
     `chraccver` against the TARGET assembly's own chromosome accession
     instead of just preferring "not superseded" — with the accession
     version stripped (`"NC_006583.2"` -> `"NC_006583"`) so a gene's
     coordinates from a slightly different accession revision would still
     match.

  3. That version-stripping was itself wrong, caught by cross-referencing
     the fix against the assembly's own GFF file: CanFam3.1's real,
     authoritative chr1 accession (per its own assembly_report.txt) is
     `NC_006583.3`. The "superseded" Gene-database record used as the
     CanFam3.1 source for ENPP1 was on `NC_006583.2` — NOT the same
     coordinate system despite matching after stripping the version, off
     by over 3 million bases (251,985-322,081 real vs. 3,251,711-3,321,650
     from the mismatched record). NCBI's Gene database's historical-record
     retention is inconsistent per-gene (PRCD/GPT/ABCB1 had NO superseded
     record at all to fall back on) and, as shown here, not reliably
     version-safe when one does exist.

  Rather than keep patching a fundamentally unreliable secondary data
  source, this now fetches the assembly's OWN official gene annotation
  (GCF_..._genomic.gff.gz, same FTP directory as the assembly_report.txt
  already used for chromosome accessions) — one file, one unambiguous
  version, matched with the EXACT accession from that same assembly's own
  report, no stripping, no guessing. Real cost: 4.6s to download and
  parse a whole dog genome's ~31,800 gene records into memory, once per
  assembly, then every lookup (by position or by name) is free.
"""
from __future__ import annotations

import gzip
import io
from dataclasses import dataclass
from urllib.parse import unquote


@dataclass
class NearbyGene:
    name: str | None
    description: str | None
    start: int | None
    stop: int | None
    variant_inside_gene_span: bool | None


@dataclass
class GeneCoordinates:
    chromosome: str
    start: int
    stop: int


@dataclass
class _GeneRecord:
    name: str | None
    description: str | None
    chromosome: str  # the label used in the input file, e.g. "1", not the accession
    start: int
    stop: int


class _AssemblyGeneAnnotation:
    """One assembly's full gene set, fetched and parsed once from its own
    official GFF3 file, then served from memory for every lookup — by
    position (triage) or by name (report_validation.py)."""

    def __init__(self) -> None:
        self._cache: dict[str, list[_GeneRecord]] = {}  # keyed by assembly.accession

    def _get(self, assembly) -> list[_GeneRecord]:
        if assembly.accession in self._cache:
            return self._cache[assembly.accession]

        from wildtype.tools.reference_genome import _chromosome_accession_map

        accession_to_chrom = {accn: chrom for chrom, accn in _chromosome_accession_map(assembly.report_url).items()}
        gff_url = assembly.report_url.replace("_assembly_report.txt", "_genomic.gff.gz")

        records: list[_GeneRecord] = []
        try:
            import requests

            resp = requests.get(gff_url, timeout=120)
            resp.raise_for_status()
            with gzip.GzipFile(fileobj=io.BytesIO(resp.content)) as f:
                for raw_line in f:
                    line = raw_line.decode("utf-8", errors="replace")
                    if line.startswith("#"):
                        continue
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 9 or cols[2] != "gene":
                        continue
                    seqid = cols[0]
                    chrom = accession_to_chrom.get(seqid)
                    if chrom is None:
                        continue  # a contig/scaffold not in the assembly report's chromosome list — skip, not a chromosome-level gene
                    try:
                        start, stop = int(cols[3]), int(cols[4])
                    except ValueError:
                        continue
                    attrs = dict(
                        pair.split("=", 1) for pair in cols[8].split(";") if "=" in pair
                    )
                    name = unquote(attrs.get("Name", "")) or None
                    description = unquote(attrs.get("description", "")) or None
                    records.append(_GeneRecord(name=name, description=description, chromosome=chrom, start=min(start, stop), stop=max(start, stop)))
        except Exception:
            records = []  # cached as empty so a bad fetch doesn't retry forever within one run; caller sees "no genes" not a crash

        self._cache[assembly.accession] = records
        return records


_annotation_cache = _AssemblyGeneAnnotation()


def find_genes_near_position(
    organism: str, chromosome: str, position: int, window_bp: int = 100000, *, assembly=None
) -> tuple[list[NearbyGene], str | None]:
    """Genes within `window_bp` of `position`, from the assembly's own
    official annotation (see module docstring). Returns (genes, error) —
    error is None on success (even if zero genes found; that's real, not
    a failure).

    `assembly` (a wildtype.tools.reference_genome.AssemblyInfo): REQUIRED
    for a correct answer on data that isn't on NCBI's current reference —
    without it, falls back to resolve_reference_assembly()'s pick.
    """
    from wildtype.tools.reference_genome import resolve_reference_assembly

    if assembly is None:
        assembly = resolve_reference_assembly(organism)
    if assembly is None:
        return [], f"No reference assembly resolvable for organism={organism!r}"

    records = _annotation_cache._get(assembly)
    if not records:
        return [], f"No gene annotation available for assembly {assembly.accession}"

    chrom_label = chromosome.lstrip("chr")
    lo, hi = max(1, position - window_bp), position + window_bp
    result: list[NearbyGene] = []
    for rec in records:
        if rec.chromosome != chrom_label and rec.chromosome != chromosome:
            continue
        if rec.stop < lo or rec.start > hi:
            continue
        inside = rec.start <= position <= rec.stop
        result.append(NearbyGene(name=rec.name, description=rec.description, start=rec.start, stop=rec.stop, variant_inside_gene_span=inside))
    return result, None


def resolve_gene_coordinates(gene_symbol: str, organism: str, assembly) -> GeneCoordinates | None:
    """The reverse of find_genes_near_position — given a gene SYMBOL (e.g.
    from an uploaded health-report PDF), find its real coordinates in a
    SPECIFIC assembly, from that assembly's own official annotation.
    Returns None if the gene isn't in this assembly's gene set at all.
    """
    records = _annotation_cache._get(assembly)
    target = gene_symbol.strip().upper()
    for rec in records:
        if rec.name and rec.name.strip().upper() == target:
            return GeneCoordinates(chromosome=rec.chromosome, start=rec.start, stop=rec.stop)
    return None
