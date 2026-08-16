"""Reference-genome resolution — the "Reference/Annotation" stage feeding
raw-genome candidate-variant discovery. Generically species-agnostic: no
species is hardcoded anywhere below. Confirmed live 2026-08-15 against dog
(Ollie's demo species) because that's the data on hand, not because
anything here is dog-specific — same three calls work for any organism
NCBI has a reference assembly for.

The chain, each step confirmed against real NCBI data:
  1. `db="assembly"` search by organism name, filtered to
     `refseq_category == "reference genome"` (or "representative genome" as
     a fallback for species without a full "reference" designation) — finds
     the canonical assembly. For dog: GCF_011100685.1 (UU_Cfam_GSD_1.0).
  2. Every NCBI assembly publishes a fixed-format `*_assembly_report.txt`
     at a predictable FTP path — column 3 is the bare chromosome number
     (matches PLINK .tped chrom values like "12"), column 7 is the RefSeq
     sequence accession (e.g. "NC_049233.1" for dog chr12). Cached per
     organism so a 230K-marker file doesn't refetch this once per marker.
  3. `run_ncbi_efetch` with that accession + seq_start/seq_stop pulls a
     real reference sequence window around any position — confirmed live
     against Ollie's actual chr12:3991 marker (real 1000bp returned).

Known real limitation: many species — especially non-model / endangered
ones, which is exactly where "species-agnostic" matters most — only have
scaffold-level assemblies, not chromosome-level ones. `resolve_chromosome`
returns None rather than guessing when a chromosome isn't in the report;
callers must treat that as "can't annotate this marker," not a bug to
route around.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache

import requests


class ReferenceGenomeError(RuntimeError):
    pass


@dataclass
class _Tile:
    start: int  # 1-based inclusive genomic coordinate of sequence[0]
    stop: int  # 1-based inclusive genomic coordinate of sequence[-1]
    sequence: str


class _ChromosomeTileCache:
    """Turns "one live NCBI efetch per marker" into "one live efetch per
    ~2Mb region" for markers that cluster together on a chromosome — the
    normal case when scanning a SNP array in file order, since PLINK .tped
    files are sorted by chromosome/position. Real measured speedup
    (2026-08-15): fetching 2,000,001 bases in one efetch call took 1.7s;
    267 separate 200bp calls covering far fewer total bases took 150.6s.

    Not a cache of interpretations or scores — every base served still
    came from a real NCBI efetch, just fetched in bigger batches and
    sliced locally afterward. Keyed by RefSeq accession (not chromosome
    label) so it can't collide across organisms/assemblies sharing a
    chromosome name. Process-lifetime only, no disk persistence — a fresh
    run always re-verifies against NCBI rather than trusting a stale file.
    """

    def __init__(self, tile_size: int = 2_000_000):
        self._tile_size = tile_size
        self._tiles: dict[str, _Tile] = {}

    def get_window(self, accession: str, position: int, window: int) -> tuple[str, int] | None:
        """Returns (sequence, window_start) for a window of `window` bases
        centered on `position`, or None if the underlying fetch fails."""
        half = window // 2
        want_start = max(1, position - half)
        want_stop = position + half

        tile = self._tiles.get(accession)
        if tile is None or want_start < tile.start or want_stop > tile.stop:
            fetch_start = max(1, position - self._tile_size // 2)
            fetch_stop = position + self._tile_size // 2
            seq = self._live_fetch(accession, fetch_start, fetch_stop)
            if seq is not None:
                tile = _Tile(start=fetch_start, stop=fetch_start + len(seq) - 1, sequence=seq)
                self._tiles[accession] = tile

        if tile is not None and tile.start <= want_start and want_stop <= tile.stop:
            offset = want_start - tile.start
            return tile.sequence[offset: offset + (want_stop - want_start + 1)], want_start

        # Tile fetch failed, or position falls outside it even after a
        # (re)fetch (e.g. right at a chromosome's start/end) — fall back to
        # one direct small fetch rather than give up on this marker.
        seq = self._live_fetch(accession, want_start, want_stop)
        if seq is None:
            return None
        return seq, want_start

    @staticmethod
    def _live_fetch(accession: str, start: int, stop: int) -> str | None:
        from proto_tools import run_ncbi_efetch
        from proto_tools.tools.database_retrieval.ncbi.efetch import NCBIEfetchConfig, NCBIEfetchInput

        try:
            out = run_ncbi_efetch(
                NCBIEfetchInput(db="nuccore", identifier=accession, return_format="fasta", seq_start=start, seq_stop=stop),
                NCBIEfetchConfig(ncbi_email="wildtype-hackathon@example.com"),
            )
        except Exception:
            # run_ncbi_efetch raises rather than returning success=False for
            # some real failure modes — live-caught 2026-08-15: NCBI answers
            # 400 Bad Request (not a graceful clamp, unlike its handling of
            # `start` near position 1) when `stop` runs past a chromosome's
            # true end, which a full-genome scan hits on every chromosome's
            # last marker. Also covers plain network hiccups (timeouts,
            # connection resets) hit live earlier this session. Either way,
            # this is "can't resolve this one window," not a reason to crash
            # an entire multi-hour run — treated the same as out.success=False.
            return None
        if not out.success or not out.fasta_records:
            return None
        record = out.fasta_records[0]
        seq = getattr(record, "sequence", None) or str(record)
        seq = re.sub(r"[^ACGTNacgtn]", "", seq)
        return seq.upper() or None


_tile_cache = _ChromosomeTileCache()


@dataclass
class AssemblyInfo:
    organism: str
    accession: str  # e.g. "GCF_011100685.1"
    name: str  # e.g. "UU_Cfam_GSD_1.0"
    report_url: str
    category: str  # "reference genome" | "representative genome"


def resolve_reference_assembly(organism: str) -> AssemblyInfo | None:
    """Find the canonical NCBI RefSeq assembly for any organism. Returns
    None (not a guess) if NCBI has no reference/representative-quality
    assembly for this species — common for species with only fragmentary
    genome data."""
    from proto_tools import run_ncbi_esearch, run_ncbi_esummary
    from proto_tools.tools.database_retrieval.ncbi.esearch import NCBIEsearchConfig, NCBIEsearchInput
    from proto_tools.tools.database_retrieval.ncbi.esummary import NCBIEsummaryConfig, NCBIEsummaryInput

    email = "wildtype-hackathon@example.com"
    search = run_ncbi_esearch(
        NCBIEsearchInput(db="assembly", search_term=f"{organism}[Organism]", max_results=100),
        NCBIEsearchConfig(ncbi_email=email),
    )
    if not search.success or not search.ids:
        return None

    summary = run_ncbi_esummary(
        NCBIEsummaryInput(db="assembly", identifier=",".join(search.ids)),
        NCBIEsummaryConfig(ncbi_email=email),
    )
    if not summary.success:
        return None

    # Prefer a full "reference genome" (one per species, NCBI's own pick);
    # fall back to "representative genome" (used for species without a
    # single designated reference, e.g. many non-model organisms).
    for wanted_category in ("reference genome", "representative genome"):
        for uid in search.ids:
            rec = summary.summary.get(uid, {})
            if rec.get("refseq_category") == wanted_category:
                accession = rec.get("assemblyaccession")
                name = rec.get("assemblyname", "")
                ftp_report = rec.get("ftppath_assembly_rpt")
                if accession and ftp_report:
                    return AssemblyInfo(
                        organism=organism,
                        accession=accession,
                        name=name,
                        report_url=ftp_report.replace("ftp://", "https://"),
                        category=wanted_category,
                    )
    return None


def resolve_assembly_by_accession(accession: str, organism: str) -> AssemblyInfo | None:
    """Construct an AssemblyInfo for a SPECIFIC, KNOWN assembly accession —
    unlike resolve_reference_assembly(), doesn't require it to hold NCBI's
    current "reference genome"/"representative genome" category. Needed
    because genotyping platforms map their coordinates to whatever
    assembly was current when the chip was designed, then keep using it
    for consistency — they don't silently follow NCBI's reference pointer
    when a newer assembly gets promoted.

    Real case this exists for (2026-08-15): Ollie's Embark TPED file uses
    CanFam3.1 (GCF_000002285.3) coordinates, confirmed two ways —
    (1) the chr2 accession's real length (85,426,708, from CanFam3.1's own
    assembly report) matches almost exactly where a full-genome scan
    against the newer UU_Cfam_GSD_1.0 started throwing "past end of
    chromosome" errors, and (2) the same numeric position (chr1:68723) on
    the two assemblies returns completely unrelated DNA sequence — not
    minor drift, a real coordinate-system mismatch. UU_Cfam_GSD_1.0 is
    what resolve_reference_assembly("Canis lupus familiaris") returns
    (it's NCBI's current pick), which is why every finding from earlier
    this session was silently comparing Ollie's genotype against the
    wrong assembly.

    Generic, not dog-specific — any caller who has independently
    determined which assembly their input data's coordinates actually use
    can pass that accession here instead of trusting NCBI's current-
    reference default.
    """
    from proto_tools import run_ncbi_esearch, run_ncbi_esummary
    from proto_tools.tools.database_retrieval.ncbi.esearch import NCBIEsearchConfig, NCBIEsearchInput
    from proto_tools.tools.database_retrieval.ncbi.esummary import NCBIEsummaryConfig, NCBIEsummaryInput

    email = "wildtype-hackathon@example.com"
    # esummary needs a numeric UID, not the "GCF_..." accession string
    # itself — confirmed live 2026-08-15 (passing the accession directly
    # to esummary returns {'uids': []}, no error, just nothing). esearch's
    # [Assembly Accession] field resolves the accession to its UID first.
    search = run_ncbi_esearch(
        NCBIEsearchInput(db="assembly", search_term=f"{accession}[Assembly Accession]", max_results=5),
        NCBIEsearchConfig(ncbi_email=email),
    )
    if not search.success or not search.ids:
        return None
    summary = run_ncbi_esummary(
        NCBIEsummaryInput(db="assembly", identifier=",".join(search.ids)),
        NCBIEsummaryConfig(ncbi_email=email),
    )
    if not summary.success:
        return None
    for rec in summary.summary.values():
        if not isinstance(rec, dict) or rec.get("assemblyaccession") != accession:
            continue
        ftp_report = rec.get("ftppath_assembly_rpt")
        if not ftp_report:
            return None
        return AssemblyInfo(
            organism=organism,
            accession=accession,
            name=rec.get("assemblyname", ""),
            report_url=ftp_report.replace("ftp://", "https://"),
            category=rec.get("refseq_category") or "explicit override",
        )
    return None


@lru_cache(maxsize=32)
def _chromosome_accession_map(report_url: str) -> dict[str, str]:
    """Parse an NCBI assembly_report.txt into {chromosome_number: refseq_accession}.
    Cached per report URL — this file covers every chromosome, fetch once
    per organism/assembly, not once per marker."""
    resp = requests.get(report_url, timeout=30)
    resp.raise_for_status()
    mapping = {}
    for line in resp.text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        _name, role, assigned_molecule, _loc_type, _genbank, _rel, refseq_accn = parts[:7]
        if role == "assembled-molecule" and refseq_accn and refseq_accn != "na":
            mapping[assigned_molecule] = refseq_accn
    return mapping


def resolve_chromosome_accession(assembly: AssemblyInfo, chromosome: str) -> str | None:
    """RefSeq accession for one chromosome of a resolved assembly, or None
    if that chromosome isn't in the report (e.g. asking for "12" against a
    species where only scaffolds exist, or a genuinely out-of-range value)."""
    mapping = _chromosome_accession_map(assembly.report_url)
    return mapping.get(chromosome) or mapping.get(chromosome.lstrip("chr"))


@dataclass
class ReferenceWindow:
    sequence: str
    start: int  # 1-based genomic coordinate of sequence[0], as returned by NCBI
    accession: str

    def base_at(self, position: int) -> str | None:
        """Reference base at a 1-based genomic `position`, or None if it
        falls outside this window (shouldn't happen for the position this
        window was fetched around, but real callers should still check)."""
        idx = position - self.start
        if 0 <= idx < len(self.sequence):
            return self.sequence[idx]
        return None


def fetch_reference_context(
    organism: str,
    chromosome: str,
    position: int,
    window: int = 500,
    assembly: AssemblyInfo | None = None,
) -> ReferenceWindow | None:
    """Real reference DNA sequence centered on `position` (1-based, as in
    PLINK .tped) on `chromosome` of `organism`'s reference genome, plus the
    genomic coordinate of the window's first base — needed to find exactly
    which base in the window corresponds to `position` (NCBI clamps `start`
    to 1 near a chromosome's beginning, so it isn't always `position - window//2`).
    Returns None (not a placeholder) if the assembly or chromosome can't be
    resolved — callers must handle that explicitly, not silently substitute
    fake sequence for a DNA foundation model to score."""
    assembly = assembly or resolve_reference_assembly(organism)
    if assembly is None:
        return None
    accession = resolve_chromosome_accession(assembly, chromosome)
    if accession is None:
        return None

    result = _tile_cache.get_window(accession, position, window)
    if result is None:
        return None
    seq, start = result
    return ReferenceWindow(sequence=seq, start=start, accession=accession)


def fetch_reference_window(
    organism: str,
    chromosome: str,
    position: int,
    window: int = 500,
    assembly: AssemblyInfo | None = None,
) -> str | None:
    """Backward-compatible string-only wrapper around fetch_reference_context."""
    ctx = fetch_reference_context(organism, chromosome, position, window=window, assembly=assembly)
    return ctx.sequence if ctx else None


if __name__ == "__main__":
    import sys

    organism = sys.argv[1] if len(sys.argv) > 1 else "Canis lupus familiaris"
    chrom = sys.argv[2] if len(sys.argv) > 2 else "12"
    pos = int(sys.argv[3]) if len(sys.argv) > 3 else 3991

    assembly = resolve_reference_assembly(organism)
    if assembly is None:
        print(f"No reference/representative assembly found for {organism!r}")
        raise SystemExit(1)
    print(f"Assembly: {assembly.accession} ({assembly.name}), category={assembly.category}")

    accession = resolve_chromosome_accession(assembly, chrom)
    print(f"Chromosome {chrom} -> {accession}")

    seq = fetch_reference_window(organism, chrom, pos, window=200, assembly=assembly)
    print(f"Window around {chrom}:{pos} ({len(seq) if seq else 0} bp): {seq}")
