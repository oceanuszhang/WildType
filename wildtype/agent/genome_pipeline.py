"""Full raw-genome pipeline — Animal Genome + Species Info -> Input
Validation -> Reference/Annotation -> Candidate Variants -> Agent
Orchestrator (iterative, LLM-driven — PROJECT_SPEC.md sections 5.1, 8) ->
Evidence Aggregation -> Risk Assessment -> JSON/CSV + Dashboard.

Supersedes the Embark-CSV path (wildtype/agent/loop.py's run_pipeline) as
the primary input: no Embark report, no pre-digested risk labels — this
starts from raw .tped/.tfam (or .vcf, via wildtype/parsers/vcf.py) files
plus species info, and does its own reference comparison to find candidate
variants. Ollie's data is the demo instance; nothing below is dog-specific
— see wildtype/tools/reference_genome.py's docstring for what's been
cross-species validated.

2026-08-15 architecture change: orchestration is now the iterative,
LLM-driven agent in wildtype/agent/iterative_agent.py — Claude decides
which tools to call and when it has enough evidence, per PROJECT_SPEC.md
section 5.1 ("The LLM is the agent, not the entire system") and section 8
("The agent should be capable of iterating"). This replaces the earlier
deterministic single-Evo2-call orchestration in this file's first version.
Deterministic code still owns parsing, candidate selection, and output
serialization — the LLM only orchestrates investigation of variants
already identified as real candidates.

Chosen scope (2026-08-15): full computational scan, no gene-focused
prefilter — "Candidate Variants" means "genuinely deviates from the
species reference genome at this position," not "falls in a known disease
gene." This is a real, principled filter (not arbitrary sampling).

Real, honest constraint this doesn't paper over: each candidate costs one
live NCBI reference fetch, and the agent's own investigation costs several
more live tool calls (literature, gene lookup, scoring) — neither instant
nor free. `max_candidates` bounds a run to what's actually executable and
verifiable today; the architecture itself doesn't cap at 230K, today's
compute/budget/time does.

2026-08-15 — added a triage stage between "candidate variants" and "full
investigation" (spec's own diagram implies this: everything gets a
computational scan, only some of that goes on to literature research).
Without it, "scan a lot" and "fully investigate each one with an 8-12
iteration LLM agent + live literature search" were the same knob, and that
knob had to stay tiny (3) to keep a run's cost/time bounded — so the scan
never actually got broad. Now they're two separate stages with two separate
budgets:
  1. find_candidate_variants — reference-deviation scan, deterministic,
     one NCBI fetch per marker. This is the "a lot of genes scanned" part;
     its cap (`max_candidates`) can be much larger now (default 150) since
     nothing here is an LLM call.
  2. triage_candidates_near_genes — deterministic positional NCBI gene
     lookup per candidate (no Claude, no literature, no foundation model
     yet), splitting candidates into "falls in/near an annotated gene"
     vs. "doesn't." Still cheap — one more NCBI call each.
  3. investigate_variant (unchanged) — the expensive, LLM-driven,
     foundation-model-scoring, literature-citing stage — now runs only on
     stage 2's survivors, bounded by `max_investigated` (default 15).
Variants that don't survive triage are NOT silently dropped — they're
still counted in TriageStats and the summary the dashboard shows, just not
expanded into full RiskFinding cards (there's genuinely nothing to
investigate about a marker that overlaps no annotated gene at all).
"""
from __future__ import annotations

import csv as csv_module
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from wildtype.agent.iterative_agent import investigate_variant
from wildtype.agent.schema import AnalysisRun, RiskFinding
from wildtype.parsers.common_variant import CommonVariant
from wildtype.parsers.tped import Individual, Marker, iter_tped, load_tfam
from wildtype.tools.reference_genome import AssemblyInfo, fetch_reference_context, resolve_reference_assembly


# --- Stage 1: Input Validation ----------------------------------------------


class GenomeInputError(RuntimeError):
    pass


@dataclass
class GenomeInput:
    tped_path: str
    tfam_path: str
    species: str
    organism: str
    individual: Individual
    marker_count: int


def validate_genome_input(tped_path: str, tfam_path: str, species: str, organism: str | None = None) -> GenomeInput:
    """Real validation, not a rubber stamp — fails loudly on empty/malformed
    input rather than letting garbage flow silently into later stages."""
    tped_p, tfam_p = Path(tped_path), Path(tfam_path)
    if not tped_p.exists():
        raise GenomeInputError(f".tped not found: {tped_path}")
    if not tfam_p.exists():
        raise GenomeInputError(f".tfam not found: {tfam_path}")

    individuals = load_tfam(tfam_path)
    if not individuals:
        raise GenomeInputError(f".tfam has no individuals: {tfam_path}")
    if len(individuals) > 1:
        raise GenomeInputError(
            f".tfam has {len(individuals)} individuals — this pipeline assumes one "
            "(matches Ollie's export shape); multi-individual .tped needs per-sample "
            "genotype columns handled explicitly, not silently taking the first."
        )

    marker_count = sum(1 for _ in iter_tped(tped_path))
    if marker_count == 0:
        raise GenomeInputError(f".tped has no markers: {tped_path}")

    return GenomeInput(
        tped_path=tped_path,
        tfam_path=tfam_path,
        species=species,
        organism=organism or species,
        individual=individuals[0],
        marker_count=marker_count,
    )


# --- Stage 2+3: Reference/Annotation -> Candidate Variants ------------------


@dataclass
class CandidateScanStats:
    markers_scanned: int = 0
    no_calls_skipped: int = 0
    reference_unresolvable_skipped: int = 0
    candidates_found: int = 0
    stopped_early: bool = False


def _marker_to_common_variant(marker: Marker, ref_base: str, sample_id: str, species: str) -> CommonVariant:
    return CommonVariant(
        sample_id=sample_id,
        species=species,
        chromosome=marker.chrom,
        position=marker.position,
        reference=ref_base,
        # "alternate" is whichever called allele differs from reference —
        # arbitrary choice between allele1/allele2 if both deviate (compound
        # het) since CommonVariant models one ref/alt pair; genotype still
        # carries both real called alleles.
        alternate=marker.allele1 if marker.allele1 != ref_base else marker.allele2,
        genotype=f"{marker.allele1}/{marker.allele2}",
    )


def find_candidate_variants(
    tped_path: str,
    species: str,
    organism: str,
    sample_id: str,
    *,
    window: int = 200,
    max_candidates: int | None = 20,
    max_markers_scanned: int | None = None,
    assembly: AssemblyInfo | None = None,
) -> tuple[list[CommonVariant], CandidateScanStats]:
    """Scan raw markers, compare each called genotype against the real
    species reference base at that position, and keep the ones that
    genuinely deviate — Reference/Annotation + Candidate Variants combined,
    since "annotation" here IS the reference comparison. Returns
    CommonVariant objects (spec section 7) ready for the iterative agent.

    `max_candidates` / `max_markers_scanned` bound a run to what's
    verifiable today (each marker costs one live NCBI fetch) — set both to
    None for a genuine full-genome scan, understanding that's a multi-hour,
    rate-limit-sensitive undertaking, not a single quick call.
    """
    assembly = assembly or resolve_reference_assembly(organism)
    if assembly is None:
        raise GenomeInputError(f"No NCBI reference/representative assembly found for organism={organism!r}")

    stats = CandidateScanStats()
    candidates: list[CommonVariant] = []

    for marker in iter_tped(tped_path):
        if max_markers_scanned is not None and stats.markers_scanned >= max_markers_scanned:
            stats.stopped_early = True
            break
        stats.markers_scanned += 1

        if marker.is_no_call:
            stats.no_calls_skipped += 1
            continue

        ctx = fetch_reference_context(organism, marker.chrom, marker.position, window=window, assembly=assembly)
        if ctx is None:
            stats.reference_unresolvable_skipped += 1
            continue
        ref_base = ctx.base_at(marker.position)
        if ref_base is None:
            stats.reference_unresolvable_skipped += 1
            continue

        if marker.allele1 != ref_base or marker.allele2 != ref_base:
            candidates.append(_marker_to_common_variant(marker, ref_base, sample_id, species))
            stats.candidates_found += 1
            if max_candidates is not None and stats.candidates_found >= max_candidates:
                stats.stopped_early = True
                break

    return candidates, stats


# --- Stage 3.5: Triage — cheap positional gene lookup, no LLM yet ----------


@dataclass
class TriageStats:
    candidates_checked: int = 0
    in_or_near_gene: int = 0
    no_gene_nearby: int = 0
    lookup_failed: int = 0
    stopped_early: bool = False


def triage_candidates_near_genes(
    candidates: list[CommonVariant],
    organism: str,
    *,
    window_bp: int = 100000,
    max_investigated: int | None = 15,
    max_checked: int | None = 2000,
) -> tuple[list[CommonVariant], TriageStats]:
    """Deterministic pre-filter between the broad reference-deviation scan
    and the expensive LLM-driven investigation. One real NCBI positional
    gene query per candidate (wildtype/tools/gene_lookup.py — same query
    the agent's own search_genes_near_position tool uses), no Claude call,
    no literature search, no foundation-model scoring yet.

    A lookup failure (NCBI hiccup) is NOT the same as "no gene nearby" —
    counted separately in TriageStats.lookup_failed and the candidate is
    still passed through to full investigation, since we genuinely don't
    know whether it's near a gene or not; better to spend an investigation
    on it than silently drop it because NCBI blinked.

    `max_investigated` caps how many gene-proximity survivors go on to
    full investigation. `max_checked` is a separate safety bound on how
    many candidates get a live lookup at all — unlike the reference scan,
    triage isn't cached (each candidate is a different genomic position, a
    different esearch), so a gene-sparse candidate list could otherwise
    make this stage run arbitrarily long searching for survivors that
    never come. Hitting `max_checked` is recorded as `stopped_early`, not
    silently swallowed.
    """
    from wildtype.tools.gene_lookup import find_genes_near_position

    stats = TriageStats()
    survivors: list[CommonVariant] = []

    for variant in candidates:
        if max_checked is not None and stats.candidates_checked >= max_checked:
            stats.stopped_early = True
            break
        stats.candidates_checked += 1
        genes, error = find_genes_near_position(organism, variant.chromosome, variant.position, window_bp)

        if error:
            stats.lookup_failed += 1
            survivors.append(variant)
        elif genes:
            stats.in_or_near_gene += 1
            survivors.append(variant)
        else:
            stats.no_gene_nearby += 1

        if max_investigated is not None and len(survivors) >= max_investigated:
            break

    return survivors, stats


# --- End-to-end run: Stage 4 (iterative agent) through JSON/CSV output -----


def run_genome_pipeline(
    tped_path: str,
    tfam_path: str,
    species: str,
    organism: str | None = None,
    *,
    window: int = 200,
    max_candidates: int = 5000,
    max_markers_scanned: int | None = 20000,
    triage_window_bp: int = 100000,
    triage_max_checked: int = 2000,
    max_investigated: int = 15,
    max_iterations_per_variant: int = 8,
    client=None,
) -> tuple[GenomeInput, AnalysisRun, CandidateScanStats, TriageStats]:
    """Defaults as of the tile-cache fix (2026-08-15): scanning 20,000
    markers now costs ~255s (was the bottleneck at ~0.56s/marker with no
    caching — 20,000 markers would've taken over 3 hours before). That's
    enough to reach multiple chromosomes, not just an early slice of chr1.
    """
    genome_input = validate_genome_input(tped_path, tfam_path, species, organism)
    assembly = resolve_reference_assembly(genome_input.organism)

    candidates, scan_stats = find_candidate_variants(
        tped_path,
        genome_input.species,
        genome_input.organism,
        genome_input.individual.individual_id,
        window=window,
        max_candidates=max_candidates,
        max_markers_scanned=max_markers_scanned,
        assembly=assembly,
    )

    survivors, triage_stats = triage_candidates_near_genes(
        candidates,
        genome_input.organism,
        window_bp=triage_window_bp,
        max_investigated=max_investigated,
        max_checked=triage_max_checked,
    )

    findings: list[RiskFinding] = []
    for variant in survivors:
        finding = investigate_variant(
            variant,
            genome_input.organism,
            assembly=assembly,
            client=client,
            max_iterations=max_iterations_per_variant,
        )
        findings.append(finding)

    run = AnalysisRun(
        sample_id=genome_input.individual.individual_id,
        species=genome_input.species,
        input_source="tped_tfam",
        variants_analyzed=scan_stats.markers_scanned,
        findings=findings,
        scan_stats=asdict(scan_stats),
        triage_stats=asdict(triage_stats),
    )
    return genome_input, run, scan_stats, triage_stats


def write_json(run: AnalysisRun, path: str) -> None:
    Path(path).write_text(json.dumps(run.to_dict(), indent=2, default=str))


def write_csv(run: AnalysisRun, path: str) -> None:
    rows = run.to_csv_rows()
    if not rows:
        Path(path).write_text("")
        return
    with open(path, "w", newline="") as f:
        writer = csv_module.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    import os
    import sys

    from dotenv import load_dotenv

    load_dotenv()
    REPO_ROOT = Path(__file__).resolve().parent.parent.parent
    tped = os.environ.get("WILDTYPE_TPED", str(REPO_ROOT / "data/ollie_embark/raw/31211050316783.tped"))
    tfam = os.environ.get("WILDTYPE_TFAM", str(REPO_ROOT / "data/ollie_embark/raw/31211050316783.tfam"))
    species = sys.argv[1] if len(sys.argv) > 1 else "dog"
    organism = sys.argv[2] if len(sys.argv) > 2 else "Canis lupus familiaris"
    max_candidates = int(sys.argv[3]) if len(sys.argv) > 3 else 150
    max_scanned = int(sys.argv[4]) if len(sys.argv) > 4 else 5000
    max_investigated = int(sys.argv[5]) if len(sys.argv) > 5 else 15

    print(f"=== raw-genome pipeline (iterative agent): {species} ({organism}) ===")
    print(f"tped={tped}\ntfam={tfam}")
    print(f"max_candidates={max_candidates} max_markers_scanned={max_scanned} max_investigated={max_investigated}\n")

    genome_input, run, scan_stats, triage_stats = run_genome_pipeline(
        tped, tfam, species, organism,
        max_candidates=max_candidates, max_markers_scanned=max_scanned, max_investigated=max_investigated,
    )
    print(f"Individual: {genome_input.individual}")
    print(f"Total markers in file: {genome_input.marker_count}")
    print(f"Scan stats: {scan_stats}")
    print(f"Triage stats: {triage_stats}\n")

    for f in run.findings:
        print("-" * 72)
        print(json.dumps(f.to_dict(), indent=2, default=str))

    out_json = "/tmp/wildtype_genome_findings.json"
    out_csv = "/tmp/wildtype_genome_findings.csv"
    write_json(run, out_json)
    write_csv(run, out_csv)
    print(f"\nWrote {out_json} and {out_csv} ({len(run.findings)} findings)")
