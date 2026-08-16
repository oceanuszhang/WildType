"""Stratified, cross-chromosome demo scan, 2026-08-15 — correct assembly
(CanFam3.1) + round-robin investigation selection across chromosomes, so
a bounded investigation budget reaches deep into the genome (chr9/13/14,
where PRCD/GPT/ABCB1 actually live) instead of being exhausted entirely
within chromosome 1. No knowledge of Embark's panel used during the run —
same blind methodology as before, just no longer chromosome-1-biased.
"""
import os
import sys
import time
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from wildtype.agent.genome_pipeline import run_genome_pipeline, write_csv, write_json
from wildtype.tools.reference_genome import resolve_assembly_by_accession

tped = os.environ.get("WILDTYPE_TPED", str(REPO_ROOT / "data/ollie_embark/raw/31211050316783.tped"))
tfam = os.environ.get("WILDTYPE_TFAM", str(REPO_ROOT / "data/ollie_embark/raw/31211050316783.tfam"))

assembly = resolve_assembly_by_accession("GCF_000002285.3", "Canis lupus familiaris")
print(f"=== Stratified scan (assembly: {assembly.accession} {assembly.name}) ===")
print(f"tped={tped}\ntfam={tfam}")
print("scope: scan through ~chr14+ (150k markers), round-robin investigation across chromosomes, 30 investigations\n")

t0 = time.time()
genome_input, run, scan_stats, triage_stats = run_genome_pipeline(
    tped, tfam, "dog", "Canis lupus familiaris",
    max_candidates=100000,
    max_markers_scanned=150000,  # covers chr1 through chr20ish, well past chr14 (ABCB1)
    triage_max_checked=100000,   # cheap now (chromosome-gene cache) — need to see all chromosomes to stratify fairly
    max_investigated=30,
    assembly_override=assembly,
    stratify_by_chromosome=True,  # round-robin across chromosomes, no gene-name knowledge involved
    checkpoint_path="/tmp/wildtype_genome_findings.json",
)
elapsed = time.time() - t0

print(f"\nTotal elapsed: {round(elapsed / 60, 1)} minutes")
print(f"Scan stats: {scan_stats}")
print(f"Triage stats: {triage_stats}")
print(f"Findings: {len(run.findings)}")
chroms_seen = sorted(set(f.variants[0].chromosome for f in run.findings if f.variants), key=lambda x: int(x) if x.isdigit() else 999)
print(f"Chromosomes represented in findings: {chroms_seen}")
for f in run.findings:
    v = f.variants[0] if f.variants else None
    locus = f"{v.chromosome}:{v.position}" if v else "?"
    print(f" - {locus} | genes={f.genes} | {f.confidence.level}")

out_json = "/tmp/wildtype_genome_findings.json"
out_csv = "/tmp/wildtype_genome_findings.csv"
write_json(run, out_json)
write_csv(run, out_csv)
print(f"\nWrote {out_json} and {out_csv} ({len(run.findings)} findings)")

from wildtype.tools.embark_panel_genes import EMBARK_PANEL_GENES

matches = [f for f in run.findings if set(f.genes) & EMBARK_PANEL_GENES]
print(f"\nPost-hoc comparison against Embark's real {len(EMBARK_PANEL_GENES)}-gene panel "
      f"(not used during the run above): {len(matches)} of {len(run.findings)} findings land on a gene Embark also tests.")
for f in matches:
    print(f"  - {f.genes} | {f.confidence.level}")
