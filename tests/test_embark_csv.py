"""Runs against Ollie's real Embark export — not a fixture — so it fails
loudly if the parser ever drifts from what Embark actually ships.
"""
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from wildtype.parsers.embark_csv import parse_report

CSV_PATH = os.environ.get(
    "WILDTYPE_EMBARK_CSV",
    str(REPO_ROOT / "data" / "ollie_embark" / "Embark-full-results-data-31211050316783.csv"),
)


def _report():
    assert Path(CSV_PATH).exists(), f"Ollie's Embark CSV not found at {CSV_PATH}"
    return parse_report(CSV_PATH)


def test_subject_metadata():
    report = _report()
    assert report.name == "Ollie"
    assert report.sex == "Female"


def test_flags_exactly_the_three_known_positives():
    report = _report()
    flagged = {v.gene_hint: v.risk_level for v in report.flagged}
    assert flagged.get("GPT") == "at_risk"
    assert flagged.get("PRCD Exon 1") == "carrier"
    assert flagged.get("FGF4 retrogene - CFA12") == "at_risk"
    assert len(report.flagged) == 3


def test_fgf4_routes_structural_not_missense():
    report = _report()
    fgf4 = next(v for v in report.flagged if "FGF4" in v.gene_hint)
    assert fgf4.variant_class == "structural"


def test_gpt_and_prcd_route_missense():
    report = _report()
    for gene in ("GPT", "PRCD Exon 1"):
        v = next(v for v in report.flagged if v.gene_hint == gene)
        assert v.variant_class == "missense_snp"


def test_mdr1_present_and_clear():
    report = _report()
    mdr1 = next(v for v in report.variants if v.gene_hint == "ABCB1")
    assert mdr1.risk_level == "clear"
