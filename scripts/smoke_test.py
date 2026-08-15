"""End-to-end pipeline run against Ollie's real Embark export.

Runs in whatever WILDTYPE_MODE is set (defaults to mock, so this works
tonight with zero external calls or keys). Swap to `local` to exercise the
real CPU ESM2 fallback, or `proto` once sponsor access is live tomorrow.

    python scripts/smoke_test.py
    WILDTYPE_MODE=local python scripts/smoke_test.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from dotenv import load_dotenv

load_dotenv()

from wildtype.agent.loop import maybe_escalate_to_af3, render_report, run_pipeline
from wildtype.parsers.embark_csv import parse_report
from wildtype.tools.base import get_toolset


def main():
    csv_path = os.environ.get(
        "WILDTYPE_EMBARK_CSV",
        str(REPO_ROOT / "data" / "ollie_embark" / "Embark-full-results-data-31211050316783.csv"),
    )
    mode = os.environ.get("WILDTYPE_MODE", "mock")
    print(f"=== WildType smoke test (mode={mode}) ===\n")

    report = parse_report(csv_path)
    print(f"Subject: {report.name} ({report.breed_name}, {report.sex})")
    print(f"Flagged variants: {len(report.flagged)}\n")

    toolset = get_toolset(mode)
    # ABCB1/MDR1 is the "confirmed clear" validation case — included even
    # though it's not in `flagged`, to exercise the true-negative path too.
    analyses = run_pipeline(report, toolset, also_include_genes=["ABCB1"])

    for analysis in analyses:
        analysis = maybe_escalate_to_af3(analysis, toolset)
        print("-" * 72)
        print(render_report(analysis))
        print()

    print("-" * 72)
    print(f"\n{len(analyses)} variant(s) processed end-to-end without error.")


if __name__ == "__main__":
    main()
