"""Parser for Embark's "full results data" CSV export.

The file is a long-format table: one row per (category, name, value) triple.
Categories seen in Ollie's export: ID, Breed Name, Breed mix, Genetic Stats,
Lineage, Health, Trait. This module cares about `Health` rows — those are
the variant calls the rest of the pipeline routes and scores.

Usage:
    from wildtype.parsers.embark_csv import load_embark_csv, extract_variants
    df = load_embark_csv("path/to/export.csv")
    variants = extract_variants(df)
    flagged = [v for v in variants if v.risk_level != "clear"]
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

RiskLevel = Literal["at_risk", "carrier", "clear", "unknown"]
VariantClass = Literal["structural", "missense_snp"]

# Names whose "value" is a number/measurement, not a genotype call — keep them
# out of the variant list so downstream code doesn't try to route/score them.
_NON_VARIANT_HEALTH_ROWS = {"inbreeding", "immune response 1", "immune response 2"}

# Keywords in the variant name that mean "this is a DNA-level structural
# change" (insertion/retrogene/etc.) rather than a scoreable missense SNP —
# routes it to the Evo1 path instead of ESM2. See docs/variant_routing.md.
_STRUCTURAL_HINTS = ("retrogene", "insertion", "deletion", "duplication", "cnv")

_GENE_PAREN_RE = re.compile(r"\(([^()]*)\)\s*$")


@dataclass
class Variant:
    name: str
    raw_value: str
    risk_level: RiskLevel
    zygosity: str
    gene_hint: str
    variant_class: VariantClass
    category: str = "Health"

    def is_flagged(self) -> bool:
        return self.risk_level in ("at_risk", "carrier")


@dataclass
class EmbarkReport:
    subject_id: str = ""
    name: str = ""
    sex: str = ""
    breed_name: str = ""
    variants: list[Variant] = field(default_factory=list)

    @property
    def flagged(self) -> list[Variant]:
        return [v for v in self.variants if v.is_flagged()]


def load_embark_csv(path: str | Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _parse_risk(value: str) -> tuple[RiskLevel, str]:
    """'at risk, heterozygote codominant' -> ('at_risk', 'heterozygote codominant')"""
    value = value.strip()
    if not value:
        return "unknown", ""
    lowered = value.lower()
    if lowered.startswith("clear"):
        return "clear", value[len("clear"):].strip(", ")
    if lowered.startswith("at risk"):
        return "at_risk", value[len("at risk"):].strip(", ")
    if lowered.startswith("carrier"):
        return "carrier", value[len("carrier"):].strip(", ")
    return "unknown", value


def _gene_hint(name: str) -> str:
    m = _GENE_PAREN_RE.search(name)
    return m.group(1).strip() if m else name.strip()


def _variant_class(name: str, gene_hint: str) -> VariantClass:
    text = f"{name} {gene_hint}".lower()
    if any(h in text for h in _STRUCTURAL_HINTS):
        return "structural"
    return "missense_snp"


def extract_variants(rows: list[dict[str, str]]) -> list[Variant]:
    variants = []
    for row in rows:
        if row.get("category") != "Health":
            continue
        name = row["name"].strip()
        if name.lower().strip() in _NON_VARIANT_HEALTH_ROWS:
            continue
        risk_level, zygosity = _parse_risk(row["value"])
        gene_hint = _gene_hint(name)
        variants.append(
            Variant(
                name=name,
                raw_value=row["value"].strip(),
                risk_level=risk_level,
                zygosity=zygosity,
                gene_hint=gene_hint,
                variant_class=_variant_class(name, gene_hint),
            )
        )
    return variants


def parse_report(path: str | Path) -> EmbarkReport:
    rows = load_embark_csv(path)
    report = EmbarkReport(variants=extract_variants(rows))
    for row in rows:
        if row.get("category") != "ID":
            continue
        n = row["name"].strip().lower()
        if n == "swab code":
            report.subject_id = row["value"].strip()
        elif n == "name":
            report.name = row["value"].strip()
        elif n == "sex":
            report.sex = row["value"].strip()
    for row in rows:
        if row.get("category") == "Breed Name" and row["name"].strip().lower() == "name":
            report.breed_name = row["value"].strip()
    return report


if __name__ == "__main__":
    import sys

    path = sys.argv[1] if len(sys.argv) > 1 else None
    if not path:
        print("usage: python -m wildtype.parsers.embark_csv <path-to-embark-csv>")
        raise SystemExit(1)
    report = parse_report(path)
    print(f"{report.name} ({report.breed_name}, {report.sex}) — {report.subject_id}")
    print(f"{len(report.variants)} health rows parsed, {len(report.flagged)} flagged\n")
    for v in report.flagged:
        print(f"  [{v.risk_level:8s}] {v.gene_hint:12s} {v.variant_class:13s} — {v.name}")
