"""Parser for PLINK .tped/.tfam pairs — the "beyond-Embark" raw marker scan.

.tfam: one row per individual — family_id, individual_id, paternal_id,
maternal_id, sex, phenotype.

.tped: one row per SNP marker — chromosome, marker_id, genetic_distance(cM,
usually 0 for array data), position(bp), then two allele columns per
individual (this file has one individual, so 6 columns total).

Ollie's file: 229,988 markers, single individual, 99.8% call rate per the
proposal. Nothing here needs sponsor compute — it's pure text parsing.

The actual "beyond-Embark" value comes from cross-referencing marker
positions against a table of known-disease-gene coordinates (OMIA, breed
panels, etc.) — that table (`KNOWN_DISEASE_POSITIONS`) is a stub tonight.
Fill it in tomorrow once Paperclip/OMIA lookups are live; the scan logic
itself is ready to run the moment the table has real entries.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, NamedTuple


class Individual(NamedTuple):
    family_id: str
    individual_id: str
    paternal_id: str
    maternal_id: str
    sex: str
    phenotype: str


@dataclass
class Marker:
    chrom: str
    marker_id: str
    genetic_dist: str
    position: int
    allele1: str
    allele2: str

    @property
    def is_no_call(self) -> bool:
        return self.allele1 == "0" or self.allele2 == "0"

    @property
    def is_heterozygous(self) -> bool:
        return not self.is_no_call and self.allele1 != self.allele2


# chrom -> [(start_bp, end_bp, gene_symbol, disease_note)]
# STUB — populate from OMIA / Paperclip lookups once live. Keeping the shape
# fixed now means the scan function below doesn't need to change tomorrow.
KNOWN_DISEASE_POSITIONS: dict[str, list[tuple[int, int, str, str]]] = {
    # "12": [(...FGF4 CFA12 coordinates..., "FGF4", "IVDD/chondrodystrophy")],
}


def load_tfam(path: str | Path) -> list[Individual]:
    individuals = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            individuals.append(Individual(*parts[:6]))
    return individuals


def iter_tped(path: str | Path) -> Iterator[Marker]:
    """Streams markers rather than loading all 230k into memory at once —
    matters once this runs against real venue laptops mid-demo."""
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            chrom, marker_id, genetic_dist, position = parts[:4]
            allele1, allele2 = parts[4], parts[5]
            yield Marker(chrom, marker_id, genetic_dist, int(position), allele1, allele2)


def call_rate(path: str | Path) -> float:
    total, called = 0, 0
    for m in iter_tped(path):
        total += 1
        called += 0 if m.is_no_call else 1
    return called / total if total else 0.0


def scan_known_positions(
    tped_path: str | Path,
    positions: dict[str, list[tuple[int, int, str, str]]] | None = None,
) -> list[dict]:
    """Cross-reference every marker against KNOWN_DISEASE_POSITIONS.
    Returns hits — markers landing inside a known disease-gene window,
    regardless of genotype, so downstream logic decides risk/carrier/clear.
    """
    table = positions if positions is not None else KNOWN_DISEASE_POSITIONS
    hits = []
    for m in iter_tped(tped_path):
        for start, end, gene, note in table.get(m.chrom, []):
            if start <= m.position <= end:
                hits.append(
                    {
                        "marker_id": m.marker_id,
                        "chrom": m.chrom,
                        "position": m.position,
                        "gene": gene,
                        "note": note,
                        "genotype": f"{m.allele1}/{m.allele2}",
                        "heterozygous": m.is_heterozygous,
                    }
                )
    return hits


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print("usage: python -m wildtype.parsers.tped <path.tped> <path.tfam>")
        raise SystemExit(1)
    tped_path, tfam_path = sys.argv[1], sys.argv[2]
    individuals = load_tfam(tfam_path)
    print(f"individuals: {individuals}")
    print("scanning markers (this streams ~230k lines, a few seconds)...")
    rate = call_rate(tped_path)
    print(f"call rate: {rate:.4%}")
    hits = scan_known_positions(tped_path)
    print(f"known-position hits: {len(hits)} (table is empty until filled in tomorrow)")
