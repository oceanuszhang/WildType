"""VCF parser — PROJECT_SPEC.md section 6.1.

Minimal, matches the spec's stated minimum fields exactly: sample ID,
chromosome, position, reference, alternate, genotype, variant ID when
available. Does not attempt QUAL/FILTER/INFO/FORMAT beyond what's needed to
locate the genotype column — spec explicitly says "do not assume every VCF
contains all optional fields," so this doesn't require them.

Single-sample VCFs only for now, matching the same one-individual
assumption wildtype/agent/genome_pipeline.py already made explicit for
TPED/TFAM (spec doesn't require multi-sample for the MVP).
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterator

from wildtype.parsers.common_variant import CommonVariant


class VCFParseError(RuntimeError):
    pass


def iter_vcf_variants(path: str | Path, species: str, sample_id: str | None = None) -> Iterator[CommonVariant]:
    """Streams variants rather than loading the whole file — same reasoning
    as tped.py's iter_tped: real VCFs can be large."""
    path = Path(path)
    sample_col_idx: int | None = None
    resolved_sample_id = sample_id

    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.lstrip("#").split("\t")
                # Standard VCF: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <sample...>
                if len(cols) < 10:
                    raise VCFParseError(
                        f"{path}: header has no sample genotype column (found {len(cols)} columns, "
                        "need at least 10 for a single-sample VCF: CHROM..FORMAT + 1 sample)"
                    )
                sample_names = cols[9:]
                if sample_id is not None:
                    if sample_id not in sample_names:
                        raise VCFParseError(f"{path}: sample_id={sample_id!r} not in VCF samples {sample_names}")
                    sample_col_idx = 9 + sample_names.index(sample_id)
                else:
                    if len(sample_names) > 1:
                        raise VCFParseError(
                            f"{path}: {len(sample_names)} samples ({sample_names}) but no sample_id given — "
                            "this parser assumes one sample per run (matches the TPED/TFAM path); pass sample_id."
                        )
                    sample_col_idx = 9
                    resolved_sample_id = sample_names[0]
                continue
            if sample_col_idx is None:
                raise VCFParseError(f"{path}: data row before #CHROM header line")

            fields = line.split("\t")
            if len(fields) <= sample_col_idx:
                continue  # malformed row — skip rather than crash the whole parse
            chrom, pos, _vid, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
            format_col = fields[8].split(":")
            sample_data = fields[sample_col_idx].split(":")
            if "GT" not in format_col:
                continue  # no genotype for this row — nothing to build a CommonVariant from
            gt = sample_data[format_col.index("GT")]

            # ALT can be comma-separated for multi-allelic sites; MVP takes
            # the first (spec doesn't ask for multi-allelic handling yet).
            alt_first = alt.split(",")[0]

            yield CommonVariant(
                sample_id=resolved_sample_id or "unknown",
                species=species,
                chromosome=chrom,
                position=int(pos),
                reference=ref,
                alternate=alt_first,
                genotype=gt,
            )


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print("usage: python -m wildtype.parsers.vcf <path.vcf> <species>")
        raise SystemExit(1)
    path, species = sys.argv[1], sys.argv[2]
    count = 0
    for v in iter_vcf_variants(path, species):
        count += 1
        if count <= 5:
            print(v)
    print(f"\n{count} variant(s) parsed")
