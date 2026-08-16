"""Common Variant Representation — PROJECT_SPEC.md section 7.

The one schema every input format (VCF, TPED/TFAM, future formats) gets
normalized into before anything downstream (the agent, foundation models,
risk assessment) touches it. Nothing past this module should know or care
whether the data came from a VCF or a PLINK TPED/TFAM pair — that's the
whole point (spec section 3: "the architecture must not be specific to
Embark" — generalized here to "not specific to any one input format").

Deliberately NOT over-engineered (spec section 7: "Do not over-engineer the
common representation in the first implementation") — six fields, matches
the spec's own example schema field-for-field.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class CommonVariant:
    sample_id: str
    species: str
    chromosome: str
    position: int
    reference: str
    alternate: str
    genotype: str  # "0/1" VCF-style, or "REF/ALT"-style allele pair — see genotype_alleles()

    def genotype_alleles(self) -> tuple[str, str]:
        """Resolve `genotype` to an (allele1, allele2) pair of actual bases,
        regardless of whether it's VCF-style ("0/1" indexing into
        reference/alternate) or already spelled out as bases ("A/G")."""
        parts = self.genotype.replace("|", "/").split("/")
        if len(parts) != 2:
            return (self.reference, self.reference)
        resolved = []
        for p in parts:
            if p == "0":
                resolved.append(self.reference)
            elif p == "1":
                resolved.append(self.alternate)
            elif p == ".":
                resolved.append("N")
            else:
                resolved.append(p)  # already a base (A/C/G/T)
        return (resolved[0], resolved[1])

    @property
    def is_no_call(self) -> bool:
        return "." in self.genotype or "N" in self.genotype_alleles()

    @property
    def deviates_from_reference(self) -> bool:
        a1, a2 = self.genotype_alleles()
        return not self.is_no_call and (a1 != self.reference or a2 != self.reference)

    @property
    def zygosity(self) -> str:
        if self.is_no_call:
            return "no_call"
        a1, a2 = self.genotype_alleles()
        r = self.reference
        if a1 == r and a2 == r:
            return "homozygous_ref"
        if a1 != r and a2 != r:
            return "homozygous_alt" if a1 == a2 else "compound_het"
        return "heterozygous"
