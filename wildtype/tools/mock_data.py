"""Canned-but-plausible responses for every tool interface, so the agent
loop, prompts, and UI can be built and demoed tonight with zero external
calls. Scores are seeded from the gene name so a given variant returns the
same numbers on every run — deterministic enough to write tests against.

None of these numbers mean anything biologically. They exist to unblock
orchestration/prompt/UI work, not to be quoted in the actual demo.
"""
from __future__ import annotations

import hashlib

from wildtype.tools.base import (
    ESM2Score,
    Evo1Score,
    FunctionalEmbedding,
    LiteratureHit,
    StructurePrediction,
    TMAlignResult,
)

# Genes we already know are "supposed to" score as disruptive, so the mock
# pipeline demos coherently against Ollie's real validation set even before
# any real model has run. See docs/variant_routing.md for why these four.
_KNOWN_DISRUPTIVE = {"gpt", "prcd", "fgf4 retrogene - cfa12", "fgf4"}


def _seed(*parts: str) -> float:
    """Stable 0.0-1.0 float from the input strings."""
    h = hashlib.sha256("|".join(parts).encode()).hexdigest()
    return int(h[:8], 16) / 0xFFFFFFFF


class MockScoringClient:
    def score_missense(self, gene: str, wt_seq: str, position: int, mut_aa: str) -> ESM2Score:
        r = _seed(gene, str(position), mut_aa)
        disruptive = gene.lower() in _KNOWN_DISRUPTIVE
        llr = -(2.0 + r * 4.0) if disruptive else -(0.1 + r * 0.8)
        return ESM2Score(
            gene=gene,
            position=position,
            wt_aa=wt_seq[position - 1] if 0 < position <= len(wt_seq) else "X",
            mut_aa=mut_aa,
            log_likelihood_ratio=round(llr, 3),
            percentile=round((1 - r) * 15 if disruptive else 40 + r * 50, 1),
        )

    def functional_embedding(self, gene: str, wt_seq: str, position: int, mut_aa: str) -> FunctionalEmbedding:
        disruptive = gene.lower() in _KNOWN_DISRUPTIVE
        return FunctionalEmbedding(
            gene=gene,
            summary=(
                f"[MOCK] embeds near a disruptive functional cluster for {gene}"
                if disruptive
                else f"[MOCK] embeds within normal functional variation for {gene}"
            ),
            disrupted_function="catalytic/structural domain (mock)" if disruptive else None,
        )


class MockStructureClient:
    def fold(self, gene: str, sequence: str, high_accuracy: bool = False) -> StructurePrediction:
        r = _seed(gene, "fold", str(high_accuracy))
        plddt = 60 + r * 30
        return StructurePrediction(
            gene=gene,
            plddt_mean=round(plddt, 1),
            plddt_per_residue=None,
            structure_ref=f"mock://structures/{gene.lower().replace(' ', '_')}.pdb",
            model="alphafold3-mock" if high_accuracy else "esmfold-mock",
        )

    def align(self, structure_a: str, structure_b: str) -> TMAlignResult:
        r = _seed(structure_a, structure_b)
        return TMAlignResult(rmsd=round(1.0 + r * 4.0, 2), tm_score=round(0.5 + (1 - r) * 0.4, 3))


class MockDNAClient:
    def score_structural(self, locus: str, context_seq: str, variant_desc: str) -> Evo1Score:
        r = _seed(locus, variant_desc)
        deviation = 2.5 + r * 3.0
        return Evo1Score(
            locus=locus,
            context_deviation=round(deviation, 2),
            interpretation=f"[MOCK] insertion context at {locus} deviates {deviation:.1f} SD from baseline",
        )


class MockLiteratureClient:
    def search(self, query: str, max_results: int = 5) -> list[LiteratureHit]:
        return [
            LiteratureHit(
                title=f"[MOCK] Literature placeholder {i + 1} for: {query}",
                doi=f"10.mock/{abs(hash((query, i))) % 100000}",
                year=2019 + (i % 6),
                summary="[MOCK] replace with a real Paperclip search result tomorrow.",
                source="paperclip-mock",
            )
            for i in range(min(max_results, 3))
        ]
