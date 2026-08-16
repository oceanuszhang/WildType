"""Tool interfaces the agent loop calls against — and the one place that
decides whether a call goes to a mock, local CPU inference, or the real
Proto/Paperclip services.

Design intent: agent/loop.py never imports MockScoringClient or
LocalESM2Client or a Proto SDK directly. It calls `tools.base.get_toolset()`
and gets back something satisfying these interfaces. Swapping WILDTYPE_MODE
in .env is the only change needed to go from tonight's offline dev to
tomorrow's live sponsor stack.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Protocol


# --- result types -----------------------------------------------------------


@dataclass
class ESM2Score:
    gene: str
    position: int
    wt_aa: str
    mut_aa: str
    log_likelihood_ratio: float  # negative = mutant less favorable than WT
    percentile: float | None = None  # filled in once we have a score distribution
    model: str = "esm2"


@dataclass
class FunctionalEmbedding:
    gene: str
    summary: str  # e.g. "embeds near catalytic-domain-disrupting cluster"
    disrupted_function: str | None
    model: str = "esm3"


@dataclass
class StructurePrediction:
    gene: str
    plddt_mean: float
    plddt_per_residue: list[float] | None
    structure_ref: str  # path/URL/id — not raw coordinates
    model: str = "esmfold"


@dataclass
class TMAlignResult:
    rmsd: float
    tm_score: float


@dataclass
class Evo1Score:
    locus: str
    context_deviation: float  # SD from local baseline distribution
    interpretation: str


@dataclass
class LiteratureHit:
    title: str
    doi: str | None
    year: int | None
    summary: str
    source: str  # "paperclip" | "uniprot" | "pdb" | ...


# --- interfaces --------------------------------------------------------------


class ScoringClient(Protocol):
    def score_missense(self, gene: str, wt_seq: str, position: int, mut_aa: str) -> ESM2Score: ...
    def functional_embedding(self, gene: str, wt_seq: str, position: int, mut_aa: str) -> FunctionalEmbedding: ...


class StructureClient(Protocol):
    def fold(self, gene: str, sequence: str, high_accuracy: bool = False) -> StructurePrediction: ...
    def align(self, structure_a: str, structure_b: str) -> TMAlignResult: ...


class DNAScoringClient(Protocol):
    def score_structural(self, locus: str, context_seq: str, variant_desc: str) -> Evo1Score: ...


class LiteratureClient(Protocol):
    def search(self, query: str, max_results: int = 5) -> list[LiteratureHit]: ...


@dataclass
class SequenceResult:
    sequence: str
    is_real: bool  # True = real UniProt fetch, False = deterministic placeholder
    source: str  # "uniprot" | "placeholder"


class SequenceClient(Protocol):
    def fetch(self, gene: str, organism: str = "Canis lupus familiaris") -> SequenceResult: ...


@dataclass
class Toolset:
    scoring: ScoringClient
    structure: StructureClient
    dna: DNAScoringClient
    literature: LiteratureClient
    sequences: SequenceClient
    mode: str


def get_toolset(mode: str | None = None) -> Toolset:
    mode = (mode or os.environ.get("WILDTYPE_MODE", "mock")).lower()

    if mode == "mock":
        from wildtype.tools.mock_data import (
            MockScoringClient,
            MockStructureClient,
            MockDNAClient,
            MockLiteratureClient,
        )
        from wildtype.tools.placeholder_sequences import PlaceholderSequenceClient

        return Toolset(
            scoring=MockScoringClient(),
            structure=MockStructureClient(),
            dna=MockDNAClient(),
            literature=MockLiteratureClient(),
            sequences=PlaceholderSequenceClient(),
            mode=mode,
        )

    if mode == "local":
        # Real ESM2 inference on CPU, no sponsor access required. Structure
        # prediction and DNA scoring have no offline-capable equivalent, so
        # those two still fall back to mocks even in "local" mode — that's
        # intentional, not a bug: ESM2 is the one model small enough to run
        # on a laptop tonight. Sequences stay placeholder too — real UniProt
        # fetch needs network access, which breaks "local" mode's offline
        # guarantee even though it needs no sponsor key.
        from wildtype.tools.esm_local import LocalESM2Client
        from wildtype.tools.mock_data import MockStructureClient, MockDNAClient, MockLiteratureClient
        from wildtype.tools.placeholder_sequences import PlaceholderSequenceClient

        return Toolset(
            scoring=LocalESM2Client(),
            structure=MockStructureClient(),
            dna=MockDNAClient(),
            literature=MockLiteratureClient(),
            sequences=PlaceholderSequenceClient(),
            mode=mode,
        )

    if mode == "proto":
        from wildtype.tools.proto_client import (
            ProtoScoringClient,
            ProtoStructureClient,
            ProtoDNAClient,
            ProtoSequenceClient,
        )
        from wildtype.tools.paperclip_client import PaperclipClient

        return Toolset(
            scoring=ProtoScoringClient(),
            structure=ProtoStructureClient(),
            dna=ProtoDNAClient(),
            literature=PaperclipClient(),
            sequences=ProtoSequenceClient(),
            mode=mode,
        )

    raise ValueError(f"unknown WILDTYPE_MODE={mode!r} (expected mock|local|proto)")
