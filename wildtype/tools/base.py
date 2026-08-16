"""Shared result dataclasses for the tool clients in this package.

2026-08-15 — this module used to also own a `Toolset`/`get_toolset()`
mode-switching factory (mock/local/proto) for the pre-pivot deterministic
Embark-CSV controller (agent/loop.py). That whole controller — loop.py,
species_mode.py, prompts.py, parsers/embark_csv.py, the Gradio UI, plus
the mock/local-only clients it alone depended on (esm_local.py,
mock_data.py, placeholder_sequences.py) — was removed once
genome_pipeline.py + iterative_agent.py (the real Claude tool-calling
agent, PROJECT_SPEC.md section 5.1) fully superseded it; those never
imported anything from Toolset/get_toolset. What's left here is just the
dataclasses proto_client.py and friends still construct directly, no
central mode-switching layer needed anymore — each tool decides for
itself how to fetch real data (Modal, NCBI, UniProt, Biohub) and reports
failure by raising/returning an error, not by falling back to a mock mode.
"""
from __future__ import annotations

from dataclasses import dataclass


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


@dataclass
class SequenceResult:
    sequence: str
    is_real: bool  # True = real UniProt fetch, False = deterministic placeholder
    source: str  # "uniprot" | "placeholder"
