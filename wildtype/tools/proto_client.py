"""Proto client — the unified MCP interface to ESMFold, AlphaFold3, Evo1,
ProteinMPNN, TM-align, UniProt.

STATUS: shape-only stub. Proto's actual connection method (MCP server URL,
auth flow, per-tool call signatures) isn't known until check-in tomorrow —
that's the whole reason mock/local modes exist tonight. What's real here is
the *interface*: these classes already satisfy tools.base's Protocols, so
the agent loop and prompts built against mock mode tonight don't change
shape when this fills in tomorrow — only the method bodies do.

TODO tomorrow, in order:
  1. Get PROTO_API_KEY / PROTO_BASE_URL (or MCP connection details) at check-in.
  2. Fill in `_call_tool()` with the real MCP client/transport.
  3. Confirm each tool's exact request/response schema against Proto's docs
     and adjust the field mapping in each method below.
  4. Flip WILDTYPE_MODE=proto in .env and re-run scripts/smoke_test.py.
"""
from __future__ import annotations

import os

from wildtype.tools.base import (
    ESM2Score,
    Evo1Score,
    FunctionalEmbedding,
    StructurePrediction,
    TMAlignResult,
)


class ProtoConfigError(RuntimeError):
    pass


class _ProtoBase:
    def __init__(self):
        self.api_key = os.environ.get("PROTO_API_KEY")
        self.base_url = os.environ.get("PROTO_BASE_URL")

    def _require_config(self):
        if not self.api_key or not self.base_url:
            raise ProtoConfigError(
                "PROTO_API_KEY / PROTO_BASE_URL not set. Fill in .env after check-in, "
                "or run with WILDTYPE_MODE=mock / local for now."
            )

    def _call_tool(self, tool_name: str, **kwargs):
        """Placeholder for the real Proto MCP call. Replace with the actual
        MCP client once we have docs — keeping this one chokepoint means
        every method below only needs the transport written once.
        """
        self._require_config()
        raise NotImplementedError(
            f"Proto MCP transport not wired yet (attempted tool={tool_name!r}, args={kwargs}). "
            "Fill in _call_tool() once Proto connection details are available."
        )


class ProtoScoringClient(_ProtoBase):
    def score_missense(self, gene: str, wt_seq: str, position: int, mut_aa: str) -> ESM2Score:
        resp = self._call_tool("esm2_score", gene=gene, sequence=wt_seq, position=position, mut_aa=mut_aa)
        return ESM2Score(gene=gene, position=position, wt_aa=wt_seq[position - 1], mut_aa=mut_aa,
                          log_likelihood_ratio=resp["llr"], percentile=resp.get("percentile"))

    def functional_embedding(self, gene: str, wt_seq: str, position: int, mut_aa: str) -> FunctionalEmbedding:
        resp = self._call_tool("esm3_embed", gene=gene, sequence=wt_seq, position=position, mut_aa=mut_aa)
        return FunctionalEmbedding(gene=gene, summary=resp["summary"], disrupted_function=resp.get("disrupted_function"))


class ProtoStructureClient(_ProtoBase):
    def fold(self, gene: str, sequence: str, high_accuracy: bool = False) -> StructurePrediction:
        tool = "alphafold3_fold" if high_accuracy else "esmfold_fold"
        resp = self._call_tool(tool, gene=gene, sequence=sequence)
        return StructurePrediction(gene=gene, plddt_mean=resp["plddt_mean"], plddt_per_residue=resp.get("plddt_per_residue"),
                                    structure_ref=resp["structure_ref"], model=tool)

    def align(self, structure_a: str, structure_b: str) -> TMAlignResult:
        resp = self._call_tool("tm_align", structure_a=structure_a, structure_b=structure_b)
        return TMAlignResult(rmsd=resp["rmsd"], tm_score=resp["tm_score"])


class ProtoDNAClient(_ProtoBase):
    def score_structural(self, locus: str, context_seq: str, variant_desc: str) -> Evo1Score:
        resp = self._call_tool("evo1_score", locus=locus, sequence=context_seq, variant=variant_desc)
        return Evo1Score(locus=locus, context_deviation=resp["deviation_sd"], interpretation=resp["interpretation"])
