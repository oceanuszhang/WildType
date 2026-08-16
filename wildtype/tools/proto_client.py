"""Proto client — programmatic access to ESMFold, AlphaFold3, Evo1,
ProteinMPNN, TM-align, UniProt via the `proto_tools` package
(https://proto.evodesign.org), Arc Institute / Stanford Laboratory of
Evolutionary Design.

STATUS as of 2026-08-15 (post Modal cheap-tier unblock):
  - No PyPI release — installed from git, see requirements.txt. No hosted
    HTTPS/MCP server yet either; this calls proto_tools functions directly
    in-process (`run_uniprot_fetch(...)`, `run_esm2_score(...)`, etc.),
    each taking a typed `*Input` + `*Config` pair from proto_tools.tools.*.
  - `fetch_uniprot_sequence()` is REAL and CONFIRMED — live-tested against
    Ollie's actual FGF4 gene (UniProt accession J9P3I7, 206 aa, real
    sequence returned), device="cpu", no GPU/Modal needed. Working
    replacement for wildtype/tools/placeholder_sequences.py.
  - GPU access was blocked on Modal requiring a payment method for
    H100/H200/A100-80GB (proto_tools' GPU_DEFAULT, its own default choice
    for a new service — not a stated hard requirement). Unblocked without
    adding a payment method: patched the *installed* proto_tools package
    (scripts/patch_proto_gpu_tier.py) so esm2/esmfold/evo1 request
    GPU_BASIC (T4/L4/A10) instead — Modal doesn't gate that tier, and
    GPU_BASIC is already a real pattern other proto_tools services use
    (esm_if1, ablang, metal3d, malinois, splice_transformer, pangolin).
    Deployed and smoke-tested 2026-08-15: `esm2-score` real call in 0.6s.
    Re-run the patch script after any `pip install`/upgrade of proto-tools
    — it doesn't survive a reinstall, it edits site-packages directly.
  - Evo2 is NOT part of this unblock and never will be via this route —
    its service file requires GPU_HOPPER (H100/H200) specifically for
    Hopper FP8 support; T4/L4/A10 are architecturally incapable of running
    it, not just untried. See wildtype/agent/genome_pipeline.py.
  - AlphaFold3 and TM-align are not yet deployed/patched — same GPU_DEFAULT
    situation ESM2/ESMFold/Evo1 were in before this fix; extend
    scripts/patch_proto_gpu_tier.py's TARGETS + `proto-tools deploy` if
    needed.
"""
from __future__ import annotations

from wildtype.tools.base import (
    ESM2Score,
    Evo1Score,
    FunctionalEmbedding,
    StructurePrediction,
    TMAlignResult,
)


class ProtoConfigError(RuntimeError):
    pass


def fetch_uniprot_sequence(gene: str, organism: str = "Canis lupus familiaris") -> str | None:
    """Real UniProt fetch — CPU-only, no GPU/Modal needed. Confirmed live
    2026-08-15 against Ollie's FGF4 (UniProt J9P3I7). Returns None if no
    confident match (caller should fall back to placeholder_sequences).

    `gene` is Embark's parsed gene_hint, which isn't always a clean gene
    symbol — e.g. Ollie's PRCD variant hint is "PRCD Exon 1", not "PRCD".
    proto_tools' run_uniprot_fetch raises a hard ValueError on no match
    rather than a success=False output (hit live 2026-08-15, crashed the
    whole pipeline run) — caught here so one ugly gene_hint degrades this
    one variant to placeholder instead of taking down every variant.
    """
    from proto_tools import run_uniprot_fetch
    from proto_tools.tools.database_retrieval.uniprot.uniprot_fetch import (
        UniProtFetchConfig,
        UniProtFetchInput,
    )

    try:
        out = run_uniprot_fetch(
            UniProtFetchInput(target_name=gene, organism=organism),
            UniProtFetchConfig(device="cpu"),
        )
    except Exception:
        return None
    if not out.success or not out.sequence:
        return None
    return out.sequence


class ProtoSequenceClient:
    """SequenceClient (wildtype/tools/base.py) backed by real UniProt fetch.
    Falls back to the deterministic placeholder if UniProt has no confident
    match for the gene/organism — agent/loop.py checks `is_real` and notes
    it in the report rather than presenting a placeholder as a real finding.
    """

    def fetch(self, gene: str, organism: str = "Canis lupus familiaris"):
        from wildtype.tools.base import SequenceResult
        from wildtype.tools.placeholder_sequences import get_placeholder_sequence

        seq = fetch_uniprot_sequence(gene, organism)
        if seq:
            return SequenceResult(sequence=seq, is_real=True, source="uniprot")
        return SequenceResult(sequence=get_placeholder_sequence(gene), is_real=False, source="placeholder")


class ProtoScoringClient:
    """ESM2 missense scoring + ESM3 functional embedding via proto_tools.

    score_missense: REAL and CONFIRMED via Modal (device="modal") — deployed
    on the cheap GPU tier (see module docstring), smoke-tested 2026-08-15:
    a real esm2-score call completed in 0.6s. Falls back to device="cpu"
    windowed scoring only if the Modal call itself fails (e.g. deploy
    drifted, network issue) — tagged distinctly in `.model` so a CPU
    fallback is never confused with the real GPU score.
    """

    # ESM2 in proto_tools caps input at 1022 residues regardless of device
    # (hit live testing this against ABCB1/MDR1, a 1321aa transporter).
    # Window around the mutated position rather than truncate from one end,
    # so the model still sees local context on both sides of the mutation.
    _MAX_ESM2_RESIDUES = 1022  # full cap now that Modal (not CPU) does the work
    _CPU_FALLBACK_MAX_RESIDUES = 220  # only applies if Modal itself fails

    def score_missense(self, gene: str, wt_seq: str, position: int, mut_aa: str) -> ESM2Score:
        from proto_tools import run_esm2_score
        from proto_tools.tools.masked_models.esm2.esm2_score import (
            ESM2ScoringConfig,
            ESM2ScoringInput,
        )

        window_seq, window_position = self._window(wt_seq, position, self._MAX_ESM2_RESIDUES)
        mut_seq = window_seq[: window_position - 1] + mut_aa + window_seq[window_position:]

        try:
            out = run_esm2_score(
                ESM2ScoringInput(sequences=[window_seq, mut_seq]),
                ESM2ScoringConfig(device="modal"),  # real weights, full checkpoint default
            )
            model_tag = "esm2-proto-modal"
        except Exception:
            # Modal deploy drifted or network hiccup — fall back to CPU with
            # a smaller window rather than fail the whole variant.
            window_seq, window_position = self._window(wt_seq, position, self._CPU_FALLBACK_MAX_RESIDUES)
            mut_seq = window_seq[: window_position - 1] + mut_aa + window_seq[window_position:]
            out = run_esm2_score(
                ESM2ScoringInput(sequences=[window_seq, mut_seq]),
                ESM2ScoringConfig(device="cpu", model_checkpoint="esm2_t12_35M_UR50D"),
            )
            model_tag = "esm2-proto-cpu-fallback"

        wt_score, mut_score = out.scores
        llr = mut_score.avg_log_likelihood - wt_score.avg_log_likelihood
        return ESM2Score(
            gene=gene,
            position=position,
            wt_aa=wt_seq[position - 1] if 0 < position <= len(wt_seq) else "X",
            mut_aa=mut_aa,
            log_likelihood_ratio=round(llr, 3),
            model=model_tag,
        )

    def _window(self, seq: str, position: int, max_residues: int) -> tuple[str, int]:
        """Slice `seq` down to `max_residues` around `position` (1-based),
        centered where possible. Returns (windowed_seq, new_1based_position).
        No-op if seq already fits."""
        if len(seq) <= max_residues:
            return seq, position
        half = max_residues // 2
        start = max(0, (position - 1) - half)
        end = min(len(seq), start + max_residues)
        start = max(0, end - max_residues)  # re-clamp if end hit the seq boundary
        return seq[start:end], position - start

    def functional_embedding(self, gene: str, wt_seq: str, position: int, mut_aa: str) -> FunctionalEmbedding:
        raise NotImplementedError(
            "proto_tools.run_esm3_score/run_esm3_embeddings exist, but which "
            "embedding-space distance counts as 'disruptive' is a design call, "
            "not wired yet — confirm approach at check-in."
        )


class ProtoStructureClient:
    """ESMFold structure prediction (+ AlphaFold3, TM-align — not yet
    wired) via proto_tools.

    fold(): REAL and CONFIRMED via Modal, deployed on the cheap GPU tier
    (see module docstring). Live-tested 2026-08-15: real esmfold-prediction
    call returned avg_plddt=0.599, ptm=0.311, avg_pae=19.5 for a 178aa test
    sequence — note avg_plddt comes back 0-1 scaled from proto_tools, but
    the rest of this codebase (MockStructureClient, prompts.py) treats
    plddt_mean as 0-100 — converted here (`* 100`) to stay consistent, not
    left for callers to silently misinterpret.

    AlphaFold3 (`high_accuracy=True`) is NOT deployed/patched yet — falls
    back to a real ESMFold call instead of AlphaFold3 when escalated,
    tagged so it's not confused with an actual AF3 result.
    """

    def fold(self, gene: str, sequence: str, high_accuracy: bool = False) -> StructurePrediction:
        if high_accuracy:
            # AlphaFold3 isn't deployed — no attempt made, straight to a
            # tagged real-ESMFold-standing-in-for-AF3 result rather than a
            # silent wrong-tool substitution.
            result = self._real_fold(gene, sequence)
            if result is not None:
                result.model = f"{result.model}+af3-not-deployed-used-esmfold"
                return result
        else:
            result = self._real_fold(gene, sequence)
            if result is not None:
                return result

        from wildtype.tools.mock_data import MockStructureClient

        result = MockStructureClient().fold(gene, sequence, high_accuracy)
        result.model = f"{result.model}+proto-call-failed"
        return result

    def _real_fold(self, gene: str, sequence: str) -> StructurePrediction | None:
        from proto_tools import run_esmfold
        from proto_tools.entities.complex import Chain, Complex
        from proto_tools.tools.structure_prediction.esmfold.esmfold import ESMFoldConfig, ESMFoldInput

        try:
            out = run_esmfold(
                ESMFoldInput(complexes=[Complex(chains=[Chain(id="A", sequence=sequence)])]),
                ESMFoldConfig(device="modal"),
            )
        except Exception:
            return None
        if not out.success or not out.structures:
            return None
        s = out.structures[0]
        plddt = getattr(s.metrics, "avg_plddt", None)
        if plddt is None:
            return None
        return StructurePrediction(
            gene=gene,
            plddt_mean=round(plddt * 100, 1),
            plddt_per_residue=None,
            structure_ref=s.source or "esmfold-proto-modal",
            model="esmfold-proto-modal",
        )

    def align(self, structure_a: str, structure_b: str) -> TMAlignResult:
        from wildtype.tools.mock_data import MockStructureClient

        return MockStructureClient().align(structure_a, structure_b)


class ProtoDNAClient:
    """DNA-level structural-variant scoring — deliberately NOT wired to
    Evo1, even though Evo1 is now cheap-tier deployable (same GPU_BASIC
    unblock as ESM2/ESMFold, see module docstring).

    Evo1 was trained on OpenGenome — bacteria, archaea, and phage only
    (its own checkpoint list is the tell: evo-1-8k-crispr, evo-1-8k-
    transposon — CRISPR systems and prokaryotic transposons, not mammalian
    biology). Scoring Ollie's canine FGF4 retrogene insertion against it
    would mean measuring deviation from a prokaryotic-genome likelihood
    model — not "less accurate," just the wrong tool, and presenting that
    number as a real finding would be actively misleading, not merely
    imprecise. (Caught by the user mid-build 2026-08-15 — the original
    workflow diagram specified Evo2 for exactly this reason; using Evo1
    here was scope drift on my part, corrected.)

    Evo2 (OpenGenome2 — spans bacteria through eukaryotes, including
    mammals) is the scientifically correct tool and is what
    wildtype/agent/genome_pipeline.py actually calls, via
    iterative_agent.py's score_dna_variant_evo2 tool (real, deployed, and
    confirmed live 2026-08-15 — the Modal H100 payment gate that blocked
    it is now cleared; deploy took 269.9s, smoke test returned real
    avg_log_likelihood scores).

    This class is a SEPARATE, still-mock-only path — the legacy
    Embark-CSV pipeline (wildtype/agent/loop.py's analyze_variant, now
    superseded by genome_pipeline.py per that module's docstring). Not
    wired to real Evo2 even though Evo2 itself is unblocked, because
    loop.py's `context_seq` here is actually a UniProt PROTEIN sequence
    (from toolset.sequences.fetch), not DNA — the same species of mistake
    that got Evo1 caught and removed from this path originally. Wiring
    real Evo2 in would need a real DNA-sequence fetch added to loop.py
    first, not just swapping in a live model call against the wrong
    sequence type. Left as MockDNAClient, tagged in `.interpretation` so
    it's never mistaken for a real score.
    """

    def score_structural(self, locus: str, context_seq: str, variant_desc: str) -> Evo1Score:
        from wildtype.tools.mock_data import MockDNAClient

        result = MockDNAClient().score_structural(locus, context_seq, variant_desc)
        result.interpretation = (
            f"[NO EUKARYOTE-TRAINED DNA MODEL AVAILABLE — mock fallback. "
            f"Evo2 is the correct tool but needs Hopper GPU (payment-gated); "
            f"Evo1 exists on the cheap tier but is prokaryote-only, not used here] "
            f"{result.interpretation}"
        )
        return result


if __name__ == "__main__":
    import sys

    gene = sys.argv[1] if len(sys.argv) > 1 else "FGF4"
    seq = fetch_uniprot_sequence(gene)
    print(f"{gene}: {len(seq) if seq else 0} aa" + (f"\n{seq}" if seq else " (no match)"))
