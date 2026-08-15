"""Real ESM2 inference on CPU — no Modal, no Proto, no sponsor keys.

This is the insurance policy: if Proto access is delayed or flaky tomorrow,
the protein-scoring half of the pipeline still produces genuine (not mocked)
numbers, using the small ESM2 checkpoint (35M params) that runs a forward
pass on a laptop CPU in well under a second for protein-length sequences.

Uses the standard "wt-marginal" variant-effect scoring: one forward pass on
the wildtype sequence, then compare the model's log-probability of the
mutant residue vs. the wildtype residue at the variant position. This is the
same scoring convention ESM2's own variant-effect papers use — swapping in
the real (larger) ESM2 via Proto tomorrow should rank variants similarly,
just with a better-calibrated model.

Setup (do this tonight — first call downloads ~150MB of weights):
    pip install fair-esm torch
    python -c "import esm; esm.pretrained.esm2_t12_35M_UR50D()"
"""
from __future__ import annotations

import os

# python.org's macOS builds don't ship a CA bundle, which makes the model
# download 401 on SSL verification (SSLCertVerificationError) the first time
# torch.hub reaches out — bit us once tonight, fixed here so it doesn't bite
# again on a fresh venv on a venue laptop.
if not os.environ.get("SSL_CERT_FILE"):
    try:
        import certifi

        os.environ["SSL_CERT_FILE"] = certifi.where()
    except ImportError:
        pass

from wildtype.tools.base import ESM2Score, FunctionalEmbedding

_DEFAULT_MODEL = "esm2_t12_35M_UR50D"


class LocalESM2Client:
    def __init__(self, model_name: str | None = None):
        self.model_name = model_name or os.environ.get("WILDTYPE_LOCAL_ESM_MODEL", _DEFAULT_MODEL)
        self._model = None
        self._alphabet = None
        self._batch_converter = None

    def _ensure_loaded(self):
        if self._model is not None:
            return
        import esm  # fair-esm
        import torch

        loader = getattr(esm.pretrained, self.model_name)
        self._model, self._alphabet = loader()
        self._model.eval()
        self._batch_converter = self._alphabet.get_batch_converter()
        self._torch = torch

    def score_missense(self, gene: str, wt_seq: str, position: int, mut_aa: str) -> ESM2Score:
        """position is 1-indexed into wt_seq."""
        self._ensure_loaded()
        if not (1 <= position <= len(wt_seq)):
            raise ValueError(f"position {position} out of range for {gene} (len {len(wt_seq)})")
        wt_aa = wt_seq[position - 1]

        data = [(gene, wt_seq)]
        _, _, tokens = self._batch_converter(data)
        with self._torch.no_grad():
            logits = self._model(tokens, repr_layers=[], return_contacts=False)["logits"]
        log_probs = self._torch.log_softmax(logits[0], dim=-1)

        token_idx = position  # index 0 is BOS, so residue `position` sits at this token index
        wt_idx = self._alphabet.get_idx(wt_aa)
        mut_idx = self._alphabet.get_idx(mut_aa)
        llr = (log_probs[token_idx, mut_idx] - log_probs[token_idx, wt_idx]).item()

        return ESM2Score(
            gene=gene,
            position=position,
            wt_aa=wt_aa,
            mut_aa=mut_aa,
            log_likelihood_ratio=round(llr, 3),
            percentile=None,  # needs a reference distribution — compute across a variant batch, not per-call
            model=self.model_name,
        )

    def functional_embedding(self, gene: str, wt_seq: str, position: int, mut_aa: str) -> FunctionalEmbedding:
        """No local ESM3 fallback — this is a documented stand-in, not a
        substitute. It reuses the ESM2 LLR as a coarse proxy so the pipeline
        has *something* to render tonight; swap for real Proto/ESM3 tomorrow.
        """
        score = self.score_missense(gene, wt_seq, position, mut_aa)
        disruptive = score.log_likelihood_ratio < -2.0
        return FunctionalEmbedding(
            gene=gene,
            summary=(
                f"[local ESM2 proxy, not ESM3] LLR={score.log_likelihood_ratio} — "
                + ("likely functionally disruptive" if disruptive else "within normal variation")
            ),
            disrupted_function=None,
            model="esm2-local-proxy",
        )
