"""Runs INSIDE the isolated `.venv-biohub` environment, never the main
project venv. Called as a subprocess by wildtype/tools/biohub_client.py.

Why this file lives in its own venv: EvolutionaryScale's official `esm` SDK
(needed to call Biohub/Forge correctly — it must pre-tokenize sequences,
raw strings get rejected) claims the same top-level `esm` import path as
`fair-esm`, which wildtype/tools/esm_local.py already depends on for the
local CPU ESM2 fallback (WILDTYPE_MODE=local). Installing both in one venv
means whichever installs last silently wins and breaks the other. Confirmed
live 2026-08-15: `.venv-biohub` has EvolutionaryScale's `esm` v3.2.1.post1;
main `.venv` keeps `fair-esm` untouched.

Protocol: reads one JSON object from stdin, writes one JSON object (or
{"error": "..."}) to stdout. No other output on stdout — logs/warnings go
to stderr only, so the parent process's json.loads() on stdout never breaks.

Input shape (score_missense — the only operation implemented so far):
    {
        "operation": "score_missense",
        "sequence": "<wild-type amino acid sequence>",
        "position": <int, 1-based>,
        "mut_aa": "<single-letter mutant amino acid>",
        "model": "esmc-300m-2024-12"   # optional, this is the default
    }

Output shape:
    {
        "wt_aa": "M", "mut_aa": "K", "position": 42,
        "log_likelihood_ratio": -3.219,
        "model": "esmc-300m-2024-12-biohub",
        "window_start": 1  # 1-based offset of `sequence` within the full
                            # protein, for callers that windowed before
                            # calling in (position in the output is already
                            # relative to the windowed `sequence` passed in)
    }

Method: wildtype-marginal scoring (single forward pass over the wild-type
sequence, compare the model's log-probability of the mutant vs wild-type
amino acid at the mutated position, conditioned on the real surrounding
sequence). This is the same "wildtype marginals" method from the ESM1v
variant-effect-prediction paper — standard practice, not something
invented for this project. Different from proto_client.py's ESM2 approach
(two full forward passes, compare average log-likelihood of the whole
mutant vs whole wild-type sequence) — both are legitimate; this one only
needs a single forward pass so it's cheaper per call.

ESMC context cap: not documented anywhere I could find, and testing a
3000-residue sequence live (2026-08-15) succeeded, so don't trust a
specific number here. biohub_client.py windows long proteins the same way
proto_client.py does (same _MAX_RESIDUES cap) as a precaution, not because
the real ESMC limit is known to be lower.
"""
from __future__ import annotations

import json
import sys


def _fail(message: str) -> None:
    print(json.dumps({"error": message}), file=sys.stdout)
    sys.exit(0)  # exit 0 — the parent distinguishes success/failure by the
    # "error" key, not the exit code, so a stray non-JSON stderr line never
    # gets misread as this process's real output.


def main() -> None:
    raw = sys.stdin.read()
    try:
        req = json.loads(raw)
    except json.JSONDecodeError as exc:
        _fail(f"invalid JSON on stdin: {exc}")
        return

    operation = req.get("operation")
    if operation != "score_missense":
        _fail(f"unsupported operation: {operation!r}")
        return

    sequence = req.get("sequence")
    position = req.get("position")
    mut_aa = req.get("mut_aa")
    model_name = req.get("model", "esmc-300m-2024-12")
    token = req.get("token")

    if not sequence or not isinstance(position, int) or not mut_aa or not token:
        _fail("missing required field: sequence, position, mut_aa, or token")
        return
    if not (1 <= position <= len(sequence)):
        _fail(f"position {position} out of range for sequence of length {len(sequence)}")
        return

    try:
        import torch
        from esm.sdk.api import ESMProtein, LogitsConfig
        from esm.sdk.forge import ESMCForgeInferenceClient
        from esm.tokenization import get_esmc_model_tokenizers
    except Exception as exc:  # pragma: no cover - environment problem, not data problem
        _fail(f"could not import EvolutionaryScale esm SDK in isolated venv: {exc}")
        return

    wt_aa = sequence[position - 1]

    try:
        client = ESMCForgeInferenceClient(model=model_name, url="https://biohub.ai", token=token)
        protein = ESMProtein(sequence=sequence)
        protein_tensor = client.encode(protein)
        out = client.logits(protein_tensor, LogitsConfig(sequence=True))
        if out.logits is None or out.logits.sequence is None:
            _fail("Biohub returned no sequence logits")
            return

        tokenizer = get_esmc_model_tokenizers()
        # +1 for the leading <cls> special token the tokenizer/encoder adds.
        token_row = out.logits.sequence[position]
        log_probs = torch.log_softmax(token_row, dim=-1)

        wt_ids = tokenizer.convert_tokens_to_ids([wt_aa])
        mut_ids = tokenizer.convert_tokens_to_ids([mut_aa])
        if wt_ids[0] is None or mut_ids[0] is None:
            _fail(f"tokenizer could not map wt_aa={wt_aa!r} or mut_aa={mut_aa!r}")
            return

        wt_logprob = log_probs[wt_ids[0]].item()
        mut_logprob = log_probs[mut_ids[0]].item()
        llr = mut_logprob - wt_logprob
    except Exception as exc:
        _fail(f"Biohub/Forge call failed: {exc}")
        return

    print(
        json.dumps(
            {
                "wt_aa": wt_aa,
                "mut_aa": mut_aa,
                "position": position,
                "log_likelihood_ratio": round(llr, 3),
                "wt_logprob": round(wt_logprob, 3),
                "mut_logprob": round(mut_logprob, 3),
                "model": f"{model_name}-biohub",
            }
        )
    )


if __name__ == "__main__":
    main()
