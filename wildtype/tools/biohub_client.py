"""Biohub client — real ESMC missense scoring via EvolutionaryScale's Forge
API (https://biohub.ai), as an ADDITIONAL path alongside Proto
(wildtype/tools/proto_client.py), not a replacement — same models Proto
already serves for free via Modal, this is a second, independent source
for cross-checking or as a fallback if Modal/proto-tools is unavailable.

Confirmed live 2026-08-15 with the user's real BIOHUB_API_KEY:
  - forge.evolutionaryscale.ai deprecates 2026-11-28; biohub.ai is the
    forward URL for the exact same API (EvolutionaryScale's SDK warns you
    to switch, so this client goes straight to biohub.ai).
  - Correct client class is ESMCForgeInferenceClient, not
    ESM3ForgeInferenceClient — using the ESM3 client against an ESMC model
    name is what produced the earlier raw-HTTP 500
    ("'SequenceStructureTokens' object has no attribute
    'secondary_structure'") — an ESM3-shaped request against an
    ESMC-shaped endpoint, not a bad API key (that 500, not 401, is what
    first confirmed the key itself was valid).
  - Requires encode() before logits() — raw ESMProtein isn't accepted,
    matching the SDK's own error message when you skip that step.

Why this shells out to a subprocess instead of importing esm directly:
EvolutionaryScale's official `esm` SDK and the already-installed `fair-esm`
package both claim the top-level `esm` import name. Installing the former
in the main project venv would silently break
wildtype/tools/esm_local.py's local CPU ESM2 fallback (WILDTYPE_MODE=local).
So the SDK lives in an isolated venv (`.venv-biohub/`, gitignored — see
setup below) and scripts/biohub_forge_runner.py runs inside it, called
here as a subprocess with JSON over stdin/stdout.

One-time setup to recreate `.venv-biohub/` on a new machine:
    python3 -m venv .venv-biohub
    .venv-biohub/bin/pip install esm httpx python-dotenv
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

from wildtype.tools.base import ESM2Score

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_BIOHUB_VENV_PYTHON = _REPO_ROOT / ".venv-biohub" / "bin" / "python"
_RUNNER_SCRIPT = _REPO_ROOT / "scripts" / "biohub_forge_runner.py"

# Same rationale as ProtoScoringClient._MAX_ESM2_RESIDUES: window rather
# than truncate from one end, so the model still sees context on both
# sides of the mutation. Biohub/ESMC's real cap is unconfirmed (a
# 3000-residue test succeeded live 2026-08-15) — this cap is a
# precaution to keep single calls fast and cheap, not a known hard limit.
_MAX_RESIDUES = 1022


class BiohubConfigError(RuntimeError):
    pass


class ProtoScoringClient_Biohub:
    """Deliberately not a subclass of proto_client.ProtoScoringClient —
    different backend, different failure modes, kept separate so a bug in
    one can't silently masquerade as the other. Same score_missense(...)
    signature/return type so callers can swap between them."""

    def score_missense(self, gene: str, wt_seq: str, position: int, mut_aa: str) -> ESM2Score:
        token = os.environ.get("BIOHUB_API_KEY", "").strip()
        if not token:
            raise BiohubConfigError("BIOHUB_API_KEY not set in environment/.env")
        if not _BIOHUB_VENV_PYTHON.exists():
            raise BiohubConfigError(
                f"{_BIOHUB_VENV_PYTHON} not found — see the setup steps in "
                "wildtype/tools/biohub_client.py's module docstring"
            )

        window_seq, window_position = self._window(wt_seq, position, _MAX_RESIDUES)

        request = {
            "operation": "score_missense",
            "sequence": window_seq,
            "position": window_position,
            "mut_aa": mut_aa,
            "model": "esmc-300m-2024-12",
            "token": token,
        }

        try:
            proc = subprocess.run(
                [str(_BIOHUB_VENV_PYTHON), str(_RUNNER_SCRIPT)],
                input=json.dumps(request),
                capture_output=True,
                text=True,
                timeout=60,
            )
        except subprocess.TimeoutExpired as exc:
            raise BiohubConfigError(f"Biohub call timed out: {exc}") from exc

        stdout = proc.stdout.strip()
        if not stdout:
            raise BiohubConfigError(
                f"biohub_forge_runner.py produced no output (exit {proc.returncode}); "
                f"stderr: {proc.stderr[-2000:]}"
            )

        try:
            out = json.loads(stdout.splitlines()[-1])
        except json.JSONDecodeError as exc:
            raise BiohubConfigError(f"could not parse runner output as JSON: {stdout[:2000]}") from exc

        if "error" in out:
            raise BiohubConfigError(f"Biohub/Forge error: {out['error']}")

        return ESM2Score(
            gene=gene,
            position=position,  # report in the caller's original coordinates
            wt_aa=out["wt_aa"],
            mut_aa=out["mut_aa"],
            log_likelihood_ratio=out["log_likelihood_ratio"],
            model=out["model"],
        )

    def _window(self, seq: str, position: int, max_residues: int) -> tuple[str, int]:
        """Identical logic to ProtoScoringClient._window — kept as its own
        copy rather than a shared import so the two clients stay fully
        independent (see class docstring)."""
        if len(seq) <= max_residues:
            return seq, position
        half = max_residues // 2
        start = max(0, (position - 1) - half)
        end = min(len(seq), start + max_residues)
        start = max(0, end - max_residues)
        return seq[start:end], position - start


def biohub_available() -> bool:
    """Cheap pre-flight check: is there a key AND is the isolated venv set
    up? Callers (e.g. the iterative agent, when deciding which scoring
    tools to expose to Claude) should check this before offering Biohub as
    an option, rather than exposing a tool that will always error."""
    return bool(os.environ.get("BIOHUB_API_KEY", "").strip()) and _BIOHUB_VENV_PYTHON.exists()


if __name__ == "__main__":
    # Manual smoke test: python -m wildtype.tools.biohub_client
    from dotenv import load_dotenv

    load_dotenv(_REPO_ROOT / ".env")
    client = ProtoScoringClient_Biohub()
    test_seq = (
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFS"
        "AGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWELVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELL"
        "SRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTL"
        "KHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"
    )
    result = client.score_missense("TEST", test_seq, 5, "P")
    print(result, file=sys.stderr)
