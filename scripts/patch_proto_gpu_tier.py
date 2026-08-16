"""Patches proto_tools' installed Modal service definitions to request the
cheap GPU tier (T4/L4/A10 — proto_tools' own GPU_BASIC) instead of the
default (H100/H200/A100-80GB — GPU_DEFAULT), for tools that have no
architectural reason to need the expensive tier.

Why this exists: Modal gates H100/H200/etc. behind a payment method on
file; T4/L4/A10 do not have that gate (confirmed 2026-08-15 — GPU_BASIC is
already the real, used tier for several other proto_tools services:
esm_if1, ablang, metal3d, malinois, splice_transformer, pangolin). Without
a payment method added, GPU_DEFAULT tools can't deploy at all; GPU_BASIC
tools can.

Scope, deliberately NOT including everything:
  - esm2, esmfold, evo1: patched. GPU_DEFAULT was proto_tools' own default
    choice for a "new service" (see gpu_profiles.py), not a stated hard
    requirement — nothing in these three models needs Hopper-class GPUs.
  - evo2: NOT patched, and never will be by this script. Its service file
    imports GPU_HOPPER specifically because Evo2 needs Hopper FP8 support
    — T4/L4/A10 are physically incapable of running it. Forcing this would
    produce a broken deploy, not a cheaper one.

This edits the installed package under site-packages directly — it does
NOT survive `pip install -r requirements.txt` re-running proto-tools' git
install, or upgrading the package. Re-run this script after either. Not
committed as a patch to proto_tools itself since we don't control that repo.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

TARGETS = {
    "masked_models/esm2_deployment/esm2_service.py",
    "structure_prediction/esmfold_deployment/esmfold_service.py",
    "causal_models/evo1_deployment/evo1_service.py",
}


def find_proto_modal_dir() -> Path:
    import proto_tools

    return Path(proto_tools.__file__).parent / "modal"


def patch_file(path: Path) -> bool:
    text = path.read_text()
    if "GPU_BASIC" in text and "GPU_DEFAULT" not in text:
        return False  # already patched
    new_text = text.replace(
        "from proto_tools.modal.gpu_profiles import GPU_DEFAULT",
        "from proto_tools.modal.gpu_profiles import GPU_BASIC as GPU_DEFAULT  # patched: scripts/patch_proto_gpu_tier.py",
    )
    if new_text == text:
        print(f"  SKIP (no GPU_DEFAULT import found): {path}")
        return False
    path.write_text(new_text)
    return True


def main() -> None:
    modal_dir = find_proto_modal_dir()
    print(f"proto_tools modal dir: {modal_dir}")
    patched = 0
    for rel in sorted(TARGETS):
        path = modal_dir / rel
        if not path.exists():
            print(f"  MISSING (proto_tools version drift?): {path}")
            continue
        if patch_file(path):
            print(f"  patched: {rel}")
            patched += 1
        else:
            print(f"  already patched or unchanged: {rel}")
    print(f"\n{patched} file(s) patched. Re-run `proto-tools deploy --apps esm2,esmfold,evo1 --test` to apply.")


if __name__ == "__main__":
    sys.exit(main())
