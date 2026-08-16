"""⚠️ PLACEHOLDER BIOLOGY — not real UniProt sequences.

Tonight's offline dev needs *some* amino acid sequence to feed the scoring
clients so the pipeline plumbing can be built and tested end-to-end. These
are deterministic pseudo-random sequences keyed by gene name — real amino
acid composition, zero biological meaning. Any ESM2 LLR computed against
them (even via the real local model in tools/esm_local.py) is a plumbing
test, not a finding.

Replace with tools/proto_client.py's real UniProt fetch tomorrow. Nothing
downstream needs to change shape — every consumer just gets a real
sequence string instead of this one.
"""
from __future__ import annotations

import hashlib

_AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def get_placeholder_sequence(gene: str, length: int = 220) -> str:
    seed = hashlib.sha256(gene.encode()).digest()
    out = []
    i = 0
    while len(out) < length:
        out.append(_AA_ALPHABET[seed[i % len(seed)] % len(_AA_ALPHABET)])
        i += 1
        if i % len(seed) == 0:
            seed = hashlib.sha256(seed).digest()
    return "".join(out)


class PlaceholderSequenceClient:
    """SequenceClient (wildtype/tools/base.py) backed by fake sequences —
    used in mock/local modes. Always reports is_real=False so callers
    (agent/loop.py) can flag placeholder biology in the report."""

    def fetch(self, gene: str, organism: str = "Canis lupus familiaris"):
        from wildtype.tools.base import SequenceResult

        # organism unused — placeholder is a pure function of gene name, but
        # accepts the param so it satisfies SequenceClient like the real one.
        return SequenceResult(sequence=get_placeholder_sequence(gene), is_real=False, source="placeholder")
