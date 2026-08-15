"""Paperclip client — live literature search across 11M+ papers, UniProt, PDB.

STATUS: shape-only stub, same situation as tools/proto_client.py — real
auth/endpoint details land at check-in. Written as a generic REST client
since that's the most likely shape; adjust `_request()` and the response
field mapping in `search()` once Paperclip's actual API docs are in hand.

TODO tomorrow:
  1. Get PAPERCLIP_API_KEY / PAPERCLIP_BASE_URL at check-in.
  2. Confirm the real request/response schema and fix the mapping below.
  3. Flip WILDTYPE_MODE=proto in .env (Paperclip ships alongside Proto mode).
"""
from __future__ import annotations

import os

import requests

from wildtype.tools.base import LiteratureHit


class PaperclipConfigError(RuntimeError):
    pass


class PaperclipClient:
    def __init__(self):
        self.api_key = os.environ.get("PAPERCLIP_API_KEY")
        self.base_url = os.environ.get("PAPERCLIP_BASE_URL", "").rstrip("/")

    def _require_config(self):
        if not self.api_key or not self.base_url:
            raise PaperclipConfigError(
                "PAPERCLIP_API_KEY / PAPERCLIP_BASE_URL not set. Fill in .env after check-in, "
                "or run with WILDTYPE_MODE=mock for now."
            )

    def search(self, query: str, max_results: int = 5) -> list[LiteratureHit]:
        self._require_config()
        resp = requests.get(
            f"{self.base_url}/search",  # TODO: confirm actual endpoint path
            headers={"Authorization": f"Bearer {self.api_key}"},
            params={"q": query, "limit": max_results},
            timeout=20,
        )
        resp.raise_for_status()
        data = resp.json()  # TODO: confirm actual response schema
        return [
            LiteratureHit(
                title=item.get("title", ""),
                doi=item.get("doi"),
                year=item.get("year"),
                summary=item.get("summary", item.get("abstract", "")),
                source=item.get("source", "paperclip"),
            )
            for item in data.get("results", [])
        ]
