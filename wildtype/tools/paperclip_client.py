"""Paperclip client — live literature search across PMC, arXiv, bioRxiv,
medRxiv, UniProt, PDB via the `gxl_paperclip` SDK (https://paperclip.gxl.ai).

STATUS: real SDK, shape confirmed against a live authenticated call on
2026-08-15 (`client.search(...)` against "FGF4 retrogene chondrodystrophy
dog IVDD" returned 3 real PMC papers). `ExecuteResult.result_data["papers"]`
is a list of dicts with document_id/source/pub_year/title/journal_title/
doi/authors/tldr/pub_date/abstract_snippet — that's the field mapping used
below.

Setup (once, per machine — not part of this repo):
    curl -fsSL https://paperclip.gxl.ai/install.sh | bash
    paperclip login          # interactive browser sign-in, caches
                              # ~/.paperclip/credentials.json
`gxl_paperclip` isn't on PyPI; the installer drops it in ~/.paperclip/lib/.
This repo's .venv picks it up via a .pth file pointing there (not committed
to requirements.txt since `pip install` can't fetch it). If you're on a
fresh machine at the venue, rerun the installer + login there too.

`PaperclipClient.from_env()` reads PAPERCLIP_API_KEY if set, otherwise
falls back to the cached OAuth credentials from `paperclip login` — so this
works with either an API key (PAPERCLIP_API_KEY in .env, better for the
demo machine so it doesn't depend on one person's login) or the personal
login used to build/test this integration.
"""
from __future__ import annotations

from wildtype.tools.base import LiteratureHit


class PaperclipConfigError(RuntimeError):
    pass


class PaperclipClient:
    def __init__(self):
        self._client = None

    def _get_client(self):
        if self._client is None:
            try:
                from gxl_paperclip import PaperclipClient as _SDKClient
            except ImportError as e:
                raise PaperclipConfigError(
                    "gxl_paperclip not importable. Run "
                    "`curl -fsSL https://paperclip.gxl.ai/install.sh | bash` "
                    "and `paperclip login`, or run with WILDTYPE_MODE=mock."
                ) from e
            self._client = _SDKClient.from_env()
        return self._client

    def search(self, query: str, max_results: int = 5) -> list[LiteratureHit]:
        client = self._get_client()

        try:
            result = client.search(query, limit=max_results, source="pmc")
        except Exception as e:
            from gxl_paperclip import AuthError

            if isinstance(e, AuthError):
                raise PaperclipConfigError(
                    "Paperclip auth rejected — run `paperclip login` again "
                    "or set PAPERCLIP_API_KEY in .env."
                ) from e
            raise

        papers = (result.result_data or {}).get("papers", [])
        return [
            LiteratureHit(
                title=p.get("title", ""),
                doi=p.get("doi"),
                year=p.get("pub_year"),
                summary=p.get("tldr") or p.get("abstract_snippet", ""),
                source=p.get("source", "paperclip"),
            )
            for p in papers[:max_results]
        ]


if __name__ == "__main__":
    import sys

    query = " ".join(sys.argv[1:]) or "FGF4 retrogene chondrodystrophy dog IVDD"
    for hit in PaperclipClient().search(query):
        print(f"- {hit.title} ({hit.year}) DOI:{hit.doi} [{hit.source}]")
        print(f"  {hit.summary}")
