"""Deterministic NCBI positional gene lookup — no Claude involved. Factored
out of wildtype/agent/iterative_agent.py's search_genes_near_position tool
so the same real NCBI query can also drive genome_pipeline.py's triage
stage (see that module's docstring: broad, cheap scan first, full
LLM-driven investigation only on survivors). Behavior is unchanged from
the original tool version — same query shape, same fields — just callable
directly instead of only through the agentic loop.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class NearbyGene:
    name: str | None
    description: str | None
    start: int | None
    stop: int | None
    variant_inside_gene_span: bool | None


def find_genes_near_position(
    organism: str, chromosome: str, position: int, window_bp: int = 100000
) -> tuple[list[NearbyGene], str | None]:
    """Real NCBI esearch (db="gene") + esummary within `window_bp` of
    `position`. Returns (genes, error) — error is None on success (even if
    zero genes found; that's a real result, not a failure)."""
    from proto_tools import run_ncbi_esearch, run_ncbi_esummary
    from proto_tools.tools.database_retrieval.ncbi.esearch import NCBIEsearchConfig, NCBIEsearchInput
    from proto_tools.tools.database_retrieval.ncbi.esummary import NCBIEsummaryConfig, NCBIEsummaryInput

    email = "wildtype-hackathon@example.com"
    lo, hi = max(1, position - window_bp), position + window_bp
    query = f"{organism}[Organism] AND {chromosome}[Chromosome] AND {lo}:{hi}[Base Position]"
    try:
        search = run_ncbi_esearch(
            NCBIEsearchInput(db="gene", search_term=query, max_results=15),
            NCBIEsearchConfig(ncbi_email=email),
        )
        if not search.success or not search.ids:
            return [], None

        summary = run_ncbi_esummary(
            NCBIEsummaryInput(db="gene", identifier=",".join(search.ids)),
            NCBIEsummaryConfig(ncbi_email=email),
        )
        genes: list[NearbyGene] = []
        for uid in search.ids:
            rec = summary.summary.get(uid, {})
            # Same discontinued-record issue as search_gene_info in
            # iterative_agent.py (see that function's comment) — matters
            # more here since a stale record's start/stop can be from a
            # completely different, retired genome assembly, which would
            # make variant_inside_gene_span compare against the wrong
            # coordinate system entirely rather than just look redundant.
            if rec.get("currentid"):
                continue
            ginfo = (rec.get("genomicinfo") or [{}])[0]
            start, stop = ginfo.get("chrstart"), ginfo.get("chrstop")
            inside = None
            if isinstance(start, int) and isinstance(stop, int):
                lo_g, hi_g = min(start, stop), max(start, stop)
                inside = lo_g <= position <= hi_g
            genes.append(
                NearbyGene(
                    name=rec.get("name"),
                    description=rec.get("description"),
                    start=start,
                    stop=stop,
                    variant_inside_gene_span=inside,
                )
            )
        return genes, None
    except Exception as e:
        return [], f"NCBI positional gene search unavailable: {type(e).__name__}: {e}"
