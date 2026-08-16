"""Species-name-only entry point — no Embark CSV, no personal genomic data
at all. Generalizes the proposal's "breed-name-only mode" beyond dogs: give
it any species (an endangered animal with no consumer genomics panel is the
real stress test for "species-agnostic"), and it discovers documented
disease-associated genes via Paperclip, then routes each through the same
analyze_variant/render_report pipeline Ollie's Embark variants use.

The one hard rule here: NEVER let Claude name a gene from general knowledge.
Extraction is constrained to genes whose symbol literally appears in the
Paperclip search results handed to it — grounded in retrieved evidence, not
recalled from training. This is slower and finds fewer genes than just
asking Claude "what genes cause disease in species X", but it's the
difference between a citation trail and a plausible-sounding guess — see
the project's own "Anthropic's ESM family" mistake earlier in this build for
why that distinction matters here specifically.
"""
from __future__ import annotations

import json

from wildtype.agent.loop import VariantAnalysis, _extract_text, _get_client, analyze_variant
from wildtype.parsers.embark_csv import Variant
from wildtype.tools.base import Toolset


class SpeciesResearchError(RuntimeError):
    pass


def discover_species_variants(
    species: str,
    toolset: Toolset,
    *,
    organism: str | None = None,
    max_genes: int = 3,
    client=None,
) -> list[Variant]:
    """Real Paperclip search + grounded extraction. Returns [] if nothing
    checks out — never fabricates a gene to avoid returning empty-handed.
    """
    client = client if client is not None else _get_client()
    if client is None:
        raise SpeciesResearchError(
            "Species-name-only mode needs Claude to extract gene candidates from "
            "literature safely — without ANTHROPIC_API_KEY this would mean guessing "
            "gene names from training data, which this pipeline deliberately doesn't do. "
            "Set the key, or use Embark-CSV mode instead."
        )

    organism = organism or species
    query = f"{species} known heritable genetic disease causal gene mutation"
    hits = toolset.literature.search(query, max_results=8)
    if not hits:
        return []

    evidence = "\n".join(
        f"{i + 1}. {h.title} ({h.year}, DOI:{h.doi}, source={h.source}): {h.summary}"
        for i, h in enumerate(hits)
    )
    prompt = (
        f'Literature search results for species "{species}":\n\n{evidence}\n\n'
        f"List up to {max_genes} genes EXPLICITLY NAMED in the titles or summaries "
        "above as associated with a heritable disease or health condition in this "
        "species. Only include a gene if its symbol literally appears in the text "
        "above — do not add genes from general knowledge, even ones you're confident "
        "are relevant. If no gene symbols appear anywhere in the text, return an "
        "empty array.\n\n"
        "Respond with ONLY a JSON array, no other text, no markdown fences:\n"
        '[{"gene": "SYMBOL", "condition": "short description", "source_index": N, '
        '"variant_class": "missense_snp"}]\n'
        '"source_index" is the 1-based number of the result above that names the gene. '
        '"variant_class" is "structural" only if the text describes an insertion, '
        'deletion, duplication, or retrogene — otherwise "missense_snp".'
    )
    resp = client.messages.create(
        model="claude-sonnet-5",
        thinking={"type": "disabled"},
        max_tokens=1500,
        messages=[{"role": "user", "content": prompt}],
    )
    text = _extract_text(resp).strip()
    if text.startswith("```"):
        text = text.strip("`")
        if text.startswith("json"):
            text = text[4:]
    try:
        candidates = json.loads(text)
    except (json.JSONDecodeError, ValueError):
        return []
    if not isinstance(candidates, list):
        return []

    variants = []
    for c in candidates[:max_genes]:
        if not isinstance(c, dict):
            continue
        gene = str(c.get("gene", "")).strip()
        if not gene:
            continue
        idx = c.get("source_index")
        source_hit = hits[idx - 1] if isinstance(idx, int) and 1 <= idx <= len(hits) else None
        condition = str(c.get("condition", "")).strip() or "documented condition"
        variant_class = c.get("variant_class") if c.get("variant_class") == "structural" else "missense_snp"
        raw_value = "documented in species literature"
        if source_hit:
            raw_value += f" — {source_hit.title} (DOI: {source_hit.doi})"
        variants.append(
            Variant(
                name=f"{condition} ({gene})",
                raw_value=raw_value,
                risk_level="unknown",
                zygosity="",
                gene_hint=gene,
                variant_class=variant_class,
            )
        )
    return variants


def run_species_pipeline(
    species: str,
    toolset: Toolset,
    *,
    organism: str | None = None,
    max_genes: int = 3,
    client=None,
) -> list[VariantAnalysis]:
    """Discover candidate genes for `species`, then run each through the
    same analyze_variant used for Embark data — one variant failing
    (bad gene symbol, no UniProt match, whatever) doesn't take down the rest.
    """
    variants = discover_species_variants(species, toolset, organism=organism, max_genes=max_genes, client=client)
    resolved_organism = organism or species

    results = []
    for v in variants:
        try:
            results.append(analyze_variant(v, toolset, species=species, organism=resolved_organism))
        except Exception as e:
            failed = VariantAnalysis(variant=v)
            failed.notes.append(f"ANALYSIS FAILED: {type(e).__name__}: {e}")
            results.append(failed)
    return results


if __name__ == "__main__":
    import sys

    from dotenv import load_dotenv

    from wildtype.agent.loop import maybe_escalate_to_af3, render_report
    from wildtype.tools.base import get_toolset

    load_dotenv()
    # Usage: python -m wildtype.agent.species_mode "horse" ["Equus caballus"]
    # Second arg (binomial name) is optional but matters — live-tested
    # 2026-08-15: UniProt fetch failed on organism="horse" (common name,
    # no match) but works on "Equus caballus". "dog"/"Canis lupus
    # familiaris" both happen to work since UniProt special-cases a few
    # major model/companion organisms; don't assume that generalizes.
    args = sys.argv[1:]
    species = args[0] if args else "cheetah"
    organism = args[1] if len(args) > 1 else species
    if organism == species:
        print(f"(no binomial name given — using {organism!r} for UniProt lookups, "
              f"which may not resolve; pass one as a 2nd arg for reliable matches)")
    print(f"=== species-name-only mode: {species!r} (organism={organism!r}) ===\n")

    toolset = get_toolset()
    analyses = run_species_pipeline(species, toolset, organism=organism)
    if not analyses:
        print("No genes found in literature — either no results, or nothing named a gene explicitly.")
    for analysis in analyses:
        analysis = maybe_escalate_to_af3(analysis, toolset)
        print("-" * 72)
        print(render_report(analysis))
        print()
