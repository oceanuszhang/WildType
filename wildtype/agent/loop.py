"""Pipeline controller.

Deliberate design choice: variant routing and tool dispatch are plain,
deterministic Python — not an autonomous Claude tool-calling loop. Two
reasons: (1) it's the same logic regardless of which mode (mock/local/proto)
is active, so it's the one thing safe to fully build and test tonight
without any sponsor key, and (2) a controller that always calls the same
tools in the same order is far easier to debug live at a demo table than an
LLM improvising tool sequences under time pressure.

Claude's job stays exactly where the proposal put it: reasoning about
tradeoffs (the ESM2/Evo1/pLDDT -> "does this variant earn AlphaFront3?"
escalation call) and writing the two-layer report. Both are isolated below
so they degrade gracefully to a deterministic fallback with no
ANTHROPIC_API_KEY at all — the rest of the pipeline still runs end to end.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from wildtype.agent.prompts import (
    ESCALATION_SYSTEM_PROMPT,
    SYNTHESIS_SYSTEM_PROMPT,
    build_escalation_prompt,
    build_synthesis_prompt,
)
from wildtype.parsers.embark_csv import EmbarkReport, Variant
from wildtype.tools.base import Evo1Score, ESM2Score, LiteratureHit, StructurePrediction, Toolset

# Starting points, not findings — calibrate against the validation set
# (GPT/PRCD as known positives, FGF4 as the known structural positive) once
# real scores are flowing. See docs/variant_routing.md.
ESM2_FLAG_LLR_THRESHOLD = -1.0
EVO1_FLAG_DEVIATION_THRESHOLD = 2.0


@dataclass
class VariantAnalysis:
    variant: Variant
    esm2: ESM2Score | None = None
    evo1: Evo1Score | None = None
    structure: StructurePrediction | None = None
    literature: list[LiteratureHit] = field(default_factory=list)
    escalated_to_af3: bool = False
    notes: list[str] = field(default_factory=list)
    sequence: str | None = None
    sequence_is_real: bool = False


def _placeholder_variant_call(gene: str, seq: str) -> tuple[int, str]:
    """Embark's export gives us a risk *label*, not the underlying
    position/allele — that level of detail comes from VEP-style annotation
    of the raw VCF/microarray probe, which is a separate ingestion path
    (see architecture diagram in the proposal). Stand-in until that's wired,
    even now that `seq` itself may be a real UniProt sequence (toolset.
    sequences) — the *position within it* is still a guess, not Ollie's
    actual mutation.
    """
    pos = (sum(map(ord, gene)) % (len(seq) - 1)) + 1
    mut_aa = "A" if seq[pos - 1] != "A" else "G"
    return pos, mut_aa


def analyze_variant(
    variant: Variant,
    toolset: Toolset,
    *,
    species: str = "dog",
    organism: str = "Canis lupus familiaris",
) -> VariantAnalysis:
    """`species`/`organism` default to dog so the existing Embark/Ollie path
    is unchanged — species-name-only mode (wildtype/agent/species_mode.py)
    passes the actual target species through both, since sequence fetch and
    literature search were previously hardcoded to canine, which defeated
    the species-agnostic claim for any non-dog input."""
    analysis = VariantAnalysis(variant=variant)
    gene = variant.gene_hint

    seq_result = toolset.sequences.fetch(gene, organism)
    analysis.sequence = seq_result.sequence
    analysis.sequence_is_real = seq_result.is_real
    if not seq_result.is_real:
        analysis.notes.append(
            f"Sequence for {gene} is a deterministic placeholder, not a real UniProt "
            "fetch (no confident match) — downstream scores are a plumbing check, not a finding."
        )

    if variant.variant_class == "structural":
        context_seq = seq_result.sequence
        analysis.evo1 = toolset.dna.score_structural(gene, context_seq, variant.name)
        if "PROTO GPU UNAVAILABLE" in analysis.evo1.interpretation:
            analysis.notes.append(
                "Evo1 score is a mock fallback, not real — proto_tools' Evo1 needs GPU/Modal "
                "(no payment method on file yet), so this is a plumbing check, not a finding."
            )
        if analysis.evo1.context_deviation >= EVO1_FLAG_DEVIATION_THRESHOLD:
            analysis.structure = toolset.structure.fold(gene, context_seq)
            if "proto-gpu-unavailable" in analysis.structure.model:
                analysis.notes.append(
                    "Structure prediction is a mock fallback, not real — same GPU/Modal gap as above."
                )
        else:
            analysis.notes.append("Evo1 deviation below flag threshold — structure prediction skipped.")
    else:
        wt_seq = seq_result.sequence
        position, mut_aa = _placeholder_variant_call(gene, wt_seq)
        # Even with a real wildtype sequence, this exact position/mutant
        # allele is still a guess, not Ollie's actual variant call — see
        # _placeholder_variant_call docstring.
        analysis.notes.append(
            f"Mutation position/allele for {gene} is still a placeholder, not Ollie's "
            "actual variant call — VEP-style VCF annotation isn't wired yet."
        )
        analysis.esm2 = toolset.scoring.score_missense(gene, wt_seq, position, mut_aa)
        if analysis.esm2.log_likelihood_ratio <= ESM2_FLAG_LLR_THRESHOLD:
            analysis.structure = toolset.structure.fold(gene, wt_seq)
            if "proto-gpu-unavailable" in analysis.structure.model:
                analysis.notes.append(
                    "Structure prediction is a mock fallback, not real — proto_tools' ESMFold "
                    "needs GPU/Modal (no payment method on file yet)."
                )
        else:
            analysis.notes.append("ESM2 score above flag threshold — structure prediction skipped.")

    # "GPT" (the gene) collides hard with "GPT" (the AI model) in a
    # literature corpus that indexes AI papers — live-tested query "GPT
    # Elevated ALT Activity (GPT) dog canine" returned ChatGPT medical-exam
    # benchmark papers, not the gene. "gene" biases the search back toward
    # genomics; doesn't fully solve it for every ambiguous gene symbol
    # (HELP, DOG1, MET, FUS all collide with common English words too),
    # just the one we caught live. `species` replaces the old hardcoded
    # "dog canine" so this isn't silently canine-only.
    query = f"{gene} gene {variant.name} {species} genetics"
    analysis.literature = toolset.literature.search(query)
    return analysis


def run_pipeline(report: EmbarkReport, toolset: Toolset, also_include_genes: list[str] | None = None) -> list[VariantAnalysis]:
    also_include_genes = also_include_genes or []
    targets = list(report.flagged)
    for v in report.variants:
        if v not in targets and any(g.lower() in v.gene_hint.lower() for g in also_include_genes):
            targets.append(v)

    results = []
    for v in targets:
        try:
            results.append(analyze_variant(v, toolset))
        except Exception as e:
            # One variant blowing up (bad gene_hint, a real-tool edge case,
            # anything) shouldn't take down the rest of the batch — hit this
            # live 2026-08-15 (PRCD's gene_hint "PRCD Exon 1" crashed
            # UniProt fetch mid-run). Surface the failure in the report
            # instead of losing the whole demo run to one variant.
            failed = VariantAnalysis(variant=v)
            failed.notes.append(f"ANALYSIS FAILED: {type(e).__name__}: {e}")
            results.append(failed)
    return results


# --- Claude: escalation + synthesis -----------------------------------------


def _get_client():
    import os

    key = os.environ.get("ANTHROPIC_API_KEY")
    if not key:
        return None
    import anthropic

    return anthropic.Anthropic(api_key=key)


def _extract_text(resp) -> str:
    """Claude Sonnet 5 returns ThinkingBlock entries ahead of the actual
    TextBlock in `resp.content` when extended thinking is active —
    `resp.content[0].text` breaks live with
    AttributeError: 'ThinkingBlock' object has no attribute 'text'.
    Find the first block that actually has text instead of assuming index 0.
    """
    for block in resp.content or []:
        if getattr(block, "type", None) == "text":
            return block.text
    return ""


def maybe_escalate_to_af3(analysis: VariantAnalysis, toolset: Toolset, client=None) -> VariantAnalysis:
    """No-op (stays SKIP) without an Anthropic key — AlphaFold3 is a
    nice-to-have per the proposal's must-have/nice-to-have split, so
    degrading silently here is the correct behavior, not a bug to fix later.
    """
    client = client if client is not None else _get_client()
    if client is None or analysis.structure is None:
        return analysis

    prompt = build_escalation_prompt(
        gene=analysis.variant.gene_hint,
        esm2_llr=analysis.esm2.log_likelihood_ratio if analysis.esm2 else None,
        evo1_deviation=analysis.evo1.context_deviation if analysis.evo1 else None,
        plddt=analysis.structure.plddt_mean,
    )
    resp = client.messages.create(
        model="claude-sonnet-5",
        # ESCALATE/SKIP + one sentence doesn't need deep reasoning — disable
        # adaptive thinking (on by default on Sonnet 5 when `thinking` is
        # omitted) so the whole max_tokens budget goes to visible text.
        # Without this, max_tokens=200 was silently eaten by thinking tokens
        # and _extract_text() got nothing — hit live 2026-08-15.
        thinking={"type": "disabled"},
        max_tokens=500,
        system=ESCALATION_SYSTEM_PROMPT,
        messages=[{"role": "user", "content": prompt}],
    )
    text = _extract_text(resp)
    if text.strip().upper().startswith("ESCALATE"):
        analysis.escalated_to_af3 = True
        analysis.structure = toolset.structure.fold(
            analysis.variant.gene_hint,
            analysis.sequence or "",
            high_accuracy=True,
        )
        analysis.notes.append(f"Escalated to AlphaFold3: {text.strip()}")
    return analysis


def _evidence_block(analysis: VariantAnalysis) -> str:
    lines = [
        f"Variant: {analysis.variant.name}",
        f"Gene: {analysis.variant.gene_hint}",
        # "Source finding" not "Embark call" — this evidence block also
        # backs species-name-only mode (wildtype/agent/species_mode.py),
        # where raw_value is a literature citation, not an Embark result.
        # The old "Embark call:" label bled into those reports and read as
        # a claim about a report that doesn't exist for that path.
        f"Source finding: {analysis.variant.raw_value}",
        f"Route: {analysis.variant.variant_class}",
    ]
    if analysis.esm2:
        lines.append(f"ESM2: LLR={analysis.esm2.log_likelihood_ratio} model={analysis.esm2.model}")
    if analysis.evo1:
        lines.append(f"Evo1: deviation={analysis.evo1.context_deviation}SD — {analysis.evo1.interpretation}")
    if analysis.structure:
        lines.append(f"Structure ({analysis.structure.model}): mean pLDDT={analysis.structure.plddt_mean}")
    for hit in analysis.literature:
        lines.append(f"Literature: {hit.title} (DOI: {hit.doi}, {hit.source})")
    for note in analysis.notes:
        lines.append(f"Note: {note}")
    return "\n".join(lines)


def _fallback_report(analysis: VariantAnalysis) -> str:
    """Deterministic template — no ANTHROPIC_API_KEY needed. Coarser than
    Claude's synthesis but keeps the full pipeline runnable tonight."""
    v = analysis.variant
    plain = [
        f"**{v.name}** — {v.raw_value}.",
    ]
    if analysis.structure:
        plain.append(f"Structural analysis (mean confidence {analysis.structure.plddt_mean}) supports this being biologically meaningful.")
    if v.risk_level == "carrier":
        plain.append("This is a carrier result: relevant for breeding decisions, no clinical action needed for this animal.")
    science = _evidence_block(analysis)
    return "Plain layer\n" + "\n".join(plain) + "\n\nScience layer\n" + science


def render_report(analysis: VariantAnalysis, client=None) -> str:
    client = client if client is not None else _get_client()
    if client is None:
        return _fallback_report(analysis)
    resp = client.messages.create(
        model="claude-sonnet-5",
        # Adaptive thinking (default on Sonnet 5) genuinely helps here —
        # weighing evidence quality is the point of this call — but
        # max_tokens is a hard cap over thinking + visible text combined.
        # 800 was sized for text-only output and left thinking nothing to
        # work with before truncating the report mid-sentence (hit live
        # 2026-08-15). Bumped with real headroom for both.
        thinking={"type": "adaptive"},
        max_tokens=4096,
        system=SYNTHESIS_SYSTEM_PROMPT,
        messages=[{"role": "user", "content": build_synthesis_prompt(_evidence_block(analysis))}],
    )
    return _extract_text(resp) or _fallback_report(analysis)
