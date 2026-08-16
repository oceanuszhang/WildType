"""Iterative, LLM-driven agent — PROJECT_SPEC.md sections 5.1, 8-10.

Replaces wildtype/agent/loop.py's deterministic controller for the new
raw-genome path: Claude decides which tools to call, evaluates whether it
has enough evidence, and iterates — rather than a fixed Python function
always calling the same tools in the same order. Bounded by
`max_iterations` to control cost (spec section 8 explicitly allows this
for the first implementation).

Deterministic code still owns everything spec section 5.1 assigns to it —
file parsing, the tool implementations themselves, output serialization.
The LLM only orchestrates; it doesn't parse VCFs or call Evo2's Python API
directly, it calls the tool wrapper, which does.
"""
from __future__ import annotations

import json
import os

import anthropic
from anthropic import beta_tool

from wildtype.agent.schema import ConfidenceAssessment, EvidenceBundle, RiskFinding, VariantRef
from wildtype.parsers.common_variant import CommonVariant
from wildtype.tools.reference_genome import AssemblyInfo, fetch_reference_context, resolve_reference_assembly

MAX_ITERATIONS = 12  # was 8 — too tight once search_genes_near_position/source_citation
# requirements gave the agent more ground to cover per variant; hit the cap on a real
# run and got a canned "no conclusion" fallback despite real evidence sitting unused
# in evidence_log. See _synthesize_from_evidence below for the actual fix.

AGENT_SYSTEM_PROMPT = """\
You are WildType's investigative agent (see PROJECT_SPEC.md). You are handed one \
candidate genomic variant — a real, confirmed deviation from the species reference \
genome, already established, not something you need to re-derive — and must decide, \
using the tools available, whether it represents a plausible genetic risk factor.

Core rules, non-negotiable:
- A foundation-model score is computational evidence, not a disease probability. \
  Never convert "Evo2 deviation = 4.3" into "90% disease risk" — state the score \
  and what it does/doesn't establish, separately.
- Ground every claim in a specific tool result. Never state a gene's function, a \
  mechanism, or a literature finding you did not actually retrieve via a tool call.
- Evidence from a different species is not equivalent to evidence from this animal's \
  own species — when you cite something, note whether it's from this species or not.
- If the evidence is weak, conflicting, or absent, the honest output is \
  confidence_level="insufficient_evidence" or "low" — not a speculative high-confidence \
  finding to seem thorough. A tool returning an error (e.g. a GPU-dependent model \
  being unavailable) is real information — say so, don't paper over it.
- Never compute or guess a protein amino-acid position from a DNA chromosome:position \
  yourself — that translation requires exon/reading-frame data no tool here provides. \
  Only call score_protein_variant_esm2 when a specific paper states the protein-level \
  change explicitly; otherwise skip protein-level scoring and say why in limitations.
- You have a limited tool-call budget. Prioritize: check whether this locus falls in or \
  near a known gene first (search_genes_near_position), then decide whether literature \
  or DNA-level scoring is worth spending the rest of the budget on.

When you've investigated enough — or are near your budget — call create_risk_finding \
EXACTLY ONCE with your conclusion. That includes concluding "insufficient evidence" \
if that's the honest read. Calling it ends the investigation.
"""


def _tool_result_text(resp) -> str:
    if isinstance(resp, str):
        return resp
    return json.dumps(resp)


def build_investigation_tools(species: str, organism: str, assembly: AssemblyInfo | None, evidence_log: list[dict]):
    """@beta_tool-decorated functions for one agent run, closed over
    species/organism (so the agent can't accidentally query the wrong one)
    and an `evidence_log` list each tool appends its call+result to — this
    is how the final RiskFinding's EvidenceBundle gets populated, rather
    than re-parsing the Tool Runner's internal message stream after the fact.
    """

    @beta_tool
    def fetch_reference_sequence(chromosome: str, position: int, window_bp: int = 200) -> str:
        """Fetch real reference-genome DNA sequence around a position.

        Args:
            chromosome: Chromosome/contig name as in the input file (e.g. "12").
            position: 1-based genomic position.
            window_bp: Total window size in base pairs, centered on position.
        """
        ctx = fetch_reference_context(organism, chromosome, position, window=window_bp, assembly=assembly)
        if ctx is None:
            result = {"error": f"No reference sequence resolvable for {organism} {chromosome}:{position}"}
        else:
            result = {"sequence": ctx.sequence, "window_start": ctx.start, "accession": ctx.accession}
        evidence_log.append({"category": "genomic_databases", "tool": "fetch_reference_sequence",
                              "args": {"chromosome": chromosome, "position": position}, "result": result})
        return _tool_result_text(result)

    @beta_tool
    def score_dna_variant_evo2(chromosome: str, position: int, alternate_allele: str, window_bp: int = 200) -> str:
        """Score this variant's effect with Evo2 — a genomic foundation model
        trained across bacteria through eukaryotes (appropriate for animal DNA,
        unlike Evo1 which is prokaryote-only). Returns a log-likelihood
        deviation, NOT a disease probability. May be unavailable (GPU-gated) —
        an error result is real information, not a failure to work around.

        Give this the variant's own chromosome/position/alternate allele —
        the tool fetches the real reference window and builds the mutant
        sequence itself (deterministic substitution at the exact indexed
        position), rather than asking you to splice sequence strings by hand.

        Args:
            chromosome: Chromosome/contig name as in the input file (e.g. "12").
            position: 1-based genomic position of the variant.
            alternate_allele: The called allele that differs from reference.
            window_bp: Total window size in base pairs, centered on position.
        """
        ctx = fetch_reference_context(organism, chromosome, position, window=window_bp, assembly=assembly)
        if ctx is None:
            result = {"error": f"No reference sequence resolvable for {organism} {chromosome}:{position}"}
            evidence_log.append({"category": "foundation_models", "tool": "score_dna_variant_evo2",
                                  "args": {"chromosome": chromosome, "position": position,
                                           "alternate_allele": alternate_allele}, "result": result})
            return _tool_result_text(result)

        idx = position - ctx.start
        reference_sequence = ctx.sequence
        mutant_sequence = reference_sequence[:idx] + alternate_allele + reference_sequence[idx + 1:]

        from proto_tools import run_evo2_score
        from proto_tools.tools.causal_models.evo2.evo2_score import Evo2ScoringConfig
        from proto_tools.tools.causal_models.shared_data_models import CausalModelScoringInput

        try:
            out = run_evo2_score(
                CausalModelScoringInput(sequences=[reference_sequence, mutant_sequence]),
                Evo2ScoringConfig(model_checkpoint="evo2_1b_base", device="modal"),
            )
            if not out.success or len(out.scores) < 2:
                result = {"error": f"Evo2 call failed: {out.errors}"}
            else:
                deviation = out.scores[1].avg_log_likelihood - out.scores[0].avg_log_likelihood
                result = {"evo2_deviation": round(deviation, 4), "model": "evo2_1b_base"}
        except Exception as e:
            result = {"error": f"Evo2 unavailable: {type(e).__name__}: {e}"}
        # Log the actual substitution, not just the request — this is what
        # makes it auditable that the right base was changed at the right
        # index, not just that a call was attempted.
        evidence_log.append({"category": "foundation_models", "tool": "score_dna_variant_evo2",
                              "args": {"chromosome": chromosome, "position": position,
                                       "reference_base_at_position": reference_sequence[idx] if 0 <= idx < len(reference_sequence) else None,
                                       "alternate_allele": alternate_allele, "window_start": ctx.start},
                              "result": result})
        return _tool_result_text(result)

    @beta_tool
    def fetch_protein_sequence(gene_symbol: str) -> str:
        """Fetch the real UniProt protein sequence for a gene in this species.

        Args:
            gene_symbol: Gene symbol, e.g. "FGF4".
        """
        from wildtype.tools.proto_client import fetch_uniprot_sequence

        seq = fetch_uniprot_sequence(gene_symbol, organism)
        result = {"error": f"No UniProt entry for {gene_symbol!r} in {organism!r}"} if seq is None else \
            {"sequence": seq, "length": len(seq)}
        evidence_log.append({"category": "genomic_databases", "tool": "fetch_protein_sequence",
                              "args": {"gene_symbol": gene_symbol}, "result": result})
        return _tool_result_text(result)

    @beta_tool
    def score_protein_variant_esm2(
        protein_sequence: str, position: int, mutant_amino_acid: str, source_citation: str
    ) -> str:
        """Score a protein missense change with ESM2. Returns a
        log-likelihood ratio (more negative = evolutionarily disfavored
        across species), NOT a disease probability.

        IMPORTANT — do not call this unless you can fill in source_citation
        honestly: there is no exon/reading-frame annotation tool available
        in this investigation, so `position` (amino acid) cannot be
        legitimately computed from the variant's DNA chromosome:position —
        that translation requires knowing exon boundaries and coding
        strand, which nothing here gives you. Only use this when a specific
        paper you retrieved via search_literature explicitly states the
        protein-level consequence (e.g. "p.Val362Met") — quote it in
        source_citation. If you don't have that, don't call this tool; note
        the missing protein-level evidence as a limitation instead of
        guessing a position.

        Reports scores from two independent ESM backends when both are
        available (Proto/Modal running full ESM2, and Biohub/Forge running
        ESMC) rather than picking one — treat rough agreement between them
        as more trustworthy than either alone, and a large disagreement as
        worth noting explicitly rather than silently preferring one.

        Args:
            protein_sequence: Wildtype amino acid sequence.
            position: 1-based amino acid position of the change, AS STATED in source_citation.
            mutant_amino_acid: Single-letter mutant amino acid code, AS STATED in source_citation.
            source_citation: The paper (title/DOI) and exact quoted notation that states this
                protein-level change. Required — not optional metadata.
        """
        from wildtype.tools.biohub_client import BiohubConfigError, ProtoScoringClient_Biohub, biohub_available
        from wildtype.tools.proto_client import ProtoScoringClient

        result: dict = {}
        try:
            score = ProtoScoringClient().score_missense("variant", protein_sequence, position, mutant_amino_acid)
            result["proto"] = {"log_likelihood_ratio": score.log_likelihood_ratio, "model": score.model}
        except Exception as e:
            result["proto"] = {"error": f"unavailable: {type(e).__name__}: {e}"}

        if biohub_available():
            try:
                score = ProtoScoringClient_Biohub().score_missense(
                    "variant", protein_sequence, position, mutant_amino_acid
                )
                result["biohub"] = {"log_likelihood_ratio": score.log_likelihood_ratio, "model": score.model}
            except BiohubConfigError as e:
                result["biohub"] = {"error": str(e)}
            except Exception as e:
                result["biohub"] = {"error": f"unavailable: {type(e).__name__}: {e}"}
        # source_citation logged alongside the result so this is auditable —
        # a missing/weak citation here is a real red flag on this finding,
        # not just bookkeeping.
        evidence_log.append({"category": "foundation_models", "tool": "score_protein_variant_esm2",
                              "args": {"position": position, "mutant_amino_acid": mutant_amino_acid,
                                       "source_citation": source_citation},
                              "result": result})
        return _tool_result_text(result)

    @beta_tool
    def search_literature(query: str, max_results: int = 5) -> str:
        """Search scientific literature (papers, preprints) for a query.

        Args:
            query: Natural-language search query.
            max_results: Max number of results to return.
        """
        from wildtype.tools.paperclip_client import PaperclipClient

        try:
            hits = PaperclipClient().search(query, max_results=max_results)
            result = [{"title": h.title, "doi": h.doi, "year": h.year, "summary": h.summary, "source": h.source}
                      for h in hits]
        except Exception as e:
            result = {"error": f"Literature search unavailable: {type(e).__name__}: {e}"}
        evidence_log.append({"category": "literature", "tool": "search_literature",
                              "args": {"query": query}, "result": result})
        return _tool_result_text(result)

    @beta_tool
    def search_gene_info(gene_symbol: str, organism_override: str | None = None) -> str:
        """Look up a gene's real NCBI record — function/description and which
        organism it's documented in. Use to check whether a locus falls near
        a known gene, and whether known biology comes from this animal's own
        species or a different one.

        Args:
            gene_symbol: Gene symbol to search for.
            organism_override: Search a different organism than this run's own
                species (e.g. check whether a gene is well-studied in humans but
                not yet in this species) — omit to use this run's species.
        """
        from proto_tools import run_ncbi_esearch, run_ncbi_esummary
        from proto_tools.tools.database_retrieval.ncbi.esearch import NCBIEsearchConfig, NCBIEsearchInput
        from proto_tools.tools.database_retrieval.ncbi.esummary import NCBIEsummaryConfig, NCBIEsummaryInput

        search_organism = organism_override or organism
        email = "wildtype-hackathon@example.com"
        try:
            search = run_ncbi_esearch(
                NCBIEsearchInput(db="gene", search_term=f"{gene_symbol}[Gene Name] AND {search_organism}[Organism]", max_results=3),
                NCBIEsearchConfig(ncbi_email=email),
            )
            if not search.success or not search.ids:
                result = {"error": f"No NCBI gene record for {gene_symbol!r} in {search_organism!r}"}
            else:
                summary = run_ncbi_esummary(
                    NCBIEsummaryInput(db="gene", identifier=",".join(search.ids)),
                    NCBIEsummaryConfig(ncbi_email=email),
                )
                result = []
                for uid in search.ids:
                    rec = summary.summary.get(uid, {})
                    # NCBI keeps retired gene records searchable across
                    # genome re-annotations, each with `currentid` pointing
                    # at the live replacement — live-caught 2026-08-15:
                    # ENPP1 came back as 3 UIDs, 2 of them discontinued
                    # duplicates of the 3rd, rendering as the same gene
                    # listed three times. Skip anything NCBI itself says is
                    # superseded rather than show it as separate evidence.
                    if rec.get("currentid"):
                        continue
                    result.append({
                        "name": rec.get("name"),
                        "description": rec.get("description"),
                        "organism": (rec.get("organism") or {}).get("scientificname"),
                        "summary": rec.get("summary", ""),
                    })
                if not result:
                    result = {"error": f"Only discontinued NCBI gene records for {gene_symbol!r} in {search_organism!r}"}
        except Exception as e:
            result = {"error": f"NCBI gene search unavailable: {type(e).__name__}: {e}"}
        evidence_log.append({"category": "genomic_databases", "tool": "search_gene_info",
                              "args": {"gene_symbol": gene_symbol, "organism": search_organism}, "result": result})
        return _tool_result_text(result)

    @beta_tool
    def search_genes_near_position(chromosome: str, position: int, window_bp: int = 100000) -> str:
        """Find real genes overlapping or near a genomic coordinate — use
        this FIRST for any variant, before search_gene_info, since you don't
        start out knowing a gene symbol for an arbitrary locus. Returns each
        gene's real NCBI coordinates so you can see how close the variant
        actually is to each one (a variant inside a gene's start/stop range
        is a very different situation from one merely near it).

        Args:
            chromosome: Chromosome/contig name as in the input file (e.g. "12").
            position: 1-based genomic position.
            window_bp: Search radius in base pairs around position (both directions).
        """
        from wildtype.tools.gene_lookup import find_genes_near_position

        genes, error = find_genes_near_position(organism, chromosome, position, window_bp, assembly=assembly)
        if error:
            result = {"error": error}
        elif not genes:
            result = {"genes_found": 0, "note": f"No NCBI gene records within {window_bp}bp of {chromosome}:{position}"}
        else:
            result = {
                "genes_found": len(genes),
                "genes": [
                    {
                        "name": g.name,
                        "description": g.description,
                        "start": g.start,
                        "stop": g.stop,
                        "variant_inside_gene_span": g.variant_inside_gene_span,
                    }
                    for g in genes
                ],
            }
        evidence_log.append({"category": "genomic_databases", "tool": "search_genes_near_position",
                              "args": {"chromosome": chromosome, "position": position, "window_bp": window_bp},
                              "result": result})
        return _tool_result_text(result)

    return [
        search_genes_near_position,
        fetch_reference_sequence,
        score_dna_variant_evo2,
        fetch_protein_sequence,
        score_protein_variant_esm2,
        search_literature,
        search_gene_info,
    ]


def _synthesize_from_evidence(variant: CommonVariant, evidence_log: list[dict], client: anthropic.Anthropic) -> dict | None:
    """Fallback for when the agent exhausts its iteration budget without
    calling create_risk_finding. The evidence it gathered along the way is
    real and shouldn't be thrown away for a canned "no conclusion" message
    — this makes one plain (non-tool) call asking Claude to synthesize a
    genuine finding from exactly what's in evidence_log, nothing more.
    Returns None if evidence_log is empty (nothing to synthesize from) or
    the call itself fails — caller falls back to the generic message then.
    """
    if not evidence_log:
        return None
    a1, a2 = variant.genotype_alleles()
    evidence_text = json.dumps(evidence_log, indent=2, default=str)
    prompt = (
        f"An investigation of this variant ran out of its tool-call budget before reaching "
        f"a conclusion. Here is every tool call it made and the real results returned:\n\n"
        f"Species: {variant.species} ({variant.chromosome}:{variant.position}, "
        f"ref {variant.reference}, genotype {a1}/{a2})\n\n{evidence_text}\n\n"
        "Synthesize the best honest conclusion from ONLY this evidence — do not invent "
        "anything not in it. If it doesn't support a real finding, say insufficient_evidence "
        "and explain what's missing, using what was actually gathered (don't just say "
        '"ran out of iterations" — say what the evidence itself shows or fails to show).\n\n'
        "Respond with ONLY a JSON object, no other text:\n"
        '{"risk_factor": "...", "confidence_level": "low|medium|high|insufficient_evidence", '
        '"confidence_score": null or 0-1, "confidence_rationale": "...", '
        '"species_relevance": "same_species|same_breed|related_species|different_species|unknown", '
        '"genes": [...], "biological_interpretation": "...", "limitations": [...], "agent_summary": "..."}'
    )
    try:
        resp = client.messages.create(
            model="claude-sonnet-5",
            thinking={"type": "disabled"},
            max_tokens=2000,
            messages=[{"role": "user", "content": prompt}],
        )
        text = next((b.text for b in resp.content if getattr(b, "type", None) == "text"), "").strip()
        if text.startswith("```"):
            text = text.strip("`")
            if text.startswith("json"):
                text = text[4:]
        return json.loads(text)
    except Exception:
        return None


def _build_evidence_bundle(evidence_log: list[dict]) -> EvidenceBundle:
    bundle = EvidenceBundle()
    for entry in evidence_log:
        item = {"tool": entry["tool"], "args": entry["args"], "result": entry["result"]}
        target = getattr(bundle, entry["category"], None)
        if target is not None:
            target.append(item)
    return bundle


def investigate_variant(
    variant: CommonVariant,
    organism: str,
    *,
    assembly: AssemblyInfo | None = None,
    client: anthropic.Anthropic | None = None,
    max_iterations: int = MAX_ITERATIONS,
) -> RiskFinding:
    """The actual iterative loop. One call per candidate variant — Claude
    decides everything from here: which tools, in what order, how many
    times, and when it has enough to conclude."""
    if client is None:
        key = os.environ.get("ANTHROPIC_API_KEY")
        if not key:
            raise RuntimeError("ANTHROPIC_API_KEY required — the iterative agent cannot run without Claude.")
        client = anthropic.Anthropic(api_key=key)
    assembly = assembly or resolve_reference_assembly(organism)

    evidence_log: list[dict] = []
    captured: dict = {}

    @beta_tool
    def create_risk_finding(
        risk_factor: str,
        confidence_level: str,
        biological_interpretation: str = "",
        agent_summary: str = "",
        genes: list[str] | None = None,
        confidence_score: float | None = None,
        confidence_rationale: str = "",
        species_relevance: str = "unknown",
        limitations: list[str] | str | None = None,
    ) -> str:
        """Report your final conclusion for this variant. Call exactly once,
        when you've investigated enough (or are near your tool-call budget)
        to state a conclusion — including "insufficient_evidence" if that's
        the honest one. Calling this ends the investigation.

        Args:
            risk_factor: Short description of the potential risk factor, or
                "No significant risk factor identified" if none found.
            confidence_level: One of "low", "medium", "high", "insufficient_evidence".
            biological_interpretation: What the evidence suggests mechanistically.
                REQUIRED IN SPIRIT even though technically optional here — see note below.
            agent_summary: 2-4 sentence plain-language summary for the animal's owner.
                REQUIRED IN SPIRIT even though technically optional here — see note below.
            genes: Gene symbol(s) involved, if any.
            confidence_score: Optional 0-1 score. NOT a clinical probability — omit if unsure.
            confidence_rationale: Why this confidence level, specifically (what evidence, what's missing).
            species_relevance: One of "same_species", "same_breed", "related_species",
                "different_species", "unknown" — is the strongest evidence from this
                animal's own species?
            limitations: Specific caveats/gaps in the evidence gathered, as a list —
                a single string is also accepted and will be treated as one item.
        """
        # biological_interpretation/agent_summary were made optional here (not
        # actually optional in the schema this was investigating against) after
        # live-catching repeated pydantic ValidationErrors on exactly these two
        # fields, 2026-08-15 — happened again even after bumping max_tokens
        # 4096->8192 and pinning thinking explicitly, so token-budget truncation
        # wasn't the (whole) cause. Rather than keep guessing at why Claude
        # sometimes omits them, this makes the tool call itself unable to
        # hard-fail on it: a blank value is recorded as a real gap (visible in
        # the finding's own limitations, not hidden), instead of raising a
        # pydantic error that burns an iteration on a retry that might fail
        # the same way again.
        if isinstance(limitations, str):
            # Claude sometimes serializes its own list to a JSON string
            # instead of passing an actual list — live-caught 2026-08-15:
            # 8 of 15 findings in one run came back with limitations as a
            # single string that was itself '["item 1", "item 2", ...]'.
            # Parse it back out rather than storing the raw JSON text as
            # one giant "limitation."
            try:
                parsed = json.loads(limitations)
                limitations_list = parsed if isinstance(parsed, list) else [limitations]
            except json.JSONDecodeError:
                limitations_list = [limitations]
        else:
            limitations_list = list(limitations or [])
        if not biological_interpretation.strip():
            limitations_list.append(
                "Agent called create_risk_finding without a biological_interpretation "
                "— recorded as a gap rather than retried indefinitely."
            )
        if not agent_summary.strip():
            agent_summary = risk_factor  # best available substitute, not blank
            limitations_list.append(
                "Agent called create_risk_finding without an agent_summary "
                "— falling back to risk_factor as the closest available summary."
            )
        captured["result"] = dict(
            risk_factor=risk_factor,
            confidence_level=confidence_level,
            confidence_score=confidence_score,
            confidence_rationale=confidence_rationale,
            species_relevance=species_relevance,
            genes=genes or [],
            biological_interpretation=biological_interpretation,
            limitations=limitations_list,
            agent_summary=agent_summary,
        )
        return "Recorded — investigation complete."

    tools = build_investigation_tools(variant.species, organism, assembly, evidence_log) + [create_risk_finding]

    a1, a2 = variant.genotype_alleles()
    prompt = (
        f"Species: {variant.species} (organism: {organism})\n"
        f"Sample: {variant.sample_id}\n"
        f"Locus: {variant.chromosome}:{variant.position}\n"
        f"Reference allele: {variant.reference}\n"
        f"Called genotype: {a1}/{a2} ({variant.zygosity})\n\n"
        "Investigate this variant using the tools available, then call "
        "create_risk_finding with your conclusion."
    )

    runner = client.beta.messages.tool_runner(
        model="claude-sonnet-5",
        # 4096 was too tight — live-caught 2026-08-15: repeated
        # create_risk_finding calls landing with biological_interpretation
        # and agent_summary (both required, both free-text) silently
        # missing, pydantic raising "Missing required argument" 3+ times in
        # a row for the same variant. Same failure mode as the earlier
        # report-truncation bug (thinking tokens eating the budget before
        # the model finishes writing a text-heavy tool call) — this call
        # never set `thinking` explicitly, so it ran on Sonnet 5's default
        # adaptive thinking with no headroom reserved for it. Bumped to
        # 8192 and thinking left adaptive (investigation genuinely
        # benefits from real reasoning across many tool calls, unlike the
        # simple yes/no decisions elsewhere that got thinking disabled).
        max_tokens=8192,
        thinking={"type": "adaptive"},
        system=AGENT_SYSTEM_PROMPT,
        tools=tools,
        messages=[{"role": "user", "content": prompt}],
    )

    iterations = 0
    for _message in runner:
        iterations += 1
        if "result" in captured or iterations >= max_iterations:
            break

    if "result" not in captured:
        synthesized = _synthesize_from_evidence(variant, evidence_log, client)
        if synthesized is not None:
            synthesized.setdefault("limitations", [])
            synthesized["limitations"].append(
                f"Investigation hit the {max_iterations}-iteration budget before calling "
                "create_risk_finding directly — this conclusion was synthesized from the "
                "evidence gathered up to that point, not stated by the agent itself."
            )
            captured["result"] = dict(
                risk_factor=synthesized.get("risk_factor", "Investigation incomplete"),
                confidence_level=synthesized.get("confidence_level", "insufficient_evidence"),
                confidence_score=synthesized.get("confidence_score"),
                confidence_rationale=synthesized.get("confidence_rationale", ""),
                species_relevance=synthesized.get("species_relevance", "unknown"),
                genes=synthesized.get("genes") or [],
                biological_interpretation=synthesized.get("biological_interpretation", ""),
                limitations=synthesized["limitations"],
                agent_summary=synthesized.get("agent_summary", ""),
            )
        else:
            # No evidence to synthesize from (budget exhausted before any
            # tool even returned), or the synthesis call itself failed —
            # this generic message is now the true last resort, not the
            # common case.
            captured["result"] = dict(
                risk_factor="Investigation incomplete",
                confidence_level="insufficient_evidence",
                confidence_score=None,
                confidence_rationale=f"Agent did not conclude within the {max_iterations}-iteration budget, "
                                      "and no usable evidence was gathered to synthesize a fallback conclusion from.",
                species_relevance="unknown",
                genes=[],
                biological_interpretation="",
                limitations=[f"Investigation stopped after {max_iterations} iterations without a create_risk_finding call."],
                agent_summary="No conclusion reached within the iteration budget.",
            )

    r = captured["result"]
    return RiskFinding(
        id=f"RF-{variant.chromosome}-{variant.position}",
        risk_factor=r["risk_factor"],
        confidence=ConfidenceAssessment(
            level=r["confidence_level"], score=r["confidence_score"], rationale=r["confidence_rationale"]
        ),
        species_relevance=r["species_relevance"],
        variants=[VariantRef(
            chromosome=variant.chromosome, position=variant.position,
            reference=variant.reference, alternate=variant.alternate, genotype=variant.genotype,
        )],
        genes=r["genes"],
        biological_interpretation=r["biological_interpretation"],
        evidence=_build_evidence_bundle(evidence_log),
        limitations=r["limitations"],
        agent_summary=r["agent_summary"],
    )
