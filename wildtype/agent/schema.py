"""Risk-finding schema — PROJECT_SPEC.md sections 11-14.

The structured output every agent run produces, replacing the old flat
VariantFinding/findings_to_dicts shape in genome_pipeline.py. Matches the
spec's field names directly (section 13) so JSON output aligns with what
James's spec/dashboard expect, rather than an independently-invented shape.

Enforces spec section 14's core rule in the type system, not just prose:
`ConfidenceLevel` and `SpeciesRelevance` are closed enums, not free text —
"low confidence" and "insufficient evidence" are first-class states, not
something the agent has to remember to say correctly every time.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Literal

ConfidenceLevel = Literal["low", "medium", "high", "insufficient_evidence"]
SpeciesRelevance = Literal["same_species", "same_breed", "related_species", "different_species", "unknown"]


@dataclass
class VariantRef:
    chromosome: str
    position: int
    reference: str
    alternate: str
    genotype: str


@dataclass
class ConfidenceAssessment:
    level: ConfidenceLevel
    score: float | None = None  # spec 12: optional, and NOT a clinical probability
    rationale: str = ""

    def to_dict(self) -> dict:
        return {"level": self.level, "score": self.score, "rationale": self.rationale}


@dataclass
class EvidenceBundle:
    foundation_models: list[dict] = field(default_factory=list)  # {"model":..., "score":..., "note":...}
    literature: list[dict] = field(default_factory=list)  # {"title":..., "doi":..., "year":..., "relevance":...}
    genomic_databases: list[dict] = field(default_factory=list)  # {"source":..., "detail":...}
    species_specific: list[dict] = field(default_factory=list)  # species-relevance notes per evidence item

    def to_dict(self) -> dict:
        return {
            "foundation_models": self.foundation_models,
            "literature": self.literature,
            "genomic_databases": self.genomic_databases,
            "species_specific": self.species_specific,
        }


@dataclass
class RiskFinding:
    id: str
    risk_factor: str
    confidence: ConfidenceAssessment
    variants: list[VariantRef]
    genes: list[str]
    biological_interpretation: str
    evidence: EvidenceBundle
    species_relevance: SpeciesRelevance = "unknown"
    limitations: list[str] = field(default_factory=list)
    agent_summary: str = ""

    def to_dict(self) -> dict:
        return {
            "id": self.id,
            "risk_factor": self.risk_factor,
            "confidence": self.confidence.to_dict(),
            "species_relevance": self.species_relevance,
            "variants": [vars(v) for v in self.variants],
            "genes": self.genes,
            "biological_interpretation": self.biological_interpretation,
            "evidence": self.evidence.to_dict(),
            "limitations": self.limitations,
            "agent_summary": self.agent_summary,
        }


@dataclass
class AnalysisRun:
    sample_id: str
    species: str
    input_source: str  # "vcf" | "tped_tfam"
    variants_analyzed: int
    findings: list[RiskFinding] = field(default_factory=list)
    generated_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    # Populated by genome_pipeline.run_genome_pipeline (None for the older
    # Embark-CSV path in agent/loop.py, which has no scan/triage stages) —
    # lets the dashboard show "N markers scanned -> M near a gene -> K
    # fully investigated" rather than just the final finding count, so a
    # small findings list reads as "most markers didn't overlap a gene",
    # not "the pipeline only looked at a few things."
    scan_stats: dict | None = None
    triage_stats: dict | None = None

    def to_dict(self) -> dict:
        return {
            "sample_id": self.sample_id,
            "species": self.species,
            "input_source": self.input_source,
            "variants_analyzed": self.variants_analyzed,
            "generated_at": self.generated_at,
            "findings": [f.to_dict() for f in self.findings],
            "scan_stats": self.scan_stats,
            "triage_stats": self.triage_stats,
            "summary": {
                "total_findings": len(self.findings),
                "confidence_counts": {
                    level: sum(1 for f in self.findings if f.confidence.level == level)
                    for level in ("high", "medium", "low", "insufficient_evidence")
                },
            },
        }

    def to_csv_rows(self) -> list[dict]:
        rows = []
        for f in self.findings:
            rows.append({
                "finding_id": f.id,
                "risk_factor": f.risk_factor,
                "confidence_level": f.confidence.level,
                "confidence_score": f.confidence.score,
                "species_relevance": f.species_relevance,
                "chromosome": f.variants[0].chromosome if f.variants else "",
                "position": f.variants[0].position if f.variants else "",
                "reference": f.variants[0].reference if f.variants else "",
                "alternate": f.variants[0].alternate if f.variants else "",
                "gene": ";".join(f.genes),
                "model_evidence_count": len(f.evidence.foundation_models),
                "literature_evidence_count": len(f.evidence.literature),
                "summary": f.agent_summary,
            })
        return rows
