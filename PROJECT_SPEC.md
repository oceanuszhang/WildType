# PROJECT_SPEC.md

# Animal Genomic Risk Analysis Agent

## 1. Project Overview

Build a minimal working prototype of an AI agent that analyzes animal genomic variant data and identifies potentially relevant genetic risk factors.

The system should:

1. Accept pre-called genomic genotype/variant data.
2. Support VCF and PLINK TPED/TFAM inputs for the MVP.
3. Use an LLM-based agent to orchestrate analysis and research.
4. Incorporate genomic/protein foundation models such as Evo 1, Evo 2, and ESM where useful.
5. Research relevant scientific literature and genomic resources.
6. Combine computational, genomic, and literature evidence.
7. Produce risk factors with an evidence-backed confidence assessment.
8. Output structured JSON/CSV results.
9. Provide a dashboard for exploring the findings.

The MVP is an interpretation and research system. It is **not** intended to diagnose disease or provide definitive clinical/veterinary predictions.

---

## 2. Core Product Goal

Given:

- an animal species,
- genomic variant/genotype data,
- and optional sample metadata,

the agent should investigate potentially relevant genetic risk factors and produce an evidence-backed report.

A risk finding should answer:

> What potentially relevant genetic risk factor was found, which variant/gene is associated with it, what evidence supports the finding, and how confident are we in the interpretation?

---

## 3. MVP Scope

### Supported Inputs

The MVP supports only:

- VCF
- TPED + TFAM

Embark TPED/TFAM files may be used as development/demo data, but the architecture must **not** be specific to Embark.

The system should treat Embark as one possible source of genotype data among many.

### Future Inputs

Explicitly out of scope for the MVP:

- FASTQ
- BAM
- CRAM
- raw sequencing reads
- de novo genome assembly
- alignment pipelines
- variant calling from sequencing reads

These can be future extensions.

---

## 4. MVP Architecture

```text
                         ┌──────────────────────┐
                         │        USER          │
                         │                      │
                         │ Species              │
                         │ VCF OR TPED/TFAM     │
                         │ Optional metadata    │
                         └──────────┬───────────┘
                                    │
                                    ▼
                         ┌──────────────────────┐
                         │    INPUT PARSER      │
                         │                      │
                         │ VCF parser           │
                         │ TPED/TFAM parser     │
                         └──────────┬───────────┘
                                    │
                                    ▼
                         ┌──────────────────────┐
                         │  COMMON VARIANT      │
                         │  REPRESENTATION      │
                         │                      │
                         │ Sample               │
                         │ Species              │
                         │ Chromosome           │
                         │ Position             │
                         │ Ref / Alt            │
                         │ Genotype             │
                         └──────────┬───────────┘
                                    │
                                    ▼
                  ┌─────────────────────────────────────┐
                  │          AGENT / LLM                │
                  │                                     │
                  │ Planning                            │
                  │ Reasoning                           │
                  │ Tool selection                      │
                  │ Evidence synthesis                  │
                  └───────────────┬─────────────────────┘
                                  │
             ┌────────────────────┼─────────────────────┐
             │                    │                     │
             ▼                    ▼                     ▼
     ┌───────────────┐    ┌───────────────┐    ┌───────────────┐
     │ Foundation    │    │ Research      │    │ Genomic       │
     │ Models        │    │ Tools         │    │ Databases     │
     │               │    │               │    │               │
     │ Evo 1         │    │ Literature    │    │ Variant info  │
     │ Evo 2         │    │ Web search    │    │ Gene info     │
     │ ESM           │    │ Papers        │    │ etc.          │
     └───────┬───────┘    └───────┬───────┘    └───────┬───────┘
             │                    │                    │
             └────────────────────┼────────────────────┘
                                  ▼
                       ┌──────────────────────┐
                       │  EVIDENCE SYNTHESIS  │
                       │                      │
                       │ Computational        │
                       │ Literature           │
                       │ Genomic              │
                       │ Species-specific     │
                       └──────────┬───────────┘
                                  │
                                  ▼
                       ┌──────────────────────┐
                       │   RISK ASSESSMENT    │
                       │                      │
                       │ Risk factors         │
                       │ Confidence           │
                       │ Supporting evidence  │
                       │ Limitations          │
                       └──────────┬───────────┘
                                  │
                         ┌────────┴────────┐
                         ▼                 ▼
                  ┌────────────┐    ┌────────────┐
                  │ JSON / CSV │    │ Dashboard  │
                  └────────────┘    └────────────┘
```

---

## 5. Important Architectural Principles

### 5.1 The LLM is the agent, not the entire system

The LLM should orchestrate specialized tools rather than perform every operation itself.

The agent can:

- decide what to investigate,
- choose tools,
- formulate research questions,
- interpret tool results,
- identify missing evidence,
- perform additional investigation,
- synthesize evidence,
- generate the final interpretation.

Deterministic code should handle:

- file parsing,
- data validation,
- structured transformations,
- model inference interfaces,
- database/API calls,
- output serialization.

### 5.2 Foundation models provide evidence, not final truth

Evo/ESM outputs should be treated as computational evidence.

Do not automatically convert:

```text
Evo score = 0.90
```

into:

```text
Disease risk = 90%
```

Instead distinguish:

- model-derived effect score,
- literature evidence,
- genomic/database evidence,
- species relevance,
- overall confidence.

### 5.3 Evidence must be traceable

Every important claim in a risk finding should ideally be linked to its source.

Sources can include:

- scientific publications,
- genomic databases,
- species-specific databases,
- foundation-model outputs,
- other authoritative resources.

### 5.4 The system should expose uncertainty

The agent should be able to say:

- evidence is strong,
- evidence is limited,
- evidence is conflicting,
- evidence exists only in another species,
- computational evidence is suggestive but not experimentally validated,
- insufficient evidence.

---

# 6. Input Design

## 6.1 VCF

The VCF adapter should extract at minimum:

- sample ID
- chromosome
- position
- reference allele
- alternate allele
- genotype
- variant ID when available

Potentially useful fields:

- QUAL
- FILTER
- INFO
- FORMAT
- allele frequency
- depth
- genotype quality

Do not assume that every VCF contains all optional fields.

## 6.2 TPED + TFAM

The TPED/TFAM adapter should extract:

### TPED

- chromosome
- marker ID
- genetic distance
- position
- allele calls
- genotype for each sample

### TFAM

- family ID
- individual ID
- paternal ID
- maternal ID
- sex
- phenotype

The implementation should validate that the genotype columns in TPED correspond correctly to the samples described by TFAM.

## 6.3 Species

Species should be an explicit input.

Example:

```json
{
  "species": "Canis lupus familiaris"
}
```

Do not infer species solely from genomic coordinates unless explicitly required.

Optional metadata can include:

- breed
- sex
- population
- phenotype
- age
- sample identifier
- reference genome/assembly

The system should preserve metadata but should not require all metadata for the MVP.

---

# 7. Common Variant Representation

After parsing VCF or TPED/TFAM, convert data into a common internal representation.

Example:

```json
{
  "sample_id": "animal_001",
  "species": "Canis lupus familiaris",
  "chromosome": "chr1",
  "position": 123456,
  "reference": "A",
  "alternate": "G",
  "genotype": "0/1"
}
```

The exact schema can evolve after testing real VCF and TPED/TFAM datasets.

Do not over-engineer the common representation in the first implementation.

---

# 8. Agent Workflow

The initial agent workflow should be:

```text
Input
  ↓
Parse
  ↓
Identify candidate variants / regions
  ↓
Agent investigates
  ↓
Foundation model analysis
  ↓
Literature research
  ↓
Genomic database research
  ↓
Evidence synthesis
  ↓
Risk assessment
  ↓
Confidence assessment
  ↓
Structured output
```

The agent should be capable of iterating.

Conceptually:

```text
Candidate variant
      ↓
     Agent
      ↓
"Is this potentially relevant?"
      ↓
    Investigate
      │
      ├── Foundation model
      ├── Literature
      └── Database
      │
      ▼
Evidence
      ↓
Agent evaluates evidence
      ↓
"Do I have enough evidence?"
      │
      ├── No → investigate further
      │
      └── Yes → produce finding
```

The first implementation can use a bounded number of iterations to control cost and complexity.

---

# 9. Foundation Model Integration

Foundation models should be exposed to the agent through well-defined tools.

## 9.1 Evo 1 / Evo 2

Use for DNA/genomic sequence analysis.

Potential use cases:

- variant effect prediction,
- sequence likelihood comparison,
- genomic sequence representations,
- potentially regulatory/noncoding sequence analysis.

The exact supported operations depend on the selected model implementation.

Conceptually:

```python
evo2_predict_variant(
    reference_sequence,
    alternate_sequence,
    context
)
```

Return structured data rather than raw model output.

Example:

```json
{
  "model": "Evo2",
  "variant": "chr1:123456:A:G",
  "effect_score": 0.87,
  "interpretation": "..."
}
```

## 9.2 ESM

Use for protein-level analysis when a variant has a relevant protein consequence.

Conceptual flow:

```text
DNA variant
  ↓
coding consequence
  ↓
amino-acid change
  ↓
protein sequence
  ↓
ESM
  ↓
protein-level evidence
```

Example interface:

```python
esm_predict_variant(
    protein_sequence,
    mutation
)
```

## 9.3 Model execution

Foundation models may be run using Modal or another GPU execution environment.

Potential architecture:

```text
Agent
  │
  ├── evo2_tool()
  │       ↓
  │     Modal
  │       ↓
  │      GPU
  │       ↓
  │     Evo 2
  │
  └── esm_tool()
          ↓
        Modal
          ↓
         GPU
          ↓
         ESM
```

The agent should not need to know the underlying infrastructure.

Model wrappers should hide:

- model loading,
- GPU configuration,
- preprocessing,
- inference,
- output formatting,
- logging.

---

# 10. Research Tools

The agent should have access to tools for external evidence gathering.

Potential tools:

```text
search_literature()
retrieve_paper()
search_species_information()
search_gene_information()
search_variant_information()
search_genomic_database()
```

The exact APIs/databases are intentionally not fixed yet.

The agent should be able to formulate targeted questions such as:

- What is the function of this gene?
- Has this variant been reported?
- Is this gene associated with a disease or phenotype?
- Is the association established in this species?
- Is the association specific to another species?
- Are there conflicting studies?
- What is the strength of the evidence?
- What biological mechanism has been proposed?

---

# 11. Evidence Model

Each risk finding should contain multiple evidence categories.

## 11.1 Computational evidence

Examples:

- Evo 1 score
- Evo 2 score
- ESM score
- other computational predictions

## 11.2 Genomic evidence

Examples:

- variant databases
- gene annotations
- allele frequencies
- conservation
- known variant annotations

## 11.3 Literature evidence

Examples:

- peer-reviewed papers
- experimental studies
- association studies
- functional studies
- review papers

## 11.4 Species relevance

Explicitly distinguish:

```text
same species
same breed
related species
different species
unknown
```

Evidence from another species should not automatically be treated as equivalent to evidence from the target species.

---

# 12. Confidence Assessment

Every risk factor should have a confidence assessment.

At minimum:

```text
Low
Medium
High
```

Optionally also provide a numerical score.

Example:

```json
{
  "confidence": {
    "level": "high",
    "score": 0.86
  }
}
```

The numerical score should not be presented as a clinical probability unless there is a validated statistical model supporting that interpretation.

Confidence should consider factors such as:

- number and quality of supporting sources,
- agreement between independent evidence types,
- strength of literature evidence,
- species relevance,
- quality of genomic evidence,
- agreement among computational models,
- conflicting evidence,
- whether the finding is directly observed or inferred.

---

# 13. Risk Finding Schema

Initial conceptual schema:

```json
{
  "id": "RF-001",

  "risk_factor": "Potential increased risk of ...",

  "confidence": {
    "level": "high",
    "score": 0.86
  },

  "variants": [
    {
      "chromosome": "chr1",
      "position": 123456,
      "reference": "A",
      "alternate": "G",
      "genotype": "0/1"
    }
  ],

  "genes": [
    "GENE1"
  ],

  "biological_interpretation": "...",

  "evidence": {
    "foundation_models": [],
    "literature": [],
    "genomic_databases": [],
    "species_specific": []
  },

  "limitations": [],

  "agent_summary": "..."
}
```

This is a starting schema and should be refined during implementation.

---

# 14. Important Terminology

Avoid conflating these concepts:

```text
Model effect score
        ≠
Functional effect
        ≠
Disease association
        ≠
Disease probability
        ≠
Overall confidence
```

For example:

> A foundation model may predict that a variant has a strong sequence-level effect. This does not by itself establish that the animal has an increased disease risk.

The agent should preserve this distinction in both machine-readable output and dashboard explanations.

---

# 15. Output

## JSON

JSON is the canonical structured output.

It should contain:

- sample information
- species
- input information
- analysis metadata
- risk findings
- evidence
- confidence
- limitations

## CSV

CSV should provide a flattened view suitable for:

- spreadsheets,
- downstream analysis,
- filtering,
- sorting.

Potential columns:

```text
finding_id
risk_factor
confidence_level
confidence_score
chromosome
position
reference
alternate
gene
model_evidence
literature_evidence_count
species_relevance
summary
```

## Dashboard

The dashboard should provide:

### Overview

- species
- sample ID
- number of variants analyzed
- number of findings
- high/medium/low confidence counts

### Risk categories

For example:

- cardiovascular
- metabolic
- neurological
- immune
- cancer
- other

Categories should be determined by the evidence rather than hardcoded to a particular animal species.

### Findings table

Show:

- risk factor
- gene
- variant
- confidence
- evidence strength

### Finding detail

Clicking a finding should show:

- variant information
- gene information
- biological interpretation
- foundation-model results
- literature evidence
- genomic database evidence
- species relevance
- conflicting evidence
- limitations

---

# 16. MVP vs Future Work

## MVP

- VCF input
- TPED/TFAM input
- basic parsing/validation
- common variant representation
- LLM agent
- research tools
- genomic database tools
- Evo 1/Evo 2 integration where practical
- ESM integration where practical
- evidence synthesis
- Low/Medium/High confidence
- JSON output
- CSV output
- basic dashboard

## Future Work

### Additional genomic inputs

- FASTQ
- BAM
- CRAM
- structural variants
- additional genotyping formats

### Bioinformatics pipeline

- sequencing QC
- alignment
- variant calling
- advanced variant normalization
- genome assembly support

### More foundation models

- additional DNA foundation models
- additional protein foundation models
- specialized phenotype/disease models
- ensemble models

### More advanced statistical modeling

- calibrated risk probabilities
- population-specific risk models
- polygenic risk scores
- haplotype analysis
- multi-variant models

### More advanced agent behavior

- autonomous hypothesis generation
- iterative literature review
- conflicting-evidence resolution
- automatic experiment/data requests
- self-evaluation
- multi-agent research

---

# 17. Development Philosophy

### Start simple

The first goal is a working end-to-end prototype.

Do not build a complete bioinformatics platform.

### Make evidence traceable

Every important conclusion should be connected to supporting evidence.

### Keep model calls modular

Foundation models should be replaceable.

For example:

```text
FoundationModel
    ├── Evo1
    ├── Evo2
    ├── ESM
    └── Future models
```

### Keep the agent model-agnostic

The orchestration layer should not depend tightly on a specific LLM provider.

Possible LLMs:

- OpenAI models
- Anthropic Claude
- other compatible foundation models

### Separate inference from interpretation

Foundation models produce computational evidence.

The agent interprets that evidence together with external evidence.

### Avoid unsupported conclusions

If evidence is insufficient, the correct output is:

```text
Insufficient evidence
```

rather than a speculative high-confidence finding.

---

# 18. Current Open Questions

These should be resolved during the next design stages.

1. What exact fields should the common variant representation contain?
2. How should VCF and TPED/TFAM variants be normalized?
3. How do we handle different reference genome assemblies?
4. Which animal species should be supported first?
5. Which genomic databases should be used?
6. Which literature/search APIs should be used?
7. Which versions of Evo 1/Evo 2 should be integrated?
8. Which ESM model should be integrated?
9. Which variants should be sent to foundation models?
10. Should the agent analyze every variant or prioritize candidates first?
11. How should computational evidence be quantitatively combined with literature evidence?
12. How should confidence be calibrated?
13. What should happen when sources disagree?
14. How should species-specific evidence be weighted?
15. Which LLM should power the agent?
16. Which agent framework should be used, if any?
17. Should Modal host foundation-model inference, and how should models be deployed?
18. What dashboard framework should be used?

---

# 19. Immediate Next Design Step

The next step is to define the **common genomic data model** and the **agent tool interface**.

Specifically:

```text
VCF / TPED+TFAM
       ↓
Variant object
       ↓
Agent
       ↓
What tools are available?
       ↓
What arguments does each tool accept?
       ↓
What does each tool return?
       ↓
How does the agent decide what to investigate?
```

This should be designed before implementing the full agent.

The first concrete deliverable after this specification should therefore be:

**A tool/API specification for the agent**, covering:

- `parse_vcf`
- `parse_tped_tfam`
- `evo1_predict`
- `evo2_predict`
- `esm_predict`
- `search_literature`
- `search_genomic_database`
- `get_gene_information`
- `get_variant_information`
- `create_risk_finding`
- `export_results`

The implementation can then be built incrementally around those interfaces.
