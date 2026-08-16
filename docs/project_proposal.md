# re:AGENT Hackathon — Project Proposal
**Team**: Computational Bio PhD + ML Compiler Engineer
**Track**: Track A — Build an AI Scientist
**Event**: re:AGENT, San Francisco · Aug 15–16, 2026

---

## Topic

A species-agnostic genomic variant interpretation agent that takes raw genomic data or consumer genomics reports (e.g. Embark, Basepaws) and produces multi-scale, mechanistically grounded health risk analysis — for any animal species. Demonstrated live on Ollie, our Mini Australian Shepherd.

---

## Research Question

> Can an agent explain the *mechanism* behind a genomic finding — not just label it — for a species that falls outside every existing interpretation tool's training and validation set?

Human-centric tools (VEP, AlphaMissense) don't transfer to non-human genomes, and the species-agnostic foundation models that could work (ESM2, Evo1) have no agentic pipeline chaining them into literature-grounded, two-audience output. That gap — not the multi-model pipeline by itself — is the actual bet we're testing. We demo it live, on our own dog.

*Scope note: "mechanistically grounded" is a claim about what the models compute (log-likelihood, structural deviation, functional embedding) and what the literature says about the gene — it is not a claim of clinical validation. The only ground truth we validate against is Embark's own genotype calls; we have no clinical outcome data to check mechanism-to-outcome accuracy against, and we should be precise about that distinction when judges ask.*

---

## Tools Used

| Tool | Role |
|---|---|
| **Paperclip** | Live literature search across 11M+ papers, UniProt, PDB — grounds every variant finding in current evidence with citations |
| **Proto** | Unified MCP interface to run ESMFold, AlphaFold3, Evo1, ProteinMPNN, TM-align, UniProt retrieval — the computational biology execution layer |
| **Claude** | Agent brain — orchestrates the pipeline, reasons about tradeoffs between model outputs, generates the two-layer report |
| ~~Biohub / ESM~~ | Dropped 2026-08-15 — `proto_tools` (Proto) already exposes `run_esm2_score`/`run_esm3_score`/`run_esmc_embeddings` directly; a second integration to the same models added no capability |
| ~~BenchFlow~~ | Dropped 2026-08-15 — `github.com/benchflow-ai` is confirmed the correct org (AI-agent evaluation/benchmarking infra: SkillsBench, ClawsBench), it's just not something this pipeline needs; descoped as out of scope, not a mismatch |
| **Benchling** | Structured storage and annotation of variant analysis results |

---

## Foundation Models

### Evo1 — Arc Institute
DNA-level genomic foundation model trained on diverse genomic sequences across species. Unlike protein models, Evo1 operates directly on nucleotide sequences — making it the right tool for structural variants like insertions and deletions that can't be meaningfully scored at the protein level. For Ollie's FGF4 retrogene insertion (the cause of her IVDD risk), Evo1 scores the DNA-level context of the insertion directly. Arc Institute is a co-host of this hackathon; this is their flagship model.

### ESM2 — Meta AI / Biohub
650M–3B parameter protein language model trained on evolutionary sequence data across millions of proteins from diverse species. Computes log-likelihood of a mutant sequence vs. wildtype — a well-validated, quantitative proxy for functional impact. Fast, cheap, species-agnostic. Used as the first-pass scoring filter for all protein-level variants.

### ESM3 — Biohub
Multimodal protein foundation model reasoning over sequence, structure, and function simultaneously. Goes beyond ESM2's single likelihood score to place the mutant in functional embedding space — revealing not just "is this variant bad" but "what function does it disrupt." Biohub is a co-host of this hackathon.

### ESMFold — Meta AI / Biohub
Fast structure prediction directly from sequence, no MSA required. Returns predicted 3D coordinates and per-residue confidence (pLDDT). Used for rapid structural analysis of surviving variants after ESM2 filtering. Much faster than AlphaFold3, sufficient for most variants.

### AlphaFold3 — Google DeepMind (via Proto)
Highest-accuracy structure prediction, including protein-ligand and protein-nucleic acid complexes. Reserved selectively for the top 1–2 highest-priority variants where structural precision matters most. Computationally expensive; used only after ESM2 and ESMFold have narrowed the candidate set.

---

## What Makes This Unique

### vs. Ensembl VEP
VEP annotates variants with consequence labels (missense, stop gained) — it's a bioinformatics annotation tool requiring command-line expertise, producing tables not explanations. No protein scoring, no structure prediction, no literature search, no natural language output. Our agent calls VEP-style annotation as one step, then layers everything VEP cannot do on top.

### vs. OMIA (Online Mendelian Inheritance in Animals)
A curated lookup database of known animal disease variants — static, manual, no computation. You search it; it doesn't reason. Our agent queries OMIA as a data source, then computes mechanistic scores and synthesizes current literature around the hits.

### vs. AlphaMissense
Human proteins only. Missense SNPs only. No explanation, no structure prediction, no literature. Cannot handle Ollie's FGF4 insertion at all.

### vs. Embark / Basepaws
Consumer genomics services that run a microarray and return risk scores and labels. They tell you *what* — our agent tells you *why*. Embark flagged Ollie as "at risk" for IVDD. Our agent explains what the FGF4 retrogene insertion does to the protein, scores it at the DNA level with Evo1, predicts structural consequences, and surfaces the supporting literature — mechanistic depth Embark has no capability to provide.

### The genuine gap this fills
No existing tool combines species-agnostic foundation model scoring + structure prediction + live literature retrieval in one agentic, natural-language pipeline. For non-human species, this gap is total: the human-centric tools don't transfer, and the bioinformatics tools require expertise and produce no explanation. This agent works for any species with a sequenced genome — companion animals today, livestock and conservation biology as the natural extension.

---

## Architecture

```
INPUT
  Embark/Basepaws CSV  →  parse flagged variants
  Breed name only      →  Paperclip: fetch known disease variants for breed
  Raw VCF / tped       →  PLINK scan → cross-reference known disease gene positions

                    ↓ per variant, run in parallel

┌─ Paperclip ──────────────→ literature: mechanism, clinical evidence, citations
├─ Proto: UniProt ─────────→ wildtype protein sequence
├─ Proto: ESM2 ────────────→ log-likelihood score (missense variants)
├─ Proto: Evo1 ────────────→ DNA-level score (insertions, deletions, structural)
├─ Proto: ESMFold ─────────→ predicted mutant structure
├─ Proto: TM-align ────────→ structural deviation vs. wildtype
└─ Proto: AlphaFold3 ──────→ high-accuracy fold (top priority variants only)

                    ↓ Claude synthesizes all signals

OUTPUT: Two-layer report
  Plain layer   → clinical significance, what to tell your vet, what to monitor
  Science layer → ESM2 score, Evo1 score, pLDDT, RMSD, citations with DOIs
```

---

## Demo: Ollie's Embark Data

Ollie (Mini Aussie, female) has three flagged variants from her Embark report:

| Variant | Gene | Embark Call | What our agent adds |
|---|---|---|---|
| IVDD / Chondrodystrophy | FGF4 retrogene, chr12 | At risk, dominant het | Evo1 DNA-level score of the retrogene insertion; ESMFold structure; activity/weight management guidance |
| Elevated ALT Activity | GPT | At risk, codominant | ESM2 score; ESM3 functional embedding; bloodwork interpretation guidance |
| Progressive Retinal Atrophy | PRCD | Carrier | ESM2 score; breeding implications; no clinical action for Ollie herself |
| Drug Sensitivity | MDR1/ABCB1 | **Clear** ✅ | Confirmed from raw tped markers; no drug restrictions needed |

**Bonus**: tped file (229,988 SNP markers, 99.8% call rate) enables a beyond-Embark scan — finding variants at known disease gene positions that Embark didn't include in their health panel.

---

## Timeline

### Day 1 — Saturday Aug 15

| Time | Task | Owner |
|---|---|---|
| 8:30–10:30am | Check-in, lightning talks, team setup, API keys | Both |
| 10:30am–12pm | Agent skeleton (Claude loop), Embark CSV parser, environment setup | Engineer |
| 10:30am–12pm | Define variant prioritization logic, validate Ollie's findings | PhD |
| 12–3pm | Paperclip integration + ESM2 scoring via Proto | Engineer |
| 12–3pm | Write Paperclip queries, validate ESM2 scores biologically | PhD |
| 3–6pm | ESMFold structure prediction + output renderer (two layers) | Engineer |
| 3–6pm | Interpret structures, write plain-language output templates | PhD |
| 7–10pm | Evo1 for FGF4 insertion; tped novel variant scan; Benchling integration | Both |

### Day 2 — Sunday Aug 16

| Time | Task | Owner |
|---|---|---|
| 8:30–9am | Breakfast, final priorities | Both |
| 9–10:30am | UI polish, breed-name-only input mode, bug fixes | Engineer |
| 9–10:30am | Validate all outputs, prepare demo script | PhD |
| 10:30–11am | End-to-end test with Ollie's data, rehearse demo | Both |
| 11am | **Submit** | Both |
| 11am–2:30pm | Project expo, live judging demo | Both |

---

## Must-Have vs. Nice-to-Have

**Must have (core demo):**
- Claude agent loop
- Embark CSV parser
- Paperclip literature search per variant
- ESM2 log-likelihood scoring via Proto
- Two-layer output (plain + scientific)
- Working UI with file upload

**Nice to have (if time):**
- Evo1 DNA-level scoring for FGF4
- ESMFold structure visualization
- tped novel variant discovery
- Breed-name-only mode (no Embark needed)
- AlphaFold3 for top variant

---

## UI Recommendation

**Gradio on Modal** — fastest to build (hours, not days), looks professional for ML demos, and runs entirely on Modal (hackathon sponsor infrastructure). Strong optics for judges.

Layout:
```
┌─────────────────────────────────────────────────────┐
│  [Upload Embark CSV]   OR   [Enter breed name]      │
│  [Species: Dog / Cat / Other]  [Run Analysis]        │
├─────────────────────────────────────────────────────┤
│  RESULTS                                            │
│  ┌──────────────────┐  ┌──────────────────────────┐ │
│  │  Pet Owner View  │  │  Researcher View         │ │
│  │  Plain English   │  │  ESM2 score, pLDDT,      │ │
│  │  + vet guidance  │  │  Evo1 score, citations   │ │
│  └──────────────────┘  └──────────────────────────┘ │
│  [3D Structure Viewer — NGL.js embedded]             │
│  [Literature citations with DOI links]               │
└─────────────────────────────────────────────────────┘
```

---

## Hosting — To Impress Judges

**Primary: Modal + Gradio**
Modal is a hackathon sponsor. Hosting the entire compute backend on Modal (which is also how Proto runs) means the demo is fully live, not running on a laptop. Judges from BenchFlow, Anthropic, and Arc will notice you're using sponsor infrastructure end-to-end. Deploy with `modal deploy` in minutes.

**Frontend option: Vercel**
If time allows, a clean Next.js frontend deployed on Vercel looks more polished than Gradio. Vercel free tier, deploys in seconds, custom domain possible.

**Fallback: Claude.ai artifact**
If anything breaks, an HTML artifact running in Claude itself still demonstrates the pipeline and is live on any browser — zero infrastructure risk.

**Recommended stack:**
```
Claude (agent) → Proto (tools on Modal) → Paperclip (API) 
      ↓
Gradio UI deployed on Modal
      ↓  
Live URL shared with judges
```

---

## One-Line Pitch

> "VEP tells you what a variant is. We built the agent that tells you what it means — computationally, structurally, and mechanistically — for any species. Demonstrated on our own dog."
