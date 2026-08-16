# WildType

Species-agnostic genomic variant interpretation agent. re:AGENT Hackathon,
Track A. Demoed live on Ollie (Mini Australian Shepherd) — her real Embark
export ships with this repo under [`data/`](data/). Full proposal (research
question, tools, architecture, timeline) is in
[`docs/project_proposal.md`](docs/project_proposal.md); the working spec is
[`PROJECT_SPEC.md`](PROJECT_SPEC.md).

**Team:** Computational Bio PhD + ML Compiler Engineer (James)

## What it does

Takes a raw genomic export — PLINK `.tped`/`.tfam` (VCF also supported) —
and an LLM tool-calling agent investigates it from scratch: no pre-digested
vendor report as input, no fixed list of pre-validated conditions. Claude
decides which tools to call and when it has gathered enough evidence
(PROJECT_SPEC.md section 5.1: "the LLM is the agent, not the entire
system").

Ground truth here is the animal's own reference genome and real published
literature, not a vendor's internal validation — "mechanistically grounded"
describes what the models compute and what the literature says about a
gene, not a claim that any finding has been clinically validated.

**Why raw data instead of a health report:** a report only tells you what
the vendor chose to test for — Embark's is a fixed, pre-validated list of
255 conditions, dog-only. Raw genotype data lets the agent investigate
*everything*, and the same pipeline works for any species with genomic data
and an NCBI-annotated reference genome — validated this project against
dog, horse, and cheetah reference resolution, not just Ollie's data.

## Pipeline

```
raw .tped/.tfam or .vcf
  -> Stage 1  Input validation                  (parsers/tped.py, vcf.py)
  -> Stage 2  Reference resolution               (tools/reference_genome.py)
  -> Stage 3  Candidate-variant scan              (agent/genome_pipeline.py)
  -> Stage 3.5 Gene-proximity triage              (tools/gene_lookup.py)
  -> Stage 4  Iterative LLM investigation         (agent/iterative_agent.py)
                foundation models (Evo2, ESM2/Biohub), live literature
                search, real NCBI gene annotation — Claude decides what
                to call and when it has enough evidence
  -> Stage 5  Structured output                   (agent/schema.py)
              JSON/CSV + published HTML dashboard
```

Stage 3 does a real, principled filter — "genuinely deviates from the
species reference genome at this position" — not gene-focused prefiltering
or arbitrary sampling. Stage 3.5 narrows that further to positions actually
worth investigating (in or near an annotated gene). Both stages are cached
per-chromosome/per-assembly so a run doesn't pay a live network cost per
marker — see the module docstrings in `reference_genome.py` and
`gene_lookup.py` for the real numbers behind those fixes.

## Repo layout

```
wildtype/
  agent/
    genome_pipeline.py    scan -> triage -> investigate -> output, the real entry point
    iterative_agent.py     the Claude tool-calling agent (per-variant investigation)
    schema.py               RiskFinding / EvidenceBundle / AnalysisRun — PROJECT_SPEC's output shape
  parsers/
    tped.py, vcf.py         raw genotype format parsers
    common_variant.py       shared internal variant representation (PROJECT_SPEC section 7)
  tools/
    reference_genome.py     NCBI assembly/reference resolution, per-chromosome sequence cache
    gene_lookup.py          gene annotation lookup, backed by an assembly's own GFF3 file
    embark_panel_genes.py   Ollie/Embark-specific: gene list for POST-HOC validation only,
                             never fed into the agent's own investigation
    proto_client.py         ESM2/ESMFold via Proto (Modal)
    biohub_client.py        ESMC via Biohub/Forge — alternative path alongside Proto
    paperclip_client.py     live literature search
data/                    demo data, checked in — Ollie's real Embark export
  ollie_embark/             CSV, PDF, and raw tped/tfam
scripts/
  run_whole_genome_scan.py  the real, full-scope run — see "Run it" below
  patch_proto_gpu_tier.py   Modal GPU-tier patch for ESM2/ESMFold/Evo1 (see comments — not Evo2)
  biohub_forge_runner.py    runs inside an isolated venv (see Biohub setup below)
docs/
  project_proposal.md     the original pitch
PROJECT_SPEC.md           working spec — schema, open questions, architecture decisions
```

## Setup

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
cp .env.template .env   # fill in keys — see inline comments for each integration's real setup process
```

Real external dependencies for a live run: `ANTHROPIC_API_KEY` (the agent's
brain), `modal setup` (GPU compute for Evo2/ESM2/ESMFold — needs a payment
method on file for Evo2 specifically, see `wildtype/tools/proto_client.py`'s
module docstring), Paperclip login (`curl -fsSL https://paperclip.gxl.ai/install.sh
| bash && paperclip login`) for literature search. NCBI E-utilities and
UniProt need no auth. Biohub is optional and needs its own isolated venv —
see `wildtype/tools/biohub_client.py`'s module docstring.

## Run it

```bash
python3 scripts/run_whole_genome_scan.py
```

Real cost, not a demo shortcut: a full run against Ollie's 229,988-marker
file takes roughly an hour (live-measured: 63 minutes for a full scan +
40 full investigations) and makes real Claude, Modal, and NCBI calls. Each
fully-investigated variant costs real API time (~50-90s) — investigating
*every* gene-adjacent candidate genome-wide (~150K+ of them) isn't
computationally feasible in one sitting; `max_investigated` bounds a run to
what's actually runnable, with `stratify_by_chromosome=True` round-robining
the investigation budget across every chromosome the scan reaches instead
of letting one chromosome (typically chr1, first in file order) consume the
whole budget.

Output: `/tmp/wildtype_genome_findings.json` / `.csv`, checkpointed after
every completed investigation (not just at the end — a long run surviving
a crash midway still keeps everything already done). Build the dashboard
from that JSON with the template + font-embedding script used for the
published demo (ask if you need the exact assembly script — not checked in
since it lives in a scratch directory during development).

## A real correctness lesson worth keeping visible

Ollie's Embark TPED coordinates turned out to be mapped to **CanFam3.1**
(the historically standard canine reference), not **UU_Cfam_GSD_1.0** (what
NCBI currently designates as dog's "reference genome"). Auto-resolving
"the current reference genome" silently produced wrong results for an
entire session before this was caught — completely different DNA at the
same numeric position between the two assemblies. `resolve_reference_
assembly()` still auto-resolves NCBI's current pick by default (the right
behavior when you don't know better), but `resolve_assembly_by_accession()`
exists specifically so a caller who *has* determined which assembly their
input data actually uses can pin it explicitly. The dashboard now always
shows which assembly a run used, for exactly this reason. Genotyping
platforms map coordinates once at chip-design time and don't silently
follow a reference pointer that moves years later — this is a real,
recurring class of bug for anyone working with raw genomic exports, not an
Ollie-specific quirk.

## Validation

There's no clinical ground truth to check this against. What exists: Ollie's
real Embark report lists her 4 actual flagged conditions (FGF4 retrogene/
IVDD, PRCD carrier, GPT/ALT activity, ABCB1 clear) against a curated
255-condition panel. `wildtype/tools/embark_panel_genes.py` holds that gene
list, used **only** as a post-hoc comparison after a run completes — never
fed into candidate selection, triage, or investigation, since doing that
would make any "match" circular rather than a real validation signal.
Current honest result: broad genome-wide runs haven't yet landed on
Embark's exact panel genes by chance (expected — 204 genes is a tiny
fraction of the genome). A real apples-to-apples test — locating Ollie's
actual TPED markers at PRCD/GPT/ABCB1's real coordinates and running
independent investigation specifically there — is designed but not yet
built (see chat history / ask before starting: FGF4's Embark-relevant
locus is a retrogene *insertion*, not the native FGF4 gene, and likely
won't resolve the same way).

## Known limitations, stated plainly

- No file-upload UI — "input" today means pointing a script at a `.tped`/
  `.tfam` path, not a web upload widget.
- Evo2 (DNA-level scoring) needs a Modal payment method on file for its
  Hopper-tier GPU requirement — real architectural constraint, not a
  workaround-able default.
- No automatic assembly detection for arbitrary input data — see the
  correctness lesson above. Works correctly once the right assembly is
  known/pinned; doesn't yet figure that out on its own.
- `max_investigated` bounds every run — full genome-wide investigation at
  real API cost is not currently feasible in one sitting for any species
  with a large genome.
