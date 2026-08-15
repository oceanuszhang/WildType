# WildType

Species-agnostic genomic variant interpretation agent. re:AGENT Hackathon,
Track A, Aug 15–16 2026. Demoed live on Ollie (Mini Australian Shepherd) —
her real Embark report ships with this repo under [`data/`](data/). Full
proposal (research question, tools, architecture, timeline) is in
[`docs/project_proposal.md`](docs/project_proposal.md).

**Team:** Computational Bio PhD + ML Compiler Engineer (James)

## What it does

Takes a raw genomic export (Embark CSV today, others later) and produces a
two-layer report per flagged variant — plain-language for the pet owner,
full scores + citations for the researcher — by routing each variant through
protein-level (ESM2/ESM3/ESMFold) or DNA-level (Evo1) scoring, escalating
the highest-priority hits to AlphaFold3, and grounding every claim in live
literature via Paperclip.

## Repo layout

```
wildtype/            source package
  agent/                pipeline controller (loop.py) + Claude synthesis prompts (prompts.py)
  parsers/               Embark CSV + PLINK tped/tfam parsers
  tools/                  scoring/structure/DNA/literature clients — mock, local (real CPU ESM2), proto backends
  ui/                     Gradio app
data/                 demo data, checked in — see data/README.md
  ollie_embark/           Ollie's real Embark export (CSV, PDF, raw tped/tfam)
  basepaws_samples/       a Basepaws (cat) sample, for the species-agnostic path
docs/
  project_proposal.md    the original pitch — research question, tools, architecture, timeline
scripts/
  smoke_test.py          full pipeline, terminal output, run this first
tests/                run against Ollie's real data, not fixtures
```

## Setup

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
cp .env.template .env   # fill in keys as they become available; mock mode needs none
```

## Run it

```bash
python3 scripts/smoke_test.py      # full pipeline, terminal output
python3 -m wildtype.ui.app         # UI at http://127.0.0.1:7860
python3 -m pytest tests/ -q        # 8 tests, run against Ollie's real CSV
```

## Mode switch

Everything routes through `WILDTYPE_MODE` in `.env`:

| Mode | What runs | Needs |
|---|---|---|
| `mock` (default) | Canned realistic scores | Nothing — fully offline |
| `local` | **Real** ESM2 (35M, CPU) for protein scoring; structure/DNA/literature still mocked | Nothing — `pip install`s the weights on first run |
| `proto` | Real Proto (ESM2/3, ESMFold, AlphaFold3, Evo1, TM-align) + real Paperclip | Sponsor API keys |

Swapping modes doesn't change calling code — `wildtype/tools/base.py` is the
one place that decides which client backs each interface. See
`wildtype/tools/proto_client.py` and `paperclip_client.py` for what's still
a shape-only stub pending real API docs (each has inline `TODO`s).

## What's still a placeholder (by design, not oversight)

- **Amino acid sequences** (`wildtype/tools/placeholder_sequences.py`) — deterministic fake sequences, clearly marked. Real UniProt fetch is a Proto call.
- **Variant position/mutant allele** — Embark's CSV gives a risk *label*, not the underlying position/allele; that comes from VEP-style annotation of the raw VCF, a separate ingestion path.
- **`KNOWN_DISEASE_POSITIONS`** in `wildtype/parsers/tped.py` — gene-coordinate table for the beyond-Embark scan, empty until populated from OMIA/Paperclip lookups.
- **Score thresholds** (`ESM2_FLAG_LLR_THRESHOLD`, `EVO1_FLAG_DEVIATION_THRESHOLD` in `wildtype/agent/loop.py`) — starting points to calibrate against the validation set, not final numbers.

## Validation set

Ollie's four Embark calls are the labeled ground truth — get all four right
before trusting the pipeline on anything novel:

| Gene | Embark call | Tests |
|---|---|---|
| FGF4 retrogene (CFA12) | At risk | DNA-level route (Evo1) |
| GPT | At risk | Protein-level route (ESM2) |
| PRCD | Carrier | Zygosity gate downgrades severity correctly |
| MDR1/ABCB1 | Clear | Pipeline doesn't false-positive |

## Architecture choice worth knowing about

`wildtype/agent/loop.py` uses a **deterministic Python controller** that
calls Claude only for two things: the AlphaFold3 escalation tradeoff and
final report synthesis — not a fully autonomous Claude tool-calling loop.
Chosen for demo reliability (predictable, debuggable under time pressure at
a live table), not a neutral default. Worth a two-minute conversation
between PhD and James before the pitch if you want the "real" agentic loop
for optics instead.

## Must-have vs. nice-to-have

**Must have:** Claude agent loop · Embark CSV parser · Paperclip literature
search per variant · ESM2 scoring via Proto · two-layer output · working UI
with file upload.

**Nice to have:** Evo1 DNA-level scoring for FGF4 · ESMFold structure
visualization · tped novel-variant discovery · breed-name-only mode ·
AlphaFold3 for the top variant.
