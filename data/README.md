# Demo data

Real consumer-genomics exports, checked into the repo so the pipeline is
runnable end-to-end by anyone who clones it — no separate data drop needed.

## `ollie_embark/`

Ollie's (Mini Australian Shepherd) full Embark report — the primary demo subject.

| File | What it is |
|---|---|
| `Embark-full-results-data-*.csv` | Full results export, long-format (`category, name, value`). Not parsed by the pipeline itself — used only as a source for `wildtype/tools/embark_panel_genes.py`'s post-hoc validation gene list (see main README's Validation section). |
| `Ollie_HealthReport.pdf` | Embark's rendered PDF report — same underlying calls, for reference/citation-checking. |
| `raw/*.tped`, `raw/*.tfam` | Raw PLINK-format genotype export, 229,988 SNP markers, 99.8% call rate. This is the real pipeline input — parsed by `wildtype.parsers.tped`, see main README's "Run it". |

Ollie's four flagged/notable calls (referenced only for post-hoc validation — see main README):

| Gene | Call |
|---|---|
| FGF4 retrogene (CFA12) | At risk — structural (drives the IVDD/chondrodystrophy demo) |
| GPT | At risk — missense |
| PRCD | Carrier — missense, tests the zygosity gate |
| MDR1/ABCB1 | Clear — true-negative check |

## `basepaws_samples/`

A Basepaws (cat) sample report, kept as a second data source to validate the
species-agnostic / breed-name-only paths against something that isn't Embark
or isn't a dog. Not yet wired into any parser — `offer_basepaw.pdf` looks
like a promotional offer rather than a results export; confirm which of the
two is the actual structured report before building a Basepaws parser
against it.

## Why this is checked in

This is pet genomic data, not human data, and the owner (Ollie's) has opted
to share it publicly as the hackathon's demo dataset. If that ever changes,
pull these files back out and re-point `.env`'s `WILDTYPE_EMBARK_CSV` /
`WILDTYPE_TPED` / `WILDTYPE_TFAM` at a local, untracked copy instead.
