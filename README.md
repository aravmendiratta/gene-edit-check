# CRISPR-Reflexion

**Find out whether a single base edit can switch a gene down in one cell state and leave another alone,
using Google DeepMind's AlphaGenome 2 as the judge.**

[![CI](https://github.com/YOUR_GITHUB_USERNAME/crispr-reflexion/actions/workflows/ci.yml/badge.svg)](https://github.com/YOUR_GITHUB_USERNAME/crispr-reflexion/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/crispr-reflexion)](https://pypi.org/project/crispr-reflexion/)
![Python](https://img.shields.io/badge/python-3.11%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Use](https://img.shields.io/badge/use-research%20only-orange)

Worked example: turn **PD-1 (*PDCD1*)** down in the exhausted T cells inside real colorectal tumours while
leaving normal-colon T cells alone. The search finds no edit that passes every safety check, and the
project reports exactly why.

![Summary figure](https://raw.githubusercontent.com/YOUR_GITHUB_USERNAME/crispr-reflexion/main/outputs/figure_crc-pd1-exhaustion.png)

> *Research use only. Everything here is a computational prediction from the AlphaGenome 2 API. It is not
> validated in a lab and is not for clinical, diagnostic or therapeutic use.*

---

## Install

```bash
pip install crispr-reflexion                 # the command-line tool and Python library
pipx install crispr-reflexion                # same, in its own isolated environment (recommended for the CLI)
pip install "crispr-reflexion[mcp]"          # also adds the server for AI coding agents
```

Needs Python 3.11 or newer. From source: `pip install git+https://github.com/YOUR_GITHUB_USERNAME/crispr-reflexion`.

## What happens when you install it

You get one command, `crispr-reflexion`, and an importable Python package. Installing does not download
data, call any API or create any files. Nothing happens until you run something.

| You have | You can do | Needs |
|---|---|---|
| **Nothing else** | Check your setup, browse AlphaGenome 2's 2,799 reference conditions ranked by how well the model actually knows each one, run the plumbing with the simulator | just the install |
| **A free API key** | Run the full pipeline on the built-in colorectal and breast cancer cohorts and get a report and figure | key, plus a one-time 396 MB atlas download |
| **Your own `.h5ad` file** | Run the same pipeline on any single-cell dataset, for any target gene | key, plus your data |

### 30 seconds, no key

```bash
crispr-reflexion doctor                                   # what is set up, and where files will go
crispr-reflexion conditions "CD8-positive"                # which reference conditions are best evidenced
```

### A real run

```bash
export ALPHAGENOME_API_KEY=...               # free, non-commercial: https://deepmind.google.com/science/alphagenome
                                             # PowerShell: $env:ALPHAGENOME_API_KEY='...'   or put it in a .env file
crispr-reflexion optimize                    # PD-1 in exhausted colorectal-cancer T cells
```

The first run downloads the Pelka et al. 2021 atlas (396 MB, once) and takes a few minutes. It makes tens of
AlphaGenome calls and caches every response, so repeating a run is fast and free. You get:

- a **terminal report**: the candidate edits, each check it passed or failed, and why
- `outputs/report_crc-pd1-exhaustion.json`: the full machine-readable record
- `outputs/figure_crc-pd1-exhaustion.png`: a four-panel summary

Without a key, a real-data run stops and tells you how to get one. It will not quietly substitute fake output.

### Your own data

```bash
crispr-reflexion run-h5ad \
  --h5ad my_atlas.h5ad --gene-id ENSG00000188389 --gene-symbol PDCD1 \
  --target-query "cell_type == 'CD8 T cell' and condition == 'tumour'" \
  --healthy-query "cell_type == 'CD8 T cell' and condition == 'normal'" \
  --cell-type "CD8-positive" --tissue "colon"
```

The file needs raw integer UMI counts (in `.raw.X`, `.layers["counts"]` or `.X`) and a column identifying the
patient or donor. Cell groups are selected with pandas queries, including `expr_GENE > 0` conditions. The
Tier 2 reference condition is chosen automatically from your `--cell-type`.

## Who this is for

- **Someone designing a state-selective regulatory edit** who wants a cheap first pass before committing to
  wet-lab work: does the model think this edit moves the gene in the cells you care about and not in the ones
  you don't, and what does it break on the way?
- **Anyone using AlphaGenome 2** who wants a better way to pick a reference condition. `crispr-reflexion
  conditions` ranks them by whether the model saw real training data or only imputed it. Only 185 of 2,799
  human conditions have observed methylation.
- **Single-cell groups** who need to condition AlphaGenome on a population but are limited to 256 cells. The
  coreset sketching and patient-disjoint hold-out are reusable on their own.
- **Researchers evaluating verifier-guided search**, with a worked and honest negative result.

It is not a guide-RNA design tool, it does not predict off-target cutting, and it does not replace experiments.

## Use it from Python

```python
from crispr_reflexion.conditions import search_conditions, evidence_profile

hits = search_conditions("CD8-positive", modality="methylation")   # strongest evidence first
print(hits[["condition", "biosample_name", "methylation"]].head())
print(evidence_profile("UBERON:0000178_CL:0000625"))               # blood CD8 T cell: all four 'observed'
```

The full pipeline is `crispr_reflexion.pipeline.run_scenario`. The scripted simulator in
`crispr_reflexion.ag_backend.simulator` lets you test code offline; its output is not a prediction.

## Use it from an AI agent

- **Skill.** `skills/alphagenome_v2_single_cell_crispr/SKILL.md` teaches an agent how to drive the CLI and how to
  report results honestly.
- **MCP server.** `crispr-reflexion-mcp` exposes `run_optimization_tool`, `list_scenarios_tool` and
  `search_conditions_tool` to any MCP-capable agent.

## What it found

Of 8 candidate edits at *PDCD1*, 7 were blocked by the CpG-methylation check, and the closest near miss fails
it at 0.146 against a limit of 0.12. Across 48 scored variants, efficacy and methylation shift are positively
associated (Spearman ρ = +0.57, reproduced across three seeds, six cohorts and six checkpoint genes). The
association is moderate, and the controls do not show it is specific to methylation. An earlier, stronger claim
(ρ = +0.94) was retracted after a larger random null did not support it.

Details, tables and the parts that did not hold up: [docs/RESULTS.md](docs/RESULTS.md).
How each stage works and why: [docs/HOW_IT_WORKS.md](docs/HOW_IT_WORKS.md).

## Where files go

| What | Where |
|---|---|
| Reports and figures | `./outputs/` in the folder you run from |
| Downloaded atlases and cached AlphaGenome responses | a per-user cache folder (`%LOCALAPPDATA%\crispr-reflexion`, `~/Library/Caches/crispr-reflexion`, or `~/.cache/crispr-reflexion`) |
| Override the cache location | set `CRISPR_REFLEXION_HOME` |

`crispr-reflexion doctor` prints the exact paths on your machine.

## Limits and terms

- AlphaGenome 2 is used under its [non-commercial terms](https://deepmind.google.com/science/alphagenome/terms),
  and you need your own API key. The code here is MIT licensed; the model and its predictions are not.
- No wet-lab validation. "Selective" means AlphaGenome 2 says so.
- The healthy comparators are imperfect (memory-like colon T cells, not naive ones), the Tier 2 hazard
  reference is blood CD8 T cells because the colon one is only imputed, and the breast cohort has 5 patients.

## Develop

```bash
git clone https://github.com/YOUR_GITHUB_USERNAME/crispr-reflexion && cd crispr-reflexion
uv venv --python 3.12 .venv && uv pip install --python .venv/Scripts/python.exe -e ".[dev,mcp]"
pytest                    # offline: no API key and no download needed
```

## Cite

Pelka et al., *Cell* 2021 (data); Wu et al., *EMBO J* 2020 (breast cohort); Avsec et al. and Google DeepMind
(AlphaGenome). See [CITATION.cff](CITATION.cff) and [references/citation.bib](references/citation.bib).
