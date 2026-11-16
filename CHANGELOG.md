# Changelog

## 0.3.0

First public release.

**Pipeline**
- Real-data pipeline on the Pelka et al. 2021 colorectal atlas, with sketches split by patient so the
  verification sketch never shares a patient with the proposal sketch.
- Edit direction is measured by an in-silico mutagenesis pre-screen instead of inferred from gradient sign,
  which agreed with the measured sign only about 70% of the time.
- Efficacy must beat a 40-edit random null; neighbouring genes in the scored window are checked; finalists are
  scored against six non-immune tissues.
- Tier 2 hazard references are chosen by the model's own evidence tier (observed, imputed-supported,
  imputed-exploratory), with a tissue-matched cross-check reported alongside rather than averaged in.
- `run-h5ad` runs the pipeline on any single-cell `.h5ad` file.

**Commands**: `optimize`, `run-h5ad`, `conditions`, `doctor`, `validate-fidelity`, `ablate-sketch`,
`sweep-coupling`, `sweep-universality`, and `crispr-reflexion-mcp` for AI agents.

**Install behaviour**
- Downloaded atlases and cached API responses go to a per-user cache folder, never the Python install
  directory. Override with `CRISPR_REFLEXION_HOME`.
- A real-data run without an API key stops with instructions instead of falling back to the simulator.
- `.env` is found in the folder you run from.
- `crispr-reflexion doctor` reports your setup.

**Corrections made before release**
- An earlier draft reported a CpG-methylation vs efficacy correlation of +0.94 at *PDCD1*. With a 40-edit random
  null it is about +0.5, the random-edit control is not significant, and splicing also correlates. The claim was
  retracted; see `docs/RESULTS.md`.
- The shipped reference-condition table had been excluded from the package by an over-broad ignore rule.
