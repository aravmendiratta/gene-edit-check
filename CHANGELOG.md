# Changelog

## Unreleased

- `scripts/make_launch_video.py` builds a 25-second silent launch video (`assets/launch.mp4`) from the project's own
  figures. It needs `pillow` and `imageio-ffmpeg`, which are not project dependencies.

- Design language, `docs/DESIGN.md`: a drawn double-helix mark with one coral rung (the edited base), a deep
  teal, cream and coral palette, serif titles, and a small helix stamp on every figure. Figures ship in light
  and dark (`fig_x.png`, `fig_x-dark.png`), with a logo, banners and a social card (`scripts/make_brand.py`).
  The terminal output uses the same two hues.

## 0.3.1

- New: check against a published wet-lab screen (`scripts/validate_enhancers.py`). Model saliency recognises two
  elements near the PDCD1 promoter but not the +14.8 and +16.1 kb enhancers; single-base edits there are
  predicted to do almost nothing. Details and caveats in `docs/RESULTS.md`.
- New: `--scan-region START:END` (repeatable) points the search at a distal element instead of the promoter, and
  `--tag` keeps the outputs of separate runs apart.
- Renamed the search loop away from LLM-agent terms: `agent/reflexion_loop.py` is now `search/search_loop.py`,
  and the report fields `critic_diagnosis` and `episodic_memory_rule_added` are now `diagnosis` and
  `constraint_added`. Reports saved by 0.3.0 need those two keys renamed to load.
- Fixed: the CLI crashed on Windows consoles and redirected output that use a legacy code page, because the
  tables contain non-ASCII characters. Output is now UTF-8 with replacement. Found by the new Windows CI job.
- CI now runs on Linux (Python 3.11 and 3.12), macOS and Windows.

- Renamed the project to `gene-edit-check`.
- Efficacy-hazard correlations now carry 95% bootstrap intervals, and the methylation-versus-splice label comes
  from a paired bootstrap of the gap instead of a fixed cutoff. For the main PDCD1 run: rho +0.57
  (CI +0.28 to +0.77); the gap to splicing is +0.41 (CI -0.02 to +0.82), so specificity to methylation is not
  established.
- New optional dose ceiling `optimize --max-effect`, which rejects edits that silence the target more than a
  given |LFC|. Off by default.
- `optimize` and `run-h5ad` now default to a 40-edit random null, matching the sweeps.
- New `docs/RELATED_WORK.md` listing prior art and what this project does and does not add.

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
`sweep-coupling`, `sweep-universality`, and `gene-edit-check-mcp` for AI agents.

**Install behaviour**
- Downloaded atlases and cached API responses go to a per-user cache folder, never the Python install
  directory. Override with `GENE_EDIT_CHECK_HOME`.
- A real-data run without an API key stops with instructions instead of falling back to the simulator.
- `.env` is found in the folder you run from.
- `gene-edit-check doctor` reports your setup.

**Corrections made before release**
- An earlier draft reported a CpG-methylation vs efficacy correlation of +0.94 at *PDCD1*. With a 40-edit random
  null it is about +0.5, the random-edit control is not significant, and splicing also correlates. The claim was
  retracted; see `docs/RESULTS.md`.
- The shipped reference-condition table had been excluded from the package by an over-broad ignore rule.
