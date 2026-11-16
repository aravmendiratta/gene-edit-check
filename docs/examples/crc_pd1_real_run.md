# Worked example: PDCD1 silencing in exhausted colorectal-tumor CD8 T cells (real data, live API)

`crispr-reflexion optimize --sketch-size 128 --pool-size 600 --max-iterations 8 --null-n 12 --seed 42`
Backend `LIVE_ALPHAGENOME_V2` (`alphagenome-2.0.0-all_folds`). Data: Pelka et al. 2021 (see README). Research use only.

**Outcome: not converged.** No candidate passed every pre-set check in 8 iterations.

## Pre-set rules (fixed before the run)
Direction: PD-1 must go **down**. Efficacy: held-out |LFC| must exceed the largest of 12 random PAM-valid edits (**0.0604**, p ≥ 0.077).
Selectivity: normal-colon shift ≤ 0.35 × target shift. Splice ≤ 0.15. CpG methylation ≤ 0.12. Zero bystander edits in window [4..8].

## Trace (LFC on held-out-patient Sketch B; healthy = normal-colon CD8)

| # | Edit | Target | Healthy | Ratio | Splice | Meth | Failed |
|---|---|---:|---:|---:|---:|---:|---|
| 1 | ABE8e `…858903 A>G` | −0.130 | −0.231 | 1.78 | 0.00 | 0.684 | selectivity, methylation, bystander (2) |
| 2 | CBE `…859997 C>T` | −0.132 | −0.031 | 0.23 | 0.00 | 0.691 | **methylation only** |
| 3 | ABE8e `…859845 A>G` | −0.077 | −0.012 | 0.15 | 0.00 | 0.146 | **methylation only (by 0.026)** |
| 4 | CBE `…858884 C>T` | −0.354 | −0.364 | 1.03 | 0.66 | 0.861 | selectivity, splice, methylation (promoter hit, non-selective) |
| 5 | CBE `…858777 G>A` | +0.041 | +0.048 | 1.20 | 0.16 | 0.059 | direction, selectivity, splice |
| 6 | ABE8e `…858890 T>C` | +0.068 | +0.071 | 1.05 | 0.15 | 0.309 | direction, selectivity, methylation |
| 7 | CBE `…858785 G>A` | +0.028 | +0.040 | 1.43 | 0.09 | 0.062 | direction, selectivity |
| 8 | ABE8e `…858912 A>G` | +0.063 | +0.073 | 1.15 | 0.00 | 0.372 | direction, selectivity, methylation |

## What this says
- **Two candidates are selective and effective** (2 and 3, ratio 0.23 and 0.15, both beat the random-edit null) and **both fail only on CpG methylation**. Candidate 2 (`protospacer CATCAAAGATGATGTCTGCG`, PAM AGG) overlaps what the project's motif table labels the PD-1 CR-C enhancer (annotation not independently verified).
- Methylation is a genuine, not a spurious, signal: random edits at this locus have median methylation shift ≈ 0.012 (11 of 12 ≤ 0.05; one outlier 0.16 that also had the largest random PD-1 effect), so 0.146 is ~12× typical and 0.69 is far beyond it. A C>T at a CpG changes the CpG itself.
- The 0.12 limit was **not loosened** to make candidate 3 pass. Whether 0.146 is acceptable is a scientific judgement for a domain expert, not something to tune away.
- Iterations 5–8 are all wrong-direction edits that the verifier correctly rejected. That is a weakness of the proposer, not a success: after the first four candidates the gradient-sign ranking prior has little left to offer, and the loop spent its remaining budget on edits that raise PD-1. A better proposer (e.g. signed ISM over more windows, or a wider PAM set such as NG) is the obvious next step.

## Caveats
Single seed, single gene, 12-edit null with one heavy-tailed outlier; Tier 2 uses a blood CD8 reference; all numbers are model predictions.
