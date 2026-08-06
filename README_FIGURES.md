# README_FIGURES — how to regenerate each figure

**Key idea: a notebook is a pipeline, not a figure.** Several paper/SI figures are not
separate analyses — they are the *same* notebook re-run with a different parameter
(usually `network_name` ∈ {email, conf, school}, or `q` ∈ {0.02, 0.04}). E.g. Fig. 2
(email) and SI Figs. S10/S11 (conference/school) are one notebook,
`fig2_eph_distributions_simple_vs_complex.ipynb`, run three times with the dataset
switched. This file is organized **by pipeline**, and each pipeline lists every
figure/panel it produces across the main text and SI, with the parameter that
distinguishes them.

Each pipeline has the same 3 stages, in order:

```
contagion_simulation/*.py   →  result_*/unweighted_{simple,complex}_<net>_....csv
        ↓
eph_computation/*.py         →  result_*/unweighted_..._EPH*.csv   (adds EPH + corr columns)
        ↓
notebooks/*.ipynb            →  reads the _EPH*.csv files, makes the plot
```

All three stages read/write relative to `../networks/` and `../result_*/` folders next
to each script — see `networks/README.md` for pointing `NETWORK_DIR` at real data, and
create the `result_*/` output folder referenced by each script's `RESULT_DIR` /
`results_dir` before running (the scripts do not auto-create it).

This is a **reconstructed pipeline order**, inferred from what each script reads/writes
(filenames, `RESULT_DIR` defaults, parameter grids) and from the paper text — it has not
been re-run end-to-end as part of this cleanup pass. Where the mapping is inferred rather
than confirmed by an explicit code comment, it says "likely".

---

## Fig. 1 (EPH concept figure — cycles A/B/C/D)

Hand-drawn/illustrative figure (toy 6/7-node cycle examples), not generated from the
empirical-network simulation pipeline. Affinity Designer source files exist under
`../affinity_outputs/` (`fig1-v1.afdesign`, `panelA_fig1.pdf`, `panelB_fig1.pdf`) — no
simulation/EPH/notebook stage applies here.

---

## Pipeline A — EPH distribution violin plots (simple vs. complex contagion)

**Produces:** Fig. 2 (email, q∈{0.02}), SI Fig. S10 (conference, q=0.02), SI Fig. S11
(school, q=0.02), SI Fig. S12 (email, q=0.04), SI Fig. S13 (conference, q=0.04), SI
Fig. S14 (school, q=0.04) — **all one notebook**, parameterized by `network_name` and the
`complex_probs` value used in the contagion stage.

1. **Contagion** — `contagion_simulation/contagion_simulation_simple_and_complex.py`.
   Set `network_list = ['email_eu_modified']` (Fig. 2 / S12) or `['conf']` (S10/S13) or
   `['school']` (S11/S14). `simple_betas = [0.02,0.03,0.04,0.05]`,
   `complex_thetas = [2,3,4,5]`, `complex_probs = [0.02]` for the main q=0.02 panels or
   `[0.04]` for the S12-S14 robustness-check panels, `stop_ratio = 0.85`, 200 sims each
   (per caption).
2. **EPH** — `eph_computation/eph_computation_v3.py` (plain, no truncation) on the
   resulting CSVs.
3. **Notebook** — `notebooks/fig2_eph_distributions_simple_vs_complex.ipynb`, re-run per
   dataset/q by changing which CSVs it loads. (`notebooks/SI_conference_school_visualize_owen.ipynb`
   is an alternate/earlier notebook that produces the same conf/school panels — see
   "Notes on redundant notebooks" below.)

**Fig. 2's right panel (truncated to 3 infection steps) is a variant within this same
figure**, not a separate pipeline:
1. **Contagion** — `contagion_simulation/contagion_simulation_stepwise_termination.py`,
   `stop_ratio = 3` (this script overloads `stop_ratio` to mean "number of infection
   loops", not fraction infected).
2. **EPH** — `eph_computation/eph_computation_v3.py` (same; no truncation flag needed
   since the simulation itself already stopped at step 3).
3. **Notebook** — same `fig2_eph_distributions_simple_vs_complex.ipynb`, which plots the
   full-dynamics and 3-step CSVs side by side. Only done for the email dataset (no SI
   equivalent found for conf/school truncated-to-3-steps).

**Fig. 2 caption's open TODO ("increase step 3 to 4,5,... maybe phase transition") is a
separate, more general mechanism** — not stepwise termination but post-hoc masking of a
fully-run simulation at a sweep of cutoffs:
1. **Contagion** — `contagion_simulation_simple_and_complex.py`, run to full 85%
   completion as usual (do **not** use `stepwise_termination` for this).
2. **EPH** — `eph_computation/eph_computation_multiprocess_with_truncation.py`, sweeps
   `up_to_step_infections = [2,3,4,5,6,7,8,9,10]` automatically, writes one
   `*_EPH-v3_step-<k>.csv` per cutoff.
3. **Notebook** — `notebooks/fig2_phase_transition_truncated_steps.ipynb`. **Status:
   likely still work-in-progress** (see main README) — treat output as unconfirmed.

---

## Pipeline B — classification (θ) + regression (θ, q) on EPH vs. corr baseline

**Produces:** Fig. 3A, 3B, 3C (email) — **and** SI Figs. S20/S21 (classification, conf/
school, q=0.02), SI Figs. S22/S23/S24 (classification, email/conf/school, q=0.04), SI
Figs. S25/S26 (θ-regression, conf/school), SI Figs. S27/S28 (q-regression, conf/school)
— **all the same notebook**, parameterized by `network_name` and `complex_probs`.

1. **Contagion** — `contagion_simulation/contagion_simulation_simple_and_complex.py` for
   the discrete-θ classification panels (Fig. 3A, S20-S24): `network_list` set to the
   target dataset, `complex_thetas = [2..7]`, `complex_probs = [0.02]` or `[0.04]`,
   800 sims/type per Fig. 3A's caption.
   For the finer continuous-β regression panels (Fig. 3C, S27, S28), use
   `contagion_simulation/contagion_simulation_regression_beta.py` instead — draws
   β ~ Uniform(0.01, 0.05) continuously (500 draws) across
   `network_list = ['email_eu_modified','conf','school']` in one run.
2. **EPH** — `eph_computation/eph_computation_v3.py` for the discrete-grid runs, or
   `eph_computation/eph_regression_beta.py` for the continuous-β regression-beta output
   (matches its `beta_regression_result_.../` directory convention).
3. **Notebook** — `notebooks/fig2fig3_classification_regression_main.ipynb`, re-run per
   dataset/q by pointing it at the relevant CSVs. It contains both the classification
   (decision tree on EPH vs. corr baseline, confusion matrix) and both regression
   sections (θ on EPH, q on EPH) in one notebook — that's why one file covers Fig. 3A/B/C
   and all of SI S8's per-dataset repeats.

---

## Pipeline C — limited-information / partial-node-observability

**Produces:** Fig. 4A main panel, Fig. 4B, Fig. 4C — one notebook family, parameterized
by dataset and by whether the target metric is classification accuracy or regression
MAPE.

1. **Contagion** — `contagion_simulation/contagion_simulation_simple_and_complex.py`,
   run to full 85% completion, per dataset (email/conf/school).
2. **EPH** — `eph_computation/eph_limited_information_partial_nodes.py` — subsamples a
   random `observable_ratio` fraction of nodes per simulation row and computes EPH/corr
   on the induced subgraph; sweep `observable_ratio` across the paper's x-axis values.
   **This script had a filtration-sign bug fixed in this cleanup pass — see "Corrections
   applied" in `README.md`. Any pre-existing Fig. 4 CSVs generated before the fix should
   be regenerated.**
3. **Notebook** — `notebooks/fig4_limited_information.ipynb` (Fig. 4A classification
   accuracy vs. % observed) and `notebooks/fig4_regression_theta_limited_information.ipynb`
   (Fig. 4B θ-MAPE and Fig. 4C q-MAPE vs. % observed — appears to cover both regressions
   in one notebook; check its sections). These are two notebooks, not per-figure copies,
   because they plot different metrics (accuracy vs. MAPE), not different datasets.

**Fig. 4A inset (mixed-dataset generalization)** is a distinct sub-analysis within the
same pipeline stage 1-2, different notebook for stage 3:
1. **Contagion** — run `contagion_simulation_simple_and_complex.py` separately on all
   three networks with the same θ/q grids.
2. **EPH** — `eph_computation/eph_computation_v3.py` per network (full-observability, not
   the partial-node version — this inset tests cross-network generalization, not partial
   observability).
3. **Notebook** — `notebooks/fig4_combining_features_mixed_dataset.ipynb` — merges the
   per-network `_EPH-v3.csv` files itself, trains on merged data, tests on a held-out
   network.

---

## SI §S2 (cycle length definitions) / SI §S3 (PH/EPH background)

Conceptual/mathematical sections, illustrated with the same toy-cycle diagrams as Fig. 1
(Affinity files) plus the algorithm box (SI Algorithm 1) that
`compute_persistent_homology()` in every `eph_computation/*.py` script implements. No
separate simulation/notebook stage.

## SI §S6 — "Topology Matters Beyond Time" (EPH vs. TDS baseline)

**Correction:** an earlier pass of this file wrongly called this a gap — TDS ("time-
discounted spread") code exists, just computed **inline inside a notebook cell** rather
than as a standalone script, which is why the first name-search (looking only for `.py`
files) missed it.

**Produces:** SI Fig. S15.

1. **Contagion** — `contagion_simulation/contagion_simulation_simple_and_complex.py`,
   run to full 85% completion, per dataset — same as Pipeline C stage 1.
2. **EPH + TDS** — `eph_computation/eph_limited_information_partial_nodes.py` (Pipeline
   C's sign-fixed script) computes EPH/corr on the partial-node-observed subgraphs, but
   the **TDS feature itself is computed inside the notebook**, not this script — see
   below.
3. **Notebook** — `notebooks/fig4_limited_information.ipynb`. Contains
   `add_tds_to_sims(sims, lambda_x=0.1)`: for each simulation row, sums
   `(count of nodes infected at step i) * (lambda_x ** i)` over observed infection steps
   to get a scalar TDS feature per row, trains a second decision-tree classifier on it
   (columns `accuracy_tds` / `ci_lower_tds` / `ci_upper_tds`), and plots it against the
   EPH-based classifier in the same figure (legend: "EPH" vs. "TDS"). **This is the same
   notebook used for Pipeline C / Fig. 4A** — SI Fig. S15 is a byproduct of that
   notebook, not a separate pipeline.

## SI §S7 — Small-world (Watts–Strogatz) benchmark

**Produces:** SI Figs. S16, S17 (absolute-threshold small-world benchmark).

1. **Graph generation** — `notebooks/SI_small_world_generation.ipynb` or
   `SI_small_world_generation_owen.ipynb` — generates the WS graph corpus at varying
   rewiring probability `p` (the "owen" one is the source of the large
   `new_samples_small_world/` `.graphml` corpus seen elsewhere in the repo, not
   duplicated into `clean_code/`).
2. **Contagion** — `contagion_simulation/contagion_simulation_basic.py` (its
   `RESULT_DIR` defaults to `result_small_world/`, and its `network_list` is
   small-world-graph-named by convention, e.g. `'smallworld_p0.1_seed0'`).
3. **EPH** — `eph_computation/eph_computation_v3.py` per small-world graph instance.
4. **Notebook** — `notebooks/SI_small_world_results.ipynb` — accuracy vs. rewiring
   probability plus average-diameter overlay.

**SI Figs. S17(second, "relative")/S18 (relative-threshold regime, EPH-distribution
side, marked "TBD" in the paper) are a separate, unfinished branch** — do not reuse
Pipeline A/B/C: `relative_threshold_exploratory/contagion_simulation_relative_threshold.py`
→ `relative_threshold_exploratory/eph_relative_threshold.py` →
`relative_threshold_exploratory/notebooks/exploring_deg_of_graphs.ipynb` /
`phase_transition_Feb14.ipynb` / `small_world_phase_transition.ipynb`. **Appears
unfinished** (see `README.md` §"Remaining / Started-but-Unfinished Analyses") — no
notebook there produces a finished θ-sweep plot matching Pipeline A/B/C's polish. Treat
as a starting point to finish, not a ready-to-run recipe.

**SI Fig. S19 (TDS vs. EPH accuracy in the relative-threshold small-world regime) is
finished, unlike S17/S18 above** — correction from an earlier pass of this file:
1. **Graph generation + contagion + EPH** — same relative-threshold branch as above
   (`contagion_simulation_relative_threshold.py` → `eph_relative_threshold.py`), run
   across the small-world graph corpus and a sweep of relative thresholds φ.
2. **Notebook** — `relative_threshold_exploratory/notebooks/small_world_results_TDS_vs_EPH.ipynb`
   (sourced from `relative_threshold/contagion_simulation/notebooks/small_world_results.ipynb`,
   12MB — this file was missed in the original `clean_code/` copy and has been added).
   Contains `compute_time_discounted_spread(df, lam)` (median-based TDS across rows),
   trains EPH-based and TDS-based classifiers per θ, and saves
   `small_world_result_TDS_EPH.pdf` — this is very likely SI Fig. S19 itself.

## SI §S9.1 — Numerical implementation details (Algorithm 1, coning trick)

Documentation of `compute_persistent_homology()` in every `eph_computation/*.py` script
— no separate figure.

## SI §S9.2 — Variance reduction (fixed RNG seeding of initial infected nodes)

Documentation of `init_nodes_var_reduction()` (present in every
`contagion_simulation/*.py` script) — no separate figure.

---

## Quick reference table

| Pipeline | Contagion script | EPH script | Notebook | Figures it produces (vary by param shown) |
|---|---|---|---|---|
| A | `contagion_simulation_simple_and_complex.py` (or `_stepwise_termination.py` for the 3-step variant) | `eph_computation_v3.py` | `fig2_eph_distributions_simple_vs_complex.ipynb` | Fig. 2 (email, q=0.02, incl. 3-step panel) · SI S10 (conf, q=0.02) · SI S11 (school, q=0.02) · SI S12 (email, q=0.04) · SI S13 (conf, q=0.04) · SI S14 (school, q=0.04) — vary `network_name`, `complex_probs` |
| A (phase-transition follow-up, WIP) | `contagion_simulation_simple_and_complex.py` | `eph_computation_multiprocess_with_truncation.py` | `fig2_phase_transition_truncated_steps.ipynb` | Fig. 2 caption's open TODO (truncation step k=2..10) |
| B | `contagion_simulation_simple_and_complex.py` (discrete θ) or `contagion_simulation_regression_beta.py` (continuous β) | `eph_computation_v3.py` or `eph_regression_beta.py` | `fig2fig3_classification_regression_main.ipynb` | Fig. 3A/3B/3C (email) · SI S20/S21 (class., conf/school, q=0.02) · SI S22-S24 (class., email/conf/school, q=0.04) · SI S25/S26 (θ-regression, conf/school) · SI S27/S28 (q-regression, conf/school) — vary `network_name`, `complex_probs` |
| C (limited info) | `contagion_simulation_simple_and_complex.py` | `eph_limited_information_partial_nodes.py` (sign-fixed) | `fig4_limited_information.ipynb` (accuracy) · `fig4_regression_theta_limited_information.ipynb` (MAPE) | Fig. 4A main, 4B, 4C — vary `network_name` and which notebook (accuracy vs. MAPE) |
| C (mixed dataset) | `contagion_simulation_simple_and_complex.py` ×3 networks | `eph_computation_v3.py` | `fig4_combining_features_mixed_dataset.ipynb` | Fig. 4A inset |
| S6 (EPH vs TDS, limited-info) | `contagion_simulation_simple_and_complex.py` | `eph_limited_information_partial_nodes.py` (sign-fixed) | `fig4_limited_information.ipynb` (TDS computed inline via `add_tds_to_sims`) | SI S15 |
| S7 (small-world, absolute) | `contagion_simulation_basic.py` | `eph_computation_v3.py` | `SI_small_world_generation*.ipynb` → `SI_small_world_results.ipynb` | SI S16, S17 |
| S7 (small-world, relative, EPH-dist — TBD) | `relative_threshold_exploratory/contagion_simulation_relative_threshold.py` | `relative_threshold_exploratory/eph_relative_threshold.py` | `exploring_deg_of_graphs.ipynb` / `phase_transition_Feb14.ipynb` / `small_world_phase_transition.ipynb` | SI S17(second)/S18 — **unfinished** |
| S7 (small-world, relative, TDS vs EPH — finished) | `relative_threshold_exploratory/contagion_simulation_relative_threshold.py` | `relative_threshold_exploratory/eph_relative_threshold.py` | `small_world_results_TDS_vs_EPH.ipynb` | SI S19 |

## Notes on redundant notebooks

`notebooks/SI_conference_school_visualize_owen.ipynb` appears to be an earlier or
parallel notebook covering some of the same conf/school EPH-distribution ground as
Pipeline A above (it produced the standalone `fig2-panelA_conf_...` /
`fig2-panelA_school_...` PDFs found in `res_SI_owen/`). It was kept alongside
`fig2_eph_distributions_simple_vs_complex.ipynb` rather than merged/deleted, since it's
unclear which one is authoritative for the *published* SI S10/S11 panels without
comparing rendered output — **if you can confirm one is strictly superseded, remove the
other and update this file.**

**Reminder:** none of this changes any algorithm logic — this file only documents *which
already-existing script/notebook to run, with which parameters, in which order, and which
figures it's reused for*. If you hit a mismatch between what a script/notebook actually
contains and what's described here, trust the code and treat this file as needing a
correction.
