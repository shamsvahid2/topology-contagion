# topology-contagion — code for "Extended Persistent Homology Distinguishes Simple and Complex Contagions with High Accuracy"

This repo holds the code behind "Extended Persistent Homology Distinguishes Simple and
Complex Contagions with High Accuracy" (Shamsaddini, Rahimian, Spencer). It was
consolidated from a much larger private working directory containing many overlapping,
iterative version folders (`code-april17-version/`, `code-final-version-jul2025/`,
`code_version_26feb/`, `relative_threshold/`, `res_SI_owen/`, `res_small_world_owen_v1/`,
etc.) — this repo keeps only the **last/most-complete version** of each script or
notebook. **No algorithm/model logic was changed** during that consolidation — files
were only copied, renamed, and (where present) had hardcoded absolute paths converted to
relative ones. No absolute `/Users/.../` or `/home/.../` paths were found in any of the
files selected here — the scripts already used paths relative to a location sitting next
to `networks/` and `result_*/` folders.

Real empirical network data files are included under `networks/` (see
`networks/README.md`). Simulation-output and EPH-feature CSVs (the `result_*/` folders)
are **not** checked in — they're regenerated locally by running the scripts; see
`README_FIGURES.md` for exactly which script to run for which figure.

One deliberate logic **correction** was made after consolidation — see "Corrections
applied" below.

### Corrections applied

- **`eph_computation/eph_limited_information_partial_nodes.py`** (drives the Fig 4
  limited-information / partial-node-observability results): the filtration passed to
  `compute_persistent_homology` was `filt` (raw infection step, un-negated), while the
  other three EPH scripts (`eph_computation_v3.py`, `eph_regression_beta.py`,
  `eph_computation_multiprocess_with_truncation.py`) all pass `-filt`. Per the paper (EPH
  filtration = reverse of infection step, i.e. earliest-infected node gets the most
  negative value) and confirmation from Vahid, `-filt` is the correct convention. Fixed
  in this repo to `filtration=-filt` on 2026-08-06. **This is a real bug fix, not just a
  cleanup rename** — the original file (in whichever version folder it came from) still
  has the wrong sign, and any Fig 4 results generated with the unfixed script should be
  treated as suspect and likely need to be regenerated with the corrected sign.

**To regenerate a specific figure, see [`README_FIGURES.md`](./README_FIGURES.md)** —
a step-by-step recipe (contagion script → EPH script → notebook, in order) for every
figure in the paper and SI.

## Folder layout

```
clean_code/
  contagion_simulation/    core SI/simple/complex(threshold+noisy) contagion simulators
  eph_computation/         EPH (extended persistent homology) + corr(order,deg) baseline computation, GUDHI-based
  baseline/                (see note below — no separate baseline script found; corr baseline lives inside eph_computation/)
  relative_threshold_exploratory/   the separate "relative threshold" branch (SI S17-19, unfinished)
  networks/                README only — points at where the real network files live (not duplicated)
  notebooks/                final analysis/plotting notebooks, renamed by the figure they (likely) produce
  README.md                 this file
```

## Traceability table: file → source → paper figure/section

| clean_code file | sourced from | role | paper mapping |
|---|---|---|---|
| `contagion_simulation/contagion_simulation_simple_and_complex.py` | `code_version_26feb/contagion_simulation/python_scripts/simulations_types_added_multiprocess.py` (newest, Feb 2026, multiprocessing) | Runs SI-style simple contagion (beta) and complex/threshold contagion (theta + prob_complex "noisy threshold", records above-/below-threshold adoption type per node) on a chosen empirical network, dumps infection-step CSVs | Underlies Figs. 2, 3 main-text simulations; the "adoption type" tracking supports the complex-contagion noise (prob_complex) sweeps referenced in Methods/S8.2 |
| `contagion_simulation/contagion_simulation_simple_and_complex_singleprocess.py` | same folder, `simulations_types_added.py` | Same logic, non-parallel version | superseded by the multiprocess version above; kept for reference |
| `contagion_simulation/contagion_simulation_basic.py` | `code_version_26feb/.../simulations.py` | Earlier/simpler simulator (no adoption-type tracking) — closest to the version used for small-world benchmark runs (`RESULT_DIR` defaults to `result_small_world`) | SI small-world benchmark (SI Figs ~16/17) |
| `contagion_simulation/contagion_simulation_regression_beta.py` | `code-april17-version/topology-contagion/python_scripts/simulations-regression-beta.py` | Weighted-network simulation sweeping beta for the weighted/regression pipeline | Fig. 3B/3C-style regression on weighted networks (beta), SI weighted-network analysis |
| `contagion_simulation/contagion_simulation_stepwise_termination.py` | `code-april17-version/.../simulations-stepwise-termination.py` | Terminates simulation after a fixed number of infection "loops"/steps rather than a stop ratio | Related to the limited-information / truncated-dynamics analyses (Fig. 4, and the Fig.2 "truncation step" follow-up) |
| `eph_computation/eph_computation_multiprocess_with_truncation.py` | `code_version_26feb/.../EPH_new_version_multiprocess.py` (newest, Feb 2026) | GUDHI coning-trick extended persistent homology per SI Algorithm 1 (dimension-1 generators, mean lifetime = EPH feature) **plus** the corr(order,degree) Spearman baseline (Cencetti et al. PRL method), **plus** an `up_to_step_infection` cutoff (2..10) that masks late infections as "never infected" before computing EPH | This is the code implementing the Fig. 2 caption's "Vahid's to do: what if we increase step 3 to 4,5,... maybe phase transition" — see `notebooks/fig2_phase_transition_truncated_steps.ipynb` |
| `eph_computation/eph_computation_v3.py` | `code_version_26feb/.../EPH-v3.py` | Same EPH+corr computation without the truncation-step sweep | Core EPH feature computation for Figs. 2-3 |
| `eph_computation/eph_limited_information_partial_nodes.py` | `code-april17-version/.../EPH-v2-limited-information-partial-node-lists.py` | EPH computed on a partial/observed subset of nodes (limited-information regime) | Fig. 4 (limited-information classification) |
| `eph_computation/eph_regression_beta.py` | `code-april17-version/.../EPH-v3-regression-beta.py` | EPH+corr computation for the beta-regression (weighted) sweep | Fig. 3B/3C-style regression, weighted networks |
| `relative_threshold_exploratory/contagion_simulation_relative_threshold.py` | `relative_threshold/contagion_simulation/python_scripts/simulations_relative_threshold.py` | Same simple/complex contagion engine, extended so `threshold_complex` can be a **relative** threshold φ∈(0,1] (adopt once `infected_neighbors ≥ ceil(φ·degree)`), not just an absolute integer count | SI §S7 "relative complex contagion regime" (caption for Supplementary Fig. 17/18/19, marked TBD in the paper) |
| `relative_threshold_exploratory/contagion_simulation_relative_causality.py` | `.../simulation_relative_casality.py` [sic, original filename typo preserved] | Similar relative-threshold contagion engine, apparently a parallel/alternate exploration ("causality" naming suggests investigating temporal mixing of simple vs. complex adoption events) | Same SI §S7 area; not clearly finished |
| `relative_threshold_exploratory/eph_relative_threshold.py` | `.../EPH-v3_relative.py` | EPH+corr computation for the relative-threshold outputs | SI §S7 |
| `relative_threshold_exploratory/notebooks/small_world_results_TDS_vs_EPH.ipynb` | `relative_threshold/contagion_simulation/notebooks/small_world_results.ipynb` (12MB, Dec 2025 — **added in this pass, was missed originally**) | Contains `compute_time_discounted_spread(df, lam)` (median-based TDS) and trains/plots EPH-vs-TDS classification accuracy across θ in the relative-threshold small-world regime, saving `small_world_result_TDS_EPH.pdf` | **Likely the actual source of SI Fig. S19** ("check if power of TDS is higher than EPH") — this part of SI §S7 is more finished than the original pass of this README reported |
| `notebooks/fig2fig3_classification_regression_main.ipynb` | `code-april17-version/topology-contagion/notebooks_new/main_notebook_reg_class.ipynb` (Feb 2026, largest/most recent "main" notebook, 8.6MB) | Contains: regression analysis (simple & complex), classification given-threshold and given-q, "Fig2-PanelA_inset confusion matrix", and final classification results | Figs. 2 (panel A inset), 3A/B/C |
| `notebooks/fig2_eph_distributions_simple_vs_complex.ipynb` | `.../notebooks_new/dist_EPH_simple_complex_together.ipynb` | EPH distribution plots for simple vs. complex contagion side-by-side | Fig. 2 (EPH distribution panels) |
| `notebooks/fig2_phase_transition_truncated_steps.ipynb` | `code_version_26feb/contagion_simulation/notebooks/analyze_results.ipynb` (newest overall notebook file, Feb 20 2026) | Explicitly titled around "more points (higher k values)" — analyzes EPH vs. truncation-step-count k, i.e., the phase-transition follow-up to Fig. 2's "Vahid's to do" caption | Fig. 2 phase-transition follow-up — **appears to be the active/most current work on this open item**, but I could not confirm from the notebook alone whether a clean phase-transition result was reached; treat as in-progress |
| `notebooks/fig4_limited_information.ipynb` | `.../notebooks_new/limited_information.ipynb` | Classification accuracy vs. % of nodes observed / limited information | Fig. 4A |
| `notebooks/fig4_regression_theta_limited_information.ipynb` | `.../notebooks_new/regression_theta_limited_information.ipynb` | Regression of theta under limited information | Fig. 4B/related SI |
| `notebooks/fig4_combining_features_mixed_dataset.ipynb` | `.../notebooks_new/combining_features.ipynb` | Mixed-dataset generalization / combined-feature classification | Fig. 4A inset (mixed-dataset generalization) |
| `notebooks/SI_regression_beta_weighted.ipynb` | `.../notebooks_new/regression_beta.ipynb` | Beta regression, weighted networks | SI weighted-network regression figures |
| `notebooks/SI_3d_regression_surface.ipynb` | `.../notebooks_new/3d_plot.ipynb` | 3D surface plot of EPH/regression across two parameters | SI 3D regression figure |
| `notebooks/SI_small_world_generation.ipynb` / `SI_small_world_generation_owen.ipynb` | `.../notebooks_new/small_world_generation.ipynb` and `res_small_world_owen_v1/generate_small_world_graphs.ipynb` | Generates the Watts-Strogatz small-world graph corpus at varying rewiring probability p | Inputs to SI small-world benchmark (SI Fig. 16/17) |
| `notebooks/SI_small_world_results.ipynb` | `.../notebooks_new/small_world_results.ipynb` (identical copies also exist in `code-final-version-jul2025`, `code_version_26feb`, `relative_threshold` — all same file, kept once here) | Small-world benchmark results/plots | SI Fig. 16/17 (small-world regime EPH vs. baseline) |
| `notebooks/SI_conference_school_visualize_owen.ipynb` | `res_SI_owen/visualize_res.ipynb` | Produces `fig2-panelA_conf_q=0.04.pdf`, `fig2-panelA_school_q=0.02/0.04.pdf` seen alongside it in `res_SI_owen/` | Fig. 2 panel A, SI robustness across conference/school networks |
| `notebooks/eph_pipeline_reference_notebook.ipynb` | `code-april17-version/topology-contagion/notebooks/EPH.ipynb` (older, but referenced explicitly by the repo's own `README.md`: "run the `EPH` notebook to compute extended persistent homology and corr (PRL method)") | Reference/tutorial notebook for the EPH+corr pipeline, matches the `.py` scripts above but interactive | Methods / SI S3 walkthrough |

### `baseline/`
**Correction (2026-08-06): TDS was found — an earlier pass of this README wrongly
reported it as a gap.** The name search only looked for standalone `.py` files, but TDS
("time-discounted spread") is computed **inline inside notebook cells**, not as a
separate module, in two places:

1. `notebooks/fig4_limited_information.ipynb` — function `add_tds_to_sims(sims,
   lambda_x=0.1)`: for each simulation row, sums `(count of nodes infected at step i) *
   (lambda_x ** i)` over observed infection steps. This feeds a second decision-tree
   classifier (`accuracy_tds` / `ci_lower_tds` / `ci_upper_tds` columns) trained
   alongside the EPH-based one, plotted against it in the same figure this notebook
   produces (labelled "EPH" vs. "TDS" in the legend). This is SI §S6's actual
   "Topology Matters Beyond Time" comparison — the code was already present in
   `clean_code/`, just not identified as TDS in the first pass.
2. `relative_threshold/contagion_simulation/notebooks/small_world_results.ipynb` (12MB,
   not originally copied into `clean_code/` — **added now** as
   `relative_threshold_exploratory/notebooks/small_world_results_TDS_vs_EPH.ipynb`) —
   a second, independent implementation, `compute_time_discounted_spread(df, lam)`
   (median-based rather than per-row), producing `small_world_result_TDS_EPH.pdf` (EPH
   vs. TDS classification accuracy across θ in the small-world/relative-threshold
   regime). This is very likely the actual source of SI Fig. S19 (captioned "To check if
   the power of TDS is higher than EPH").

The corr(order, degree) Spearman baseline from Cencetti et al. remains a separate thing,
computed inline inside every `eph_computation/*.py` script (the `corr` column) — so
`baseline/` still isn't a real folder; both baselines (corr and TDS) live inline in
notebooks/scripts rather than as standalone modules.

## What in the paper has NO code found in this repo
- ~~No dedicated TDS baseline~~ — **retracted, see correction above; TDS code exists in
  two notebooks.**
- SI Fig. S17(second)/S18 (q, θ effect on EPH under the **relative** threshold regime) —
  the `relative_threshold/` branch has simulation + EPH code
  (`relative_threshold_exploratory/`) and, per the correction above, a finished
  TDS-vs-EPH accuracy comparison (`small_world_results_TDS_vs_EPH.ipynb`) exists for SI
  Fig. S19. However `phase_transition_Feb14.ipynb` and `small_world_phase_transition.ipynb`
  (which look like they'd cover the q/θ EPH-distribution side, SI Figs. S17/S18) are
  still exploratory (EMD/"temporal mixing index" metrics, no markdown cells describing
  final conclusions) — **those two still appear unfinished/TBD**, even though the TDS
  comparison itself (S19) is more complete than previously reported.

## Code with no clear corresponding paper output (orphaned / exploratory)
- `relative_threshold_exploratory/contagion_simulation_relative_causality.py` — naming
  ("casality"/causality) and content suggest an exploration of temporal separation
  between simple vs. complex adoption events; not obviously tied to any specific
  figure — likely exploratory groundwork for SI §S7.
- `relative_threshold_exploratory/notebooks/exploring_deg_of_graphs.ipynb`,
  `phase_transition_Feb14.ipynb`, `small_world_phase_transition.ipynb` — exploratory,
  no clear finished figure; `small_world_phase_transition.ipynb` ends on "## Cover
  letter" as its last markdown cell, suggesting it was left mid-thought.
- `notebooks/fig2_phase_transition_truncated_steps.ipynb` (`analyze_results.ipynb`) —
  is the most recently touched notebook in the whole repo (Feb 20, 2026) and directly
  addresses the Fig. 2 "Vahid's to do" item, but its markdown only says "more points
  (higher k values)" — I could not confirm a finished/conclusive phase-transition
  result from its structure alone. **Likely still active work-in-progress.**
- `code-april17-version/topology-contagion/notebooks/bias_analysis.ipynb`,
  `concat_df.ipynb`, `testcodes.ipynb`, `global_clustering_coeffecient.ipynb` — smaller
  utility/scratch notebooks not obviously tied to a specific figure; not copied into
  `clean_code/` since they look like scratch work, but noted here for completeness.

## Remaining / Started-but-Unfinished Analyses

1. **Fig. 2 phase-transition follow-up** ("Vahid's to do: what if we increase step 3 to
   4,5,... maybe phase transition"): the code to do this now exists and looks fairly
   mature — `eph_computation/eph_computation_multiprocess_with_truncation.py` sweeps
   `up_to_step_infections = [2,3,4,5,6,7,8,9,10]`, and
   `notebooks/fig2_phase_transition_truncated_steps.ipynb` is the newest file in the
   entire repository (Feb 20, 2026), so this looks like the most recently active thread
   of work. Whether a clean phase-transition signal was actually found/written up could
   not be confirmed without deeper notebook-output inspection (not done at this pass).

2. **SI Fig. S17(second)/S18 "relative complex contagion regime" (marked TBD in the
   paper)**: the `relative_threshold/` folder is a separate exploratory branch with its
   own simulation engine (relative threshold φ, `ceil(φ·degree)`), EPH computation, and
   notebooks computing metrics like Earth Mover's Distance and a "Temporal Mixing Index"
   between simple/complex adoption timelines — these look like intermediate diagnostics,
   not a finished q/θ sweep with classification accuracy plotted against the main-text
   Figs. 3/4 style. **Update:** SI Fig. S19 specifically ("check if power of TDS is
   higher than EPH") is **not** part of this unfinished set — its code is finished, see
   the TDS correction above (`small_world_results_TDS_vs_EPH.ipynb`). It's the q/θ →
   EPH-distribution side (S17-second/S18) that still **appears unfinished**, consistent
   with those specific "TBD" captions. This branch was copied into
   `clean_code/relative_threshold_exploratory/` separately from the main pipeline rather
   than merged in, to keep the unfinished/experimental status visible. Main-text figures
   do not use relative threshold — only the SI mentions it.

3. **refine_feed_back_sep26/refine-feedback-priority-sorted.txt**: this is an automated
   AI-review ("refine.ink") of the paper's **prose/text**, not of code — all 10 preview
   comments are about wording/notation clarity in the Methods and SI text (e.g., the
   cycle-length off-by-one in the 7-cycle example in SI §S2, the L(C) vs. EPH lifetime
   sign convention, the coning-trick algorithm description in SI Algorithm 1, the
   Fig. 3A/4A axis-label mismatch with the Methods text in S8.2). None of these describe
   unfinished code — they are writing fixes for the manuscript. Worth relaying to the
   authors as text-only TODOs; only 10 of an estimated 24 comments were unlocked (the
   file itself says a full paid review would surface 14 more).

4. Numerous near-duplicate `result_*` output folders across `code-april17-version/`
   (e.g. `result_3loopsstop_10%init_nodes-v3`, `result_5loopsstop_10%init_nodes-v3`,
   `result_50stop_10init_nodes`, `result_85stop_10%init_nodes-v3[-archived/-limited-v4]`,
   etc.) suggest many parameter-sweep iterations (stop ratio, seed count, loop-count
   cutoff) were tried and abandoned/superseded before settling on the `stop_ratio=0.85`,
   `10%`-initial-nodes convention used by the final scripts copied here. These raw
   output CSVs were not copied into `clean_code/` (they are regenerable via the scripts
   here) but are left in place in the original folders for audit if needed.

## Note on scope of this pass
Given the very large number of files (170+ Python/notebook files across 9+ version
folders, many multi-MB notebooks with embedded output), this pass prioritized: (a) file
mtimes/sizes to identify the most recent/complete version of each logical script, (b)
reading full source for all Python scripts ultimately included here, and (c) reading
notebook markdown-cell titles (not full cell-by-cell code/output) to infer notebook
purpose and figure mapping. Where inference rather than direct confirmation was used,
this is flagged with "likely"/"appears" above. A deeper pass re-executing notebooks or
diffing cell-by-cell against the paper's actual saved figure PDFs would increase
confidence in the traceability table.
