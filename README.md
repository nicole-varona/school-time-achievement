# Extended School Time and Student Achievement in Argentine Primary Education

Replication code for an MSc dissertation (MY499, MSc Social Research Methods,
Department of Methodology, London School of Economics).

**Question.** Does the expansion of extended- and full-day schooling in Argentine
primary education improve grade-6 achievement in mathematics and language?

**Design.** School-level staggered difference-in-differences, linking the annual
school census (*Relevamiento Anual*, 2011–2025) to the census-based *Aprender*
assessment (waves 2016, 2018, 2021, 2022, 2023, 2025). The primary estimator is
chosen by **design fit** rather than by result — an argument the reader can check
against the data instead of having to take on trust: because adoption is
non-absorbing (47% of ever-treated schools later contract on the annual record,
and 32% adopt, contract and adopt again) and Callaway–Sant'Anna assumes absorbing adoption, 
the primary estimator is de Chaisemartin–D'Haultfœuille `did_multiplegt_dyn`, with
Callaway–Sant'Anna and FEct as alternative identification strategies and two-way
fixed effects as a benchmark.

**Treatment.** The main specification is a binary indicator, `main_treatment`:
a school is treated in a year if at least 50% of its primary enrolment is in
extended- or full-day schooling that year. It states the observed annual
state of the policy, which is what makes the treatment path the non-absorbing
0-1-0-1 series the estimator is built for. Sensitivities vary the cut-off
(>0, 40%, 60%, 75%) and the arm (combined, full-day), holding the rule fixed.
A continuous dose-response on the share was removed: the share measures 
how many pupils are under the regime, not how long their day is, 
so it is a dose in coverage where the mechanism at issue is a dose in time. 
Intensity enters instead through the hours-verified arms of `08_hora_mas.R`. 
Every treatment definition is built in one place, `02_build_estimation_sample.R`, 
and every other script reads it: an estimator robustness varies its own dimension 
while holding `main_treatment` fixed, so a `treat_*` variable appears only in 
the five scripts whose purpose is to vary the definition itself.

---

## Data availability statement

**The microdata are not included in this repository and cannot be redistributed.**

The analysis uses restricted administrative microdata provided under a data-access
agreement by the Subsecretaría de Información y Evaluación Educativa, within the
Secretaría de Educación of Argentina's Ministerio de Capital Humano, comprising:

- *Relevamiento Anual* (annual school census), 2011–2025;
- *Aprender* standardised assessment, grade-6 primary, six waves 2016–2025.

Records are linked through an anonymised fictitious school identifier assigned by
the data provider. The files contain no direct identifiers, but school-level
administrative records are treated as **restricted** under the LSE information
classification: they are not redistributable, and combining them with external
sources could in principle raise re-identification risk. Accordingly:
- `data/` is excluded from version control in full, including derived panels in
  `data/processed/`;
- only aggregate results (estimation output, summary tables, figures) are
  published, in `output/`;
- researchers wishing to reproduce the analysis from raw data should approach the
  Subsecretaría de Información y Evaluación Educativa directly.

**Ethics.** The project uses secondary, pre-anonymised administrative data and
involves no contact with human participants. The relevant ethics consideration is
the handling of restricted microdata, mitigated as described above.

---

## Repository structure

```
.
├── MY499-62628.qmd                   # Manuscript source (renders to PDF)
├── MY499-62628.pdf                   # Compiled manuscript
├── references.bib                    # Bibliography
├── scripts/                          # Analysis pipeline (see below)
├── output/
│   ├── tables/                       # Aggregate results (.csv)
│   └── figures/                      # Figures (.png)
├── .here                             # Project-root marker (here package)
└── README.md
```

`data/` exists locally but is not tracked (see the statement above).

---

## Pipeline

All scripts assume the working directory is the project root and are intended to
run in numerical order. Each sources `00_setup.R`, which holds configuration,
paths, column maps and shared helpers.

| Script | Purpose |
|---|---|
| `run_all.R` | Runs the whole reproducible pipeline below, in order |
| `00_setup.R` | Configuration, paths, readers, column maps, shared helpers |
| `01_build_panel.R` | Builds the unified school × year panel (2011–2025), restricted to primary schools. Derives no treatment indicator: it produces `share_expanded`, which the definitions in 02 are built from |
| `02_build_estimation_sample.R` | Single source of truth for the estimation sample and every treatment definition (loaded by all estimation scripts) |
| `03_descriptives.R` | Sample overview, coverage, score distributions |
| `04_threshold_sensitivity.R` | Treatment definition: the primary estimator re-run at alternative contemporaneous cut-offs (>0 / 40 / 50 / 60 / 75%), rather than confounded |
| `05_reversion_check.R` | Reversal diagnostics on the annual census series: how often schools expand and later contract |
| `06_estimation.R` | Main estimation: CDH (primary), Callaway–Sant'Anna, TWFE benchmark. Writes every placebo each specification admits, not only the first |
| `06c_covid_robustness.R` | Wave composition: COVID waves, the 2022 sample wave, and the 2023-cohort diagnostic |
| `06d_scale_robustness.R` | Outcome-scale robustness (anchored vs within-wave) |
| `06e_grade_alignment.R` | Whole-primary treatment vs grade-6 outcome alignment: the pigeonhole bound on how often grade-6 exposure is forced, and the estimate on the schools where it is. Also splits by 6- vs 7-year primary, which is **not reported**: that split is a provincial split and cannot separate alignment from the five-hour contamination (see the note in the script) |
| `06f_clustering_sensitivity.R` | Provincial structure: do provincial shocks confound (province-specific year trends) or only induce dependence (clustering level, effective clusters) |
| `06g_horizon_sensitivity.R` | Aggregation window; every placebo the main specification admits (asserted to equal the event-study leads); and the CDH options the design allows (`controls`, `trends_lin`, `same_switchers`) |
| `06h_composition_check.R` | Does the treatment change school composition, and what the SES covariate contributes |
| `06i_switcher_populations.R` | Who supports each horizon, and how those populations differ before treatment |
| `06j_exposure_partition.R` | Event time in years rather than waves: the primary estimator partitioned by calendar time since first adoption |
| `06k_leave_one_province.R` | Is the estimate carried by a single jurisdiction? |
| `06l_exposure_validity.R` | Calendar time since the switch vs years actually treated: under a non-absorbing treatment these are different quantities, and only the second is exposure |
| `06m_calendar_support.R` | What population each candidate outcome calendar identifies: support, cohorts, and the calendar time since adoption inside each horizon |
| `06n_calendar_estimates.R` | The same estimand under each candidate outcome calendar, by calendar time since adoption, with the placebo-cancellation and leave-one-cohort-out diagnostics |
| `06p_spell_length_classes.R` | Does the effect depend on how long the first treated spell lasts? Switchers partitioned by first-spell length (1 / 2 / 3+ years), estimated within class. Descriptive across classes, not a dose-response: spell length is realised after adoption. Not cited in the manuscript |
| `06q_cdh_bootstrap.R` | Province-clustered inference for the primary estimator: analytic errors at the school and at the province against a nonparametric pairs cluster bootstrap over whole provinces, with the failed draws and the provincial support each draw retains. The slowest step in the pipeline |
| `07_pretrends.R` | Event studies, placebo tests, HonestDiD sensitivity |
| `07b_pretrend_decomposition.R` | Decomposition of the Mathematics placebo by cohort and region |
| `07c_dynamic_profile.R` | Full post-treatment dynamic profile |
| `08_hora_mas.R` | Separates genuine 6h+ extension from the federal 1h+ programme |
| `08b_recording_rates.R` | How well the census records each length of school day: recording rates by directors' reported hours and by jurisdiction, and the size of the five-hour margin left in the comparison group |
| `09_trajectories.R` | Secondary outcome: grade progression (composition check) |
| `10_heterogeneity.R` | Pre-registered heterogeneity (policy extension) |
| `11_reversals.R` | FEct estimator; carryover test; non-reverting subsample |
| `11b_like_for_like.R` | All estimators on one common sample; matched event time |
| `11c_pretrend_forensics.R` | Why the pre-trend test rejects: the 2023-cohort diagnosis |
| `11d_identifying_variation.R` | Is the post-exclusion convergence mechanical? |
| `90_style.R` | Shared figure/table style: palette, theme, fonts (a module) |
| `91_figures.R` | Renders every manuscript figure from saved results, one style |
| `92_estimate_registry.R` | One registry of every ATT-level estimate the pipeline produces, on a fixed schema, with the role each plays. Runs last, reading the committed tables |
| `92b_palette_check.R` | Colour-vision validation of the figure palette (design-time check; writes nothing, run by hand when `90_style.R` changes, not in `run_all.R`) |

---

## Reproducing the analysis

Requires R (developed under 4.6.0). Every step runs in R.

```bash
# From the project root — runs the whole pipeline in order:
Rscript scripts/run_all.R
```

`run_all.R` runs each numbered script in a fresh R process, in dependency order,
starting with `01_build_panel.R` and `02_build_estimation_sample.R` (which every
estimation script then loads). `06q_cdh_bootstrap.R` is much the slowest step —
part A alone takes about an hour — and it is in the pipeline rather than run by
hand because the manuscript quotes its output.

**Package versions** are pinned in `renv.lock` (R 4.6.0; key packages
`data.table`, `did`, `DIDmultiplegtDYN`, `fect`, `fixest`, `HonestDiD`, `haven`,
`ggplot2`, `polars`, `here`, `ragg`, `systemfonts`, `kableExtra`). To recreate the
environment, install `renv` and run `renv::restore()`. Note that `polars` is
installed from R-universe (`https://rpolars.r-universe.dev`), not CRAN.

Paths are anchored to the project root with the `here` package, so scripts can be
run from any working directory inside the project (the `.here` marker locates the
root); there is no need to `setwd()`.

**Reproducible standard errors.** Bootstrapped standard errors are seeded per
call (via `set_call_seed()` in `00_setup.R`), not once per script, so a result
depends only on its own inputs and never on the order in which specifications
run. Point estimates are deterministic and every standard error reproduces
exactly under the pinned package versions.

The primary estimator `did_multiplegt_dyn` is called with analytic clustered 
standard errors throughout: the one specification that used its `bootstrap` option 
— the continuous dose — was removed, and with it the `continuous` and `bootstrap` 
arguments of the local `run_cdh()` wrapper.

---

## Outputs

Tables are written to `output/tables/` and figures to `output/figures/`,
named by the script that produces them (for example, `11b_like_for_like.csv` is
produced by `11b_like_for_like.R`). The manuscript reads these files directly, so
the figures and tables in the PDF are generated from the committed outputs rather
than transcribed by hand. Most figures are rendered by `91_figures.R` from the
plot-data CSV of the same stem (e.g. `07_event_study_cdh.png` from
`07_event_study_cdh.csv`), so those figures have a companion data table. Three
exceptions plot the raw panel rather than saved estimates and are therefore drawn
where the panel is read: `03_provincial_timing.png` in `03_descriptives.R` (whose
data are committed separately as `03_provincial_timing.csv`), and
`05_school_trajectories.png` and `05_national_totals.png` in `05_reversion_check.R`,
which have no companion table. They go through the shared theme and the shared
`save_fig()` of `90_style.R` all the same, so every figure in the repository is
written by one routine at one resolution.

### Output → exhibit map

Traceability runs both ways: every exhibit in the manuscript is generated from a
committed output (or built inline from documentary facts, where the table is a
definition rather than an estimate), and every committed output is classified
below as a body exhibit, an appendix exhibit, a plot-data companion, an output
quoted in prose, or a supporting output not cited.

**Body exhibits** — 20 (14 tables, 6 figures), in manuscript order. Section
numbers refer to the manuscript.

| Exhibit | Section | Source output(s) |
|---|---|---|
| `tbl-framework` | §3 Context | *(built inline: the five legal and programmatic instruments)* |
| `tbl-jurisdictions` | §3 Context | `03_jurisdiction_profile.csv` |
| `fig-provincial-timing` | §3 Context | `03_provincial_timing.png` (drawn in `03`; data in `03_provincial_timing.csv`) |
| `tbl-match` | §4.1 Sources and linkage | `01_match_rate_by_wave.csv` |
| `tbl-treatdef` | §4.2 Treatment construction | *(built inline: every construction of treatment and its role)* |
| `fig-trajectories` | §5.2 Estimator choice | `05_school_trajectories.png` (drawn in `05`; no companion table) |
| `tbl-time` | §5.3 Event time | *(built inline: event time, calendar time, exposure)* |
| `tbl-threats` | §5 Empirical strategy | *(built inline: threats and the check that addresses each)* |
| `fig-transitions` | §6.1 What the treatment is | `05_transitions.png` ← `05_transitions.csv` |
| `tbl-provincias` | §6.1 What the treatment is | `03_provincial_timing.csv`, `08b_recording_by_province.csv` (+ documentary columns) |
| `fig-contamination` | §6.1 What the treatment is | `08_contamination_2023.png` ← `08_contamination_2023.csv` |
| `tbl-arms` | §6.1 What the treatment is | `08_arm_results.csv`, `06_estimation_main.csv` |
| `tbl-main` | §6.2 The effect on achievement | `06_estimation_main.csv`, `06j_exposure_partition.csv`, `11b_like_for_like.csv` |
| `tbl-timing-placebo` | §6.2 The effect on achievement | `06j_exposure_partition.csv` |
| `fig-eventstudy` | §6.3 Diagnostics and estimator disagreement | `07_event_study_cdh.png` ← `07_event_study_cdh.csv` (+ `11d_switchers_by_horizon.csv` for the support counts) |
| `tbl-decomp` | §6.3 Diagnostics and estimator disagreement | `07b_cohort_decomposition.csv` |
| `fig-hetero` | §6.5 Heterogeneity | `10_heterogeneity.png` ← `10_heterogeneity_results.csv` |
| `tbl-rob-treatdef` | §7.1 Robustness | `04_threshold_sensitivity_cdh.csv` |
| `tbl-rob-specs` | §7.1 Robustness | `06_estimation_main.csv` |
| `tbl-rob-covid` | §7.1 Robustness | `06c_covid_robustness.csv` |

**Appendix exhibits** — 21 (17 tables, 4 figures), in manuscript order.

| Exhibit | Appendix | Source output(s) |
|---|---|---|
| `tbl-cohorts` | A | `03_adoption_cohorts.csv` |
| `fig-app-genuine` | A | `08_event_study_genuine.png` ← `08_event_study_genuine.csv` |
| `tbl-balance` | A | `03_baseline_balance.csv` |
| `fig-intensity` | A | `06_adoption_intensity.png` ← `06_adoption_intensity.csv` |
| `tbl-reversion` | A | `11_reversal_accounting.csv` |
| `tbl-switching` | A | `03_treatment_transitions.csv` + `03_switches_per_school.csv` |
| `tbl-contamination` | A | `08_contamination_2023.csv` |
| `fig-app-fectdyn` | B | `11_fect_dynamic.png` ← `11_fect_dynamic.csv` |
| `tbl-app-matched` | B | `11b_dynamic_matched.csv` |
| `tbl-app-carryover` | B | `11_fect_carryover.csv` |
| `tbl-app-leads` | B | `11c_pretrend_leads.csv` |
| `tbl-app-idvar` | B | `11d_identifying_variation.csv` |
| `tbl-app-cdhboot` | B | `06q_cdh_bootstrap.csv` |
| `tbl-app-resampling` | B | `06f_resampling.csv` |
| `tbl-fect` | B | `11_fect_results.csv` |
| `tbl-lfl` | B | `11b_like_for_like.csv` |
| `tbl-excl2023` | B | `11c_sample_grid.csv` |
| `tbl-app-traj` | C | `09_trajectory_results.csv` |
| `tbl-app-hetero` | C | `10_heterogeneity_results.csv` |
| `tbl-app-interaction` | C | `10_interaction_test.csv` |
| `fig-traj` | C | `09_event_study_trajectories.png` ← `09_event_study_trajectories.csv` |

**Plot-data companions of a cited figure.** Read by `91_figures.R`, not by the
manuscript, and committed so that every rendered figure has the table behind it:
`05_transitions.csv`, `06_adoption_intensity.csv`,
`08_event_study_genuine.csv`, `09_event_study_trajectories.csv`, and
`11_fect_dynamic.csv`.

**Quoted in prose without an exhibit of their own.** These outputs supply numbers
the manuscript states inline: `03_assessment_coverage.csv`, `03_extended_day_by_year.csv`,
`03_full_day_by_year.csv`, `06_placebos.csv`, `06d_scale_robustness.csv`,
`06e_g6_certain_bound.csv`, `06e_grade_alignment.csv`, `06f_clustering.csv`,
`06f_effective_clusters.csv`, `06f_inference_diagnostics.csv`,
`06g_cdh_options.csv`, `06g_horizon_sensitivity.csv`, `06g_max_placebos.csv`,
`06j_exposure_classes.csv`, `06k_leave_one_province.csv`, `06l_exposure_validity.csv`,
`07_event_study_cdh.csv`, `07_event_study_cdh_controls.csv`,
`07_event_study_cs_covariates.csv`, `07b_kminus2_by_region.csv`,
`08b_five_hour_margin.csv`, `08b_hours_by_sector.csv`, `08b_labelling_rule.csv`,
`08b_recording_by_hours.csv`, `11_fect_tests.csv`, `11_nonreverter_results.csv`,
`11d_cohort_sizes.csv`, `11d_switchers_by_horizon.csv`.

**Supporting and descriptive outputs, not cited** (kept for the working document):

- *Descriptives (`03_descriptives.R`):* `03_sample_overview.csv`,
  `03_missingness_by_wave.csv`, `03_summary_stats_pooled.csv`,
  `03_scores_by_wave.csv`, `03_threshold_sensitivity.csv`,
  `03_extended_vs_full.csv`, `03_reversions_by_province.csv`,
  `03_treated_by_wave.csv`, and the descriptive figures
  `03_treatment_intensity.png`, `03_enrolment_dist.png`,
  `03_score_lang_density.png`, `03_score_math_density.png`.
- *Reversal (`05`):* `05_national_totals.png`.
- *Calendar and definition work (`06m`–`06p`):* `06m_*`, `06n_*`, `06p_*`.
  These support the decision, still open, of which outcome calendar the headline
  should use; the sign of the ATT depends on it.
- *Provincial structure and inference (`06f`, `06g`, `06q`):*
  `06f_province_trends.csv`, `06g_ses_trends_feasibility.csv`, and
  `06q_bootstrap_support.csv`, which records zero failed draws and the
  provincial support each retained.

  Two working files of `06q` are deliberately **not** committed, and the
  `.gitignore` says why: `06q_bootstrap_draws.csv`, the 299 replications per
  subject, an intermediate of one estimate that nothing reads; and
  `06q_cdh_bootstrap_partA.csv`, which is a cache rather than a result — `06q`
  re-reads it only under `PARTS="B"`, to skip the hour part A costs, so shipping
  it would hand a reader a way to skip the step they came to reproduce. Both are
  rebuilt by re-running `scripts/06q_cdh_bootstrap.R`, and what they were
  produced for is published in `06q_cdh_bootstrap.csv`.
- *Switcher populations and composition (`06h`, `06i`):* `06h_*`, `06i_*`.
- *Produced by the pipeline, not shown in the manuscript:*
  `07c_dynamic_profile.png` with its companion `07c_dynamic_profile.csv`,
  `07_event_study_cdh_controls.png`, `07_event_study_cs_covariates.png`, and
  `07_honestdid_panel.png` with the four
  `07_honestdid_{language,mathematics}_{uncontrolled,baselinecovariates}.csv`
  behind it. The two event-study **tables** are a different matter and stay under
  *quoted in prose* above: the manuscript uses those numbers without the figure.
  Also here, superseded by the combined panel, the four `07_honestdid_*.png`
  drawn one per specification in `07_pretrends.R`.
- *Other uncited estimation outputs:* `09_event_study_cdh.csv`,
  `10_quartile_direction.csv`, `10_switchers_per_cell.csv`,
  `11_reverter_balance.csv`.
- *Not reported, and deliberately:* `06e_by_primary_structure.csv` — the 6- vs
  7-year primary split is a provincial split and cannot separate grade-6
  alignment from the five-hour contamination; see the note in the script.
- *Registry:* `92_estimate_registry.csv` indexes every estimate with its role.
