# ==============================================================================
# run_all.R — Reproduce the full analysis pipeline, in order
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Runs every numbered script that is part of the reproducible pipeline, in the
# order their outputs depend on one another, in a FRESH R session each (so no
# object leaks between scripts and each is verified to stand on its own).
#
# Usage, from the project root:
#     Rscript scripts/run_all.R
#
# Requirements:
#   - the restricted microdata under data/ (see README; not in the repo);
#   - the R packages listed in the README;
#
# Not run as steps here:
#   - 00_setup.R, 90_style.R              : modules, sourced by other scripts;
#   - 92b_palette_check.R                 : design-time check, run by hand when the
#                                           palette in 90_style.R changes;
# 06q is the slowest step by a wide margin (part A alone is about an hour). It is
# in the pipeline rather than run by hand because the manuscript quotes it: the
# province bootstrap of §7.2 and the appendix table of inference procedures both
# read 06q_cdh_bootstrap.csv, so a run that skips it cannot rebuild the PDF.
# 91_figures.R runs LAST: it reads the saved results and renders every figure in
# one shared style (90_style.R), so figures never re-estimate and always match.
#
# Bootstrapped standard errors are reproducible: randomness is seeded per call
# via set_call_seed() (00_setup.R), so results do not depend on run order.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")

t_start <- Sys.time()

# The pipeline, in dependency order. 02 must run before any estimation script,
# since it builds the estimation sample they all load.
PIPELINE <- c(
  "01_build_panel.R",              # raw microdata -> unified annual panel
  "02_build_estimation_sample.R",  # panel -> estimation sample + treatment defs
  "03_descriptives.R",             # sample overview, coverage, score distributions
  "04_threshold_sensitivity.R",    # treatment-threshold sensitivity (benchmark)
  "05_reversion_check.R",          # non-absorbing-treatment diagnostics
  "06_estimation.R",               # MAIN estimation (CDH, CS, TWFE)
  "06c_covid_robustness.R",        # pandemic-wave / 2023-cohort robustness
  "06d_scale_robustness.R",        # outcome-scale robustness
  "06e_grade_alignment.R",         # treatment-grade alignment
  "06f_clustering_sensitivity.R",  # clustering level and effective clusters
  "06g_horizon_sensitivity.R",     # aggregation window and three CDH options
  "06h_composition_check.R",       # does treatment change school composition?
  "06i_switcher_populations.R",    # who supports each horizon, and how they differ
  "06j_exposure_partition.R",      # event time in years: partition by exposure
  "06l_exposure_validity.R",       # elapsed time vs years actually treated
  "06m_calendar_support.R",        # what each candidate outcome calendar identifies
  "06n_calendar_estimates.R",      # the same estimand under each calendar x definition
  "06p_spell_length_classes.R",    # does the effect depend on how long the first spell lasts
  "06k_leave_one_province.R",      # is the estimate carried by one jurisdiction?
  "06q_cdh_bootstrap.R",           # province-clustered bootstrap of the primary estimator
  "07_pretrends.R",                # event studies, placebos, HonestDiD
  "07b_pretrend_decomposition.R",  # Maths placebo by cohort and region
  "07c_dynamic_profile.R",         # full dynamic profile
  "08_hora_mas.R",                 # genuine 6h+ vs the 1h+ programme
  "08b_recording_rates.R",         # how well the census records each length of day
  "09_trajectories.R",             # secondary outcome: grade progression
  "10_heterogeneity.R",            # pre-registered heterogeneity
  "11_reversals.R",                # FEct; carryover; non-reverting subsample
  "11b_like_for_like.R",           # all estimators on one common sample
  "11c_pretrend_forensics.R",      # the 2023-cohort diagnosis
  "11d_identifying_variation.R",   # is the convergence mechanical?
  "91_figures.R",                  # all manuscript figures, one shared style
  "92_estimate_registry.R"         # one registry of every estimate produced above
)

# --- Run each script in its own process; stop on the first failure ----------
run_one <- function(script) {
  path <- here::here("scripts", script)
  if (!file.exists(path)) stop("missing script: ", path)
  cat(sprintf("─── %-30s %s ", script, format(Sys.time(), "%H:%M:%S")))
  t0 <- Sys.time()
  # A fresh Rscript process: guarantees no state leaks between steps.
  status <- system2("Rscript", path, stdout = "", stderr = "")
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  if (status != 0)
    stop(sprintf("\n%s exited with status %d after %.1f min. Pipeline halted.",
                 script, status, dt))
  cat(sprintf("done %.1f min\n", dt))
  data.table::data.table(script = script, minutes = dt)
}

if (!dir.exists(here::here("data")))
  stop("Project root not found: expected a 'data/' folder next to 'scripts/'.")

log <- data.table::rbindlist(lapply(PIPELINE, run_one))

total <- round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1)
cat(sprintf("\n─── PIPELINE COMPLETE: %d scripts in %.1f min.\n", nrow(log), total))
print(log[order(-minutes)])
