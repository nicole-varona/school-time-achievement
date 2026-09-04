# ==============================================================================
# 11d_identifying_variation.R — How much identifying variation does the
#                               2023-cohort exclusion actually remove?
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# 11c showed that excluding the 2023 adoption cohort makes the estimators
# converge on small, insignificant estimates and clears every pre-trend
# diagnostic. The obvious objection, raised in review, is that convergence has
# two possible causes:
#     (1) the estimators all become correct once a contaminated comparison goes, or
#     (2) so little identifying variation is left that there is nothing to
#         disagree about.
# A standard-error comparison does not settle this, because SEs respond to sample
# size while identification responds to TREATMENT TIMING VARIATION. Those are not
# the same thing: dropping the single largest adoption cohort could leave the
# number of observations almost intact while gutting the variation that actually
# identifies the effect.
#
# This script measures the second quantity directly:
#   - ever-treated schools, and how many are removed
#   - distinct adoption cohorts and their sizes (the timing variation itself)
#   - 0 -> 1 switches observed within the outcome-wave panel
#   - switchers supporting each event-time horizon (CDH reports these, which is
#     precisely why: they are the units doing the identifying work at each k)
#
# READING THE TWO SWITCHER COUNTS: they are NOT the same quantity, and the
# horizon counts can exceed the 0 -> 1 count without contradiction. `observed_0to1`
# counts adoptions only; the per-horizon counts include BOTH directions, because
# with non-absorbing treatment a reversal also identifies an effect.
#
# If the switcher counts survive in proportion, the convergence is informative.
# If they collapse, the convergence is mechanical and the interpretation in 11c
# has to be weakened accordingly. The check is set up so that it can fail.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(polars) })
# Randomness: per-call via set_call_seed() (00_setup.R), not a global seed.

# ---- Load the estimation sample ---------------------------------------------
# Panel, treatment definitions and baseline covariates are built ONCE in
# 02_build_estimation_sample.R. Rebuilding them here (as this script used to do)
# is what allowed the same quantity to be measured differently across scripts.
#   d                : school x outcome wave, ready to estimate
#   schools          : one row per school (reversal, adoption cohort, subsamples)
#   treatment_annual : treatment on the full annual census grid
d                <- readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds"))
schools          <- readRDS(file.path(DIR_PROCESSED, "school_attributes.rds"))
treatment_annual <- readRDS(file.path(DIR_PROCESSED, "treatment_annual.rds"))

ann <- treatment_annual[!is.na(main_treatment)][order(school_id, year)]
coh23_ids  <- schools[cohort_2023 == TRUE, school_id]        # 2023 adopters

SAMPLES <- list("full" = d, "excl. 2023 cohort" = d[!school_id %in% coh23_ids])

# ---- Structural measures of identifying variation ---------------------------
describe <- function(dat, nm) {
  x <- copy(dat); setorder(x, gid, wave)
  # First treated WAVE per school = its adoption cohort (0 = never treated here).
  coh <- x[, .(firstG = if (any(main_treatment == 1)) min(wave[main_treatment == 1]) else 0L,
               ever    = any(main_treatment == 1)), by = gid]
  # 0 -> 1 switches actually OBSERVED between consecutive waves: the transitions
  # a DiD design can use. Always-treated schools contribute none.
  sw <- x[, .(n_sw = sum(diff(main_treatment) > 0)), by = gid]
  data.table(
    sample              = nm,
    schools             = uniqueN(x$gid),
    school_waves        = nrow(x),
    ever_treated        = coh[ever == TRUE, .N],
    treated_school_waves= x[main_treatment == 1, .N],
    adoption_cohorts    = coh[ever == TRUE & firstG > 1, uniqueN(firstG)],
    observed_0to1       = sw[, sum(n_sw)],
    schools_with_switch = sw[n_sw > 0, .N])
}
struct <- rbindlist(lapply(names(SAMPLES), function(n) describe(SAMPLES[[n]], n)))

ret <- struct[2, ] ; base <- struct[1, ]
pct <- data.table(metric = c("schools", "school_waves", "ever_treated",
                             "treated_school_waves", "observed_0to1", "schools_with_switch"),
                  full = unlist(base[, .(schools, school_waves, ever_treated,
                                         treated_school_waves, observed_0to1, schools_with_switch)]),
                  excl_2023 = unlist(ret[, .(schools, school_waves, ever_treated,
                                             treated_school_waves, observed_0to1, schools_with_switch)]))
pct[, pct_retained := round(100 * excl_2023 / full, 1)]
fwrite(pct, file.path(DIR_TABLES, "11d_identifying_variation.csv"))

# Cohort sizes: where the timing variation actually lives.
cohort_tab <- rbindlist(lapply(names(SAMPLES), function(n) {
  x <- SAMPLES[[n]]; setorder(x, gid, wave)
  c1 <- x[, .(firstG = if (any(main_treatment == 1)) min(wave[main_treatment == 1]) else 0L), by = gid]
  # `by` is what makes this a cohort table: without it .N counts the whole
  # subset and the same total is emitted once per school (the bug this had
  # until 07/08 — 12,181 rows, 12 distinct, `schools` constant across cohorts).
  c1[firstG > 0, .(sample = n, cohort_year = WAVES[.BY[[1]]], schools = .N),
     by = .(cohort_wave = firstG)][order(cohort_wave)]
}), fill = TRUE)
fwrite(cohort_tab, file.path(DIR_TABLES, "11d_cohort_sizes.csv"))

# ---- Switchers per horizon: the variation CDH actually uses -----------------
# did_multiplegt_dyn reports the number of switchers supporting each dynamic
# effect, which is the cleanest available measure of identifying variation at
# each horizon. Compared before and after the exclusion.
run_cdh <- function(dat, yname, label) {
  set_call_seed(label)
  tryCatch(DIDmultiplegtDYN::did_multiplegt_dyn(
    df = as.data.frame(dat), outcome = yname, group = "gid", time = "wave",
    treatment = "main_treatment", effects = 3, placebo = 2, cluster = "gid",
    graph_off = TRUE), error = function(e) { NULL })
}
sw_tab <- rbindlist(lapply(names(SAMPLES), function(n) {
  o <- run_cdh(SAMPLES[[n]], "y_lang", sprintf("%s -> Language", n))
  if (is.null(o)) return(NULL)
  e <- as.data.frame(o$results$Effects)
  # The command pads its row names to a fixed width, so "Effect_1" arrives with
  # trailing spaces and every consumer has to trim before matching on it. Trimmed
  # once here instead: the manuscript and 91_figures.R both key on this column.
  data.table(sample = n, horizon = trimws(rownames(e)),
             estimate = round(e[["Estimate"]], 4),
             switchers = e[["Switchers"]])
}), fill = TRUE)
fwrite(sw_tab, file.path(DIR_TABLES, "11d_switchers_by_horizon.csv"))

