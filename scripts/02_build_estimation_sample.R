# ==============================================================================
# 02_build_estimation_sample.R — Single source of truth for the estimation
#                                sample and the treatment definition
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# WHY THIS SCRIPT EXISTS
#
# The estimation panel and the treatment variables used to be rebuilt inside six
# separate scripts (06, 10, 11, 11b, 11c, 11d), each with its own copy of the
# treatment rule. Duplicated construction is a silent-divergence risk: changing
# the treatment definition in one place leaves the other five estimating
# something subtly different, with no error raised. This script builds all of it
# once, saves it, and every downstream script loads it instead of rebuilding.
#
# It also resolves a real inconsistency that the duplication had hidden.
# Reversal (schools that expand and later contract) was measured two ways:
#   - on the ANNUAL census history  -> 32.0% of ever-treated schools
#   - on the six outcome waves only -> 22.1%
# Reversal is a property of the annual treatment history, so the annual measure
# is the correct one; measuring it on six irregularly spaced waves misses
# reversals that happen between waves. Both are exported here under explicit
# names (`reverts_annual`, `reverts_waves`) so the distinction is visible rather
# than accidental, and `reverts_annual` is the one downstream code should use.
#
# INPUT   data/processed/panel.rds     (unified annual school x year panel)
#         data/processed/ra_panel.rds  (annual census, treatment source)
# OUTPUT  data/processed/estimation_sample.rds   school x wave, ready to estimate
#         data/processed/school_attributes.rds   one row per school
#         data/processed/treatment_annual.rds    treatment on the annual grid
#         data/processed/baseline_covariates.rds one row per school in the panel
#
# One table per unit of analysis, and nothing downstream derives any of it: the
# annual panel is 01's, the three views above are this script's, and every other
# script only reads. 07, 08 and 09 used to rebuild the panel, the outcome scaling,
# the wave index and `gid` themselves -- so `gid` indexed a different set of
# schools in each -- and 09 derived the baseline covariates from a different
# source argument than this script did.
#
# The construction below deliberately reproduces the previous in-script code
# line for line (including the order in which `gid` is assigned) so that this
# refactor is output-identical. See the validation block at the end.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))

panel    <- readRDS(file.path(DIR_PROCESSED, "panel.rds"))
ra_panel <- readRDS(file.path(DIR_PROCESSED, "ra_panel.rds"))
setorder(ra_panel, school_id, year)          # shift() below depends on this order

# ---- Treatment parameters (single place to change them) ---------------------
THR_MAIN <- 0.50      # headline binary: >=50% of primary enrolment
THR_HIGH <- 0.75      # alternative threshold (the antimode of the share)
THR_ANY  <- 1e-9      # "any expansion" (>0), the Nistal-style benchmark

# ---- The binary treatment definition ------------------------------------------
# A school is treated in year t if its share is at or above the threshold in t,
# full stop. The indicator states the observed annual state of the policy and
# nothing else, and that is what makes the treatment path the non-absorbing
# 0-1-0-1 the estimator is built for.
#
# A two-year persistence ("sustained") rule was the primary until August 2026 and
# was retired entirely on 12 August 2026, together with the script that compared
# the two (`06o`, retired from the project and not part of this repository).
# Three properties decided it: the
# rule used future information (t depended on t+1); it deleted one-year episodes
# that did happen; and it could not be satisfied in the last census year, so
# schools crossing the threshold in 2025 were recorded as untreated for want of a
# 2026 to confirm them -- a property of where the data end, not of the treatment.
#
# The definition is applied to the FULL annual history and then carried to the
# outcome waves, never computed on the waves alone.
#
# MISSING SHARES ARE NA, which is a fix rather than a convention: a school-year of
# UNKNOWN treatment must not enter the comparison pool as an untreated one. That
# is 7,305 school-years annually and 202 school-waves in the estimation sample.
make_treatment <- function(share_col, thr = THR_MAIN) {
  ra_panel[, .(school_id, year,
               treated_s = fifelse(is.na(get(share_col)), NA_integer_,
                                   as.integer(get(share_col) >= thr)))]
}

# Full-day share is derived here so every downstream script sees the same one.
ra_panel[, share_full := fifelse(!is.na(enroll_primary) & enroll_primary > 0,
                                 n_full_day / enroll_primary, NA_real_)]

# PRIMARY, and its full-day arm. The name is deliberate: `main_treatment` is the
# variable of the main specification, and every estimator robustness varies its
# own dimension while holding THIS variable fixed. One grep therefore answers
# "does this script estimate the main specification's treatment".
main_treatment      <- make_treatment("share_expanded", THR_MAIN); setnames(main_treatment,      "treated_s", "main_treatment")
main_treatment_full <- make_treatment("share_full",     THR_MAIN); setnames(main_treatment_full, "treated_s", "main_treatment_full")
# Alternative cut-offs are NOT built here. 04 constructs them at the point of
# use, so that varying the cut-off varies the cut-off alone.
# ---- Annual treatment history ------------------------------------------------
# Some analyses need the treatment on the FULL annual grid rather than on the six
# outcome waves: adoption timing for the intensity figure (06), reversal
# accounting (11), and the adoption-cohort definitions (11c, 11d). Exported here
# so that no downstream script has to re-derive the treatment.
treatment_annual <- merge(main_treatment, main_treatment_full,
                          by = c("school_id", "year"), all = TRUE)
setorder(treatment_annual, school_id, year)

# ---- Estimation sample (school x outcome wave) ------------------------------
# NOTE ON `gid`: it is assigned over the merged ANNUAL panel BEFORE the rows
# without outcomes are dropped, exactly as the previous in-script code did. The
# values therefore index all schools in the census, not only those that end up
# in the estimation sample. Keeping this order is what makes the refactor
# output-identical; do not move the assignment below the filter.
d <- merge(panel[, .(school_id, year, score_lang, score_math, share_expanded)],
           main_treatment, by = c("school_id", "year"), all.x = TRUE)
d <- merge(d, main_treatment_full, by = c("school_id", "year"), all.x = TRUE)
d[, `:=`(y_lang  = score_lang / 100,
         y_math  = score_math / 100,
         share_c = pmin(pmax(share_expanded, 0), 1),   # continuous dose, capped
         wave    = match(year, WAVES),
         gid     = as.integer(factor(school_id)))]
# The primary treatment decides the sample. Both definitions are NA on exactly the
# same school-years (wherever the share is unknown), so this filter selects the
# same rows either way; it names the primary so it stays true if that changes.
d <- d[!is.na(main_treatment) & !is.na(y_lang)]

# Baseline covariates: one time-invariant vector per school, fixed at 2016.
# Computed ONCE here and saved as baseline_covariates.rds, because 09 works on the
# annual panel (a different unit of analysis) and needs the same vectors. It used
# to call make_baseline_covars() itself, which is how the same quantity ended up
# being derived twice from two different sources.
baseline_covars <- make_baseline_covars(panel, ra_panel)
d <- merge(d, baseline_covars, by = "school_id", all.x = TRUE)

# Province, carried so the estimators can cluster on it. Treatment is assigned by
# the province (national programme with provincial adhesion, provincial funding
# decisions switching it on and off), so that is the level at which residuals are
# plausibly correlated and the level the headline clusters on. It is time-
# invariant per school and 01 normalises it, so one value per school is exact
# rather than a choice; stop() if that ever ceases to hold, because a split
# province silently becomes a spurious extra cluster.
prov <- unique(ra_panel[!is.na(province), .(school_id, province)], by = NULL)
if (anyDuplicated(prov$school_id))
  stop("province is not unique within school; clustering key would be split")
d <- merge(d, prov, by = "school_id", all.x = TRUE)
if (anyNA(d$province)) stop("province missing for some schools in the estimation sample")

setorder(d, gid, wave)

# ---- School-level attributes ------------------------------------------------
# Reversal, adoption cohort and the derived subsamples used across the analysis.
reverts_seq <- function(x) { x <- x[!is.na(x)]; length(x) >= 2 && any(diff(x) < 0) }

# Adoption date, reversal and re-adoption are properties of the annual path, so
# they are derived from the treatment the project estimates: the primary
# indicator, and nothing else. A script that reads
# `first_treated_year` or `cohort_2023` is therefore describing the treatment it
# is estimating, which is what makes the diagnostics and the estimand agree.
path_attrs <- function(tv) {
  a <- treatment_annual[!is.na(get(tv))][order(school_id, year)]
  a[, .(ever_treated_annual = any(get(tv) == 1),
        reverts_annual      = reverts_seq(get(tv)),
        first_treated_year  = if (any(get(tv) == 1)) min(year[get(tv) == 1]) else NA_integer_,
        n_census_years      = .N), by = school_id]
}
attr_annual <- path_attrs("main_treatment")

attr_waves <- d[, .(
  ever_treated_waves = any(main_treatment == 1),
  reverts_waves      = reverts_seq(main_treatment),            # <- understates; kept for comparison
  first_treated_wave = if (any(main_treatment == 1)) min(wave[main_treatment == 1]) else NA_integer_,
  n_waves            = .N), by = .(school_id, gid)]

schools <- merge(attr_waves, attr_annual, by = "school_id", all.x = TRUE)
for (v in c("reverts_annual", "ever_treated_annual"))
  schools[is.na(get(v)), (v) := FALSE]
schools[, `:=`(
  cohort_2023   = !is.na(first_treated_year) & first_treated_year == 2023,
  non_reverting = reverts_annual == FALSE)]

# ---- Predetermined socioeconomic covariate ----------------------------------
# Mother's secondary education averaged over the waves that precede the school's
# OWN first adoption, so the covariate cannot have been affected by the
# treatment. Two construction choices carry this and both matter.
#
#   (1) "Pre" is defined on the ANNUAL treatment history, not on the first treated
#       WAVE. A school that adopts in 2017 and is assessed in 2018 has a
#       post-treatment 2018 wave even if 2018 is not one of its treated waves
#       (the share may have dipped below the threshold). Dating "pre" by the wave
#       would count that 2018 observation as pre-treatment. It is not.
#
#   (2) The measure is WAVE-DEMEANED before averaging. Maternal education rises
#       steeply across waves -- about +5 percentage points per wave-pair even
#       among schools that never change treatment -- so a raw mean over each
#       school's own pre-window would encode WHEN that window fell. Schools that
#       adopt early would get systematically lower values than late adopters and
#       never-treated schools, turning the covariate into a proxy for adoption
#       timing, which is the opposite of what a predetermined covariate is for.
#       Demeaning expresses each school as a deviation from the national mean of
#       its own wave, which is comparable across windows of different length and
#       date. Never-treated schools use every wave they appear in.
#
# Schools with no pre-treatment wave get NA. They are left-censored: no
# pre-treatment information about them exists, and returning NA is the honest
# outcome rather than silently using a post-treatment value. It is also the point
# of the construction -- it makes that contamination visible instead of absorbing
# it. `n_pre_waves` is carried so the cost is reportable.
mom_w <- panel[!is.na(pct_mother_secondary) & year %in% WAVES,
               .(school_id, year, mom = pct_mother_secondary)]
mom_w[, mom_dm := mom - mean(mom), by = year]
mom_w <- merge(mom_w, schools[, .(school_id, first_treated_year)], by = "school_id")
pre_ses <- mom_w[is.na(first_treated_year) | year < first_treated_year,
                 .(bl_mom_pre = mean(mom_dm), n_pre_waves = .N), by = school_id]
schools         <- merge(schools, pre_ses, by = "school_id", all.x = TRUE)
d               <- merge(d, pre_ses, by = "school_id", all.x = TRUE)
baseline_covars <- merge(baseline_covars, pre_ses, by = "school_id", all.x = TRUE)
setorder(d, gid, wave)

# ---- Validation -------------------------------------------------------------
# A refactor that silently changes the sample is caught here.
#
# Values updated on 10 August 2026, when the primary treatment became the
# CONTEMPORANEOUS indicator and missing shares became NA under both definitions
# (see the definitions block above). Two things moved and both are intended:
#   - the sample lost the 202 school-waves (6 schools) whose treatment share is
#     unknown and which the sustained rule used to record as untreated;
#   - every path attribute is now measured on the primary rather than on the
#     sustained rule, so ever-treated goes 7,801 -> 10,749, reverters 2,496 ->
#     5,031 and the 2023 cohort 2,119 -> 3,695.
#
# The persistence rule itself was removed on 12 August 2026, with the script that
# compared the two definitions (`06o`) retired from the project. Every path
# attribute below is measured on `main_treatment` and there is no second rule.
# The previous constants were 99,899 / 22,941 / 7,801 / 2,496 / 2,119 / 20,822 /
# 20,464, set on 30 July when the panel was restricted to primary schools.
expect <- function(label, got, want, tol = 0) {
  ok <- abs(got - want) <= tol
  ok
}
ok <- c(
  expect("school-waves in estimation sample", nrow(d), 99697),
  expect("schools in estimation sample",      uniqueN(d$gid), 22935),
  expect("outcome waves",                     uniqueN(d$wave), 6),
  expect("ever-treated (annual, all schools)", attr_annual[ever_treated_annual == TRUE, .N], 10749),
  expect("reverters (annual, all schools)",    attr_annual[ever_treated_annual == TRUE & reverts_annual, .N], 5031),
  # Scope matters here: the count is over schools that sit at least one Aprender
  # wave, not over the whole census, and it is the one to quote when describing
  # what the diagnostic exclusion removes.
  expect("2023 cohort (in estimation sample)", schools[cohort_2023 == TRUE, .N], 3695),
  expect("schools left excluding 2023 cohort", schools[cohort_2023 == FALSE, .N], 19240),
  expect("non-reverting (est. sample)",        schools[non_reverting == TRUE, .N], 17979))
if (!all(ok)) stop("Validation mismatch: the estimation sample changed. ",
                      "Investigate before re-running any estimation script.")

# ---- Save -------------------------------------------------------------------
saveRDS(d,                file.path(DIR_PROCESSED, "estimation_sample.rds"))
saveRDS(schools,          file.path(DIR_PROCESSED, "school_attributes.rds"))
saveRDS(treatment_annual, file.path(DIR_PROCESSED, "treatment_annual.rds"))
saveRDS(baseline_covars,  file.path(DIR_PROCESSED, "baseline_covariates.rds"))
