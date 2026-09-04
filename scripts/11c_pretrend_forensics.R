# ==============================================================================
# 11c_pretrend_forensics.R — Why does FEct reject its own pre-trend model?
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# 11/11b left one question open, and it is the one that decides how much weight
# the FEct results can carry: FEct's joint pre-period F-test rejects for BOTH
# subjects, including Language, which is otherwise well behaved on every other
# diagnostic (CDH placebos clean at -0.001 and -0.013, both n.s.). That
# inconsistency has to be explained before anything is rewritten.
#
# Four questions, each answered by a specific comparison:
#   Q1 Is the rejection driven by a single lead?
#        -> per-lead estimates with confidence intervals.
#   Q2 Is it economically large, or just a very large sample with high power?
#        -> lead magnitudes against the 0.02-0.06 SD effects being claimed,
#           reported next to the number of treated observations per lead.
#   Q3 Does it survive on the sample CDH actually uses?
#        -> refit on the non-reverting sample used in 11b.
#   Q4 Does it come from the 2023 adoption cohort, i.e. the SAME cohort that
#      generates the CDH Mathematics placebo (07b, 06c)?
#        -> refit excluding schools whose first adoption year is 2023.
#
# Q4 is the one that could unify the whole story. If dropping the 2023 cohort
# removes FEct's Language rejection AND narrows the estimator disagreement, then
# there is ONE identification problem in these data (a single cohort adopting
# across the pandemic wave) rather than two unrelated ones, and the write-up
# collapses into a single narrative instead of a list of caveats.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(fect); library(polars) })
# Randomness: per-call via set_call_seed() (00_setup.R), not a global seed.
# NBOOTS comes from 00_setup.R.

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

# ---- Sample definitions ------------------------------------------------------
ann <- treatment_annual[!is.na(main_treatment)][order(school_id, year)]

nonrev_ids <- schools[non_reverting == TRUE, school_id]      # the 11b sample
coh23_ids  <- schools[cohort_2023 == TRUE, school_id]        # 2023 adopters

SAMPLES <- list(
  "full"                = d,
  "non-reverting (CDH)" = d[school_id %in% nonrev_ids],
  "excl. 2023 cohort"   = d[!school_id %in% coh23_ids])

fect_att <- function(o) if (is.null(o) || is.null(o$est.avg)) c(NA_real_, NA_real_) else as.numeric(o$est.avg)[1:2]

# Fit FEct (FE imputation) with leave-one-out and the joint pre-period F-test.
# Wrapped in tryCatch so a single failed fit does not abort the whole grid.
fit_fect <- function(dat, yname) {
  set_call_seed(paste(yname, nrow(dat)))
  args <- list(formula = as.formula(sprintf("%s ~ main_treatment", yname)),
               data = as.data.frame(dat[, .(gid, wave, y_lang, y_math, main_treatment)]),
               index = c("gid", "wave"), method = "fe", force = "two-way",
               se = TRUE, nboots = NBOOTS, parallel = TRUE, min.T0 = 1,
               seed = SEED, loo = TRUE)
  tryCatch(do.call(fect::fect, args),
           error = function(e) { NULL })
}
fp <- function(o) if (is.null(o) || is.null(o$test.out$f.p)) NA_real_ else as.numeric(o$test.out$f.p)[1]

# ---- Q1/Q2: per-lead estimates and magnitudes (full sample) -----------------
leads <- rbindlist(lapply(c("y_lang", "y_math"), function(yn) {
  o <- fit_fect(d, yn)
  if (is.null(o) || is.null(o$est.att)) return(NULL)
  x <- as.data.table(as.data.frame(o$est.att), keep.rownames = "k")
  x[, k := as.integer(as.character(k))]
  x[k < 0, .(subject = if (yn == "y_lang") "language" else "mathematics", k,
             att = round(ATT, 4), se = round(S.E., 4),
             ci_lo = round(CI.lower, 4), ci_hi = round(CI.upper, 4),
             p = signif(p.value, 3), n_obs = count,
             sig = fifelse(p.value < 0.05, "*", ""))]
}), fill = TRUE)
fwrite(leads, file.path(DIR_TABLES, "11c_pretrend_leads.csv"))

# ---- Q1: is the rejection driven by a single lead? --------------------------
# Answered directly by the per-lead table above: Language has TWO adjacent
# significant leads of opposite sign (k = -3 at +0.038, k = -2 at -0.036), i.e.
# oscillation rather than a single outlier.
#
# A leave-one-lead-out F-test was attempted here but removed: fect's `pre.periods`
# argument is interpreted as a [min, max] RANGE, not a set, so it cannot exclude
# an interior lead and every "drop one" run silently tested the same full range.
# The per-lead estimates answer the question without it.

# ---- Q3/Q4: does the rejection survive on other samples? --------------------
grid <- rbindlist(lapply(names(SAMPLES), function(sn) {
  rbindlist(lapply(c("y_lang", "y_math"), function(yn) {
    o <- fit_fect(SAMPLES[[sn]], yn)
    a <- fect_att(o)
    data.table(sample = sn, subject = if (yn == "y_lang") "language" else "mathematics",
               n_schools = uniqueN(SAMPLES[[sn]]$school_id),
               fect_att = round(a[1], 4), fect_se = round(a[2], 4),
               pretrend_f_p = signif(fp(o), 3))
  }))
}), fill = TRUE)
fwrite(grid, file.path(DIR_TABLES, "11c_sample_grid.csv"))

