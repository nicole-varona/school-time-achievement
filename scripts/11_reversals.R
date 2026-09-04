# ==============================================================================
# 11_reversals.R — Responses to the NON-ABSORBING treatment
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Two robustness blocks, both prompted by supervisor feedback on the fact that
# the treatment switches on and off (schools expand, then contract):
#
#   PART A — FEct (Liu, Wang & Xu 2024, AJPS): COUNTERFACTUAL IMPUTATION.
#            A different identification logic from de Chaisemartin-D'Haultfoeuille:
#            fit the two-way model on UNTREATED observations, impute Y(0) for the
#            treated ones, average the residuals. It accommodates treatment
#            reversal natively, so it triangulates the headline CDH estimate
#            across estimator FAMILIES, not merely across variants of one family.
#            Its real payoff here is diagnostic:
#              - PLACEBO TEST : hides pre-treatment periods and tests whether the
#                               imputed "effect" there is zero.
#              - EQUIVALENCE  : inverts the burden of proof. A conventional
#                TEST (TOST)    placebo that fails to reject supports parallel
#                               trends only weakly (a low-power test rejects
#                               nothing; cf. Roth 2022). The equivalence test asks
#                               whether pre-treatment effects lie INSIDE a
#                               pre-specified bound, so passing it is positive
#                               evidence FOR parallel trends rather than mere
#                               absence of evidence against it. This speaks
#                               directly to the maths k=-2 placebo (see 07/07b).
#
#   PART B — NON-REVERTING SUBSAMPLE. The supervisor suggested dropping reverting
#            schools "if these reversals only happen in a few cases". They do not:
#            the accounting below quantifies exactly how many. Reported as a
#            ROBUSTNESS, never as the main specification, for two reasons stated
#            in the write-up: (i) it removes a large share of treated schools and
#            so changes the target population and the estimand; (ii) it conditions
#            on POST-treatment behaviour (a school's future treatment path), i.e.
#            selection, which is why the balance table at the end is reported.
#            It does buy one clean argument: on this subsample treatment is
#            ABSORBING BY CONSTRUCTION, so the Callaway-Sant'Anna irreversibility
#            assumption is satisfied rather than merely assumed.
#
# NOTE — where reversal is measured. 06_estimation.R flags reverters using the 6
#   OUTCOME WAVES only. Reversal is a property of the ANNUAL treatment history,
#   so a school can revert between waves and go undetected there. This script
#   measures reversal on the FULL ANNUAL RA history (2011-2025) and reports both,
#   which is why its non-reverting subsample is stricter than 06's.
#
# Event time is in OUTCOME WAVES, not equal calendar years (see 06_estimation.R).
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
for (p in c("fect", "DIDmultiplegtDYN", "did", "fixest"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(fect); library(did); library(fixest); library(polars); library(ggplot2)
})
# Randomness: per-call via set_call_seed() (00_setup.R), not a global seed.
# NBOOTS comes from 00_setup.R.
IFE_R  <- 2            # number of interactive factors, FIXED (no cross-validation)

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
stopifnot(!any(duplicated(d[, .(gid, wave)])))   # fect requires a unique id x time key

# ============================================================================
# PART A — FEct: counterfactual imputation with treatment reversal
# ============================================================================
# `min.T0 = 1` keeps schools with a single pre-treatment wave: the panel has only
# 6 waves and a large share of treated schools have few pre-periods (left
# censoring, see 07_pretrends.R), so the package default would discard much of
# the treated sample. The number of units actually used is reported for each fit.
D_FECT <- as.data.frame(d[, .(gid, wave, y_lang, y_math, main_treatment)])

# Defensive extractor: fect's slot names vary across versions, so pull the ATT
# from whichever of the known containers is present rather than assuming one.
fect_att <- function(o) {
  if (is.null(o)) return(c(NA_real_, NA_real_, NA_real_))
  if (!is.null(o$est.avg)) {
    v <- as.numeric(o$est.avg)                    # ATT, S.E., CI.lower, CI.upper, p
    return(c(v[1], v[2], if (length(v) >= 5) v[5] else NA_real_))
  }
  c(if (!is.null(o$att.avg)) as.numeric(o$att.avg)[1] else NA_real_, NA_real_, NA_real_)
}

# Dynamic (event-time) effects, as a tidy table.
fect_dynamic <- function(o, label) {
  if (is.null(o)) return(NULL)
  m <- o$est.att
  if (is.null(m)) return(NULL)
  dt <- as.data.table(as.data.frame(m), keep.rownames = "k")
  dt[, `:=`(spec = label, k = as.integer(as.character(k)))]
  dt[]
}

run_fect <- function(yname, tname = "main_treatment", label = "", method = "fe",
                     placebo = FALSE, placebo.period = NULL, loo = FALSE,
                     carryover = FALSE, carryover.period = NULL, carryover.rm = NULL,
                     tost.threshold = NULL, dat = D_FECT, nboots = NBOOTS) {
  set_call_seed(paste(label, method, yname))
  args <- list(formula = as.formula(sprintf("%s ~ %s", yname, tname)),
               data = dat, index = c("gid", "wave"), method = method,
               force = "two-way", se = TRUE, nboots = nboots, parallel = TRUE,
               min.T0 = 1, seed = SEED)
  # Fix the factor count instead of cross-validating it: with 6 waves CV is both
  # expensive and unstable, and r is reported as a stated choice.
  if (method == "ife") { args$r <- IFE_R; args$CV <- FALSE }
  if (placebo) { args$placeboTest <- TRUE; args$placebo.period <- placebo.period }
  if (loo)     args$loo <- TRUE
  if (carryover) { args$carryoverTest <- TRUE; args$carryover.period <- carryover.period }
  if (!is.null(carryover.rm))   args$carryover.rm   <- carryover.rm
  if (!is.null(tost.threshold)) args$tost.threshold <- tost.threshold
  out <- tryCatch(do.call(fect::fect, args),
                  error = function(e) { NULL })
  if (!is.null(out)) {
    a <- fect_att(out)
    # First run: expose the object's structure so the test extractors below can be
    # matched to this package version.
  }
  out
}

fect_res <- list()
fect_res$lang_fe <- run_fect("y_lang", label = "Language  (FE imputation)")
fect_res$math_fe <- run_fect("y_math", label = "Maths     (FE imputation)")

# Interactive fixed effects: relaxes parallel trends to a factor structure, i.e.
# unobserved time-varying confounders with heterogeneous loadings. A stricter test
# of whether the estimate survives a weaker identifying assumption.
fect_res$lang_ife <- run_fect("y_lang", label = "Language  (interactive FE)", method = "ife")
fect_res$math_ife <- run_fect("y_math", label = "Maths     (interactive FE)", method = "ife")

# ---- Diagnostics: placebo, equivalence, carryover ---------------------------
# placebo.period = c(-2,-1): hide the two waves before adoption and test whether
# the imputed effect there is zero. Deliberately centred on the maths k=-2
# placebo that 07/07b identified as the study's main vulnerability.
fect_res$lang_placebo <- run_fect("y_lang", label = "Language  (placebo -2 to -1)",
                                  placebo = TRUE, placebo.period = c(-2, -1))
fect_res$math_placebo <- run_fect("y_math", label = "Maths     (placebo -2 to -1)",
                                  placebo = TRUE, placebo.period = c(-2, -1))

# Leave-one-out pre-period estimates: what the TOST equivalence test is built on.
# Run TWICE. The package default picks the equivalence bound from the data, and
# here it lands around 0.16 SD -- three to eight times the effects this study
# estimates (0.02-0.06 SD), so "passing" at that bound is uninformative about
# whether pre-trends are small RELATIVE TO THE CLAIM. TOST_STRICT re-runs the test
# against a bound equal to the largest effect the study reports, which is the
# comparison that actually bears on the argument. Both are reported.
TOST_STRICT <- 0.05
fect_res$lang_equiv <- run_fect("y_lang", label = "Language  (equivalence/loo)", loo = TRUE)
fect_res$math_equiv <- run_fect("y_math", label = "Maths     (equivalence/loo)", loo = TRUE)
fect_res$lang_equiv_strict <- run_fect("y_lang", label = "Language  (equivalence, bound 0.05)",
                                       loo = TRUE, tost.threshold = TOST_STRICT)
fect_res$math_equiv_strict <- run_fect("y_math", label = "Maths     (equivalence, bound 0.05)",
                                       loo = TRUE, tost.threshold = TOST_STRICT)

# CARRYOVER — the decisive diagnostic for the FEct/CDH disagreement.
# FEct imputes Y(0) from the pool of UNTREATED observations. With a reverting
# treatment that pool contains post-reversal school-years: schools that expanded
# the day and later contracted. If any benefit persists after contraction
# (trained teachers, infrastructure, materials), those observations sit ABOVE
# their true never-treated counterfactual. Fitting Y(0) on them pushes the imputed
# counterfactual up, and therefore the treated residual -- the ATT -- DOWN. That
# is exactly the direction in which FEct departs from CDH here, and CDH is not
# exposed to it because it compares switchers with stayers around the switch
# instead of pooling all untreated periods.
#   carryoverTest : is the outcome still elevated 1-2 waves AFTER treatment exits?
#   carryover.rm  : drop those post-exit waves from the untreated pool and re-fit.
# If the ATT moves toward the CDH estimate once they are removed, carryover
# contamination is the explanation. If it does not, the disagreement is real and
# must be reported as such -- this test is set up to be able to fail.
fect_res$lang_carry <- run_fect("y_lang", label = "Language  (carryover test)",
                                carryover = TRUE, carryover.period = c(1, 2))
fect_res$math_carry <- run_fect("y_math", label = "Maths     (carryover test)",
                                carryover = TRUE, carryover.period = c(1, 2))
fect_res$lang_carry_rm <- run_fect("y_lang", label = "Language  (post-exit waves removed)",
                                   carryover.rm = 2)
fect_res$math_carry_rm <- run_fect("y_math", label = "Maths     (post-exit waves removed)",
                                   carryover.rm = 2)

# Test statistics live in $test.out; names verified against this fect version.
grab_test <- function(o, spec) {
  if (is.null(o)) return(NULL)
  t <- o$test.out
  g <- function(nm) if (!is.null(t) && !is.null(t[[nm]])) as.numeric(t[[nm]])[1] else NA_real_
  a <- fect_att(o)
  data.table(spec = spec, att = a[1], se = a[2],
             placebo_p         = g("placebo.p"),          # H0: pre-treatment effect = 0
             placebo_equiv_p   = g("placebo.equiv.p"),    # TOST: effect inside the bound
             f_stat            = g("f.stat"),             # joint test of all pre-periods
             f_p               = g("f.p"),
             tost_equiv_p      = g("tost.equiv.p"),
             tost_threshold    = g("tost.threshold"),
             carryover_p       = g("carryover.p"),        # H0: no effect after exit
             carryover_equiv_p = g("carryover.equiv.p"))
}
tests_tab <- rbindlist(list(
  grab_test(fect_res$lang_placebo,      "language placebo(-2 to -1)"),
  grab_test(fect_res$math_placebo,      "maths placebo(-2 to -1)"),
  grab_test(fect_res$lang_equiv,        "language equivalence (default bound)"),
  grab_test(fect_res$math_equiv,        "maths equivalence (default bound)"),
  grab_test(fect_res$lang_equiv_strict, "language equivalence (bound 0.05)"),
  grab_test(fect_res$math_equiv_strict, "maths equivalence (bound 0.05)"),
  grab_test(fect_res$lang_carry,        "language carryover(1-2 after exit)"),
  grab_test(fect_res$math_carry,        "maths carryover(1-2 after exit)")), fill = TRUE)
fwrite(tests_tab, file.path(DIR_TABLES, "11_fect_tests.csv"))

# Does removing the post-exit waves reconcile FEct with CDH?
carry_tab <- rbindlist(list(
  data.table(subject = "language",    spec = "FEct baseline",
             att = fect_att(fect_res$lang_fe)[1],       se = fect_att(fect_res$lang_fe)[2]),
  data.table(subject = "language",    spec = "FEct, post-exit waves removed",
             att = fect_att(fect_res$lang_carry_rm)[1], se = fect_att(fect_res$lang_carry_rm)[2]),
  data.table(subject = "mathematics", spec = "FEct baseline",
             att = fect_att(fect_res$math_fe)[1],       se = fect_att(fect_res$math_fe)[2]),
  data.table(subject = "mathematics", spec = "FEct, post-exit waves removed",
             att = fect_att(fect_res$math_carry_rm)[1], se = fect_att(fect_res$math_carry_rm)[2])),
  fill = TRUE)
carry_tab[, `:=`(att = round(att, 4), se = round(se, 4))]
fwrite(carry_tab, file.path(DIR_TABLES, "11_fect_carryover.csv"))

# ---- FEct summary + dynamic path --------------------------------------------
fect_row <- function(o, subject, method) {
  a <- fect_att(o)
  data.table(estimator = sprintf("FEct (%s)", method), subject = subject,
             att = a[1], se = a[2], p = a[3],
             n_units = if (!is.null(o$N)) as.integer(o$N) else NA_integer_)
}
fect_tab <- rbindlist(list(
  fect_row(fect_res$lang_fe,  "language",    "FE imputation"),
  fect_row(fect_res$math_fe,  "mathematics", "FE imputation"),
  fect_row(fect_res$lang_ife, "language",    "interactive FE"),
  fect_row(fect_res$math_ife, "mathematics", "interactive FE")), fill = TRUE)
fect_tab[, `:=`(att = round(att, 4), se = round(se, 4), p = round(p, 4))]
fwrite(fect_tab, file.path(DIR_TABLES, "11_fect_results.csv"))

dyn_tab <- rbindlist(list(
  fect_dynamic(fect_res$lang_fe, "language"),
  fect_dynamic(fect_res$math_fe, "mathematics")), fill = TRUE)
if (nrow(dyn_tab)) fwrite(dyn_tab, file.path(DIR_TABLES, "11_fect_dynamic.csv"))
# Figure 11_fect_dynamic.png is rendered from 11_fect_dynamic.csv by 91_figures.R
# in the shared style.

# ============================================================================
# PART B — Non-reverting subsample
# ============================================================================
# B.1 Accounting: how common are reversals, really? And of those that revert, how
#     many come back (treated -> untreated -> treated again) rather than exit for
#     good? Measured on the FULL ANNUAL RA history (the definition that makes the
#     treatment non-absorbing), and, for contrast, on the 6 outcome waves only
#     (06_estimation.R's view).

# Both measures are computed once in 02_build_estimation_sample.R and carried on
# `schools`, so this script reports them rather than re-deriving them:
#   reverts_annual -> the correct measure (annual census history)
#   reverts_waves  -> the six-wave view, which understates reversal
ann_all <- treatment_annual[!is.na(main_treatment)][order(school_id, year)]
n_ever_ann <- ann_all[, .(ever = any(main_treatment == 1)), by = school_id][ever == TRUE, .N]
n_rev_ann  <- ann_all[, .(rev = .N >= 2 && any(diff(main_treatment) < 0),
                          ever = any(main_treatment == 1)), by = school_id][ever & rev, .N]

# Re-adoption: the ordered treatment sequence (NA dropped) has two or more separate
# "on" blocks, i.e. the school is treated, drops out, and is treated again at least
# once. A strict subset of reverters -- how many toggle back on instead of exiting
# for good. Note the primary indicator does NOT filter these: it records the annual
# state, so a one-year
# reporting blips. Counted per school as the number of runs of 1s (>= 2 = re-adopter).
on_blocks <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(0L)
                           r <- rle(x); sum(r$values == 1) }
blk_ann  <- treatment_annual[order(school_id, year), .(nb = on_blocks(main_treatment)), by = school_id]
blk_wave <- d[order(school_id, wave),                .(nb = on_blocks(main_treatment)), by = school_id]
est_ids  <- unique(d$school_id)

acct <- data.table(
  scope = c("RA annual history (all schools)",
            "RA annual history (estimation sample)",
            "6 outcome waves only (estimation sample)"),
  ever_treated = c(n_ever_ann,
                   schools[ever_treated_annual == TRUE, .N],
                   schools[ever_treated_waves  == TRUE, .N]),
  reverters    = c(n_rev_ann,
                   schools[ever_treated_annual == TRUE & reverts_annual, .N],
                   schools[ever_treated_waves  == TRUE & reverts_waves,  .N]),
  re_adopters  = c(blk_ann[nb >= 2, .N],
                   blk_ann[school_id %in% est_ids][nb >= 2, .N],
                   blk_wave[nb >= 2, .N]))
acct[, `:=`(pct_reverting  = round(100 * reverters   / ever_treated, 1),
            pct_readopting = round(100 * re_adopters / ever_treated, 1))]
setcolorder(acct, c("scope", "ever_treated", "reverters", "pct_reverting",
                    "re_adopters", "pct_readopting"))
fwrite(acct, file.path(DIR_TABLES, "11_reversal_accounting.csv"))

# B.2 Build the subsample: never-treated schools + schools that adopt and never
#     revert anywhere in the annual history.
keep_ids <- schools[non_reverting == TRUE, school_id]
d_nr <- d[school_id %in% keep_ids]

# B.3 Re-estimate. CDH is the headline estimator, so the comparison that matters
#     is full sample vs non-reverting subsample under the SAME estimator.
run_cdh <- function(dat, yname, tname, label, effects = 2, placebo = 1) {
  set_call_seed(label)
  out <- tryCatch(DIDmultiplegtDYN::did_multiplegt_dyn(
    df = as.data.frame(dat), outcome = yname, group = "gid", time = "wave",
    treatment = tname, effects = effects, placebo = placebo, cluster = "gid",
    graph_off = TRUE),
    error = function(e) { NULL })
  if (!is.null(out)) {
    ate <- out$results$ATE
  }
  out
}
nr <- list()
nr$cdh_lang <- run_cdh(d_nr, "y_lang", "main_treatment", "Non-reverting -> Language")
nr$cdh_math <- run_cdh(d_nr, "y_math", "main_treatment", "Non-reverting -> Mathematics")

# CS on this subsample: treatment is absorbing BY CONSTRUCTION here, so the
# irreversibility assumption of Callaway-Sant'Anna is satisfied rather than
# assumed. Always-treated schools (adopted before the first wave) are dropped:
# they have no pre-period. Cohort G = first treated wave.
# Already non-reverting by construction, so the reverter filter is declared off
# rather than relied on to be a no-op.
d_nr_cs <- cs_sample(d_nr, drop_reverters = FALSE)

run_cs <- function(dat, yname, label, xformla = ~1, keep_rows = NULL) {
  cvars <- all.vars(xformla)
  set_call_seed(label)
  dd <- dat[!is.na(get(yname))]
  if (!is.null(keep_rows)) dd <- dd[eval(keep_rows, dd)]
  if (length(cvars)) dd <- dd[complete.cases(dd[, ..cvars])]
  fit <- tryCatch(
    att_gt(yname = yname, tname = "wave", idname = "gid", gname = "G",
           xformla = xformla, data = dd, control_group = "notyettreated",
           allow_unbalanced_panel = TRUE, bstrap = TRUE, cband = FALSE, biters = BITERS, est_method = "dr"),
    error = function(e) { NULL })
  if (is.null(fit)) return(NULL)
  agg <- suppressWarnings(aggte(fit, type = "simple", na.rm = TRUE))
  list(fit = fit, simple = agg, n_schools = uniqueN(dd$gid))
}
comp_cc <- quote(bl_comp_complete == TRUE)
nr$cs_lang <- run_cs(d_nr_cs, "y_lang", "Non-reverting -> Language",    XF_COMP, comp_cc)
nr$cs_math <- run_cs(d_nr_cs, "y_math", "Non-reverting -> Mathematics", XF_COMP, comp_cc)

# TWFE on the same subsample, to keep the estimator-vs-design contrast visible.
tw_l <- feols(y_lang ~ main_treatment | gid + wave, data = d_nr, cluster = ~gid)
tw_m <- feols(y_math ~ main_treatment | gid + wave, data = d_nr, cluster = ~gid)

grab_cdh <- function(r) if (is.null(r)) c(NA_real_, NA_real_) else r$results$ATE[1:2]
nr_tab <- rbindlist(list(
  data.table(estimator = "CDH", sample = "non-reverting", subject = "language",
             att = grab_cdh(nr$cdh_lang)[1], se = grab_cdh(nr$cdh_lang)[2]),
  data.table(estimator = "CDH", sample = "non-reverting", subject = "mathematics",
             att = grab_cdh(nr$cdh_math)[1], se = grab_cdh(nr$cdh_math)[2]),
  data.table(estimator = "CS (absorbing by construction)", sample = "non-reverting",
             subject = "language",
             att = if (is.null(nr$cs_lang)) NA_real_ else nr$cs_lang$simple$overall.att,
             se  = if (is.null(nr$cs_lang)) NA_real_ else nr$cs_lang$simple$overall.se),
  data.table(estimator = "CS (absorbing by construction)", sample = "non-reverting",
             subject = "mathematics",
             att = if (is.null(nr$cs_math)) NA_real_ else nr$cs_math$simple$overall.att,
             se  = if (is.null(nr$cs_math)) NA_real_ else nr$cs_math$simple$overall.se),
  data.table(estimator = "TWFE", sample = "non-reverting", subject = "language",
             att = coef(tw_l)[["main_treatment"]], se = se(tw_l)[["main_treatment"]]),
  data.table(estimator = "TWFE", sample = "non-reverting", subject = "mathematics",
             att = coef(tw_m)[["main_treatment"]], se = se(tw_m)[["main_treatment"]])), fill = TRUE)
nr_tab[, `:=`(att = round(att, 4), se = round(se, 4))]
fwrite(nr_tab, file.path(DIR_TABLES, "11_nonreverter_results.csv"))

# B.4 Selection diagnostic. Dropping reverters conditions on POST-treatment
#     behaviour, so the two groups need not be comparable. Reported so the
#     robustness is read as what it is: a different population, not a cleaner
#     version of the same one.
bal_src <- unique(merge(d[, .(school_id, gid, sector, area, log_enroll,
                              bl_score_lang, bl_score_math, bl_pct_repeater)],
                        schools[, .(school_id, ever_treated = ever_treated_annual, reverts_annual)],
                        by = "school_id"), by = "school_id")
bal <- bal_src[ever_treated == TRUE, .(
  n              = .N,
  pct_public     = round(100 * mean(grepl("stat|p.blic|publi", sector, ignore.case = TRUE), na.rm = TRUE), 1),
  pct_rural      = round(100 * mean(grepl("rural", area, ignore.case = TRUE), na.rm = TRUE), 1),
  mean_log_enrol = round(mean(log_enroll, na.rm = TRUE), 3),
  mean_bl_lang   = round(mean(bl_score_lang, na.rm = TRUE), 3),
  mean_bl_math   = round(mean(bl_score_math, na.rm = TRUE), 3),
  mean_repeater  = round(mean(bl_pct_repeater, na.rm = TRUE), 3)),
  by = .(group = fifelse(reverts_annual, "reverters", "non-reverters"))]
fwrite(bal, file.path(DIR_TABLES, "11_reverter_balance.csv"))

# Fitted model objects are NOT saved. They were, and nothing ever read them:
# 855 MB of caches sat in data/processed alongside the three real tables and
# made the data flow hard to follow. Every number this script produces is in
# its CSV outputs, and the registry (92) indexes them.
