# ==============================================================================
# 11b_like_for_like.R — Why do the estimators disagree?
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# 11_reversals.R produced a sign conflict: CDH and Callaway-Sant'Anna give small
# POSITIVE estimates, TWFE and FEct give NEGATIVE ones. Two candidate
# explanations were tested there and BOTH FAILED:
#   - treatment reversal      : on the non-reverting subsample, where treatment is
#                               absorbing by construction, TWFE is still strongly
#                               negative and CDH/CS still positive.
#   - carryover contamination : the carryover test does not reject (p = 0.27 / 0.38),
#                               and removing post-exit waves makes FEct MORE
#                               negative, i.e. the opposite of that hypothesis.
#
# What remains is the counterfactual TREND MODEL. CDH and CS identify from LOCAL
# comparisons against not-yet-treated schools. TWFE and FEct (FE imputation)
# instead impose ONE common set of time effects on every school, so if treated and
# comparison schools follow different trends -- which 07b documented for the 2023
# cohort and the Northeast -- the imposed common trend is wrong and the estimate
# inherits that error. The prediction this script tests: relaxing the common-trend
# restriction (interactive fixed effects, which let schools load differently on
# common factors) should move the FEct estimate toward the CDH/CS estimates.
#
# Everything below is run on ONE identical sample -- the non-reverting subsample --
# so that sample composition cannot explain any remaining gap, and dynamic effects
# are compared at MATCHED event times (FEct k = 0 is the first treated wave, which
# is CDH's Effect_1; FEct k = 1 is CDH's Effect_2).
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({
  library(fect); library(did); library(fixest); library(polars)
})
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

# ---- The single common sample: schools that never revert in the ANNUAL history
keep_ids <- schools[non_reverting == TRUE, school_id]
d_nr <- d[school_id %in% keep_ids]

D_F <- as.data.frame(d_nr[, .(gid, wave, y_lang, y_math, main_treatment)])
fect_att <- function(o) if (is.null(o) || is.null(o$est.avg)) c(NA_real_, NA_real_) else as.numeric(o$est.avg)[1:2]

run_fect <- function(yname, method = "fe", label = "") {
  set_call_seed(paste(label, method, yname))
  args <- list(formula = as.formula(sprintf("%s ~ main_treatment", yname)), data = D_F,
               index = c("gid", "wave"), method = method, force = "two-way",
               se = TRUE, nboots = NBOOTS, parallel = TRUE, min.T0 = 1, seed = SEED)
  if (method == "ife") { args$r <- 2; args$CV <- FALSE }
  o <- tryCatch(do.call(fect::fect, args),
                error = function(e) { NULL })
  o
}
run_cdh <- function(yname, label, effects = 3) {
  set_call_seed(label)
  o <- tryCatch(DIDmultiplegtDYN::did_multiplegt_dyn(
    df = as.data.frame(d_nr), outcome = yname, group = "gid", time = "wave",
    treatment = "main_treatment", effects = effects, placebo = 2, cluster = "gid",
    graph_off = TRUE), error = function(e) { NULL })
  o
}

res <- list()
res$fe_lang  <- run_fect("y_lang", "fe",  "Language  ")
res$fe_math  <- run_fect("y_math", "fe",  "Maths     ")
res$ife_lang <- run_fect("y_lang", "ife", "Language  ")
res$ife_math <- run_fect("y_math", "ife", "Maths     ")
res$cdh_lang <- run_cdh("y_lang", "Language")
res$cdh_math <- run_cdh("y_math", "Mathematics")

# CS on the same sample (irreversibility satisfied by construction).
d_cs <- cs_sample(d_nr, drop_reverters = FALSE)
run_cs <- function(yname, label) {
  set_call_seed(label)
  dd <- d_cs[!is.na(get(yname))]
  fit <- tryCatch(att_gt(yname = yname, tname = "wave", idname = "gid", gname = "G",
                         xformla = ~1, data = dd, control_group = "notyettreated",
                         allow_unbalanced_panel = TRUE, bstrap = TRUE, cband = FALSE, biters = BITERS, est_method = "dr"),
                  error = function(e) { NULL })
  if (is.null(fit)) return(NULL)
  a <- suppressWarnings(aggte(fit, type = "simple", na.rm = TRUE))
  a
}
res$cs_lang <- run_cs("y_lang", "Language")
res$cs_math <- run_cs("y_math", "Mathematics")

tw_l <- feols(y_lang ~ main_treatment | gid + wave, data = d_nr, cluster = ~gid)
tw_m <- feols(y_math ~ main_treatment | gid + wave, data = d_nr, cluster = ~gid)

# ---- One table, one sample --------------------------------------------------
cdh_ate <- function(o) if (is.null(o)) c(NA_real_, NA_real_) else o$results$ATE[1:2]
cs_ate  <- function(a) if (is.null(a)) c(NA_real_, NA_real_) else c(a$overall.att, a$overall.se)
row_ <- function(est, assumption, subj, v)
  data.table(estimator = est, counterfactual_trend = assumption, subject = subj,
             att = v[1], se = v[2])
tab <- rbindlist(list(
  row_("CDH (did_multiplegt_dyn)", "local: not-yet-treated", "language",    cdh_ate(res$cdh_lang)),
  row_("CDH (did_multiplegt_dyn)", "local: not-yet-treated", "mathematics", cdh_ate(res$cdh_math)),
  row_("Callaway-Sant'Anna",       "local: not-yet-treated", "language",    cs_ate(res$cs_lang)),
  row_("Callaway-Sant'Anna",       "local: not-yet-treated", "mathematics", cs_ate(res$cs_math)),
  row_("FEct (interactive FE)",    "factor structure",       "language",    fect_att(res$ife_lang)),
  row_("FEct (interactive FE)",    "factor structure",       "mathematics", fect_att(res$ife_math)),
  row_("FEct (FE imputation)",     "common two-way trend",   "language",    fect_att(res$fe_lang)),
  row_("FEct (FE imputation)",     "common two-way trend",   "mathematics", fect_att(res$fe_math)),
  row_("TWFE",                     "common two-way trend",   "language",    c(coef(tw_l)[["main_treatment"]], se(tw_l)[["main_treatment"]])),
  row_("TWFE",                     "common two-way trend",   "mathematics", c(coef(tw_m)[["main_treatment"]], se(tw_m)[["main_treatment"]]))))
tab[, `:=`(att = round(att, 4), se = round(se, 4))]
fwrite(tab, file.path(DIR_TABLES, "11b_like_for_like.csv"))

# ---- Matched event times: FEct k vs CDH Effect_(k+1) ------------------------
dyn_fect <- function(o, subj) {
  if (is.null(o) || is.null(o$est.att)) return(NULL)
  x <- as.data.table(as.data.frame(o$est.att), keep.rownames = "k")
  x[, .(subject = subj, source = "FEct (FE imputation)", k = as.integer(as.character(k)),
        att = round(ATT, 4), se = round(S.E., 4))]
}
dyn_cdh <- function(o, subj) {
  if (is.null(o)) return(NULL)
  e <- as.data.table(as.data.frame(o$results$Effects), keep.rownames = "lab")
  nm <- names(e)
  x <- e[, .(subject = subj, source = "CDH", k = seq_len(.N) - 1L,
             att = round(get(nm[2]), 4), se = round(get(nm[3]), 4))]
  x[]
}
dyn <- rbindlist(list(dyn_fect(res$fe_lang, "language"), dyn_fect(res$fe_math, "mathematics"),
                      dyn_cdh(res$cdh_lang, "language"), dyn_cdh(res$cdh_math, "mathematics")),
                 fill = TRUE)
fwrite(dyn, file.path(DIR_TABLES, "11b_dynamic_matched.csv"))

# Fitted model objects are NOT saved. They were, and nothing ever read them:
# 855 MB of caches sat in data/processed alongside the three real tables and
# made the data flow hard to follow. Every number this script produces is in
# its CSV outputs, and the registry (92) indexes them.
