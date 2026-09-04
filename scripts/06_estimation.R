# ==============================================================================
# 06_estimation.R — Main estimation
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Three-layer strategy (per methodological review):
#   MAIN   — de Chaisemartin-D'Haultfoeuille (did_multiplegt_dyn), BINARY
#            CONTEMPORANEOUS treatment (>=50% share in the year, no persistence
#            rule; decision of 10/08). Chosen by DESIGN FIT: adoption is
#            non-absorbing and Callaway-Sant'Anna assumes absorbing adoption.
#   ROB 3  — Callaway-Sant'Anna, doubly robust, on the absorbing subsample:
#            conditions on baseline covariates (XF_COMP, see 00_setup.R).
#            Threshold sensitivity lives in 04_threshold_sensitivity.R, run with
#            the PRIMARY estimator rather than with CS.
#   BENCH  — static TWFE (Nistal & Edo) as a (biased) reference.
#
# Treatment: primary-enrolment share in extended/full-day, annual state --
# the indicator states the observed annual state of the policy, which is what
# makes the treatment path the non-absorbing 0-1-0-1 series the estimator is
# built for. 04 varies the CUT-OFF of that indicator; the two-year persistence
# rule that preceded it was retired on 12/08 and is no longer a specification.
# Outcome  : Aprender language/maths SCHOOL MEANS (repeated cross-sections of
#            different grade-6 students each wave — a school-level exposure/ITT
#            estimand, not student value-added), in SD units of the anchored 2016
#            scale (/100). Time = the 6 outcome waves, indexed 1-6.
# NOTE on event time: the index is in OUTCOME WAVES, not equal calendar years
#   (2018->2021 = 3y, 2021->2022 = 1y, 2023->2025 = 2y). Dynamic "Effect_k" must
#   be read in wave units and translated to calendar time (mapping printed below).
# SEs: analytical clustered-at-school where the estimator provides them (CDH);
#   CS uses the multiplier bootstrap (its standard for uniform bands).
#
# COVID: the headline is the full six-wave panel; wave exclusions are sensitivity
#   only, in 06c_covid_robustness.R, because dropping waves changes the estimand.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
for (p in c("DIDmultiplegtDYN", "did", "fixest"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
# DIDmultiplegtDYN uses the 'polars' backend and requires it to be attached.
suppressPackageStartupMessages({library(did); library(fixest); library(polars); library(ggplot2)})
# Randomness: seeded per call via set_call_seed() (see 00_setup.R), not globally,
# so each bootstrap is independent of the order in which specifications run.

# ---- Load the estimation sample ---------------------------------------------
# The panel, every treatment definition and the baseline covariates are built
# ONCE in 02_build_estimation_sample.R and loaded here. They used to be rebuilt
# inside this script (and five others), which is how the same quantity ended up
# being measured two different ways across scripts without any error being
# raised. Run 02 first if these files are missing.
#   d                : school x outcome wave, ready to estimate
#   treatment_annual : treatment on the full annual census grid (for the
#                      adoption-intensity figure at the end of this script)
#   ra_panel         : raw census, only for `share_expanded` in that figure
d                <- readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds"))
treatment_annual <- readRDS(file.path(DIR_PROCESSED, "treatment_annual.rds"))
ra_panel         <- readRDS(file.path(DIR_PROCESSED, "ra_panel.rds"))

# Event time is in OUTCOME WAVES, not equal calendar years: translate accordingly.

# ---- Baseline covariates -----------------------------------------------------
# Merged into `d` by script 02 (one time-invariant vector per school, fixed at
# 2016). The covariate SETS (XF_STRUCT, XF_COMP_NOSES, XF_COMP) are defined in
# 00_setup.R so that the five scripts that condition on them cannot diverge.
# XF_COMP is complete-case: schools first tested after 2016 lack composition.

results <- list()

# ============================================================================
# 1. Reference: de Chaisemartin-D'Haultfoeuille, BINARY (non-absorbing)
# ============================================================================
# Analytic clustered SEs throughout, which is the package default and is
# deterministic. The `continuous` and `bootstrap` arguments were removed with the
# dose-response specification (11/08): nothing here draws on the RNG any more.
run_cdh <- function(dat, yname, tname, label, effects = 2,
                    placebo = max_placebo(effects, length(WAVES)),
                    extra = NULL) {
  args <- list(df = as.data.frame(dat), outcome = yname, group = "gid",
               time = "wave", treatment = tname, effects = effects,
               placebo = placebo, cluster = "gid", graph_off = TRUE)
  if (!is.null(extra))     args <- modifyList(args, extra)
  set_call_seed(label)
  out <- tryCatch(do.call(DIDmultiplegtDYN::did_multiplegt_dyn, args),
                  error = function(e) NULL)
  # The command reports N per horizon but not for the aggregate, so the
  # school-waves it is given are recorded here; both estimators then report N on
  # the same basis.
  if (!is.null(out)) out$n_obs <- nrow(dat[!is.na(get(yname)) & !is.na(get(tname))])
  out
}
# Baseline covariates are time-invariant, and CDH's within-school first-difference
# design absorbs them (ΔX = 0), so CDH needs no controls.
results$cdh_comb_lang <- run_cdh(d, "y_lang", "main_treatment", "Combined -> Language")
results$cdh_comb_math <- run_cdh(d, "y_math", "main_treatment", "Combined -> Mathematics")
results$cdh_full_lang <- run_cdh(d[!is.na(main_treatment_full)], "y_lang", "main_treatment_full", "Full-day -> Language")
results$cdh_full_math <- run_cdh(d[!is.na(main_treatment_full)], "y_math", "main_treatment_full", "Full-day -> Mathematics")

# Two design choices re-estimated with the PRIMARY estimator, which is what the
# rule requires: varying a design choice and measuring it with the robustness
# estimator answers "what does CS depend on", which is not the question. Both
# were previously available only through Callaway-Sant'Anna, on the mistaken
# belief that did_multiplegt_dyn takes neither weights nor a choice of
# comparison group. It takes both (`weight`, `only_never_switchers`).
#   ENROLMENT-WEIGHTED re-targets the estimand from the average SCHOOL to the
#   average PUPIL, which is the quantity policy cares about. `weight` enters the
#   cell weights N(g,t), which appear both in the within-horizon average and in
#   the aggregation across horizons, so this re-targets the estimand and is not
#   a precision device. Two weights are run because they answer to different
#   populations: total primary enrolment, and grade-6 enrolment, which is the
#   population the outcome is measured on. Both are baseline (pre-treatment) and
#   fixed within school; a contemporaneous weight could move with treatment.
#   NEVER-SWITCHERS restricts the comparison to schools that never change
#   treatment, removing later adopters -- and the 2023 cohort with them -- from
#   the control pool. Note this changes the comparison group and therefore the
#   population the effect is identified for; it is not a free robustness check.
results$cdh_wt_lang <- run_cdh(d, "y_lang", "main_treatment", "Enrolment-weighted -> Language",
                               extra = list(weight = "enroll_primary"))
results$cdh_wt_math <- run_cdh(d, "y_math", "main_treatment", "Enrolment-weighted -> Mathematics",
                               extra = list(weight = "enroll_primary"))
results$cdh_wt6_lang <- run_cdh(d, "y_lang", "main_treatment", "Grade-6-weighted -> Language",
                                extra = list(weight = "enroll_g6"))
results$cdh_wt6_math <- run_cdh(d, "y_math", "main_treatment", "Grade-6-weighted -> Mathematics",
                                extra = list(weight = "enroll_g6"))
results$cdh_nsw_lang <- run_cdh(d, "y_lang", "main_treatment", "Never-switchers -> Language",
                                extra = list(only_never_switchers = TRUE))
results$cdh_nsw_math <- run_cdh(d, "y_math", "main_treatment", "Never-switchers -> Mathematics",
                                extra = list(only_never_switchers = TRUE))

# ============================================================================
# 2. ROB 3: Callaway-Sant'Anna, binary >=50%, doubly robust
# ============================================================================
# Same treatment as everywhere else: the contemporaneous binary, >=50% of primary
# enrolment in the year itself. What changes is the SAMPLE, not the definition --
# Callaway-Sant'Anna requires absorbing adoption, so `cs_sample()` keeps the
# never-treated plus the schools that adopt once and never revert, and drops the
# always-treated (treated in wave 1, no pre-period). Cohort G = first treated wave.
# BASELINE COVARIATES enter here: the doubly-robust estimator (xformla) conditions
# on pre-treatment X, so parallel trends need only hold given X.
setorder(d, gid, wave)
d_cs <- cs_sample(d)   # absorbing subsample + cohort G (00_setup.R)

# Doubly-robust att_gt with optional baseline covariates. `keep_rows` restricts
# to a fixed sample (e.g. composition-complete) so controlled vs uncontrolled
# fits on the same sample are directly comparable; covariate complete-cases are
# then enforced automatically.
run_cs <- function(dat, yname, label, xformla = ~1, keep_rows = NULL, wname = NULL,
                   control = "notyettreated") {
  cvars <- all.vars(xformla)
  set_call_seed(label)
  dd <- dat[!is.na(get(yname))]
  if (!is.null(keep_rows)) dd <- dd[eval(keep_rows, dd)]
  if (length(cvars)) dd <- dd[complete.cases(dd[, ..cvars])]
  set_call_seed(label)
  # The reason a fit failed is KEPT, not discarded. Returning NULL made a
  # specification that cannot be estimated indistinguishable from one nobody ran,
  # and the table below then wrote a row of blank cells with nothing to say why.
  fit <- tryCatch(
    att_gt(yname = yname, tname = "wave", idname = "gid", gname = "G",
           xformla = xformla, data = dd, control_group = control,
           weightsname = wname, allow_unbalanced_panel = TRUE,
           bstrap = TRUE, cband = FALSE, biters = BITERS, est_method = "dr"),
    error = function(e) list(cs_error = conditionMessage(e)))
  if (!is.null(fit$cs_error)) return(list(error = fit$cs_error))
  list(fit = fit, simple = suppressWarnings(aggte(fit, type = "simple", na.rm = TRUE)),
       n_schools = uniqueN(dd$gid), n_obs = nrow(dd))
}

# For each subject: (i) uncontrolled on the full absorbing sample; (ii) structural
# covariates on the full sample; (iii) uncontrolled on the composition-complete
# sample [isolates the sample change]; (iv) MAIN: + composition on that sample.
comp_cc <- quote(bl_comp_complete == TRUE)
results$cs_comb_lang         <- run_cs(d_cs, "y_lang", "Combined -> Language")
results$cs_comb_math         <- run_cs(d_cs, "y_math", "Combined -> Mathematics")
results$cs_comb_lang_struct  <- run_cs(d_cs, "y_lang", "Combined -> Language",     XF_STRUCT)
results$cs_comb_math_struct  <- run_cs(d_cs, "y_math", "Combined -> Mathematics",  XF_STRUCT)
results$cs_comb_lang_ncc     <- run_cs(d_cs, "y_lang", "Combined -> Language",     ~1,      comp_cc)
results$cs_comb_math_ncc     <- run_cs(d_cs, "y_math", "Combined -> Mathematics",  ~1,      comp_cc)
results$cs_comb_lang_comp    <- run_cs(d_cs, "y_lang", "Combined -> Language",     XF_COMP, comp_cc)
results$cs_comb_math_comp    <- run_cs(d_cs, "y_math", "Combined -> Mathematics",  XF_COMP, comp_cc)

# ============================================================================
# ROB 1: de Chaisemartin-D'Haultfoeuille, CONTINUOUS treatment (dose-response)
# ============================================================================
# ROB 0a: never-treated comparison group
# ============================================================================
# The headline uses not-yet-treated, which is more efficient but puts later
# adopters in the comparison group. That is exactly where the 2023 cohort does
# its damage: its untreated stretch runs 2021->2023, across the post-pandemic
# recovery. Restricting the comparison to schools that NEVER adopt removes them
# from the control pool, so this is a falsifiable prediction of the diagnosis in
# 11c: if the cohort is the problem, the disagreement should attenuate here.
# REMOVED 07/08: this check now runs on the PRIMARY estimator instead
# (`only_never_switchers` above). Keeping both answered the same question twice,
# once on the estimator that does not carry the headline.

# ============================================================================
# ROB 0b: what the SES covariate contributes
# ============================================================================
# Mother's secondary education is now IN the main covariate set (XF_COMP, see
# 00_setup.R): social vulnerability was an explicit school-selection criterion
# of the provincial programmes, so SES is a confounder rather than a precision
# covariate. This spec refits on the PREVIOUS set (no SES) so its contribution
# is reported rather than assumed.
#
# IT CANNOT BE FITTED, and the reason is the one already recorded for the 75%
# cut-off further down: att_gt fails with a singular design matrix in the
# (g,t) = (4,5) cell, the cohort whose first treated wave is the 2022 SAMPLE
# wave. Diagnosed 13/08; `faster_mode = FALSE`, which the package's own error
# suggests, does not recover it either. The irony is worth stating because it
# looks backwards: the LARGER set (XF_COMP) succeeds where this smaller one
# fails, because requiring `bl_mom_pre` drops the schools with no pre-treatment
# wave and with them enough of the degenerate cell.
#
# The rows are still written, now carrying that reason in `status` rather than
# four blank cells. What the SES covariate contributes is therefore not reported
# from here; 06h_ses_covariate.csv holds what can be estimated.
results$cs_noses_lang <- run_cs(d_cs, "y_lang", "No SES covariate -> Language",    XF_COMP_NOSES, comp_cc)
results$cs_noses_math <- run_cs(d_cs, "y_math", "No SES covariate -> Mathematics", XF_COMP_NOSES, comp_cc)

# ============================================================================
# REMOVED 11/08: continuous (dose-response) specification
# ============================================================================
# The continuous share is no longer estimated. It answered a different question
# from the rest of the grid -- the effect of INCREASING exposure rather than of
# crossing a threshold -- and it was the only specification whose standard errors
# came from a bootstrap the package documents as not backed by a proven
# asymptotic normality result. The reason the treatment is not modelled as
# continuous is argued in the manuscript rather than shown as a row.

# ============================================================================
# NOTE: threshold sensitivity moved out of this script (04_threshold_sensitivity.R)
# ============================================================================
# Alternative binary cut-offs are now re-estimated with the PRIMARY estimator
# (CDH) rather than with CS. Two reasons: the claim should be about the
# estimator the results rest on, and the CS grid could not be completed -- at
# the 75% cut-off att_gt fails with a singular design matrix in the (g,t)=(4,5)
# cell, a 36-school cohort whose first treated wave is the 2022 sample wave.

# ============================================================================
# ROB 4 (REMOVED 11/08): quasi value-added — condition on the 2016 baseline score
# ============================================================================
# Dropped from the dissertation, not merely unreported. Conditioning on a lagged
# outcome does not check the identifying assumption of a levels DiD, it replaces
# it, so the row was answering a different question from the grid it sat in. It
# also could not be estimated: it runs only on the doubly-robust estimator, whose
# 2022 cohort holds fewer than twenty schools per wave and yields a singular
# design matrix once the covariate set changes.

# ============================================================================
# ROB 5: enrollment-weighted CS (student-average estimand)
# ============================================================================
# The MAIN CS spec weights each school equally (school-average: the effect on the
# typical SCHOOL, coherent with the school-level ITT). Weighting by baseline
# primary enrolment re-targets the estimand to the typical STUDENT (larger schools
# count more) -- a decision-relevant robustness. Weights are baseline enrolment
# (time-invariant, pre-treatment); att_gt's `weightsname` propagates them through
# the influence functions, so the re-weighting is a clean change of estimand.
# The CDH counterpart is estimated above: `weight` enters the cell weights N(g,t),
# which appear both in the within-horizon average and in the aggregation across
# horizons, so it re-targets the estimand rather than only re-weighting within a
# horizon. The two weighted numbers are still not the same estimand -- different
# populations and dynamic counterfactuals.
d_cs[, w_enroll := pmax(enroll_primary, 1)]
results$cs_comb_lang_wtd <- run_cs(d_cs, "y_lang", "Combined -> Language (enrol-wtd)",    XF_COMP, comp_cc, wname = "w_enroll")
results$cs_comb_math_wtd <- run_cs(d_cs, "y_math", "Combined -> Mathematics (enrol-wtd)", XF_COMP, comp_cc, wname = "w_enroll")
# Same weight as the CDH counterpart, so the 2x2 in the working document compares
# estimators at a fixed weight instead of mixing estimator with weight definition.
d_cs[, w_g6 := pmax(enroll_g6, 1)]
results$cs_comb_lang_wt6 <- run_cs(d_cs, "y_lang", "Combined -> Language (g6-wtd)",    XF_COMP, comp_cc, wname = "w_g6")
results$cs_comb_math_wt6 <- run_cs(d_cs, "y_math", "Combined -> Mathematics (g6-wtd)", XF_COMP, comp_cc, wname = "w_g6")

# ============================================================================
# 3. Benchmark: static TWFE (Nistal & Edo specification)
# ============================================================================
tw_l <- feols(y_lang ~ main_treatment | gid + wave, data = d, cluster = ~gid)
tw_m <- feols(y_math ~ main_treatment | gid + wave, data = d, cluster = ~gid)

# ============================================================================
# 4. Summary table
# ============================================================================
grab_cdh <- function(r) if (is.null(r)) c(NA_real_, NA_real_) else r$results$ATE[1:2]
# `status` is empty on every row that carries an estimate and states the reason on
# every row that does not, so a blank att/se in the committed table is never left
# for the reader to interpret. Two rows use it: the no-SES covariate set cannot be
# fitted on the absorbing subsample (see the note where it is estimated).
cdh_row <- function(r, regime, subject, estimator = "CDH (headline)")
  data.table(estimator = estimator, regime = regime, subject = subject,
             covariates = "none (absorbed)", n_schools = NA_integer_,
             n_obs = if (is.null(r)) NA_integer_ else as.integer(r$n_obs),
             att = grab_cdh(r)[1], se = grab_cdh(r)[2],
             status = if (is.null(r)) "not estimated" else "")
cs_row <- function(r, subject, covariates, estimator = "CS (absorbing)", regime = "combined") {
  failed <- is.null(r) || !is.null(r$error)
  data.table(estimator = estimator, regime = regime, subject = subject,
             covariates = covariates,
             n_schools = if (failed) NA_integer_ else as.integer(r$n_schools),
             n_obs = if (failed) NA_integer_ else as.integer(r$n_obs),
             att = if (failed) NA_real_ else r$simple$overall.att,
             se  = if (failed) NA_real_ else r$simple$overall.se,
             status = if (!failed) "" else if (is.null(r)) "not estimated" else r$error)
}

summary_tab <- rbindlist(list(
  cdh_row(results$cdh_comb_lang, "combined", "language"),
  cdh_row(results$cdh_comb_math, "combined", "mathematics"),
  cdh_row(results$cdh_full_lang, "full-day", "language"),
  cdh_row(results$cdh_full_math, "full-day", "mathematics"),
  cs_row(results$cs_comb_lang,        "language",    "none"),
  cs_row(results$cs_comb_math,        "mathematics", "none"),
  cs_row(results$cs_comb_lang_struct, "language",    "structural"),
  cs_row(results$cs_comb_math_struct, "mathematics", "structural"),
  cs_row(results$cs_comb_lang_ncc,    "language",    "none (comp-sample)"),
  cs_row(results$cs_comb_math_ncc,    "mathematics", "none (comp-sample)"),
  cs_row(results$cs_comb_lang_comp,   "language",    "composition"),
  cs_row(results$cs_comb_math_comp,   "mathematics", "composition"),
  cs_row(results$cs_noses_lang,       "language",    "composition, no SES"),
  cs_row(results$cs_noses_math,       "mathematics", "composition, no SES"),
  cdh_row(results$cdh_wt_lang,  "combined", "language",    "CDH (enrol-weighted)"),
  cdh_row(results$cdh_wt_math,  "combined", "mathematics", "CDH (enrol-weighted)"),
  cdh_row(results$cdh_wt6_lang, "combined", "language",    "CDH (grade-6-weighted)"),
  cdh_row(results$cdh_wt6_math, "combined", "mathematics", "CDH (grade-6-weighted)"),
  cdh_row(results$cdh_nsw_lang, "combined", "language",    "CDH (never-switchers)"),
  cdh_row(results$cdh_nsw_math, "combined", "mathematics", "CDH (never-switchers)"),
  cs_row(results$cs_comb_lang_wtd, "language",    "composition", "CS (enrol-weighted)", "combined"),
  cs_row(results$cs_comb_math_wtd, "mathematics", "composition", "CS (enrol-weighted)", "combined"),
  cs_row(results$cs_comb_lang_wt6, "language",    "composition", "CS (grade-6-weighted)", "combined"),
  cs_row(results$cs_comb_math_wt6, "mathematics", "composition", "CS (grade-6-weighted)", "combined"),
  data.table(estimator="TWFE (benchmark)", regime="combined", subject="language",    covariates="none (absorbed)", n_schools=NA_integer_, n_obs=as.integer(nobs(tw_l)), att=coef(tw_l)[["main_treatment"]], se=se(tw_l)[["main_treatment"]], status=""),
  data.table(estimator="TWFE (benchmark)", regime="combined", subject="mathematics", covariates="none (absorbed)", n_schools=NA_integer_, n_obs=as.integer(nobs(tw_m)), att=coef(tw_m)[["main_treatment"]], se=se(tw_m)[["main_treatment"]], status="")
), fill = TRUE)
summary_tab[, `:=`(att = round(att, 4), se = round(se, 4))]
fwrite(summary_tab, file.path(DIR_TABLES, "06_estimation_main.csv"))

# ---- Placebos of every CDH call ---------------------------------------------
# cdh_row() extracts only the ATE, so the placebos never reached a committed
# table and could not be checked against the text. They are the diagnostic that
# decides the full-day arm, so they get their own file rather than a new column
# in the main table, which other scripts read on a fixed schema.
#
# EVERY placebo the specification returns, not the first. The extractor used to
# take row 1 while `run_cdh` asked for one placebo, so a two-horizon window --
# which admits two -- was being judged on half its pre-treatment evidence. The
# headline's SECOND placebo is the one that rejects, so the eleven other
# specifications were being cleared on the wrong half. `run_cdh` now asks for
# max_placebo() and this extractor returns the full table.
pl_row <- function(r, regime, subject, estimator) {
  pl <- if (is.null(r)) NULL else r$results$Placebos
  if (is.null(pl) || !nrow(pl))
    return(data.table(estimator = estimator, regime = regime, subject = subject,
                      placebo = "Placebo_1", att = NA_real_, se = NA_real_))
  data.table(estimator = estimator, regime = regime, subject = subject,
             placebo = paste0("Placebo_", seq_len(nrow(pl))),
             att = as.numeric(pl[, "Estimate"]), se = as.numeric(pl[, "SE"]))
}
placebo_tab <- rbindlist(list(
  pl_row(results$cdh_comb_lang,  "combined", "language",    "CDH (headline)"),
  pl_row(results$cdh_comb_math,  "combined", "mathematics", "CDH (headline)"),
  pl_row(results$cdh_full_lang,  "full-day", "language",    "CDH (headline)"),
  pl_row(results$cdh_full_math,  "full-day", "mathematics", "CDH (headline)"),
  pl_row(results$cdh_wt_lang,    "combined", "language",    "CDH (enrol-weighted)"),
  pl_row(results$cdh_wt_math,    "combined", "mathematics", "CDH (enrol-weighted)"),
  pl_row(results$cdh_wt6_lang,   "combined", "language",    "CDH (grade-6-weighted)"),
  pl_row(results$cdh_wt6_math,   "combined", "mathematics", "CDH (grade-6-weighted)"),
  pl_row(results$cdh_nsw_lang,   "combined", "language",    "CDH (never-switchers)"),
  pl_row(results$cdh_nsw_math,   "combined", "mathematics", "CDH (never-switchers)")))
fwrite(placebo_tab, file.path(DIR_TABLES, "06_placebos.csv"))
# Fitted model objects are NOT saved. They were, and nothing ever read them:
# 855 MB of caches sat in data/processed alongside the three real tables and
# made the data flow hard to follow. Every number this script produces is in
# its CSV outputs, and the registry (92) indexes them.
#
# The cost of that decision is real and was paid once: correcting the placebo
# extractor above required re-estimating everything, because no fit survived the
# run that produced it. The substitute is a rule about the CSVs rather than a
# cache -- EXTRACT EVERYTHING THE SPECIFICATION RETURNS, not the subset today's
# text happens to quote. `pl_row()` now follows it; the old one did not, which is
# why the second placebo went missing for eleven specifications.

# ============================================================================
# 5. Descriptive: treatment-intensity event study (validates the binary indicator)
# ============================================================================
# Mean expanded-day SHARE relative to the first treated RA year, pooled over
# adopting schools. Flat and low before, a sharp jump at k = 0 and a high level
# after is what the binary indicator is meant to be tracking.
adopt <- treatment_annual[main_treatment == 1, .(g = min(year)), by = school_id]
ai <- merge(ra_panel[, .(school_id, year, share_expanded)], adopt, by = "school_id")
ai[, k := year - g]
ai_agg <- ai[!is.na(share_expanded) & k >= -6 & k <= 6,
             .(mean_share = mean(share_expanded), n = .N), by = k][order(k)]
fwrite(ai_agg, file.path(DIR_TABLES, "06_adoption_intensity.csv"))
# Figure: rendered from this table (06_adoption_intensity.csv) by 91_figures.R
# in the shared style. Estimation scripts write data; figures live in 91.
