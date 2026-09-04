# ==============================================================================
# 08_hora_mas.R — Separating genuine extended day from the "Hora+ / 1h+" programme
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# The Relevamiento Anual records the federal "Hora+/1h+" programme (one hour
# added to a simple 4h day -> 5h) as "extended day", conflating it with a
# genuine 6h+ extension. The Aprender director questionnaire reports grade-6
# daily hours (censal in 2021/2022/2023: items dp23/dp18/dp47), which lets us
# separate the two. Diagnostic finding: 73% of the schools that ADOPT in 2023
# (the flow that identifies the effect) report only 5h -> Hora+, not genuine
# extension, against 5% of the pre-2023 stock. (2,282 of 3,114 in the flow;
# read off 08_contamination_2023.csv, which this script writes. The figure said
# ~61% until 13/08, from an earlier treatment definition.) This script:
#   (1) parses director hours and classifies each school-wave into
#       none / horamas (5h) / genuine (6h+);
#   (2) documents the contamination (stock vs flow; 2024 explicit-flag check);
#   (3) estimates a MULTI-ARM design (Decision: "Option D"):
#         - GENUINE arm: genuine 6h+ vs untreated, EXCLUDING ever-Hora+ schools;
#         - HORA+  arm: 5h vs untreated, EXCLUDING ever-genuine schools;
#   (4) re-checks the mathematics pre-trend on the genuine arm (does removing
#       Hora+ clean the k=-2 placebo flagged in 07_pretrends.R?).
#
# Hours are observed only in 2021/2022/2023. 2016/2018 predate Hora+ (launched
# 2022) -> RA-extended is genuine there. 2025 has no hours item -> RA fallback,
# flagged, with a drop-2025 robustness.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
for (p in c("DIDmultiplegtDYN", "did", "fixest"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({library(did); library(fixest); library(ggplot2); library(polars)})
theme_set(theme_bw(base_size = 11))
# Seeded per call inside run_cs(), not once here (00_setup.R). The `label`
# argument the estimators already take is what identifies each stream; it used to
# be passed and then ignored.

panel    <- readRDS(file.path(DIR_PROCESSED, "panel.rds"))
ra_panel <- readRDS(file.path(DIR_PROCESSED, "ra_panel.rds"))
HOURS_WAVES  <- c(2021L, 2022L, 2023L)     # director grade-6 hours available
PRE_HM_WAVES <- c(2016L, 2018L)            # predate Hora+ (2022) -> RA = genuine

# ---- (1) Classify director-reported grade-6 daily hours ----------------------
# `parse_hours` is in 00_setup.R (08b classifies on the same item). Bins:
#   <=4h -> simple (untreated by hours) | 5h -> Hora+ | >=6h -> genuine extension.
hours_class <- function(h)
  fifelse(is.na(h), NA_character_,
   fifelse(h <= 4, "none", fifelse(h == 5, "horamas", "genuine")))

dh <- panel[year %in% HOURS_WAVES & !is.na(dir_hours_g6),
            .(school_id, year, hours = parse_hours(dir_hours_g6))]
dh[, hclass := hours_class(hours)]

# ---- (2) Contamination diagnostic --------------------------------------------

# RA extended treatment: share >= 0.50 in the year itself.
# Treatment comes from 02 (the single source of truth) instead of being rebuilt
# here. The persistence rule used to live in nine copies, and shift() silently
# depends on the panel being sorted -- an invariant no script declared.
ra_panel <- merge(ra_panel,
               readRDS(file.path(DIR_PROCESSED, "treatment_annual.rds"))[
                 , .(school_id, year, main_treatment)],
               by = c("school_id", "year"), all.x = TRUE)

# Adoption cohort on the primary treatment indicator (over outcome waves) for stock/flow.
tw      <- ra_panel[year %in% WAVES, .(school_id, year, main_treatment)]
first_g <- tw[main_treatment == 1, .(firstwave = min(year)), by = school_id]

p23 <- merge(panel[year == 2023, .(school_id, share_expanded)],
             dh[year == 2023, .(school_id, hclass)], by = "school_id")
p23 <- merge(p23, first_g, by = "school_id", all.x = TRUE)
p23[, ra_treated := as.integer(!is.na(share_expanded) & share_expanded >= 0.5)]
p23[, when := fifelse(is.na(firstwave), "never treated in a wave",
              fifelse(firstwave == 2023, "adopted 2023 (flow)", "adopted <2023 (stock)"))]
contam <- dcast(p23[ra_treated == 1 & !is.na(hclass) & when != "never treated in a wave"],
                when ~ hclass, value.var = "school_id", fun.aggregate = length)
fwrite(contam, file.path(DIR_TABLES, "08_contamination_2023.csv"))

# Explicit-flag validation (2024 sample, dp21): does self-reported Hora+ line up
# with 5h director hours (2023) and RA treatment?
# REMOVED 13/08: the 2024 explicit-flag cross-check was computed into `val` and
# then never read -- no output, no assertion, no number in the manuscript. The
# item it validated against (dp21, the 2024 sample) is described in the data
# dictionary; if the check is wanted back it needs an output to land in.

# ---- (3) Build arm treatment states (school-wave) ----------------------------
# The estimation panel comes from 02 (outcome scaling, wave index, `gid`,
# treatment definitions and baseline covariates all included); only the director-
# hours classification, which is specific to this script, is merged on.
d <- readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds"))
d <- merge(d, dh[, .(school_id, year, hours, hclass)], by = c("school_id", "year"), all.x = TRUE)

# Two treatment definitions:
# (a) PER-WAVE (hours-refined): each wave classified by director hours where
#     observed; RA fallback in pre-Hora+ (2016/18) and unverifiable 2025 waves.
#     Kept only as a SENSITIVITY: hours flip 5<->6+ year to year, so treatment
#     switches on/off spuriously and induces artefactual event-study dynamics
#     (a spurious negative k=-1; see 08_event_study_genuine.png).
d[, state := NA_character_]
d[year %in% HOURS_WAVES & !is.na(hclass),               state := hclass]
d[year %in% HOURS_WAVES & is.na(hclass),                state := fifelse(main_treatment == 1, "genuine", "none")]
d[year %in% PRE_HM_WAVES,                               state := fifelse(main_treatment == 1, "genuine", "none")]
d[year == 2025L,                                        state := fifelse(main_treatment == 1, "genuine", "none")]
d[, `:=`(treat_gen_pw = as.integer(state == "genuine"),
         treat_hm_pw  = as.integer(state == "horamas"),
         hours_unverified = as.integer(year == 2025L))]

# (b) STABLE (school-level, PRIMARY): label each school time-invariantly by its
#     PEAK grade-6 hours in 2021-2023, then take adoption TIMING from the stable
#     RA treatment indicator (main_treatment). This removes the year-to-year switching
#     while still separating genuine (>=6h) extenders from Hora+ (5h) schools.
sch_label <- dh[, .(max_h = suppressWarnings(max(hours, na.rm = TRUE))), by = school_id]
sch_label[!is.finite(max_h), max_h := NA_integer_]
sch_label[, sch_class := fifelse(is.na(max_h), "unobs",
                         fifelse(max_h >= 6, "genuine",
                         fifelse(max_h == 5, "horamas", "simple")))]
d <- merge(d, sch_label[, .(school_id, sch_class)], by = "school_id", all.x = TRUE)
d[is.na(sch_class), sch_class := "unobs"]
d[, ever_ra := any(main_treatment == 1), by = school_id]
d[, treat_gen := as.integer(sch_class == "genuine" & main_treatment == 1)]  # RA timing, verified genuine
d[, treat_hm  := as.integer(sch_class == "horamas" & main_treatment == 1)]  # RA timing, verified Hora+

# clean never-treated controls: never RA-treated and NOT a Hora+ school
ctrl_ok   <- quote(!ever_ra & sch_class %in% c("unobs", "simple"))

# ---- Estimators --------------------------------------------------------------
run_cdh <- function(dat, yname, tname, label, effects = 2, placebo = 1) {
  out <- tryCatch(DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(dat), outcome = yname, group = "gid", time = "wave",
      treatment = tname, effects = effects, placebo = placebo, cluster = "gid",
      graph_off = TRUE),
    error = function(e) { NULL })
  # School-waves the estimator actually uses, recorded alongside the fit so the
  # arms table can report N on the same basis for both estimators.
  if (!is.null(out)) out$n_obs <- nrow(dat[!is.na(get(yname)) & !is.na(get(tname))])
  out
}
run_cs <- function(dat, yname, tname, label, xformla = XF_COMP) {
  x <- copy(dat); setorder(x, gid, wave)
  xcs <- cs_sample(x, tname)
  cvars <- all.vars(xformla); xcs <- xcs[!is.na(get(yname))]
  if (length(cvars)) xcs <- xcs[complete.cases(xcs[, ..cvars])]
  set_call_seed(label)
  fit <- tryCatch(att_gt(yname = yname, tname = "wave", idname = "gid", gname = "G",
           xformla = xformla, data = xcs, control_group = "notyettreated",
           allow_unbalanced_panel = TRUE, bstrap = TRUE, cband = FALSE, biters = BITERS, est_method = "dr"),
         error = function(e) { NULL })
  if (is.null(fit)) return(NULL)
  agg <- suppressWarnings(aggte(fit, type = "simple", na.rm = TRUE))
  list(fit = fit, simple = agg, n_schools = uniqueN(xcs$gid), n_obs = nrow(xcs))
}

# ---- (3a) GENUINE arm (PRIMARY, stable): genuine 6h+ schools (RA timing) vs
#          clean never-RA-treated controls (Hora+ schools excluded) ------------
d_gen <- d[(sch_class == "genuine" & ever_ra) | eval(ctrl_ok)]
res <- list()
res$gen_cdh_lang <- run_cdh(d_gen, "y_lang", "treat_gen", "Genuine -> Language")
res$gen_cdh_math <- run_cdh(d_gen, "y_math", "treat_gen", "Genuine -> Mathematics")
res$gen_cs_lang  <- run_cs (d_gen, "y_lang", "treat_gen", "Genuine -> Language (CS+covars)")
res$gen_cs_math  <- run_cs (d_gen, "y_math", "treat_gen", "Genuine -> Mathematics (CS+covars)")

# ---- (3b) HORA+ arm (stable): Hora+ 5h schools (RA timing) vs clean controls -
d_hm <- d[(sch_class == "horamas" & ever_ra) | eval(ctrl_ok)]
res$hm_cdh_lang <- run_cdh(d_hm, "y_lang", "treat_hm", "Hora+ -> Language")
res$hm_cdh_math <- run_cdh(d_hm, "y_math", "treat_hm", "Hora+ -> Mathematics")
# NOT estimated with Callaway-Sant'Anna on this arm: att_gt kills the R session
# (SIGSEGV, exit 139), which no tryCatch can intercept, so the call cannot be made
# at all rather than returning an NA row.
#
# Cause: the arm is one cohort with stragglers. Of its 1,720 treated schools in
# the CS frame, 1,304 first appear treated at wave 5; cohorts 2, 3 and 4 hold 3,
# 15 and 1. CS fits a doubly-robust regression inside each (g,t) cell, which those
# cannot support. That follows from the measurement finding itself -- three
# quarters of the 2023 adoption flow is five-hour.
#
# Two things this comment claimed until 13/08, both checked and both false: the
# arm is NOT the thinner of the two (1,720 treated against the genuine arm's
# 1,096), and structural covariates only do NOT run either -- XF_STRUCT aborts the
# same way. Dropping the three small cohorts does run (+0.006, s.e. 0.038) but
# estimates something else, so it is not reported. The arm is reported with the
# primary estimator, which handles it.

# ---- (4) Pre-trend re-check: does removing Hora+ clean the maths k=-2 placebo?
# Primary = stable definition. We also plot the PER-WAVE definition to show that
# its year-to-year switching produces a SPURIOUS negative k=-1 (an artefact, not
# a real pre-trend). Key finding: under the stable definition the maths k=-2
# placebo does NOT shrink -> the pre-trend is inherent to genuine extenders, not
# Hora+ contamination.
es_frame <- function(r, subject, spec) {
  if (is.null(r)) return(NULL)
  parts <- list()
  if (!is.null(r$results$Placebos)) {
    p <- as.data.table(r$results$Placebos)
    p[, k := -as.integer(sub("Placebo_", "", rownames(r$results$Placebos)))]; parts$p <- p }
  e <- as.data.table(r$results$Effects); e[, k := as.integer(sub("Effect_", "", rownames(r$results$Effects)))]
  parts$e <- e
  out <- rbindlist(parts, fill = TRUE)
  setnames(out, c("Estimate","LB CI","UB CI"), c("est","lo","hi"), skip_absent = TRUE)
  rbind(out[, .(k, est, lo, hi)], data.table(k = 0, est = 0, lo = NA, hi = NA))[
    , `:=`(subject = subject, spec = spec)][order(k)]
}
ever_hm_pw <- unique(d[treat_hm_pw == 1]$school_id)
d_gen_pw   <- d[!(school_id %in% ever_hm_pw)]
es_gen <- rbindlist(list(
  es_frame(run_cdh(d_gen,    "y_lang", "treat_gen",    "ES stable Lang",   3, 2), "Language",    "stable (school-level)"),
  es_frame(run_cdh(d_gen,    "y_math", "treat_gen",    "ES stable Math",   3, 2), "Mathematics", "stable (school-level)"),
  es_frame(run_cdh(d_gen_pw, "y_lang", "treat_gen_pw", "ES per-wave Lang", 3, 2), "Language",    "per-wave (hours switching)"),
  es_frame(run_cdh(d_gen_pw, "y_math", "treat_gen_pw", "ES per-wave Math", 3, 2), "Mathematics", "per-wave (hours switching)")), fill = TRUE)
if (nrow(es_gen)) {
  fwrite(es_gen, file.path(DIR_TABLES, "08_event_study_genuine.csv"))
  # Figure 08_event_study_genuine.png is rendered from 08_event_study_genuine.csv
  # by 91_figures.R in the shared style.
}

# ---- Summary table -----------------------------------------------------------
grab <- function(r, kind) {
  if (is.null(r)) return(c(NA_real_, NA_real_, NA_integer_, NA_integer_))
  if (kind == "CDH") c(r$results$ATE[1], r$results$ATE[2], NA_integer_, r$n_obs)
  else c(r$simple$overall.att, r$simple$overall.se,
         as.integer(r$n_schools), as.integer(r$n_obs))
}
arm_row <- function(r, kind, arm, subject) {
  g <- grab(r, kind)
  data.table(estimator = kind, arm = arm, subject = subject,
             att = g[1], se = g[2], n_schools = as.integer(g[3]),
             n_obs = as.integer(g[4]))
}
summary_tab <- rbindlist(list(
  arm_row(res$gen_cdh_lang, "CDH", "genuine 6h+ (stable)", "language"),
  arm_row(res$gen_cdh_math, "CDH", "genuine 6h+ (stable)", "mathematics"),
  arm_row(res$gen_cs_lang,  "CS",  "genuine 6h+ (stable)", "language"),
  arm_row(res$gen_cs_math,  "CS",  "genuine 6h+ (stable)", "mathematics"),
  arm_row(res$hm_cdh_lang,  "CDH", "1h+ 5h (stable)",    "language"),
  arm_row(res$hm_cdh_math,  "CDH", "1h+ 5h (stable)",    "mathematics")
), fill = TRUE)
summary_tab[, `:=`(att = round(att, 4), se = round(se, 4))]
fwrite(summary_tab, file.path(DIR_TABLES, "08_arm_results.csv"))
# Fitted model objects are NOT saved. They were, and nothing ever read them:
# 855 MB of caches sat in data/processed alongside the three real tables and
# made the data flow hard to follow. Every number this script produces is in
# its CSV outputs, and the registry (92) indexes them.
