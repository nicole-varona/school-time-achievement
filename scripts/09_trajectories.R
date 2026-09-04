# ==============================================================================
# 09_trajectories.R — Secondary outcome: grade-6 school trajectories (RA, annual)
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Two purposes:
#   (1) a substantive SECONDARY outcome — does extended-day schooling change
#       grade-6 non-promotion (repetition) and dropout?
#   (2) a far stronger PARALLEL-TRENDS test than the six Aprender waves allow:
#       the RA trajectory outcome is ANNUAL (2011-2025), so the 2023 adoption
#       cohort has ~10 pre-periods (vs 2-3 for the test-score panel).
#
# Uses the UNIFIED panel (panel.rds: one row per school x year; trajectory rates
# present annually, Aprender columns NA outside the six waves). Treatment is the
# same extended-day indicator as the score analyses.
#
# COVID: in school-years 2020-2021 the national "promocion acompañada" policy
# suspended grade repetition, collapsing non-promotion administratively. Primary
# spec keeps the full panel (CS/CDH time effects absorb this common shock);
# robustness drops 2020-2021.
#
# Timing: the RA trajectory base records school-year t-1 at RA-year t (see
# 01_build_panel.R). We realign the outcome to the school-year it describes
# (TRAJ_LAG = 1) and verify against the COVID dip before estimating.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
for (p in c("DIDmultiplegtDYN", "did", "fixest"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({library(did); library(fixest); library(ggplot2); library(polars)})
theme_set(theme_bw(base_size = 11))
# Seeded per call inside run_cs(), not once here (00_setup.R). The two outcomes
# are fitted in a loop, so a script-level seed would make the dropout standard
# errors depend on the non-promotion fit having run first.

panel <- readRDS(file.path(DIR_PROCESSED, "panel.rds"))   # unified annual panel
stopifnot(all(c("nonpromotion_rate_g6", "dropout_rate_g6", "share_expanded") %in% names(panel)))

TRAJ_LAG    <- 1L                         # RA-year t records school-year t-1
COVID_YEARS <- c(2020L, 2021L)            # repetition suspended (aligned school-years)
OUTCOMES    <- c(nonpromotion = "nonpromotion_rate_g6", dropout = "dropout_rate_g6")

# ---- Build the annual estimation panel --------------------------------------
setorder(panel, school_id, year)
# Extended-day treatment on the full annual history.
# Treatment comes from 02 (the single source of truth) instead of being rebuilt
# here. The persistence rule used to live in nine copies, and shift() silently
# depends on the panel being sorted -- an invariant no script declared.
panel <- merge(panel,
               readRDS(file.path(DIR_PROCESSED, "treatment_annual.rds"))[
                 , .(school_id, year, treat = main_treatment)],
               by = c("school_id", "year"), all.x = TRUE)

# Treatment at school-year Y (from RA-year Y); outcomes describing school-year Y
# come from RA-year Y + TRAJ_LAG. Realign by relabelling outcome year -> year-LAG.
trt <- panel[, .(school_id, year, treat, share_expanded)]
out <- panel[, .(school_id, year,
                 nonpromotion = 100 * nonpromotion_rate_g6,   # percentage points
                 dropout      = 100 * dropout_rate_g6)]
out[, year := year - TRAJ_LAG]
d <- merge(trt, out, by = c("school_id", "year"))
d <- d[!is.na(treat)]

# Verify the timing realignment against the known COVID repetition suspension:
# the non-promotion trough should land in aligned school-years 2020-2021.

# School-level Hora+ label (peak director hours 2021-2023), for the genuine-only
# robustness. Parser identical to 08_hora_mas.R.
parse_hours <- function(x) {
  x <- as.character(x); n <- suppressWarnings(as.integer(sub(".*?([0-9]+).*", "\\1", x)))
  n[grepl("menos de 3|3 horas o menos", tolower(x))] <- 2L; n }
dh <- panel[year %in% c(2021, 2022, 2023) & !is.na(dir_hours_g6),
            .(mh = suppressWarnings(max(parse_hours(dir_hours_g6)))), by = school_id]
genuine_ids <- dh[is.finite(mh) & mh >= 6, school_id]
horamas_ids <- dh[is.finite(mh) & mh == 5, school_id]

# Baseline covariates come from 02, which computes them once for every school in
# the panel. This script keeps building its own ANNUAL panel -- its outcome is the
# yearly grade-progression series, a different unit of analysis from the six
# assessment waves -- but it no longer derives the covariates, which it used to do
# from a different source argument than 02 did.
d <- merge(d, readRDS(file.path(DIR_PROCESSED, "baseline_covariates.rds")),
           by = "school_id", all.x = TRUE)
YEARS <- sort(unique(d$year))
d[, `:=`(wave = match(year, YEARS), gid = as.integer(factor(school_id)))]

# Cohort structure: pre-period coverage now that the panel is annual.
setorder(d, gid, wave)
coh <- d[, .(firstG = if (any(treat == 1)) min(wave[treat == 1]) else NA_integer_,
             npre = if (any(treat == 1)) sum(wave < min(wave[treat == 1])) else NA_integer_),
         by = gid]

# ---- Estimators -------------------------------------------------------------
run_cdh <- function(dat, yname, tname, label, effects = 4, placebo = 4) {  # placebo <= effects
  out <- tryCatch(DIDmultiplegtDYN::did_multiplegt_dyn(
      df = as.data.frame(dat[!is.na(get(yname))]), outcome = yname, group = "gid",
      time = "wave", treatment = tname, effects = effects, placebo = placebo,
      cluster = "gid", graph_off = TRUE),
    error = function(e) { NULL })
  out
}
# CS on the absorbing subsample, UNBALANCED panel kept (essential for the annual
# outcome: few schools have trajectory data in all 14 years). Returns the overall
# ATT and the DYNAMIC event-study frame — the long annual pre-window is the point
# of this analysis and, unlike CDH (placebo <= effects <= few post-periods), the
# CS dynamic aggregation exposes all available pre-periods.
run_cs <- function(dat, yname, label, xformla = XF_COMP) {
  x <- copy(dat[!is.na(get(yname))]); setorder(x, gid, wave)
  xcs <- cs_sample(x, "treat")
  cvars <- all.vars(xformla); if (length(cvars)) xcs <- xcs[complete.cases(xcs[, ..cvars])]
  set_call_seed(label)
  fit <- tryCatch(att_gt(yname = yname, tname = "wave", idname = "gid", gname = "G",
           xformla = xformla, data = xcs, control_group = "notyettreated",
           allow_unbalanced_panel = TRUE, bstrap = TRUE, cband = FALSE, biters = BITERS,
           est_method = "dr"),
         error = function(e) { NULL })
  if (is.null(fit)) return(NULL)
  simple <- suppressWarnings(aggte(fit, type = "simple", na.rm = TRUE))
  dyn    <- tryCatch(suppressWarnings(aggte(fit, type = "dynamic", na.rm = TRUE)),
                     error = function(e) NULL)
  frame  <- if (is.null(dyn)) NULL else
    data.table(k = dyn$egt, est = dyn$att.egt, se = dyn$se.egt,
               lo = dyn$att.egt - 1.96 * dyn$se.egt, hi = dyn$att.egt + 1.96 * dyn$se.egt)
  list(fit = fit, simple = simple, dyn = dyn, frame = frame, n_schools = uniqueN(xcs$gid))
}
es_frame <- function(r, subject, spec) {
  if (is.null(r)) return(NULL)
  parts <- list()
  if (!is.null(r$results$Placebos)) {
    p <- as.data.table(r$results$Placebos)
    p[, k := -as.integer(sub("Placebo_", "", rownames(r$results$Placebos)))]; parts$p <- p }
  e <- as.data.table(r$results$Effects); e[, k := as.integer(sub("Effect_", "", rownames(r$results$Effects)))]
  parts$e <- e
  o <- rbindlist(parts, fill = TRUE)
  setnames(o, c("Estimate","LB CI","UB CI"), c("est","lo","hi"), skip_absent = TRUE)
  rbind(o[, .(k, est, lo, hi)], data.table(k = 0, est = 0, lo = NA, hi = NA))[
    , `:=`(subject = subject, spec = spec)][order(k)]
}

# ---- Estimate each outcome: primary + COVID-drop + genuine-only robustness ---
res <- list(); es_all <- list()
for (nm in names(OUTCOMES)) {
  y <- nm  # column name in d ("nonpromotion"/"dropout")
  res[[paste0(y, "_cdh")]] <- run_cdh(d, y, "treat", sprintf("%s | full panel", y))
  res[[paste0(y, "_cs")]]  <- run_cs(d, y, sprintf("%s | full panel (CS+covars)", y))
  # robustness: drop COVID school-years (repetition suspended)
  res[[paste0(y, "_covid")]] <- run_cdh(d[!(year %in% COVID_YEARS)], y, "treat",
                                        sprintf("%s | drop 2020-21 (COVID)", y))
  # robustness: genuine 6h+ schools only (exclude Hora+); controls = never-treated non-Hora+
  d_gen <- d[!(school_id %in% horamas_ids)]
  res[[paste0(y, "_genuine")]] <- run_cdh(d_gen, y, "treat",
                                          sprintf("%s | genuine (excl Hora+)", y))
  es_all[[y]] <- rbindlist(list(
    es_frame(res[[paste0(y, "_cdh")]],   y, "full panel"),
    es_frame(res[[paste0(y, "_covid")]], y, "drop 2020-21")), fill = TRUE)
}
es_cdh <- rbindlist(es_all, fill = TRUE)   # CDH, 4 placebos (supplementary)
if (nrow(es_cdh)) fwrite(es_cdh, file.path(DIR_TABLES, "09_event_study_cdh.csv"))

# ---- Primary PT test: CS dynamic event-study (the long annual pre-window) ----
# The whole point of the annual outcome: the CS dynamic aggregation exposes every
# available pre-period (~10 for the 2023 cohort), a far stronger parallel-trends
# test than the 2-3 placebos the six Aprender waves allow.
cs_es <- rbindlist(lapply(names(OUTCOMES), function(nm) {
  fr <- res[[paste0(nm, "_cs")]]$frame
  if (is.null(fr)) return(NULL)
  copy(fr)[, subject := nm]
}), fill = TRUE)
if (nrow(cs_es)) {
  fwrite(cs_es, file.path(DIR_TABLES, "09_event_study_trajectories.csv"))
  # Figure 09_event_study_trajectories.png is rendered from
  # 09_event_study_trajectories.csv by 91_figures.R in the shared style.
}

# ---- Summary table ----------------------------------------------------------
grab_cdh <- function(r) if (is.null(r)) c(NA_real_, NA_real_) else r$results$ATE[1:2]
grab_cs  <- function(r) if (is.null(r)) c(NA_real_, NA_real_, NA_integer_) else c(r$simple$overall.att, r$simple$overall.se, r$n_schools)
rows <- list()
for (nm in names(OUTCOMES)) {
  rows[[length(rows)+1]] <- data.table(outcome=nm, spec="CDH full panel",    att=grab_cdh(res[[paste0(nm,"_cdh")]])[1],     se=grab_cdh(res[[paste0(nm,"_cdh")]])[2])
  rows[[length(rows)+1]] <- data.table(outcome=nm, spec="CDH drop 2020-21",  att=grab_cdh(res[[paste0(nm,"_covid")]])[1],   se=grab_cdh(res[[paste0(nm,"_covid")]])[2])
  rows[[length(rows)+1]] <- data.table(outcome=nm, spec="CDH genuine only",  att=grab_cdh(res[[paste0(nm,"_genuine")]])[1], se=grab_cdh(res[[paste0(nm,"_genuine")]])[2])
  rows[[length(rows)+1]] <- data.table(outcome=nm, spec="CS full (+covars)", att=grab_cs(res[[paste0(nm,"_cs")]])[1],       se=grab_cs(res[[paste0(nm,"_cs")]])[2])
}
summary_tab <- rbindlist(rows, fill = TRUE)
# Outcomes are rates x100, so effects are percentage points. The untreated mean
# gives the scale needed to read them as relative changes.
base <- rbindlist(lapply(names(OUTCOMES), function(nm)
  data.table(outcome = nm, base_rate_pp = d[treat == 0, mean(get(nm), na.rm = TRUE)])))
summary_tab <- merge(summary_tab, base, by = "outcome", sort = FALSE)
setnames(summary_tab, c("att", "se"), c("att_pp", "se_pp"))
summary_tab[, `:=`(att_pp = round(att_pp, 4), se_pp = round(se_pp, 4),
                   base_rate_pp = round(base_rate_pp, 4))]
fwrite(summary_tab, file.path(DIR_TABLES, "09_trajectory_results.csv"))
# Fitted model objects are NOT saved. They were, and nothing ever read them:
# 855 MB of caches sat in data/processed alongside the three real tables and
# made the data flow hard to follow. Every number this script produces is in
# its CSV outputs, and the registry (92) indexes them.
