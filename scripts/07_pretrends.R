# ==============================================================================
# 07_pretrends.R — Parallel-trends evidence: event-study + Honest-DiD
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# (a) Event-study from de Chaisemartin-D'Haultfoeuille (dynamic effects +
#     placebos): the estimator's own tool for assessing parallel trends and
#     no-anticipation JOINTLY. TWO is the design's ceiling, not a choice: the
#     command caps placebos at the number of dynamic effects requested, and with
#     three effects the data support two. Only cohorts adopting
#     after 2016 identify pre-trends (57% of treated schools have >=1 pre-wave).
# (b) Honest-DiD (Rambachan & Roth 2023) sensitivity on the Callaway-Sant'Anna
#     event-study (the standard bridge): how large a pre-trend violation would
#     have to be to overturn the result. This is the tool for a short pre-window.
# Inputs  : data/processed/panel.rds, ra_panel.rds
# Outputs : output/figures/07_event_study_cdh.png
#           output/figures/07_honestdid_*.png ; output/tables/07_*.csv
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
for (p in c("DIDmultiplegtDYN", "did", "HonestDiD", "ggplot2"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({library(did); library(HonestDiD); library(ggplot2); library(polars)})
# Event-study figures are rendered by 91_figures.R from the saved CSVs. Only the
# Honest-DiD sensitivity plot is drawn here, because it needs the fitted model
# object (the original confidence set), which cannot be rebuilt from a table. It
# uses the shared style so it matches the rest.
source(here::here("scripts", "90_style.R"))
# NO SEED IS SET HERE, and none is needed: nothing in this script draws on the
# random number generator. `did_multiplegt_dyn` is called with its analytic
# clustered variance, and `att_gt` runs with bstrap = FALSE because HonestDiD
# needs the influence function rather than bootstrap draws. The script used to
# carry a script-level set.seed() with a note that it was "not yet on per-call
# seeding"; it was seeding nothing.

# ---- Load the estimation sample ---------------------------------------------
# Read from 02, not rebuilt. This script used to reconstruct the panel, the
# outcome scaling, the wave index, `gid` and the baseline covariates on its own,
# which meant `gid` indexed a different set of schools here than everywhere else
# and a change to the covariate set silently broke it.
# The contemporaneous composition columns are the one thing 02 does not carry
# (they vary within school), so they are merged from the annual panel.
d     <- readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds"))
panel <- readRDS(file.path(DIR_PROCESSED, "panel.rds"))
CDH_CTRL_SRC <- c("pct_female", "pct_repeater", "pct_preschool")
d <- merge(d, panel[, c("school_id", "year", CDH_CTRL_SRC), with = FALSE],
           by = c("school_id", "year"), all.x = TRUE)
CDH_CTRL <- CDH_CTRL_SRC   # contemporaneous composition

# ============================================================================
# (a) Event-study — de Chaisemartin-D'Haultfoeuille (uncontrolled headline)
# ============================================================================
run_cdh_es <- function(dat, yname, controls = NULL) {
  DIDmultiplegtDYN::did_multiplegt_dyn(
    df = as.data.frame(dat), outcome = yname, group = "gid", time = "wave",
    treatment = "main_treatment", effects = 3, placebo = max_placebo(3L, length(WAVES)),
    cluster = "gid",
    controls = controls, graph_off = TRUE)
}
cdh_lang <- tryCatch(run_cdh_es(d, "y_lang"), error = function(e) { NULL })
cdh_math <- tryCatch(run_cdh_es(d, "y_math"), error = function(e) { NULL })

# Assemble the event-study frame: placebos at negative event-time, a normalised
# reference at 0, and dynamic effects at positive event-time.
es_frame <- function(r, subject, spec = "uncontrolled") {
  if (is.null(r)) return(NULL)
  parts <- list()
  if (!is.null(r$results$Placebos)) {
    p <- as.data.table(r$results$Placebos)
    p[, k := -as.integer(sub("Placebo_", "", rownames(r$results$Placebos)))]
    parts$p <- p
  }
  e <- as.data.table(r$results$Effects)
  e[, k := as.integer(sub("Effect_", "", rownames(r$results$Effects)))]
  parts$e <- e
  out <- rbindlist(parts, fill = TRUE)
  setnames(out, c("Estimate", "LB CI", "UB CI"), c("est", "lo", "hi"), skip_absent = TRUE)
  out <- rbind(out[, .(k, est, lo, hi)],
               data.table(k = 0, est = 0, lo = NA, hi = NA))   # normalised reference
  out[, `:=`(subject = subject, spec = spec)][order(k)]
}
es <- rbindlist(list(es_frame(cdh_lang, "Language"), es_frame(cdh_math, "Mathematics")))
if (nrow(es)) {
  fwrite(es, file.path(DIR_TABLES, "07_event_study_cdh.csv"))
  # Figure 07_event_study_cdh.png -> 91_figures.R (from 07_event_study_cdh.csv).

  # The k = -1 and k = -2 leads ARE the two placebos of the main specification --
  # the same estimates, not two diagnostics that agree -- and the manuscript says
  # so in print. Asserted here, where the second of the two artefacts is produced,
  # because 06g writes the first and runs earlier: checking it there compared a
  # fresh number against the previous run's file and halted the pipeline.
  mpf <- file.path(DIR_TABLES, "06g_max_placebos.csv")
  if (file.exists(mpf)) {
    lead <- es[k %in% c(-1L, -2L), .(subject = tolower(subject),
                                     parameter = paste0("Placebo_", -k),
                                     lead = round(est, 4))]
    cmp <- merge(lead, fread(mpf)[regime == "combined" & parameter %like% "Placebo",
                                  .(subject, parameter, estimate)],
                 by = c("subject", "parameter"))
    if (nrow(cmp) != 4L || !isTRUE(all.equal(cmp$lead, cmp$estimate)))
      stop("the event-study leads no longer equal the placebos in ",
           "06g_max_placebos.csv; the manuscript states that they are the same ",
           "estimates. Got ", paste(cmp$lead, collapse = "/"), " against ",
           paste(cmp$estimate, collapse = "/"))
  }
}

# ---- (a2) Does conditioning on compositional trends shrink the pre-trend? -----
# On the composition-complete sample, compare the uncontrolled event-study with
# one that adds time-varying composition CONTROLS. In pre-periods these controls
# are pre-treatment (treatment has not switched on yet) -> legitimate, not bad
# controls. Focus: the mathematics placebo at k=-2 (the flagged pre-trend).
d_cc <- d[complete.cases(d[, ..CDH_CTRL])]
es_cmp <- rbindlist(list(
  es_frame(tryCatch(run_cdh_es(d_cc, "y_math"),          error=function(e) NULL), "Mathematics", "uncontrolled"),
  es_frame(tryCatch(run_cdh_es(d_cc, "y_math", CDH_CTRL),error=function(e) NULL), "Mathematics", "composition-controlled"),
  es_frame(tryCatch(run_cdh_es(d_cc, "y_lang"),          error=function(e) NULL), "Language",    "uncontrolled"),
  es_frame(tryCatch(run_cdh_es(d_cc, "y_lang", CDH_CTRL),error=function(e) NULL), "Language",    "composition-controlled")
), fill = TRUE)
if (nrow(es_cmp)) {
  fwrite(es_cmp, file.path(DIR_TABLES, "07_event_study_cdh_controls.csv"))
  # Figure 07_event_study_cdh_controls.png -> 91_figures.R (from the CSV above).
}

# ============================================================================
# (b) CS event-study + Honest-DiD, uncontrolled vs baseline-covariate-adjusted
# ============================================================================
# The core conditional pre-trend test: the doubly-robust CS event-study (with
# the notyettreated comparison) is re-estimated conditioning on BASELINE
# covariates (xformla). If the pre-period placebos move toward zero once we
# condition on X, baseline differences explain (part of) the pre-trend, and the
# conditional parallel-trends assumption is more credible. Honest-DiD then
# quantifies how far the post-treatment result survives a residual violation.
setorder(d, gid, wave)
d_cs <- cs_sample(d)

# CS event-study frame from a dynamic aggregation (pre-periods = placebos).
cs_es_frame <- function(dyn, subject, spec) {
  if (is.null(dyn)) return(NULL)
  data.table(k = dyn$egt, est = dyn$att.egt, se = dyn$se.egt,
             lo = dyn$att.egt - 1.96 * dyn$se.egt, hi = dyn$att.egt + 1.96 * dyn$se.egt,
             subject = subject, spec = spec)
}

honest_for <- function(yname, subject, xformla = ~1, spec = "uncontrolled") {
  cvars <- all.vars(xformla)
  dd <- d_cs[!is.na(get(yname))]
  if (length(cvars)) dd <- dd[complete.cases(dd[, ..cvars])]
  # bstrap=FALSE so the influence function is available; this did version does
  # not populate V_analytical, so we build the event-study vcov from the
  # influence function (the standard HonestDiD-CS bridge).
  fit <- tryCatch(att_gt(yname = yname, tname = "wave", idname = "gid", gname = "G",
                         xformla = xformla, data = dd, control_group = "notyettreated",
                         allow_unbalanced_panel = TRUE, bstrap = FALSE, est_method = "dr", base_period = "universal"),
                  error = function(e) { NULL })
  if (is.null(fit)) return(NULL)
  dyn <- tryCatch(suppressWarnings(aggte(fit, type = "dynamic", na.rm = TRUE)), error = function(e) NULL)
  if (is.null(dyn)) return(NULL)
  frame <- cs_es_frame(dyn, subject, spec)
  inf <- dyn$inf.function$dynamic.inf.func.e
  V   <- t(inf) %*% inf / nrow(inf)^2
  npre <- sum(dyn$egt < 0); npost <- sum(dyn$egt >= 0)
  if (npre < 1 || npost < 1) { return(list(dyn = dyn, frame = frame, sens = NULL)) }
  sens <- tryCatch(HonestDiD::createSensitivityResults_relativeMagnitudes(
            betahat = dyn$att.egt, sigma = V,
            numPrePeriods = npre, numPostPeriods = npost, Mbarvec = c(0.5, 1, 1.5, 2)),
          error = function(e) { NULL })
  if (!is.null(sens)) {
    fwrite(as.data.table(sens), file.path(DIR_TABLES,
           sprintf("07_honestdid_%s_%s.csv", tolower(subject), gsub("[^a-z]", "", spec))))
    orig <- HonestDiD::constructOriginalCS(betahat = dyn$att.egt, sigma = V,
                                           numPrePeriods = npre, numPostPeriods = npost)
    pl <- HonestDiD::createSensitivityPlot_relativeMagnitudes(sens, orig) +
      labs(title = sprintf("Honest-DiD sensitivity: %s (%s)", subject, spec)) +
      # Override the package's red/cyan with the thesis palette (blue/orange).
      scale_colour_manual(values = PAL_CAT) +
      theme_thesis()
    save_fig(pl, sprintf("07_honestdid_%s_%s.png", tolower(subject), gsub("[^a-z]", "", spec)),
             width = 7, height = 5)
  }
  list(dyn = dyn, frame = frame, sens = sens)
}

h_lang_u <- honest_for("y_lang", "Language",    ~1,      "uncontrolled")
h_lang_c <- honest_for("y_lang", "Language",    XF_COMP, "baseline-covariates")
h_math_u <- honest_for("y_math", "Mathematics", ~1,      "uncontrolled")
h_math_c <- honest_for("y_math", "Mathematics", XF_COMP, "baseline-covariates")

# ---- CS event-study: uncontrolled vs baseline-covariate-adjusted -------------
cs_es <- rbindlist(lapply(list(h_lang_u, h_lang_c, h_math_u, h_math_c),
                          function(h) if (is.null(h)) NULL else h$frame), fill = TRUE)
if (nrow(cs_es)) {
  fwrite(cs_es, file.path(DIR_TABLES, "07_event_study_cs_covariates.csv"))
  # Figure 07_event_study_cs_covariates.png -> 91_figures.R (from the CSV above).
}

