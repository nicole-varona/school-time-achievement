# ==============================================================================
# 06f_clustering_sensitivity.R — Provincial structure: dependence and confounding
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Inputs  : data/processed/estimation_sample.rds
# Outputs : output/tables/06f_clustering.csv        SEs by clustering level
#           output/tables/06f_effective_clusters.csv  how concentrated they are
#           output/tables/06f_province_trends.csv     province-specific year trends
#           output/tables/06f_inference_diagnostics.csv  evidence behind the two
#                                                        inference decisions
#           output/tables/06f_resampling.csv        province-clustered resampling
#                                                   inference, where it is available
#
# Treatment is assigned mainly by the province, which raises TWO separate
# questions that are easy to conflate. Whether provincial shocks CONFOUND the
# estimate is a question about the point estimate, answered by letting each
# province have its own year effects (`trends_nonparam`). Whether they induce
# DEPENDENCE is a question about the standard error, answered by clustering.
# This script does both, in that order. The province matters here because the
# programme is national with provincial adhesion and provincial funding decisions
# switch it on and off. The PROVINCE is therefore the substantively relevant
# level of dependence: adoption is coordinated there, the funding that reverses
# it is decided there, and the measurement error in the treatment is structured
# by jurisdiction (06m). That case is made by the design, not by which level
# returns the wider interval.
#
# Two findings below are easy to over-read, so what each does and does not
# license is stated here rather than inferred from the numbers:
#   (1) clustering at the province makes the PRIMARY estimator's SEs SMALLER, at
#       every horizon and in both subjects, and smaller still at region. This is
#       NOT by itself evidence that the analytic formula is wrong at that level:
#       CDH sums group-level contributions WITHIN a cluster before squaring them,
#       so offsetting contributions can legitimately lower the estimate, and no
#       theorem makes coarser clusters give larger SEs.
#   (2) the provinces are few and very unequal — 24 jurisdictions against an
#       effective count reported in 06f_effective_clusters.csv — so what IS in
#       doubt is whether the analytic approximation is reliable here, not whether
#       the level is the right one.
# Neither the analytic formula nor the pairs bootstrap of part 5 has strong
# finite-sample guarantees with this many effective clusters. They are reported
# together, and the conclusions are not allowed to turn on which of the two puts
# a given estimate under a conventional threshold. The school-level analytic
# variance remains the command's default and is reported for comparability.
#
# The static TWFE benchmark is included for one narrow purpose: it shows, on the
# same panel, that intra-provincial dependence EXISTS (its SEs widen the
# conventional way). It does NOT calibrate how much CDH's variance should widen —
# the two have different estimands, weights and influence functions.
#
# Part 5 asks the question the analytic province-clustered SEs cannot answer on
# their own: with this few effective clusters the cluster-robust asymptotics
# carry no strong finite-sample guarantee, and the remedy is to resample. WHICH
# resampling scheme is available depends on the estimator, and the difference is
# not cosmetic — the wild cluster bootstrap is defined on regression residuals
# and is the better-behaved scheme with few clusters, so only TWFE gets it.
# Part 5 covers the two regression-based estimators; the primary estimator's own
# province-clustered bootstrap is in 06q_cdh_bootstrap.R, which is separate only
# because it costs hours to run.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(did); library(fixest); library(fect); library(polars)})

d <- readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds"))

# Region: the five standard Argentine regions, used only as a third, coarser
# clustering level to show the direction the SEs move in as clusters are removed.
region_of <- function(p)
  fifelse(grepl("Jujuy|Salta|Tucum|Catamarca|Santiago", p), "NOA",
  fifelse(grepl("Formosa|Chaco|Corrientes|Misiones", p),    "NEA",
  fifelse(grepl("Mendoza|San Juan|San Luis|La Rioja", p),   "Cuyo",
  fifelse(grepl("Neuqu|R.o Negro|Chubut|Santa Cruz|Tierra", p), "Patagonia", "Centro"))))

d[, `:=`(prov = as.integer(factor(province)),
         reg  = as.integer(factor(region_of(province))))]

LEVELS <- list(school = "gid", province = "prov", region = "reg")

# ---- 1. CDH (primary) and CS (robustness) at each clustering level -----------
cdh <- function(y, clus, lab) {
  set_call_seed(lab)
  o <- tryCatch(DIDmultiplegtDYN::did_multiplegt_dyn(
         df = as.data.frame(d), outcome = y, group = "gid", time = "wave",
         treatment = "main_treatment", effects = 2, placebo = 1, cluster = clus,
         graph_off = TRUE), error = function(e) NULL)
  if (is.null(o)) return(NULL)
  data.table(estimator = "CDH (primary)", att = o$results$ATE[1], se = o$results$ATE[2],
             se_eff1 = o$results$Effects[1, 2], se_eff2 = o$results$Effects[2, 2],
             se_placebo1 = o$results$Placebos[1, 2])
}

# CS needs the absorbing subsample, built here exactly as in 06_estimation.R.
setorder(d, gid, wave)
dcs <- cs_sample(d)

cs <- function(y, clusvar, lab) {
  set_call_seed(lab)
  fit <- tryCatch(att_gt(yname = y, tname = "wave", idname = "gid", gname = "G",
           xformla = XF_COMP, data = dcs, control_group = "notyettreated",
           allow_unbalanced_panel = TRUE, bstrap = TRUE, biters = BITERS,
           est_method = "dr", clustervars = clusvar), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  a <- suppressWarnings(aggte(fit, type = "simple", na.rm = TRUE))
  data.table(estimator = "CS (robustness)", att = a$overall.att, se = a$overall.se,
             se_eff1 = NA_real_, se_eff2 = NA_real_, se_placebo1 = NA_real_)
}

twfe <- function(y, clusvar) {
  m <- feols(as.formula(paste(y, "~ main_treatment | gid + wave")), data = d,
             cluster = as.formula(paste0("~", clusvar)))
  data.table(estimator = "TWFE (benchmark)", att = coef(m)[[1]], se = se(m)[[1]],
             se_eff1 = NA_real_, se_eff2 = NA_real_, se_placebo1 = NA_real_)
}

res <- rbindlist(lapply(c("y_lang", "y_math"), function(y) {
  subj <- ifelse(y == "y_lang", "language", "mathematics")
  rbindlist(lapply(names(LEVELS), function(lv) {
    v <- LEVELS[[lv]]
    # CS clusters on the id when asked for the school level (its own default).
    csvar <- if (lv == "school") NULL else v
    rbindlist(list(cdh(y, v, paste("cdh", y, lv)),
                   cs(y, csvar, paste("cs", y, lv)),
                   twfe(y, if (lv == "school") "gid" else v)), fill = TRUE)[
      , `:=`(subject = subj, cluster_level = lv)]
  }), fill = TRUE)
}), fill = TRUE)

setcolorder(res, c("estimator", "subject", "cluster_level", "att", "se"))
num <- c("att", "se", "se_eff1", "se_eff2", "se_placebo1")
res[, (num) := lapply(.SD, function(x) round(x, 4)), .SDcols = num]
fwrite(res, file.path(DIR_TABLES, "06f_clustering.csv"))

# ---- 2. How many clusters actually carry information ------------------------
# G is 24, but the provinces are very unequal, so the count overstates how much
# independent information they carry. The effective count is the inverse
# Herfindahl of each province's share, on three bases: all schools, unique
# switching schools, and switcher-waves at each horizon. Identification comes
# from the treatment CHANGES, so the switcher-based figures are the informative
# ones -- a province with 10,000 schools and 3 switchers carries almost nothing.
g_eff <- function(x) sum(x)^2 / sum(x^2)
d[, dD := main_treatment - shift(main_treatment), by = gid]

base_counts <- list(
  schools  = d[, .(n = uniqueN(gid)), by = province]$n,
  switchers = d[, .(n = uniqueN(gid[dD != 0 & !is.na(dD)])), by = province][n > 0]$n
)
eff <- rbindlist(lapply(names(base_counts), function(b)
  data.table(base = b, n_clusters = length(base_counts[[b]]),
             effective_clusters = round(g_eff(base_counts[[b]]), 1))))

# The basis closest to the ideal is the Herfindahl over the estimator's OWN
# effective weights, which the package does not expose. Per horizon is as close as
# it gets: the command reports how many switchers support each effect but not
# their provincial composition, so that is recovered from the estimation sample it
# saves -- `did_effect` records which event-study effect each (g,t) cell is used
# for. Reported so the claim that the switcher-based approximation is stable can
# be checked rather than asserted; this is one analytic call, no bootstrap.
hz <- tryCatch({
  set_call_seed("effective clusters by horizon")
  o <- DIDmultiplegtDYN::did_multiplegt_dyn(
         df = as.data.frame(d), outcome = "y_lang", group = "gid", time = "wave",
         treatment = "main_treatment", effects = 3, placebo = 0, cluster = "gid",
         graph_off = TRUE, save_sample = TRUE)
  as.data.table(o$save_sample)[!is.na(did_effect), .(n = .N), by = .(province, did_effect)]
}, error = function(e) NULL)
if (!is.null(hz))
  eff <- rbind(eff, hz[order(did_effect), .(base = paste0("switcher-waves, Effect_", did_effect),
                                            n_clusters = .N,
                                            effective_clusters = round(g_eff(n), 1)),
                       by = did_effect][, did_effect := NULL])
fwrite(eff, file.path(DIR_TABLES, "06f_effective_clusters.csv"))

# ---- 3. Do provincial shocks confound the estimate? -------------------------
# Province-specific year effects, via the estimator's own `trends_nonparam`.
# This is the check a reader asks for once told that adoption is provincial:
# the 2024 funding cut, the staggered programme launches and the differing
# provincial recoveries are all absorbed. It is feasible here, which is not
# obvious and was measured before running it: demeaning the treatment
# transitions by province-year retains 83.5% of their variance (79.0% of the
# level's), because schools within a province do NOT all switch in the same
# year — only 7 of 120 province-year cells have no internal variation. The
# exception is the one province that executed its programme in a single year.
# Reported as SPECIFICATION robustness on the PRIMARY estimator: it varies an
# identifying assumption (parallel trends now holds within province-year), so
# it is not a free check, but the conclusion should survive it if provincial
# assignment is exogenous to school trajectories as the design argues.
trends <- rbindlist(lapply(c("y_lang", "y_math"), function(y) {
  subj <- ifelse(y == "y_lang", "language", "mathematics")
  base <- cdh(y, "gid", paste("trends base", y))
  tnp  <- tryCatch({
    set_call_seed(paste("trends nonparam", y))
    o <- DIDmultiplegtDYN::did_multiplegt_dyn(
           df = as.data.frame(d), outcome = y, group = "gid", time = "wave",
           treatment = "main_treatment", effects = 2, placebo = 1, cluster = "gid",
           trends_nonparam = "prov", graph_off = TRUE)
    data.table(att = o$results$ATE[1], se = o$results$ATE[2],
               switchers = o$N_switchers_effect_1)
  }, error = function(e) NULL)
  rbindlist(list(
    data.table(spec = "main", subject = subj, att = base$att, se = base$se),
    if (!is.null(tnp)) data.table(spec = "province-specific year trends",
                                  subject = subj, att = tnp$att, se = tnp$se)
  ), fill = TRUE)
}), fill = TRUE)
trends[, (c("att", "se")) := lapply(.SD, function(x) round(x, 4)), .SDcols = c("att", "se")]
fwrite(trends, file.path(DIR_TABLES, "06f_province_trends.csv"))

# ---- 4. The evidence behind the two decisions -------------------------------
# Three diagnostics that the working document cites as facts, committed here so
# they are regenerable rather than computed by hand.
#   (a) SWITCHER DIRECTION. Reported as a DESCRIPTIVE feature of the design: how
#       much of the switching runs in each direction. It is deliberately NOT
#       offered as the explanation of why province clustering shrinks the SEs.
#       The formula-level mechanism is enough for that — contributions are summed
#       within a cluster before being squared, so offsetting terms can lower the
#       estimate, and no theorem makes coarser clusters give larger SEs. Tying
#       the shrinkage to THIS count would require showing that these are the
#       contributions doing the offsetting, and the package does not expose the
#       influence functions, so the link is not claimed.
#   (b) PROVINCE-YEAR VARIANCE. Province-specific year trends are only feasible
#       if treatment variation survives demeaning by province-year. It does,
#       because schools within a province do not all switch in the same year.
#   (c) EIGENVALUES of the primary estimator's joint effects-and-placebos
#       covariance matrix. HonestDiD needs a valid covariance matrix; this one
#       is not positive semi-definite, which is why the sensitivity analysis is
#       not run on the primary estimator.
diag <- list()

sm <- tryCatch({
  set_call_seed("switcher direction")
  o <- DIDmultiplegtDYN::did_multiplegt_dyn(
         df = as.data.frame(d), outcome = "y_lang", group = "gid", time = "wave",
         treatment = "main_treatment", effects = 1, placebo = 0, cluster = "gid",
         graph_off = TRUE, save_sample = TRUE)
  s <- as.data.table(o$save_sample)
  tag <- grep("did_sample", names(s), value = TRUE)[1]
  s[, .(gid, tag = get(tag))][tag %like% "Switcher"]
}, error = function(e) NULL)
if (!is.null(sm)) {
  ins <- sm[tag %like% "in", .N]; outs <- sm[tag %like% "out", .N]
  diag[[length(diag) + 1]] <- data.table(
    diagnostic = "switcher direction",
    quantity = c("switcher-waves in", "switcher-waves out", "% in minority direction"),
    value = c(ins, outs, round(100 * min(ins, outs) / (ins + outs), 1)))
}

# (b) counting only, no estimation
dd <- d[!is.na(dD)]
dd[, `:=`(dD_dm = dD - mean(dD),
          tc_dm = main_treatment - mean(main_treatment)), by = .(province, year)]
cells <- d[!is.na(dD), .(n = .N, ch = sum(dD != 0), up = sum(dD > 0), dn = sum(dD < 0)),
           by = .(province, year)]
diag[[length(diag) + 1]] <- data.table(
  diagnostic = "province-year trends: feasibility",
  quantity = c("% of transition variance retained", "% of level variance retained",
               "province-year cells", "cells with no internal variation"),
  value = c(round(100 * var(dd$dD_dm) / var(dd$dD), 1),
            round(100 * var(dd$tc_dm) / var(dd$main_treatment), 1),
            nrow(cells),
            cells[ch == 0 | (ch == n & (up == 0 | dn == 0)), .N]))

# (c) is the joint covariance matrix usable for HonestDiD?
ev <- tryCatch({
  set_call_seed("vcov eigenvalues")
  o <- DIDmultiplegtDYN::did_multiplegt_dyn(
         df = as.data.frame(d), outcome = "y_math", group = "gid", time = "wave",
         treatment = "main_treatment", effects = 2, placebo = 2, cluster = "gid",
         graph_off = TRUE)
  eigen(o$coef$vcov, only.values = TRUE)$values
}, error = function(e) NULL)
if (!is.null(ev))
  diag[[length(diag) + 1]] <- data.table(
    diagnostic = "vcov of effects and placebos (mathematics)",
    quantity = c(paste("eigenvalue", seq_along(ev)), "positive semi-definite"),
    value = c(round(ev, 6), as.integer(all(ev > -1e-10))))

fwrite(rbindlist(diag), file.path(DIR_TABLES, "06f_inference_diagnostics.csv"))

# ---- 5. Resampling inference at the province, where the estimator allows it --
# With this few effective clusters the analytic province-clustered SEs of part 1
# carry no strong finite-sample guarantee, and the remedy is to resample. Which
# scheme is available is a property of the estimator, not a choice:
#   TWFE is one coefficient in a linear regression, so the WILD cluster bootstrap
#     applies — the scheme designed for exactly this case.
#   FEct's ATT is an average of imputed counterfactuals rather than a
#     coefficient, so the wild bootstrap does not apply; the package's own
#     nonparametric cluster bootstrap does, and takes the province as the block.
#   CDH is in the same position as FEct and is handled in 06q_cdh_bootstrap.R.
# Both are reported against their school-clustered counterparts. The question is
# whether the negative estimates persist once provincial dependence is allowed —
# NOT which estimator has the smaller p-value, since the two receive different
# resampling schemes and their p-values are not comparable as evidence.
WCB_B <- 999L   # (B+1)*alpha integer at the conventional levels

# Wild cluster bootstrap-t, restricted null, Rademacher weights [@cameron2008;
# @mackinnon2018]. The fixed effects are swept by Frisch-Waugh once, so a
# replication costs one demeaning and one univariate fit instead of a refit:
# under H0 the bootstrap outcome is the restricted fit plus reweighted residuals,
# and the restricted fit lies in the span of the fixed effects, so it vanishes on
# demeaning. The degrees-of-freedom factor is constant across replications and
# cancels in the comparison |t*| >= |t|, so it is omitted from both sides;
# `wcb_check` below verifies against fixest what the omission costs.
wcb_p <- function(dat, y, D, FE, gvec, b0 = 0, B = WCB_B, lab = "", chunk = 200L) {
  yv <- dat[[y]] - b0 * dat[[D]]
  Dt <- as.numeric(demean(dat[[D]], FE))
  u  <- as.numeric(demean(yv, FE))          # restricted residuals
  DD <- sum(Dt^2)
  gf <- factor(gvec); G <- nlevels(gf)
  tstat <- function(Y) {                    # Y: n x m demeaned outcomes
    b <- as.numeric(crossprod(Dt, Y) / DD)
    S <- rowsum(Dt * (Y - Dt %o% b), gf)    # G x m cluster score sums
    b / sqrt(colSums(S^2) / DD^2)
  }
  t0 <- tstat(matrix(u, ncol = 1))
  set_call_seed(lab)
  hit <- 0L
  for (s in seq(1L, B, by = chunk)) {
    m <- min(chunk, B - s + 1L)
    V <- matrix(sample(c(-1, 1), G * m, replace = TRUE), nrow = G)[as.integer(gf), , drop = FALSE]
    hit <- hit + sum(abs(tstat(demean(V * u, FE))) >= abs(t0))
  }
  list(t = t0, p = hit / B, G = G)
}

# The bootstrap CI inverts the test: the values of beta0 the WCB does not reject
# at 5%. A grid rather than a bisection, because a Monte Carlo p-value is noisy
# and bisecting on it would report a precision the procedure does not have. The
# grid step is a fifth of the analytic standard error, and is reported.
wcb_ci <- function(dat, y, D, FE, gvec, bhat, se, lab, half = 4, step = 0.2) {
  grid <- bhat + seq(-half, half, by = step) * se
  keep <- vapply(seq_along(grid), function(i)
    wcb_p(dat, y, D, FE, gvec, b0 = grid[i], lab = paste(lab, i))$p, 0) > 0.05
  if (!any(keep)) return(c(NA_real_, NA_real_, step * se))
  c(min(grid[keep]), max(grid[keep]), step * se)
}

res_tw <- rbindlist(lapply(c("y_lang", "y_math"), function(y) {
  subj <- ifelse(y == "y_lang", "language", "mathematics")
  dt <- d[!is.na(get(y))]
  FE <- dt[, .(gid, wave)]
  m_s <- feols(as.formula(paste(y, "~ main_treatment | gid + wave")), data = dt, cluster = ~gid)
  m_p <- feols(as.formula(paste(y, "~ main_treatment | gid + wave")), data = dt, cluster = ~prov)
  bhat <- coef(m_p)[["main_treatment"]]; se_p <- se(m_p)[["main_treatment"]]
  w  <- wcb_p(dt, y, "main_treatment", FE, dt$prov, lab = paste("wcb", y))
  ci <- wcb_ci(dt, y, "main_treatment", FE, dt$prov, bhat, se_p, paste("wcbci", y))
  rbindlist(list(
    data.table(estimator = "TWFE (benchmark)", subject = subj, cluster_level = "school",
               inference = "analytic cluster-robust", att = coef(m_s)[[1]],
               se = se(m_s)[[1]], p = pvalue(m_s)[[1]]),
    data.table(estimator = "TWFE (benchmark)", subject = subj, cluster_level = "province",
               inference = "analytic cluster-robust", att = bhat, se = se_p,
               p = pvalue(m_p)[[1]]),
    data.table(estimator = "TWFE (benchmark)", subject = subj, cluster_level = "province",
               inference = "wild cluster bootstrap", att = bhat, se = NA_real_, p = w$p,
               ci_lo = ci[1], ci_hi = ci[2], ci_step = ci[3], n_boot = WCB_B,
               # audit: the omitted dof factor is exactly sqrt(G/(G-1)), because
               # the school effects are nested in the province and fixest drops
               # them from K. Equality of these two columns is the check that
               # the hand-rolled t reproduces the package's.
               t_ratio = w$t / (bhat / se_p), dof_factor = sqrt(w$G / (w$G - 1)))
  ), fill = TRUE)
}), fill = TRUE)

# FEct: same fit as 11_reversals.R, varying only the bootstrap block.
fect_boot <- function(y, cl, lab) {
  set_call_seed(lab)
  o <- tryCatch(fect::fect(as.formula(paste(y, "~ main_treatment")),
         data = as.data.frame(d[, .(gid, wave, y_lang, y_math, main_treatment, prov)]),
         index = c("gid", "wave"), method = "fe", force = "two-way", se = TRUE,
         nboots = NBOOTS, parallel = TRUE, min.T0 = 1, seed = SEED, cl = cl),
       error = function(e) NULL)
  if (is.null(o)) return(NULL)
  a <- o$est.avg
  data.table(estimator = "FEct (FE imputation)", subject = ifelse(y == "y_lang", "language", "mathematics"),
             cluster_level = if (is.null(cl)) "school" else "province",
             inference = "nonparametric cluster bootstrap",
             att = a[1, 1], se = a[1, 2], p = a[1, 5],
             ci_lo = a[1, 3], ci_hi = a[1, 4], n_boot = NBOOTS)
}

res_fe <- rbindlist(lapply(c("y_lang", "y_math"), function(y)
  rbindlist(list(fect_boot(y, NULL,   paste("fect school", y)),
                 fect_boot(y, "prov", paste("fect prov", y))), fill = TRUE)), fill = TRUE)

resamp <- rbindlist(list(res_tw, res_fe), fill = TRUE)
nm <- c("att", "se", "p", "ci_lo", "ci_hi", "ci_step", "t_ratio", "dof_factor")
resamp[, (nm) := lapply(.SD, function(x) round(x, 4)), .SDcols = nm]
fwrite(resamp, file.path(DIR_TABLES, "06f_resampling.csv"))
