# ==============================================================================
# 06q_cdh_bootstrap.R — Province-clustered bootstrap of the PRIMARY estimator
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Inputs  : data/processed/estimation_sample.rds
# Outputs : output/tables/06q_cdh_bootstrap.csv       estimate, analytic and
#                                                     bootstrapped SEs, per parameter
#           output/tables/06q_bootstrap_draws.csv     the bootstrap distribution
#           output/tables/06q_bootstrap_support.csv   failed draws and how many
#                                                     distinct provinces each holds
#
# 06f establishes that clustering at the province SHRINKS the primary estimator's
# analytic standard errors. That is not by itself a defect: the CDH analytic
# formula sums group-level contributions WITHIN a cluster before squaring them,
# so contributions of opposite sign inside a province can lower the estimated
# variance. There is no theorem making coarser clusters give larger SEs. What the
# shrinkage does not settle is whether the analytic approximation is reliable
# when the effective number of provinces is near 8, and that is what resampling
# answers — not by being correct where the formula is wrong, but by showing
# whether the two procedures agree.
#
# The command's `bootstrap` argument resamples at whatever level `cluster` names
# (package source: `bs_group <- if (!is.null(cluster)) cluster`), drawing whole
# clusters with replacement. It is a nonparametric pairs cluster bootstrap. The
# wild cluster bootstrap, which is the better-behaved scheme with few clusters,
# is not available here: it is defined on the residuals of a linear regression
# and this estimator's ATT is not a regression coefficient.
#
# Part B repeats the same resampling by hand. Not out of distrust: the package
# returns only the standard deviation of the draws, and the quantities needed to
# judge whether the resampling is informative at all — how many draws fail, how
# much provincial support each draw retains, and the shape of the distribution
# whose SD is being reported — are not exposed. Agreement between the two SEs is
# the check that the hand-rolled resampling matches the package's.
# ==============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages(library(polars))

# In-process replication; the package default spawns one R subprocess per draw,
# which costs more than the estimation itself.
options(DID_BOOTSTRAP_METHOD = "row_replication")
BOOT_B <- as.integer(Sys.getenv("BOOT_B", "999"))   # package bootstrap replications
DRAW_B <- as.integer(Sys.getenv("DRAW_B", "299"))   # hand-rolled diagnostic draws
# Part A costs about an hour and its results are deterministic given BOOT_B, so
# the two parts can be run separately; PARTS = "B" reuses A from its checkpoint.
PARTS  <- Sys.getenv("PARTS", "AB")

d <- readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds"))
d[, prov := as.integer(factor(province))]

# The headline specification, and only it: the bootstrap is a statement about
# the estimate the dissertation reports, not a second robustness grid.
EFFECTS <- 2L
PLACEBO <- max_placebo(EFFECTS, length(WAVES))

cdh_fit <- function(dat, y, cluster = "prov", extra = NULL) {
  args <- list(df = as.data.frame(dat), outcome = y, group = "gid", time = "wave",
               treatment = "main_treatment", effects = EFFECTS, placebo = PLACEBO,
               cluster = cluster, graph_off = TRUE)
  if (!is.null(extra)) args <- modifyList(args, extra)
  do.call(DIDmultiplegtDYN::did_multiplegt_dyn, args)
}

# Estimate and SE of every reported parameter, from a fitted object.
params_of <- function(o) rbindlist(list(
  data.table(parameter = "ATE", att = o$results$ATE[1], se = o$results$ATE[2]),
  data.table(parameter = paste0("Effect_", seq_len(nrow(o$results$Effects))),
             att = o$results$Effects[, 1], se = o$results$Effects[, 2]),
  if (!is.null(o$results$Placebos))
    data.table(parameter = paste0("Placebo_", seq_len(nrow(o$results$Placebos))),
               att = o$results$Placebos[, 1], se = o$results$Placebos[, 2])))

# ---- A. Analytic SEs at both levels, and the package's own bootstrap ---------
main <- if (!grepl("A", PARTS)) {
  f <- file.path(DIR_TABLES, "06q_cdh_bootstrap_partA.csv")
  if (!file.exists(f)) stop("PARTS excludes A but no part-A checkpoint exists")
  fread(f)
} else rbindlist(lapply(c("y_lang", "y_math"), function(y) {
  subj <- ifelse(y == "y_lang", "language", "mathematics")
  a_sch <- params_of(cdh_fit(d, y, cluster = "gid"))[, .(parameter, att, se_school = se)]
  a_prv <- params_of(cdh_fit(d, y, cluster = "prov"))[, .(parameter, se_province = se)]
  set_call_seed(paste("cdh boot", y))
  bt <- params_of(cdh_fit(d, y, cluster = "prov",
                          extra = list(bootstrap = c(reps = BOOT_B, seed = SEED))))
  merge(merge(a_sch, a_prv, by = "parameter"),
        bt[, .(parameter, se_boot_package = se)], by = "parameter")[, subject := subj]
}), fill = TRUE)

# Checkpointed before part B: part A costs the better part of an hour and part B
# costs more, so a failure in the second should not discard the first.
if (grepl("A", PARTS)) fwrite(main, file.path(DIR_TABLES, "06q_cdh_bootstrap_partA.csv"))

# ---- B. The same resampling by hand, to expose what the package absorbs ------
# A province drawn twice must enter as two independent copies, so the school
# identifier is rebuilt per copy; leaving it unchanged would merge the copies
# into one school and silently shrink the panel.
provs <- sort(unique(d$prov))
G <- length(provs)

one_draw <- function(y, b, offset) {
  # Deliberately NOT set_call_seed(): that helper hashes its label as the SUM of
  # the character codes, so any two labels that are anagrams share a seed — and
  # every permutation of the digits of `b` is one. Indexing draws by label
  # collapsed 999 draws onto 54 distinct resamples when this was first run.
  # Draws are therefore seeded arithmetically, which cannot collide.
  set.seed(SEED + offset + b)
  sel <- sample(provs, G, replace = TRUE)
  dat <- rbindlist(lapply(seq_along(sel), function(j) {
    z <- d[prov == sel[j]]
    z[, gid := gid * 1000L + j]
    z
  }))
  o <- tryCatch(cdh_fit(dat, y, cluster = "prov"), error = function(e) NULL)
  list(ok = !is.null(o), n_prov = uniqueN(sel),
       par = if (is.null(o)) NULL else params_of(o)[, .(parameter, att)])
}

draws <- rbindlist(lapply(c("y_lang", "y_math"), function(y) {
  subj <- ifelse(y == "y_lang", "language", "mathematics")
  offset <- if (y == "y_lang") 0L else 1000000L
  rbindlist(lapply(seq_len(DRAW_B), function(b) {
    r <- one_draw(y, b, offset)
    if (!r$ok) return(data.table(subject = subj, draw = b, n_prov = r$n_prov,
                                 parameter = NA_character_, att = NA_real_))
    r$par[, `:=`(subject = subj, draw = b, n_prov = r$n_prov)]
  }), fill = TRUE)
}), fill = TRUE)

fwrite(draws, file.path(DIR_TABLES, "06q_bootstrap_draws.csv"))

# Support: draws that failed outright, and how much of the provincial structure
# survives resampling with replacement. A draw holding 15 of 24 provinces is the
# expected behaviour of the scheme, not a defect, but it is what the standard
# error is averaging over and the reader should see it.
support <- draws[, .(draws = uniqueN(draw),
                     failed = uniqueN(draw[is.na(parameter)]),
                     min_distinct_provinces = min(n_prov),
                     median_distinct_provinces = as.numeric(median(n_prov)),
                     max_distinct_provinces = max(n_prov)), by = subject]
fwrite(support, file.path(DIR_TABLES, "06q_bootstrap_support.csv"))

own <- draws[!is.na(parameter), .(se_boot_own = sd(att, na.rm = TRUE),
                                  q025 = quantile(att, 0.025, na.rm = TRUE),
                                  q975 = quantile(att, 0.975, na.rm = TRUE),
                                  n_used = .N), by = .(subject, parameter)]

out <- merge(main, own, by = c("subject", "parameter"), all.x = TRUE)
setcolorder(out, c("subject", "parameter", "att", "se_school", "se_province",
                   "se_boot_package", "se_boot_own", "q025", "q975", "n_used"))
nm <- setdiff(names(out), c("subject", "parameter", "n_used"))
out[, (nm) := lapply(.SD, function(x) round(x, 4)), .SDcols = nm]
setorder(out, subject, parameter)
fwrite(out, file.path(DIR_TABLES, "06q_cdh_bootstrap.csv"))
