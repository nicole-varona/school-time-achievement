# ==============================================================================
# 06g_horizon_sensitivity.R — Horizon window and four estimator options
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Inputs  : data/processed/estimation_sample.rds
#           data/processed/panel.rds   (contemporaneous composition, for controls)
# Outputs : output/tables/06g_horizon_sensitivity.csv  ATT by aggregation window
#           output/tables/06g_cdh_options.csv          same_switchers, trends_lin,
#                                                      controls, SES-stratum trends
#           output/tables/06g_ses_trends_feasibility.csv  is that last one feasible
#
# Five checks the PRIMARY estimator provides and the pipeline did not run. Each
# varies one implementation choice and is estimated with CDH, per the rule that
# specification robustness belongs on the estimator the headline rests on.
#
#   (1) AGGREGATION WINDOW. The headline ATT averages two horizons. That is a
#       choice, and a live one: on the non-reverting subsample the same estimator
#       gives +0.069 over two horizons and +0.049 over three (11_nonreverter vs
#       11b_like_for_like), a 29% gap from the window alone. Re-estimated here
#       over one, two and three horizons on the headline sample.
#   (2) same_switchers. Switchers collapse across horizons (2,297 -> 1,205 ->
#       790), so each Effect_k is estimated on a different set of schools and the
#       dynamic profile mixes real dynamics with compositional change. This
#       option restricts every effect to the switchers that support all of them,
#       which is the estimator's own answer to that problem.
#   (3) trends_lin. Group-specific linear trends: estimation moves to the
#       outcome's first difference and the event-study effects are summed back
#       up. It relaxes parallel trends to parallel SECOND differences. The
#       average cumulative effect per treatment unit is not defined under this
#       option, so the comparison runs on the per-horizon effects.
#   (4) controls. Time-varying composition, which is what the option is for
#       (baseline covariates are differenced out). Run on the complete-case
#       subsample together with its own uncontrolled counterpart, so the two
#       differ by the controls and not by the sample.
#   (5) trends_nonparam over SES strata. The theoretical worry about socioeconomic
#       status is that it predicts differential TRENDS, which no covariate can
#       express and which a time-invariant one cannot even enter. This is the
#       estimator's own way to allow it: each stratum gets its own year effects.
#
# SEs are analytic and clustered at the school throughout, as in the headline.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(polars)})

d <- readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds"))
if (!"bl_mom_pre" %in% names(d))
  stop("bl_mom_pre missing from the estimation sample; run 02 first")

SUBJECTS <- c(y_lang = "language", y_math = "mathematics")

# One CDH call, returned long: one row per reported parameter. `switchers` is
# the count supporting each effect, which is the identifying variation at that
# horizon. ATE is absent under trends_lin, hence the guard.
cdh <- function(dat, yname, spec, tname = "main_treatment", ...) {
  set_call_seed(paste(spec, yname, tname))
  o <- tryCatch(DIDmultiplegtDYN::did_multiplegt_dyn(
         df = as.data.frame(dat), outcome = yname, group = "gid", time = "wave",
         treatment = tname, cluster = "gid", graph_off = TRUE, ...),
       error = function(e) NULL)
  if (is.null(o)) return(NULL)
  grab <- function(m) {
    if (is.null(m)) return(NULL)
    m <- as.data.frame(m)
    data.table(parameter = trimws(rownames(m)), estimate = m[[1]], se = m[[2]],
               switchers = if ("Switchers" %in% names(m)) m[["Switchers"]] else NA_real_)
  }
  out <- rbindlist(list(grab(o$results$Effects), grab(o$results$Placebos)), fill = TRUE)
  if (!is.null(o$results$ATE))
    out <- rbind(out, data.table(parameter = "ATE", estimate = o$results$ATE[1],
                                 se = o$results$ATE[2], switchers = NA_real_), fill = TRUE)
  out[, `:=`(spec = spec, subject = unname(SUBJECTS[yname]))]
}

run_both <- function(dat, spec, ...)
  rbindlist(lapply(names(SUBJECTS), function(y) cdh(dat, y, spec, ...)), fill = TRUE)

tidy_out <- function(x, first = NULL) {
  setcolorder(x, c(first, "spec", "subject", "parameter", "estimate", "se", "switchers"))
  x[, `:=`(estimate = round(estimate, 4), se = round(se, 4))][]
}

# ---- 1. Aggregation window ---------------------------------------------------
# placebo = 1 throughout so the pre-trend evidence is held fixed and the only
# thing varying is how many post-treatment horizons the ATT averages over.
# Run on BOTH treatment regimes: the full-day arm carries the largest positive of
# the primary estimator, so the window has to be varied there too or the check
# covers the smaller estimate only.
REGIMES <- list("combined" = list(dat = d, tname = "main_treatment"),
                "full-day" = list(dat = d[!is.na(main_treatment_full)], tname = "main_treatment_full"))
windows <- rbindlist(lapply(names(REGIMES), function(rg) {
  a <- REGIMES[[rg]]
  rbindlist(lapply(1:3, function(k)
    run_both(a$dat, sprintf("%d horizon(s)", k), tname = a$tname,
             effects = k, placebo = 1)), fill = TRUE)[, regime := rg]
}), fill = TRUE)
fwrite(tidy_out(windows, "regime"), file.path(DIR_TABLES, "06g_horizon_sensitivity.csv"))

# ---- 1b. Every placebo the design admits -------------------------------------
# The table above holds the placebo at one so that only the number of horizons
# varies, which is what that comparison needs. The cost is that it never reports a
# SECOND placebo, and the main specification admits one. The package allows as
# many placebos as effects, bounded by the periods left before the switch, so two
# horizons over six waves admit two (min(2, 6 - 2 - 1)). Requesting one was a
# floor set by the call and not by the design, so the pre-treatment evidence of
# the primary specification was being read off half of what it has. It is also
# silent: the package truncates an over-request with a message rather than
# failing, so nothing ever flagged the shortfall.
#
# THE SECOND PLACEBO IS NOT A NEW ESTIMATE. Placebo_1 and Placebo_2 are the k = -1
# and k = -2 leads of the same event study, identical to the last decimal against
# 07_event_study_cdh.csv, which is why the check below is an equality and not a
# comparison. This table exists so the prose can cite them AS the specification's
# own placebos: the significant mathematics lead at k = -2 was being discussed as
# an event-study result separate from "the placebo of the primary estimator",
# and they are one object. Two exhibits, one diagnostic -- so they never
# corroborate each other.
PL_MAIN <- max_placebo(2L, length(WAVES))
max_pl <- rbindlist(lapply(names(REGIMES), function(rg) {
  a <- REGIMES[[rg]]
  run_both(a$dat, "main (2 horizons), all placebos", tname = a$tname,
           effects = 2, placebo = PL_MAIN)[, regime := rg]
}), fill = TRUE)
if (!"Placebo_2" %in% max_pl$parameter)
  stop("the second placebo was not returned; the design no longer supports it and ",
       "the manuscript reads it")

# The identity with the event-study leads is ASSERTED IN 07, not here. It was here
# first, and that was wrong: 07 runs later in the pipeline, so on a clean run this
# script compared its fresh placebos against the PREVIOUS run's event study and
# halted. A check must not depend on a file produced after it.
fwrite(tidy_out(max_pl, "regime"), file.path(DIR_TABLES, "06g_max_placebos.csv"))

opts <- list()

# ---- 2. same_switchers -------------------------------------------------------
# Three effects and two placebos, with and without the restriction: the gap is
# how much of the dynamic profile is composition rather than dynamics.
opts$a <- run_both(d, "3 horizons, all switchers",  effects = 3, placebo = 2)
opts$b <- run_both(d, "3 horizons, same switchers", effects = 3, placebo = 2,
                   same_switchers = TRUE, same_switchers_pl = TRUE)
# The three-horizon restriction keeps only 337 schools and answers a different
# question from the headline, which aggregates two horizons. This pair fixes the
# switcher set at exactly the window the headline uses, so the comparison speaks
# to `Av_tot_eff` directly instead of by analogy.
opts$b2 <- run_both(d, "2 horizons, same switchers", effects = 2, placebo = 1,
                    same_switchers = TRUE, same_switchers_pl = TRUE)
# The fourth cell of the 2x2. `same_switchers` fixes WHICH schools identify the
# parameter; `weight` fixes what each of them counts for inside that population.
# The two operations are independent, and nothing requires them to move the
# estimate in the same direction -- which is the whole point of running the cross.
opts$b3 <- run_both(d, "2 horizons, same switchers, pupil-weighted", effects = 2,
                    placebo = 1, same_switchers = TRUE, same_switchers_pl = TRUE,
                    weight = "enroll_g6")

# ---- 3. trends_lin -----------------------------------------------------------
# Reported for what it is: with 43% of treated schools left-censored and two
# usable placebos, group-specific linear trends are thinly identified here. The
# check is whether the effects survive the weaker trend assumption at all.
opts$c <- run_both(d, "main (2 horizons)", effects = 2, placebo = 1)
opts$d <- run_both(d, "group-specific linear trends", effects = 2, placebo = 1,
                   trends_lin = TRUE)

# ---- 4. controls -------------------------------------------------------------
# Contemporaneous composition, merged from the annual panel (the estimation
# sample carries baseline covariates only). In pre-periods these are
# pre-treatment, so using them is legitimate rather than a bad control.
CTRL <- c("pct_female", "pct_repeater", "pct_preschool")
panel <- readRDS(file.path(DIR_PROCESSED, "panel.rds"))
dc <- merge(d, panel[, c("school_id", "year", CTRL), with = FALSE],
            by = c("school_id", "year"), all.x = TRUE)
if (nrow(dc) != nrow(d)) stop("composition merge changed the number of rows")
dc <- dc[complete.cases(dc[, ..CTRL])]

opts$e <- run_both(dc, "complete-case, uncontrolled", effects = 2, placebo = 1)
opts$f <- run_both(dc, "complete-case, composition controls", effects = 2,
                   placebo = 1, controls = CTRL)

# ---- 5. SES-stratum-specific year trends ------------------------------------
# The theoretically motivated worry is not that socioeconomic status shifts the
# LEVEL of achievement -- the within-school differencing absorbs that -- but that
# it predicts differential TRENDS. A covariate cannot express that, and in this
# estimator a time-invariant one cannot enter at all (delta-X = 0). The native way
# to allow it is `trends_nonparam`: let each socioeconomic stratum carry its own
# year effects, so parallel trends is required only WITHIN stratum-year. This is
# the same tool applied to the province in 06f, on the dimension the selection
# argument points to instead (Veleda documents social vulnerability as an explicit
# school-selection criterion of the provincial programmes).
#
# Strata: terciles of the PREDETERMINED SES measure (`bl_mom_pre`, built in 02 from
# the waves preceding each school's own adoption), plus a fourth stratum for the
# schools that have no pre-treatment wave. Those are kept rather than dropped: they
# are largely the left-censored schools, a distinct population that carries more
# switchers than any tercile, and giving them their own year effects is more
# honest than letting the estimator discard them.
#
# Feasibility is measured before estimating, as it was for the province, because
# the option is only meaningful if treatment variation survives demeaning by
# stratum-year.
sch <- unique(d[, .(gid, bl_mom_pre)])
br  <- quantile(sch$bl_mom_pre, probs = seq(0, 1, length.out = 4), na.rm = TRUE)
sch[, ses_st := fifelse(is.na(bl_mom_pre), 4L,
                        as.integer(cut(bl_mom_pre, br, include.lowest = TRUE)))]
ds <- merge(d, sch[, .(gid, ses_st)], by = "gid")
setorder(ds, gid, wave)
ds[, dD := main_treatment - shift(main_treatment), by = gid]
dd <- ds[!is.na(dD)]
dd[, dD_dm := dD - mean(dD), by = .(ses_st, year)]
cells <- dd[, .(n = .N, ch = sum(dD != 0), up = sum(dD > 0), dn = sum(dD < 0)),
            by = .(ses_st, year)]
feas <- rbind(
  data.table(quantity = "% of transition variance retained",
             value = round(100 * var(dd$dD_dm) / var(dd$dD), 1)),
  data.table(quantity = "stratum-year cells", value = as.numeric(nrow(cells))),
  data.table(quantity = "cells with no internal variation",
             value = as.numeric(cells[ch == 0 | (ch == n & (up == 0 | dn == 0)), .N])),
  ds[, .(quantity = paste("switchers, stratum", ses_st),
         value = as.numeric(uniqueN(gid[dD != 0 & !is.na(dD)]))), by = ses_st][
           order(ses_st), .(quantity, value)])
fwrite(feas, file.path(DIR_TABLES, "06g_ses_trends_feasibility.csv"))

opts$g <- run_both(ds, "SES-stratum-specific year trends", effects = 2, placebo = 1,
                   trends_nonparam = "ses_st")

# ---- 8. more_granular_demeaning ----------------------------------------------
# The package documents its default variance as possibly conservative when the
# treatment changes more than once. This option defines finer switcher cohorts,
# which can reduce that conservatism when those cohorts hold enough switchers.
# Run ONCE, as an appendix sensitivity, and not adopted for the headline: with a
# borderline primary result and conclusions that rest on robustness rather than
# precision, optimising the variance downward adds no substantive information.
opts$h <- run_both(d, "more granular demeaning", effects = 2, placebo = 1,
                   more_granular_demeaning = TRUE)

fwrite(tidy_out(rbindlist(opts, fill = TRUE)),
       file.path(DIR_TABLES, "06g_cdh_options.csv"))
