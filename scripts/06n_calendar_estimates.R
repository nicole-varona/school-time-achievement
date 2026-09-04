# ==============================================================================
# 06n_calendar_estimates.R — The same estimand under three outcome calendars and
#                            the annual treatment definition
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Inputs  : data/processed/estimation_sample.rds
#           data/processed/calendar_designs.rds      (built by 06m)
#           output/tables/06m_calendar_support.csv
# Outputs : output/tables/06n_calendar_estimates.csv     aggregate, per cell
#           output/tables/06n_exposure_estimates.csv     effect + placebo by class
#           output/tables/06n_placebo_cancellation.csv   is the aggregate placebo
#                                                        zero, or cancellation?
#           output/tables/06n_cohort_diagnostic.csv      leave-one-cohort-out
#
# 06m established what each cell identifies. This script estimates the same
# quantity in each of them.
#
#   CALENDAR   six waves / five censuses / three censuses. This is the dimension
#              that moves the answer: the sign of the ATT depends on which waves
#              enter, and on neither the estimator nor the treatment definition.
#   TREATMENT  the main treatment (>=50% share in that year, full stop) and
#              nothing else. The two-year persistence rule was crossed with the
#              calendar here until it was retired on 12/08, because the two
#              interacted -- the rule was enforced on the annual history, so which
#              annual years a calendar could see decided what the rule removed
#              from the estimator's view. With one rule left the grid is a single
#              row, and the `definition` column keeps its single value, "main".
#
# Four things are run in every cell:
#   1. the aggregate effect at one and two horizons, with `same_switchers` and
#      with grade-6 enrolment weights (the pupil-average estimand rather than the
#      school-average one);
#   2. the same effect separately by CALENDAR TIME SINCE ADOPTION, so that no
#      reported number pools switchers measured at different distances from their
#      own switch. That is not accumulated exposure, which 06m reports alongside;
#   3. the placebo for every one of those classes, because an aggregate placebo
#      near zero is only reassuring if it is zero in each class rather than
#      cancellation between classes of opposite sign;
#   4. the estimate leaving out the 2023 and the 2024 adoption cohorts in turn,
#      as a diagnostic and not as an alternative estimator. 2023 carries both the
#      mathematics positive and the pre-trend problem; 2024 is the largest cohort
#      left in the three-census calendar once 2023 is removed, and that calendar
#      is the one whose negative does not dissolve with 2023.
#
# The classification is read from 06m and never re-derived here, so the two
# scripts cannot describe different schools. Placebos are requested at the maximum
# the calendar can carry, which for three periods and two horizons is none; a
# specification missing from the output is one the cell cannot support.
#
# A RULE, stated before the numbers: neither the calendar nor the treatment
# definition is chosen from these coefficients. Picking the cell that returns the
# largest estimate, after seeing them, is specification search. The criterion is
# clarity of the estimand, support, and the identification diagnostics, which is
# what 06m measures.
# ==============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(DIDmultiplegtDYN); library(polars)})

SUBJECTS <- c(y_lang = "language", y_math = "mathematics")
MIN_SW   <- 150   # same floor as 06j: below it a class is reported, not estimated
COHORTS  <- c(2023L, 2024L)   # the leave-one-out diagnostics
SPECS    <- list(
  "1 horizon"                  = list(effects = 1L, extra = NULL),
  "2 horizons"                 = list(effects = 2L, extra = NULL),
  "2 horizons, same switchers" = list(effects = 2L, extra = list(same_switchers = TRUE,
                                                                same_switchers_pl = TRUE)),
  "2 horizons, pupil-weighted" = list(effects = 2L, extra = list(weight = "enroll_g6")))

d   <- as.data.table(readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds")))
DES <- readRDS(file.path(DIR_PROCESSED, "calendar_designs.rds"))
sup <- fread(file.path(DIR_TABLES, "06m_calendar_support.csv"))

CELLS <- CJ(definition = names(DES$defs), design = names(DES$designs), sorted = FALSE)

# ---- One estimation ----------------------------------------------------------
# Waves are re-indexed 1..K inside the calendar, which is what makes the estimator
# treat them as consecutive periods -- the very convention under test.
frame <- function(dn, design, ids = NULL) {
  yrs <- DES$designs[[design]]
  x <- d[year %in% yrs & !is.na(get(DES$defs[[dn]]))]
  if (!is.null(ids)) x <- x[school_id %in% ids]
  x[, wv := match(year, yrs)]
  as.data.frame(x)
}

# The design bound on placebos: the periods a calendar leaves before the switch,
# and the package's own rule that there be no more placebos than effects. It is
# computed rather than discovered because the package truncates an over-request
# with a message instead of failing, and a truncation nobody reads is how a table
# ends up claiming more pre-treatment evidence than it has. The package can still
# tighten this further when the switch timing leaves no pre-period to compare, and
# that message is deliberately left visible.

# `cluster` is not passed: the package discards it when it equals `group`, and the
# group is the school, so passing it would document a choice that has no effect.
fit_cdh <- function(dat, yname, tv, effects, extra, label, n_periods) {
  set_call_seed(label)
  args <- list(df = dat, outcome = yname, group = "gid", time = "wv",
               treatment = tv, effects = effects, graph_off = TRUE,
               placebo = max_placebo(effects, n_periods))
  if (!is.null(extra)) args <- modifyList(args, extra)
  tryCatch(do.call(did_multiplegt_dyn, args), error = function(e) NULL)
}

tidy <- function(fit) {
  m2dt <- function(m) if (is.null(m)) NULL else
    data.table(parameter = trimws(rownames(m)), estimate = m[, "Estimate"],
               se = m[, "SE"], switchers = m[, "Switchers"])
  rbind(m2dt(fit$results$Effects), m2dt(fit$results$Placebos),
        data.table(parameter = "ATE", estimate = fit$results$ATE[1],
                   se = fit$results$ATE[2], switchers = NA_real_))
}

# Runs one specification and returns it tidied, tagged with `tags`.
run <- function(dat, dn, design, yname, spec, tags) {
  label <- paste(c(unlist(tags), yname, spec), collapse = " ")
  fit <- fit_cdh(dat, yname, DES$defs[[dn]], SPECS[[spec]]$effects, SPECS[[spec]]$extra,
                 label, length(DES$designs[[design]]))
  if (is.null(fit)) return(NULL)
  cbind(as.data.table(tags), subject = unname(SUBJECTS[yname]), spec = spec, tidy(fit))
}

over_cells <- function(f) rbindlist(lapply(seq_len(nrow(CELLS)), function(i)
  f(CELLS$definition[i], CELLS$design[i])), fill = TRUE)

# ---- 1. Aggregate estimate in each cell --------------------------------------
agg <- over_cells(function(dn, nm) {
  dat <- frame(dn, nm)
  rbindlist(lapply(names(SUBJECTS), function(y)
    rbindlist(lapply(names(SPECS), function(sp)
      run(dat, dn, nm, y, sp, list(definition = dn, design = nm))))))
})
if (!nrow(agg)) stop("every calendar-level estimation failed")

# The switcher count the estimator reports at the first horizon must match the
# reconstruction in 06m. Checked on language, whose outcome has no missing values;
# mathematics has 306, so the estimator drops those rows and its counts run lower.
chk <- merge(agg[spec == "1 horizon" & subject == "language" & parameter == "Effect_1",
                 .(definition, design, got = as.integer(switchers))],
             sup[, .(definition, design, want = supports_k1)],
             by = c("definition", "design"))
if (nrow(chk) != nrow(CELLS) || !identical(chk$got, chk$want))
  stop("switcher support disagrees with 06m: ",
       paste(chk$definition, chk$design, chk$got, "vs", chk$want, collapse = "; "))

# The six-wave cell of the main treatment IS the headline specification, re-indexed, so
# it must reproduce what 06 reports. Checked against that file rather than against
# fixed numbers, so the guard survives a re-estimation and still catches the thing
# it is for: re-indexing the waves silently moving the estimate, which would make
# every comparison across cells a comparison of two changes at once.
hf <- file.path(DIR_TABLES, "06_estimation_main.csv")
if (file.exists(hf)) {
  hl <- agg[definition == "main" & design == "six waves" &
            spec == "2 horizons" & parameter == "ATE"][order(subject), round(estimate, 4)]
  want <- fread(hf)[estimator == "CDH (headline)" & regime == "combined"][
            order(subject), round(att, 4)]
  if (!isTRUE(all.equal(hl, want)))
    stop("the main six-wave cell no longer reproduces the headline of ",
         "06_estimation_main.csv: got ", paste(hl, collapse = "/"), ", expected ",
         paste(want, collapse = "/"))
}

setorder(agg, definition, design, subject, spec, parameter)
fwrite(agg[, .(definition, design, subject, spec, parameter,
               estimate = round(estimate, 4), se = round(se, 4), switchers)],
       file.path(DIR_TABLES, "06n_calendar_estimates.csv"))

# ---- 2-3. Effect and placebo by calendar time since adoption -----------------
# Each class is estimated on the never-switchers of that cell plus the switchers
# of that class only, and effects are never pooled across classes. Switchers that
# support no horizon (their step wave is unobserved) belong to no class and are
# excluded rather than assigned to one; 06m counts them.
classes <- DES$schools[!is.na(years_since_adoption), .(switchers = .N),
                       by = .(definition, design, years_since_adoption)]
classes <- classes[switchers >= MIN_SW][order(definition, design, years_since_adoption)]

expo <- rbindlist(lapply(seq_len(nrow(classes)), function(i) {
  dn <- classes$definition[i]; nm <- classes$design[i]
  e  <- classes$years_since_adoption[i]
  ids <- DES$schools[definition == dn & design == nm &
                     (is.na(first_switch_wave) | years_since_adoption == e), school_id]
  dat <- frame(dn, nm, ids)
  rbindlist(lapply(names(SUBJECTS), function(y)
    rbindlist(lapply(c("1 horizon", "2 horizons"), function(sp)
      run(dat, dn, nm, y, sp, list(definition = dn, design = nm, years_since_adoption = e))))))
}), fill = TRUE)
if (!nrow(expo)) stop("every exposure-class estimation failed")

setorder(expo, definition, design, years_since_adoption, subject, spec, parameter)
fwrite(expo[, .(definition, design, years_since_adoption, subject, spec, parameter,
                estimate = round(estimate, 4), se = round(se, 4), switchers)],
       file.path(DIR_TABLES, "06n_exposure_estimates.csv"))

# ---- Is the aggregate placebo zero, or cancellation? -------------------------
# The check 06j made mandatory, automated: an aggregate placebo near zero can be
# an average of class placebos of opposite sign. `spread` is what distinguishes
# the two -- a genuine zero has class placebos near zero as well.
#
# The first placebo is the same pre-period difference whichever window is
# aggregated, so one row per cell and subject suffices -- which is also what lets
# the three-period calendar, that carries no placebo at two horizons, appear here
# at all. Checked rather than assumed.
KEY <- c("definition", "design", "subject")
p1 <- function(x, sp) x[parameter == "Placebo_1" & spec == sp]
inv <- merge(p1(agg, "1 horizon")[, c(KEY, "estimate"), with = FALSE],
             p1(agg, "2 horizons")[, c(KEY, "estimate"), with = FALSE],
             by = KEY, suffixes = c("_1h", "_2h"))
if (!isTRUE(all.equal(inv$estimate_1h, inv$estimate_2h)))
  stop("the first placebo moved with the aggregation window; it can no longer be ",
       "reported once per cell")

cancel <- p1(expo, "1 horizon")[
  , .(classes = .N,
      min_switchers = min(switchers),
      placebo_min = round(min(estimate), 4),
      placebo_max = round(max(estimate), 4),
      spread = round(max(estimate) - min(estimate), 4),
      placebo_weighted = round(weighted.mean(estimate, switchers), 4),
      opposite_signs = any(estimate > 0) & any(estimate < 0)),
  by = KEY]
cancel <- merge(p1(agg, "1 horizon")[, .(definition, design, subject,
                                         placebo_aggregate = round(estimate, 4),
                                         se_aggregate = round(se, 4),
                                         switchers_aggregate = switchers)],
                cancel, by = KEY)
setorder(cancel, definition, design, subject)
fwrite(cancel, file.path(DIR_TABLES, "06n_placebo_cancellation.csv"))

# ---- 4. Leave-one-cohort-out diagnostic --------------------------------------
# Reported as a set per cell, with the weight each cohort carries at each horizon.
# Excluding a cohort removes treatment variation and changes the estimand, so this
# locates where an estimate comes from; it is not a preferred specification. The
# cohorts are dated on the annual path of the SAME definition, so under the
# main treatment they are larger.
coh <- over_cells(function(dn, nm) {
  ids_all <- DES$schools[definition == dn & design == nm, school_id]
  ct <- DES$cohorts[definition == dn]
  rbindlist(lapply(c(NA_integer_, COHORTS), function(g) {
    ids <- if (is.na(g)) NULL else
      setdiff(ids_all, ct[first_treated_year == g, school_id])
    s <- if (is.na(g)) "all schools" else paste("excluding the", g, "cohort")
    dat <- frame(dn, nm, ids)
    rbindlist(lapply(names(SUBJECTS), function(y)
      run(dat, dn, nm, y, "2 horizons",
          list(definition = dn, design = nm, sample = s))))
  }))
})
if (!nrow(coh)) stop("every cohort diagnostic failed")

# The share of each horizon's switchers that belongs to each cohort: the weight
# whose change is the alternative explanation for any change in the estimate.
hc <- merge(DES$horizons[horizon <= 2L], DES$cohorts,
            by = c("definition", "school_id"), all.x = TRUE)
share <- rbindlist(lapply(COHORTS, function(g)
  hc[, .(cohort = g, pct = round(100 * mean(first_treated_year == g, na.rm = TRUE), 1)),
     by = .(definition, design, horizon)]))
share[, horizon := paste0("pct_k", horizon)]
share <- dcast(share, definition + design + cohort ~ horizon, value.var = "pct")

setorder(coh, definition, design, subject, sample, parameter)
fwrite(coh[, .(definition, design, subject, sample, parameter,
               estimate = round(estimate, 4), se = round(se, 4), switchers)],
       file.path(DIR_TABLES, "06n_cohort_diagnostic.csv"))
fwrite(share[order(definition, design, cohort)],
       file.path(DIR_TABLES, "06n_cohort_weight.csv"))
