# ==============================================================================
# 06j_exposure_partition.R — Event time in years, not in waves
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Inputs  : data/processed/estimation_sample.rds
#           data/processed/treatment_annual.rds   (treatment on the annual grid)
#           data/processed/ra_panel.rds           (raw annual share)
# Outputs : output/tables/06j_exposure_classes.csv    switchers per class
#           output/tables/06j_exposure_partition.csv  ATT by class
#
# `did_multiplegt_dyn` assumes evenly spaced periods. The Aprender calendar is not
# (2016, 2018, 2021, 2022, 2023, 2025: gaps of 2, 3, 1, 1, 2 years), while the
# treatment is known ANNUALLY. Indexing event time by wave therefore lets a single
# `Effect_k` average materially different exposure lengths: measured here, the
# first effect mixes 0 to 8 years, with 69% at the mode. (The share is written to
# 06j_exposure_classes.csv and the manuscript reads it from there rather than
# from this comment, which said 58% until 13/08 -- a figure left over from the
# treatment definition retired on 12/08.)
#
# The package documents the fix for outcomes observed less often than treatment:
# partition switchers by where the first change falls relative to the outcome
# dates, run the command once per partition, and use the annual information ONLY
# to select the subsamples -- the estimation stays at the outcome frequency. This
# script applies that logic to an irregular calendar, where the partition is by
# EXPOSURE LENGTH rather than by odd/even period. Each partition is estimated on
# all never-switcher school-waves plus the switchers of that class, and effects
# are not pooled across partitions.
#
# ONE TREATMENT DEFINITION is run, the project's primary one. Two were, until the
# two-year persistence rule was retired on 12/08: the classes are built from the
# date a school first adopts, and the two rules dated adoption differently, so the
# comparison was worth making while both existed. It no longer is.
#
# The `definition` column survives in both outputs with a single value, "main".
# It is kept rather than dropped because the manuscript filters on it in three
# places, and because a column that names which definition produced a row is the
# right shape for this table if a second one ever returns.
# ==============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(DIDmultiplegtDYN); library(polars)})

SUBJECTS <- c(y_lang = "language", y_math = "mathematics")

d  <- as.data.table(readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds")))
ta <- as.data.table(readRDS(file.path(DIR_PROCESSED, "treatment_annual.rds")))
ra <- as.data.table(readRDS(file.path(DIR_PROCESSED, "ra_panel.rds")))
setorder(d, school_id, wave)

# The treatment comes from 02. This script used to rebuild the main
# one from the share, which is the duplication 02 exists to prevent: the same
# quantity was being derived here, in 04 and in 05, under a third name (`treat_raw`, now `main_treatment`)
# and with its own missing-data handling.
DEFS <- list(
  main = list(wave = "main_treatment",
                         annual = ta[!is.na(main_treatment), .(school_id, year, tr = main_treatment)]))

# ---- Exposure classes for one treatment definition ---------------------------
# F is the first OBSERVED wave whose treatment differs from the school's first
# observed wave; the switch year is the first ANNUAL year, from that wave onwards,
# at which the annual series differs. Exposure at horizon 1 is the gap between.
build_classes <- function(tv, annual) {
  x <- d[!is.na(get(tv))]
  fs <- x[, {
    ch <- which(get(tv) != get(tv)[1L])
    .(F = if (length(ch)) wave[ch[1L]] else NA_integer_,
      d1 = get(tv)[1L], y0 = WAVES[wave[1L]])
  }, by = school_id]
  sw <- fs[!is.na(F)]
  yr <- merge(annual[!is.na(tr)], sw[, .(school_id, d1, y0)], by = "school_id")
  chg <- yr[year >= y0][, {
    ch <- which(tr != d1[1L])
    .(switch_year = if (length(ch)) year[ch[1L]] else NA_integer_)
  }, by = school_id]
  sw <- merge(sw, chg, by = "school_id")[!is.na(switch_year)]
  sw[, exposure := WAVES[F] - switch_year]
  list(switchers = sw, never = fs[is.na(F), school_id])
}

MIN_SW <- 150
SPECS <- list("1 horizon"                  = list(effects = 1, same = FALSE),
              "2 horizons"                 = list(effects = 2, same = FALSE),
              "3 horizons"                 = list(effects = 3, same = FALSE),
              "2 horizons, same switchers" = list(effects = 2, same = TRUE))

run_class <- function(ids, tv, e, yname, spec, def) {
  dd <- as.data.frame(d[school_id %in% ids & !is.na(get(tv))])
  set_call_seed(paste(def, e, yname, spec))
  # Every placebo the window admits, not the first. These classes are the exhibit
  # that pairs each effect with its own pre-treatment coefficient, and for the
  # headline it is the SECOND placebo that rejects, so asking for one would judge
  # each class on half its pre-treatment evidence. max_placebo() also keeps the
  # request inside what the command will honour instead of letting it truncate.
  eff <- SPECS[[spec]]$effects
  args <- list(df = dd, outcome = yname, group = "gid", time = "wave",
               treatment = tv, effects = eff,
               placebo = max_placebo(eff, length(WAVES)),
               cluster = "gid", graph_off = TRUE)
  if (isTRUE(SPECS[[spec]]$same))
    args <- c(args, list(same_switchers = TRUE, same_switchers_pl = TRUE))
  fit <- tryCatch(do.call(did_multiplegt_dyn, args), error = function(err) NULL)
  if (is.null(fit)) return(NULL)
  rows <- rbind(
    data.table(parameter = rownames(fit$results$Effects),
               estimate = fit$results$Effects[, "Estimate"],
               se = fit$results$Effects[, "SE"],
               switchers = fit$results$Effects[, "Switchers"]),
    data.table(parameter = "ATE", estimate = fit$results$ATE[1],
               se = fit$results$ATE[2], switchers = NA_real_))
  pl <- fit$results$Placebos
  if (!is.null(pl) && nrow(pl))
    rows <- rbind(rows, data.table(parameter = paste0("Placebo_", seq_len(nrow(pl))),
                                   estimate  = pl[, "Estimate"],
                                   se        = pl[, "SE"],
                                   switchers = pl[, "Switchers"]))
  rows[, `:=`(definition = def, exposure_years = e,
              subject = unname(SUBJECTS[yname]), spec = spec)]
  rows[]
}

cls <- list(); res <- list()
for (def in names(DEFS)) {
  b  <- build_classes(DEFS[[def]]$wave, DEFS[[def]]$annual)
  tb <- b$switchers[, .(switchers = .N), by = exposure][order(exposure)]
  tb[, `:=`(definition = def, pct = round(100 * switchers / sum(switchers), 1))]
  cls[[def]] <- tb
  for (e in tb[switchers >= MIN_SW, exposure]) {
    ids <- c(b$never, b$switchers[exposure == e, school_id])
    for (yname in names(SUBJECTS)) for (sp in names(SPECS))
      res[[length(res) + 1]] <- run_class(ids, DEFS[[def]]$wave, e, yname, sp, def)
  }
}

classes <- rbindlist(cls)
setcolorder(classes, c("definition", "exposure", "switchers", "pct"))
fwrite(classes, file.path(DIR_TABLES, "06j_exposure_classes.csv"))

out <- rbindlist(res, fill = TRUE)
# A silent tryCatch returning nothing is indistinguishable from a run that found
# nothing, so the failure is made loud.
if (!nrow(out)) stop("every partitioned estimation failed")
out[, parameter := trimws(parameter)]
setcolorder(out, c("definition", "exposure_years", "spec", "subject", "parameter",
                   "estimate", "se", "switchers"))
out[, `:=`(estimate = round(estimate, 4), se = round(se, 4))]
setorder(out, definition, exposure_years, spec, subject, parameter)
fwrite(out, file.path(DIR_TABLES, "06j_exposure_partition.csv"))
