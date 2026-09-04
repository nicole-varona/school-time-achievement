# ==============================================================================
# 06p_spell_length_classes.R — Does the effect depend on how long the first
#                              treated spell lasts?
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Inputs  : data/processed/estimation_sample.rds
#           data/processed/treatment_annual.rds
#           data/processed/calendar_designs.rds     (built by 06m)
# Outputs : output/tables/06p_spell_length_support.csv   who is in each class
#           output/tables/06p_spell_length_estimates.csv effect + placebo by class
#
# The treatment is defined on coverage alone: a school is treated in a year if the
# share clears the threshold that year. Duration is therefore not inside the
# variable, and the obvious question -- do effects appear only when the expansion
# LASTS? -- has to be asked separately rather than assumed away. That is what this
# script does: switchers are partitioned by the length of their FIRST treated
# spell on the annual path, and the effect is estimated within each class.
#
# The two-year rule that was the project's definition until August 2026 answered
# that question by decree instead: it deleted the schools whose expansion lasted a
# single year. Estimating within the one-year class is what that rule made
# unaskable, and it is the honest replacement for it.
#
# The class column is `first_spell`, not `persistence`: it records the length a
# spell TURNED OUT to have, which is a realised outcome, and naming it after the
# retired rule invited reading it as that rule surviving in the outputs. Renamed
# 13/08; nothing downstream reads this file.
#
# WHAT THIS IS NOT. Episode length is realised AFTER adoption, so these classes are
# not assigned: a school that keeps its extended day for three years may differ
# from one that drops it after a year in ways that also move achievement, and
# dropping out may itself respond to how the first year went. Within a class the
# estimand is well defined for that class's switchers; ACROSS classes this is a
# description of who the effect comes from, not a dose-response curve, and it
# cannot be read as the effect of making a programme last longer.
#
# NOT IN THE MANUSCRIPT. The results below are produced and committed but no
# number from them is cited: see the pending note in the working log.
# ==============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(DIDmultiplegtDYN); library(polars)})

SUBJECTS <- c(y_lang = "language", y_math = "mathematics")
DEFN     <- "main"
MIN_SW   <- 150
SPECS    <- list("1 horizon" = 1L, "2 horizons" = 2L)

d   <- as.data.table(readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds")))
DES <- readRDS(file.path(DIR_PROCESSED, "calendar_designs.rds"))
ta  <- as.data.table(readRDS(file.path(DIR_PROCESSED, "treatment_annual.rds")))

# Spell structure, derived here rather than read from a saved object: this is the
# only script that partitions on it, so the run-length encoding belongs with the
# code that uses it. Runs are computed over the years with a RECORDED treatment,
# so a year of unknown share joins rather than breaks a spell.
setorder(ta, school_id, year)
ep  <- ta[!is.na(main_treatment)][order(school_id, year), {
  r   <- rle(main_treatment)
  len <- r$lengths[r$values == 1L]
  .(episodes      = length(len),
    first_len     = if (length(len)) len[1L] else NA_integer_,
    longest       = if (length(len)) max(len) else NA_integer_,
    treated_years = sum(main_treatment))
}, by = school_id][episodes > 0]
TV  <- DES$defs[[DEFN]]

# ---- The classes -------------------------------------------------------------
# Length of the first treated episode, top-coded at three years because beyond
# that the classes thin out and the question is one year, two, or more.
ep[, first_spell := fifelse(first_len >= 3L, "3 or more years",
                            paste(first_len, fifelse(first_len == 1L, "year", "years")))]
CLASSES <- c("1 year", "2 years", "3 or more years")
if (!all(ep$first_spell %in% CLASSES))
  stop("an episode length fell outside the three classes")

# Support, per calendar: how many of each class actually reach the first horizon,
# and how their paths differ afterwards. `pct_returns` matters for reading the
# one-year class: those schools are not untreated ever after, they come back.
sup <- rbindlist(lapply(names(DES$designs), function(nm) {
  h1 <- DES$horizons[definition == DEFN & design == nm & horizon == 1L, school_id]
  z  <- ep[school_id %in% h1]
  z[, .(design = nm, switchers = .N,
        mean_episodes = round(mean(episodes), 2),
        mean_treated_years = round(mean(treated_years), 2),
        pct_returns = round(100 * mean(episodes > 1), 1)), by = first_spell]
}))
sup[, pct := round(100 * switchers / sum(switchers), 1), by = design]
setorder(sup, design, first_spell)
fwrite(sup, file.path(DIR_TABLES, "06p_spell_length_support.csv"))

# ---- Estimation --------------------------------------------------------------
# Each class on the never-switchers of that calendar plus the switchers of that
# class only, exactly as 06n partitions by exposure length, and never pooled.
frame <- function(design, ids) {
  yrs <- DES$designs[[design]]
  x <- d[year %in% yrs & !is.na(get(TV)) & school_id %in% ids]
  x[, wv := match(year, yrs)]
  as.data.frame(x)
}

tidy <- function(fit) {
  m2dt <- function(m) if (is.null(m)) NULL else
    data.table(parameter = trimws(rownames(m)), estimate = m[, "Estimate"],
               se = m[, "SE"], switchers = m[, "Switchers"])
  rbind(m2dt(fit$results$Effects), m2dt(fit$results$Placebos),
        data.table(parameter = "ATE", estimate = fit$results$ATE[1],
                   se = fit$results$ATE[2], switchers = NA_real_))
}

grid <- sup[switchers >= MIN_SW, .(design, first_spell)]
res <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  nm <- grid$design[i]; cl <- grid$first_spell[i]
  keep <- ep[first_spell == cl, school_id]
  ids  <- DES$schools[definition == DEFN & design == nm &
                      (is.na(first_switch_wave) | school_id %in% keep), school_id]
  dat  <- frame(nm, ids); k <- length(DES$designs[[nm]])
  rbindlist(lapply(names(SUBJECTS), function(y)
    rbindlist(lapply(names(SPECS), function(sp) {
      label <- paste(nm, cl, y, sp)
      set_call_seed(label)
      fit <- tryCatch(did_multiplegt_dyn(
        df = dat, outcome = y, group = "gid", time = "wv", treatment = TV,
        effects = SPECS[[sp]], placebo = max_placebo(SPECS[[sp]], k),
        graph_off = TRUE), error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      cbind(design = nm, first_spell = cl, subject = unname(SUBJECTS[y]), spec = sp,
            tidy(fit))
    }))))
}), fill = TRUE)
if (!nrow(res)) stop("every first_spell-class estimation failed")

setorder(res, design, first_spell, subject, spec, parameter)
fwrite(res[, .(design, first_spell, subject, spec, parameter,
               estimate = round(estimate, 4), se = round(se, 4), switchers)],
       file.path(DIR_TABLES, "06p_spell_length_estimates.csv"))
