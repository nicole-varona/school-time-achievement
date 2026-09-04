# ==============================================================================
# 04_threshold_sensitivity.R — Does the estimated effect depend on the threshold?
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# The binary treatment needs a cut-off on the expanded share, and the headline
# uses >=50%. This script re-runs the PRIMARY estimator under alternative
# cut-offs, changing nothing else, so that the estimate can be shown not to
# depend on where the line was drawn.
#
# The estimator matters here. A static TWFE has no event time, so it cannot
# detect the threat that actually applies to a staggered design: that moving the
# cut-off moves a school's ADOPTION DATE, and hence the comparisons the
# estimator is built on. The grid is therefore run with de Chaisemartin-
# D'Haultfoeuille, which is the specification the results rest on.
#
# Thresholds: >0 (the definition used by the closest antecedent), >=50%
# (headline, a majority of primary enrolment, aligned with the official
# four-hour definition of a simple day) and >=75% (the antimode of the share
# distribution). All are built like the headline, varying the cut-off alone.
#
# INPUT   data/processed/estimation_sample.rds
# OUTPUT  output/tables/04_threshold_sensitivity_cdh.csv
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
if (!requireNamespace("DIDmultiplegtDYN", quietly = TRUE))
  install.packages("DIDmultiplegtDYN", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(library(polars))   # DIDmultiplegtDYN backend

d  <- readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds"))
ra <- readRDS(file.path(DIR_PROCESSED, "ra_panel.rds"))

# The grid varies the CUT-OFF and nothing else. 0.40 and 0.60 bracket the headline
# closely enough to say whether the estimate turns on where the line sits, which
# the wider 0-and-0.75 pair cannot: they move the treated population so much that
# they answer a different question.
#
# Four cut-offs, not five: 0.50 is NOT rebuilt here, because that row IS the main
# specification and uses `main_treatment` from 02. Rebuilding it would put a
# fourth copy of the same definition in the project and let the headline row of
# this grid drift away from the headline itself without any error being raised.
#
# This block used to describe a fifth row that varied the RULE rather than the
# cut-off, separating coverage from persistence. That row went with the two-year
# persistence rule on 12/08; the grid has varied coverage alone since. `06o`, the
# script the note pointed to, is retired from the project and not in this
# repository.
# 0.50 uses
# `main_treatment` from 02. Rebuilding it would put a fourth copy of the same
# definition in the project and let the headline row of this grid drift away from
# the headline itself without any error being raised.
CUTS <- c("any (>0)" = 1e-9, "0.40" = 0.40, "0.60" = 0.60, "0.75 (data valley)" = 0.75)
for (i in seq_along(CUTS))
  d[, (paste0("cut", i)) := fifelse(is.na(share_expanded), NA_integer_,
                                    as.integer(share_expanded >= CUTS[i]))]

THRESHOLDS <- c("any (>0)"                   = "cut1",
                "0.40"                       = "cut2",
                "0.50 (headline)"            = "main_treatment",
                "0.60"                       = "cut3",
                "0.75 (data valley)"         = "cut4")
SUBJECTS   <- c(language = "y_lang", mathematics = "y_math")

# Same call as the headline in 06: two post-treatment horizons, EVERY placebo the
# design admits, clustered at the school.
#
# The placebo count is computed, not typed. This grid used to ask for one where
# the design carries two, which is the defect 06_estimation.R records and fixed
# for itself: for the headline it is the SECOND placebo that rejects, so a table
# that reports only the first clears each cut-off on half its pre-treatment
# evidence -- and this table's whole point is to pair every effect with its own
# placebo. Placebo_1 is unaffected by the request (the command computes each on
# its own support); the second is new information.
EFFECTS <- 2L
PLACEBO <- max_placebo(EFFECTS, length(WAVES))
est_one <- function(tname, yname, label) {
  set_call_seed(label)
  fit <- DIDmultiplegtDYN::did_multiplegt_dyn(
    df = as.data.frame(d[!is.na(get(tname))]), outcome = yname, group = "gid",
    time = "wave", treatment = tname, effects = EFFECTS, placebo = PLACEBO,
    cluster = "gid", graph_off = TRUE)
  # The command truncates an over-request with a message rather than failing, so
  # a placebo that was asked for and not returned must read as NA, not as row 1.
  pl   <- fit$results$Placebos
  pget <- function(i, j) if (!is.null(pl) && nrow(pl) >= i) pl[i, j] else NA_real_
  data.table(att = fit$results$ATE[1], se = fit$results$ATE[2],
             eff1 = fit$results$Effects[1, 1], eff2 = fit$results$Effects[2, 1],
             switchers1 = fit$results$Effects[1, "Switchers"],
             placebo1 = pget(1, 1), placebo1_se = pget(1, 2),
             placebo2 = pget(2, 1), placebo2_se = pget(2, 2))
}

res <- rbindlist(lapply(seq_along(THRESHOLDS), function(i)
  rbindlist(lapply(seq_along(SUBJECTS), function(j)
    cbind(threshold = names(THRESHOLDS)[i], subject = names(SUBJECTS)[j],
          est_one(THRESHOLDS[i], SUBJECTS[j],
                  paste(THRESHOLDS[i], SUBJECTS[j])))))))
num <- setdiff(names(res), c("threshold", "subject"))
res[, (num) := lapply(.SD, round, 4), .SDcols = num]
fwrite(res, file.path(DIR_TABLES, "04_threshold_sensitivity_cdh.csv"))
