# ==============================================================================
# 06k_leave_one_province.R — Is the headline carried by one jurisdiction?
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Inputs  : data/processed/estimation_sample.rds
# Outputs : output/tables/06k_leave_one_province.csv
#
# Adoption is coordinated provincially and the 2023 expansion is concentrated in
# the Northeast, so the natural question is whether the headline rests on a single
# jurisdiction. This drops one province at a time and re-estimates the primary
# specification.
#
# These are POINT ESTIMATES, not an inference procedure: with 24 jurisdictions and
# roughly seven effective clusters, a fourth set of p-values would add nothing.
# What the exercise answers is a leverage question -- how far the estimate moves
# when any one province is removed -- and that reads off the points alone.
# ==============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(DIDmultiplegtDYN); library(polars)})

SUBJECTS <- c(y_lang = "language", y_math = "mathematics")

d <- as.data.table(readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds")))
provinces <- sort(unique(d$province[!is.na(d$province)]))

run_drop <- function(prov, yname) {
  dd <- if (is.na(prov)) d else d[province != prov]
  set_call_seed(paste("lopo", prov, yname))
  fit <- tryCatch(did_multiplegt_dyn(
    df = as.data.frame(dd), outcome = yname, group = "gid", time = "wave",
    treatment = "main_treatment", effects = 2, placebo = 1,
    cluster = "gid", graph_off = TRUE), error = function(e) NULL)
  if (is.null(fit)) return(data.table(dropped = prov, subject = unname(SUBJECTS[yname]),
                                      att = NA_real_, se = NA_real_, schools = NA_integer_))
  data.table(dropped = prov, subject = unname(SUBJECTS[yname]),
             att = fit$results$ATE[1], se = fit$results$ATE[2],
             schools = uniqueN(dd$school_id))
}

grid <- CJ(prov = c(NA_character_, provinces), yname = names(SUBJECTS), sorted = FALSE)
out <- rbindlist(lapply(seq_len(nrow(grid)), function(i)
  run_drop(grid$prov[i], grid$yname[i])))

# A silent tryCatch turning every cell into NA is indistinguishable from a run
# that produced nothing, so the failure is made loud here.
if (all(is.na(out$att))) stop("every leave-one-province estimation failed")

out[, dropped := fifelse(is.na(dropped), "none (full sample)", dropped)]
ref <- out[dropped == "none (full sample)", .(subject, att_full = att)]
out <- merge(out, ref, by = "subject")
out[, shift := round(att - att_full, 4)]
out[, `:=`(att = round(att, 4), se = round(se, 4), att_full = NULL)]
setorder(out, subject, -shift)
fwrite(out, file.path(DIR_TABLES, "06k_leave_one_province.csv"))
