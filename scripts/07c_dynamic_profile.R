# ==============================================================================
# 07c_dynamic_profile.R — Shape of the dynamic treatment response
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# The static WAS is ~0 while the dynamic effects are positive, suggesting the
# effect emerges with exposure. Before leaning on that "exposure dynamics" story,
# we must look at the SHAPE of the post-treatment profile: does it rise (roughly)
# monotonically / plateau (a credible accumulation), or is it jumpy (be cautious)?
# We plot ATT(e) for all event-times with CIs, from Callaway-Sant'Anna (dynamic
# aggregation) for Language and Mathematics. NOTE: event-time is in OUTCOME WAVES,
# not equal calendar years (irregular spacing) — interpret accordingly.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(did); library(ggplot2); library(polars)})
# Seeded per call inside dyn_profile(), not once here: see 00_setup.R.
theme_set(theme_bw(base_size = 11))

panel <- readRDS(file.path(DIR_PROCESSED, "panel.rds"))
# Treatment comes from 02 (the single source of truth) instead of being rebuilt
# here. The persistence rule used to live in nine copies, and shift() silently
# depends on the panel being sorted -- an invariant no script declared.
panel <- merge(panel,
               readRDS(file.path(DIR_PROCESSED, "treatment_annual.rds"))[
                 , .(school_id, year, main_treatment)],
               by = c("school_id", "year"), all.x = TRUE)
d <- panel[year %in% WAVES & !is.na(score_lang),
           .(school_id, year, y_lang = score_lang / 100, y_math = score_math / 100, main_treatment)]
d <- d[!is.na(main_treatment)]
d[, `:=`(wave = match(year, WAVES), gid = as.integer(factor(school_id)))]
d_cs <- cs_sample(d)

dyn_profile <- function(yname, subject) {
  set_call_seed(paste("dynamic profile", subject))
  fit <- att_gt(yname = yname, tname = "wave", idname = "gid", gname = "G",
                data = d_cs[!is.na(get(yname))], control_group = "notyettreated",
                allow_unbalanced_panel = TRUE, bstrap = TRUE, biters = BITERS, est_method = "dr")
  dyn <- aggte(fit, type = "dynamic", na.rm = TRUE)
  data.table(k = dyn$egt, att = dyn$att.egt, se = dyn$se.egt,
             lo = dyn$att.egt - 1.96 * dyn$se.egt, hi = dyn$att.egt + 1.96 * dyn$se.egt,
             subject = subject)
}
prof <- rbindlist(list(dyn_profile("y_lang", "Language"), dyn_profile("y_math", "Mathematics")))
prof_r <- copy(prof)[, c("att","se","lo","hi") := lapply(.SD, round, 4), .SDcols = c("att","se","lo","hi")]
fwrite(prof_r, file.path(DIR_TABLES, "07c_dynamic_profile.csv"))

# Figure 07c_dynamic_profile.png is rendered from 07c_dynamic_profile.csv by
# 91_figures.R in the shared style.
