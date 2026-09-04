# ==============================================================================
# 05_reversion_check.R — Is the treatment really non-absorbing? (visual check)
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Purpose : Show, in a table and figures, how the raw treatment counts
#           (n_extended_day, n_full_day) behave over 2011-2025, so the
#           non-absorbing pattern can be seen rather than asserted.
# Inputs  : data/processed/ra_panel.rds
# Outputs : output/tables/05_transitions.csv
#           output/figures/05_national_totals.png
#           output/figures/05_school_trajectories.png
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
ra <- readRDS(file.path(DIR_PROCESSED, "ra_panel.rds"))
setorder(ra, school_id, year)
# Diagnostic figures drawn here use the shared style (90_style.R) so they match.
source(here::here("scripts", "90_style.R"))
theme_set(theme_thesis())

# ---- Year-over-year transitions ---------------------------------------------
# On the PRIMARY treatment, read from 02 rather than rebuilt here: deriving it
# twice is how two scripts end up measuring the same thing differently. 01 used to
# carry a `treated` column of its own, which collided by name with this merge;
# that derivation was removed rather than worked around.
#
# This comment used to say the opposite -- that the indicator here was
# "deliberately the raw one, NOT the sustained definition the estimation uses" --
# which stopped being true when the two-year persistence rule was retired on
# 12/08 and the merge below started reading `main_treatment`. There is one
# indicator now, and this is it. The 2024 exit spike the old note discounted as
# an artefact of the raw series is therefore a feature of the treatment the
# project estimates: it is driven by schools recorded as extended in 2023 and
# 2025 but not 2024, four-fifths of them in one province running a five-hour
# programme.
ra <- merge(ra, readRDS(file.path(DIR_PROCESSED, "treatment_annual.rds"))[
              , .(school_id, year, treated = main_treatment)],
            by = c("school_id", "year"), all.x = TRUE)
setorder(ra, school_id, year)   # merge reorders, and shift() below depends on it
ra[, treated_prev := shift(treated), by = school_id]
ra[, year_prev := shift(year), by = school_id]
transitions <- ra[!is.na(treated) & !is.na(treated_prev) & year_prev == year - 1,
                  .(adopt_0to1    = sum(treated_prev == 0 & treated == 1),
                    revert_1to0   = sum(treated_prev == 1 & treated == 0),
                    stay_treated  = sum(treated_prev == 1 & treated == 1),
                    base_treated_prev = sum(treated_prev == 1)),
                  by = year][order(year)]
transitions[, pct_revert := round(100 * revert_1to0 / base_treated_prev, 1)]
fwrite(transitions, file.path(DIR_TABLES, "05_transitions.csv"))

# ---- Figure A: national totals of extended / full-day students by year ------
nat <- ra[!is.na(enroll_primary) & enroll_primary > 0,
          .(extended = sum(n_extended_day, na.rm = TRUE),
            full     = sum(n_full_day,     na.rm = TRUE)), by = year][order(year)]
natL <- melt(nat, id.vars = "year", variable.name = "regime", value.name = "students")
natL[, regime := factor(regime, labels = c("extended day", "full day"))]
pA <- ggplot(natL, aes(year, students / 1000, colour = regime)) +
  geom_vline(xintercept = 2022.5, linetype = "dashed", colour = "grey60") +
  annotate("text", x = 2022.5, y = Inf, label = "1h+ (2023)", vjust = 1.5,
           hjust = 1.05, size = 3, colour = "grey40") +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  labs(caption = paste("National totals: primary students in extended vs full-day, RA 2011-2025.",
                       "Full-day grows smoothly; extended day is stable then jumps in 2023."),
       x = NULL, y = "students (thousands)", colour = NULL)
save_fig(pA, "05_national_totals.png", width = 9, height = 5)

# ---- Figure B: raw trajectories of example schools, by archetype ------------
traj <- ra[!is.na(treated), .(
  years_obs   = .N,
  ever        = any(treated == 1),
  first_treat = if (any(treated == 1)) min(year[treated == 1]) else NA_integer_,
  n_switches  = sum(abs(diff(treated))),
  # Where the school ENDS. A switch count alone cannot tell an exit from a
  # re-entry: `n_switches == 2` is satisfied both by 1-0-1 and by 0-1-0, and the
  # archetypes below turn on which of the two it is. `ra` is ordered by
  # school_id, year above, so the last element is the last census year observed.
  last_treated = treated[.N]), by = school_id]

# The only randomness in the project that is not a bootstrap: which schools get
# drawn as illustrations. Seeded through the same helper as everything else so
# there is one idiom and no loose constant -- it used to be a bare set.seed(42),
# the one seed in the pipeline that did not derive from SEED.
set_call_seed("trajectory archetypes")
pick <- function(cond, k = 3) {
  ids <- traj[eval(cond), school_id]
  if (length(ids) == 0) return(character(0))
  sample(ids, min(k, length(ids)))
}
# Each archetype is defined by the SHAPE of the path, not by the switch count
# alone. The exit and the adopter both have exactly one switch and are told apart
# by where they end; `n_switches == 2` was previously labelled "one-way exit" and
# drew re-entrants (1-0-1) alongside genuine exits, which contradicted the
# manuscript's own caption -- it says these schools "cross the line and fall
# back", and a school that ends treated does not.
sel <- rbind(
  data.table(school_id = pick(quote(ever & n_switches == 0)),               type = "1. Always treated"),
  data.table(school_id = pick(quote(ever & n_switches == 1 & last_treated == 1 &
                                     first_treat %between% c(2014, 2019))),  type = "2. Clean adopter"),
  data.table(school_id = pick(quote(ever & n_switches == 1 & last_treated == 0)),
                                                                             type = "3. One-way exit"),
  data.table(school_id = pick(quote(n_switches >= 4)),                       type = "4. Oscillator")
)
if (nrow(sel) != 12L)
  stop("expected three illustrative schools per archetype, got ", nrow(sel))
# The fictitious identifier is not informative to a reader and competes with the
# archetype for attention, so the panels are labelled by what they illustrate.
# One panel per school, three per archetype, so each row of the grid is one
# archetype and no two lines share a panel. The fictitious identifier is dropped:
# it is not informative to a reader and competes with the archetype for attention.
sel[, label := paste0(type, "\nillustrative school ", seq_len(.N)), by = type]

# The SHARE, not the counts. Counts put every panel on its own axis, which makes
# the archetypes incomparable and hides the treatment: the threshold of @eq-treat
# is defined on the share, so only the share shows a school crossing it and
# falling back. One line per school -- D_gt is defined on extended and full-day
# enrolment combined, which `share_expanded` already is (built once in 02).
traj_data <- merge(ra[school_id %in% sel$school_id,
                      .(school_id, year, share_expanded)],
                   sel, by = "school_id")

pB <- ggplot(traj_data, aes(year, share_expanded, group = school_id)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 0.8, colour = PAL$ink) +
  geom_point(size = 1.2, colour = PAL$ink) +
  facet_wrap(~ label, ncol = 3) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  # No title or subtitle: this figure enters the manuscript, where the caption and
  # the note below the panel carry what these said. theme_thesis() blanked them.
  labs(x = NULL, y = "expanded share of primary enrolment") +
  theme(strip.text = element_text(size = 8))
save_fig(pB, "05_school_trajectories.png", width = 11, height = 7)

