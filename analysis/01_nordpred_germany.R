# 01_nordpred_germany.R
# Project lip and oral cavity cancer burden in Germany to 2050 using Nordpred
# (Møller et al. 2003; Chen et al. Cancer Medicine 2025 used same approach)
#
# Inputs:  GBD 2023 Results extract for Germany, 1990-2023, 5-year age groups, all sexes
# Output:  Forecast tables + Figure 1 (historical + projection)
#
# Periods (5-year): 1990-94, 1995-99, 2000-04, 2005-09, 2010-14, 2015-19
# Forecast horizon: 2020-24, 2025-29, 2030-34, 2035-39, 2040-44, 2045-49, 2050-54
# Observed 2020-23 data is held back as a visual validation overlay.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(nordpred)
})

# ---- Paths ----------------------------------------------------------------
DATA_DIR <- "."  # repository root
CAUSE_CSV <- file.path(DATA_DIR, "data/gbd_cause",
                      "IHME-GBD_2023_DATA-684738bc-1.csv")
OUT_DIR  <- file.path(DATA_DIR, "output")
FIG_DIR  <- file.path(OUT_DIR, "figures")
TBL_DIR  <- file.path(OUT_DIR, "tables")
DATA_OUT <- file.path(OUT_DIR, "data")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT, recursive = TRUE, showWarnings = FALSE)

# ---- Constants ------------------------------------------------------------
AGE_5YR_RAW <- c("<5 years","5-9 years","10-14 years","15-19 years","20-24 years",
                 "25-29 years","30-34 years","35-39 years","40-44 years","45-49 years",
                 "50-54 years","55-59 years","60-64 years","65-69 years","70-74 years",
                 "75-79 years","80-84 years","85-89 years","90-94 years","95+ years")
# Nordpred expects exactly 18 5-year age groups with 85+ pooled as final
AGE_5YR <- c("<5 years","5-9 years","10-14 years","15-19 years","20-24 years",
             "25-29 years","30-34 years","35-39 years","40-44 years","45-49 years",
             "50-54 years","55-59 years","60-64 years","65-69 years","70-74 years",
             "75-79 years","80-84 years","85+ years")
# Helper: collapse 85-89/90-94/95+ -> 85+
collapse_85plus <- function(x_named) {
  out <- x_named[setdiff(names(x_named), c("85-89 years","90-94 years","95+ years"))]
  out["85+ years"] <- sum(x_named[c("85-89 years","90-94 years","95+ years")],
                          na.rm = TRUE)
  out[AGE_5YR]
}

# Historical periods (six exact 5-year periods ending 2023; first 4 years of
# GBD data 1990-1993 dropped to keep equal-length periods, which is the Nordpred
# convention).
HIST_PERIODS <- list(c(1994,1998), c(1999,2003), c(2004,2008),
                     c(2009,2013), c(2014,2018), c(2019,2023))
HIST_PERIOD_LBL <- c("1994-1998","1999-2003","2004-2008",
                     "2009-2013","2014-2018","2019-2023")
HIST_PERIOD_MID <- c(1996, 2001, 2006, 2011, 2016, 2021)

# Project 4 future 5-year periods to align with Chen et al. Cancer Medicine
# 2025 methodology (projections to 2040 / 2041 period midpoint). This avoids
# extrapolating beyond the established literature horizon for Nordpred-based
# cancer projections.
FCST_PERIOD_LBL <- c("2024-2028","2029-2033","2034-2038","2039-2043")
FCST_PERIOD_MID <- c(2026, 2031, 2036, 2041)

NPROJ <- length(FCST_PERIOD_LBL)

# Starting age group for fitting (index): 25-29 = position 6 in AGE_5YR (oral
# cancer cases below ~25 years are negligible)
START_EST_AGE <- 6
START_USE_AGE <- 6

# ---- Load data ------------------------------------------------------------
cat("Loading cause data...\n")
d <- fread(CAUSE_CSV)
cat(sprintf("  rows: %d\n", nrow(d)))

# Filter to Germany, 5-year age groups (raw 20)
dg <- d[location_name == "Germany" & age_name %in% AGE_5YR_RAW]
cat(sprintf("  Germany 5-year-age rows: %d\n", nrow(dg)))

# Collapse 85-89, 90-94, 95+ into a single "85+ years" bin for Nordpred
dg[, age_bin := fifelse(age_name %in% c("85-89 years","90-94 years","95+ years"),
                       "85+ years", age_name)]
dg_n <- dg[metric_name == "Number",
           .(val = sum(val, na.rm = TRUE)),
           by = .(measure_name, sex_name, year, age_bin)]
dg_n[, metric_name := "Number"]
# For rates, we cannot simply sum — recompute by aggregating cases and
# back-deriving population from the original (count / rate * 1e5) for each
# raw 5-year bin, then summing populations and dividing.
dg_n_raw <- dg[metric_name == "Number", .(measure_name, sex_name, year, age_name, n = val)]
dg_r_raw <- dg[metric_name == "Rate",   .(measure_name, sex_name, year, age_name, r = val)]
dg_pop   <- merge(dg_n_raw, dg_r_raw, by = c("measure_name","sex_name","year","age_name"))
dg_pop[, pop := ifelse(r > 0, n / r * 1e5, NA_real_)]
dg_pop[, age_bin := fifelse(age_name %in% c("85-89 years","90-94 years","95+ years"),
                            "85+ years", age_name)]
dg_pop_agg <- dg_pop[, .(n = sum(n, na.rm = TRUE),
                         pop = sum(pop, na.rm = TRUE)),
                     by = .(measure_name, sex_name, year, age_bin)]
dg_pop_agg[, rate_per_100k := ifelse(pop > 0, n / pop * 1e5, 0)]

# Reshape so build_matrices() can use a uniform interface:
# We replace dg with the aggregated frame containing n, pop, rate.
dg <- dg_pop_agg[, .(measure_name, sex_name, age_name = age_bin, year,
                     cases_n = n, population = pop, rate_per_100k)]
dg[, age_name := factor(age_name, levels = AGE_5YR)]

# ---- Helper: build period-aggregated cases & pyr matrices -----------------
# Operates on the already-aggregated `dg` frame (with cases_n, population).
build_matrices <- function(dat, sex, measure) {
  df <- dat[sex_name == sex & measure_name == measure,
            .(year, age_name, cases = cases_n, pyr = population)]

  df[, period := NA_character_]
  for (i in seq_along(HIST_PERIODS)) {
    rng <- HIST_PERIODS[[i]]
    df[year >= rng[1] & year <= rng[2], period := HIST_PERIOD_LBL[i]]
  }
  agg <- df[!is.na(period),
            .(cases = sum(cases, na.rm = TRUE),
              pyr   = sum(pyr,   na.rm = TRUE)),
            by = .(age_name, period)]

  cases_mat <- dcast(agg, age_name ~ period, value.var = "cases", drop = FALSE)
  cases_mat <- as.matrix(cases_mat[, HIST_PERIOD_LBL, with = FALSE])
  rownames(cases_mat) <- AGE_5YR

  pyr_mat <- dcast(agg, age_name ~ period, value.var = "pyr", drop = FALSE)
  pyr_mat <- as.matrix(pyr_mat[, HIST_PERIOD_LBL, with = FALSE])
  rownames(pyr_mat) <- AGE_5YR

  cases_mat[is.na(cases_mat)] <- 0
  pyr_mat[is.na(pyr_mat) | pyr_mat <= 0] <- 1

  # GBD provides continuous mean estimates; Nordpred's Poisson GLM expects
  # integer counts.
  cases_mat <- round(cases_mat)
  storage.mode(cases_mat) <- "integer"

  # Extend pyr with future periods (static population assumption — refinable)
  last_pop <- pyr_mat[, ncol(pyr_mat), drop = FALSE]
  future_pyr <- matrix(rep(last_pop, NPROJ), nrow = nrow(pyr_mat), ncol = NPROJ)
  colnames(future_pyr) <- FCST_PERIOD_LBL
  pyr_mat_full <- cbind(pyr_mat, future_pyr)

  # CRITICAL: Nordpred's internal code (nordpred.prediction) uses
  # data.frame-style indexing (`x[, 1]`); matrices drop dimensions on single-
  # row subset, which breaks the function. Convert both to data.frames.
  cases_df <- as.data.frame(cases_mat)
  pyr_df   <- as.data.frame(pyr_mat_full)
  rownames(cases_df) <- AGE_5YR
  rownames(pyr_df)   <- AGE_5YR

  list(cases = cases_df, pyr = pyr_df)
}

# ---- Run nordpred for one configuration -----------------------------------
run_nordpred <- function(mats, label) {
  cat(sprintf("\n  nordpred for %s ... ", label))
  fit <- nordpred(cases       = mats$cases,
                  pyr         = mats$pyr,
                  startestage = START_EST_AGE,
                  startuseage = START_USE_AGE,
                  noperiods   = NPROJ,
                  cuttrend    = c(0, 0.25, 0.5, 0.75),
                  linkfunc    = "power5")
  cat("done\n")
  fit
}

# ---- Pull projected age-specific rates and crude/total counts -------------
extract_projection <- function(fit, mats, label) {
  # nordpred.getpred(fit, incidence=TRUE, standpop=NULL) returns rates per 100k
  rates_all <- nordpred.getpred(fit, incidence = TRUE, standpop = NULL)
  # rates_all is matrix age × (hist+future periods)
  # The columns are HIST_PERIOD_LBL followed by future period midpoints

  # Get counts via age-specific rate × age-specific pyr from the LAST historical
  # period (assume stable population structure forward). nordpred internally
  # uses the most recent person-years for projection.
  # Simpler: use nordpred.getpred(..., incidence=FALSE) to get counts directly.
  counts_all <- nordpred.getpred(fit, incidence = FALSE, standpop = NULL)

  list(rates = rates_all, counts = counts_all)
}

# ---- Compute age-standardized rate (ASR) using GBD world standard --------
# GBD 2023 world standard age weights (proportion of population in each 5-year
# age group). Source: GBD methods appendix. These approximate values are used
# to back-out ASRs from age-specific rates.
# 18-group ASW (85+ pooled): 0.0044 + 0.0015 + 0.0001 = 0.0060
GBD_ASW <- c(`<5 years`=0.0886, `5-9 years`=0.0869, `10-14 years`=0.0860,
             `15-19 years`=0.0847, `20-24 years`=0.0822, `25-29 years`=0.0793,
             `30-34 years`=0.0761, `35-39 years`=0.0715, `40-44 years`=0.0659,
             `45-49 years`=0.0604, `50-54 years`=0.0537, `55-59 years`=0.0455,
             `60-64 years`=0.0372, `65-69 years`=0.0296, `70-74 years`=0.0221,
             `75-79 years`=0.0152, `80-84 years`=0.0091, `85+ years`=0.0060)
GBD_ASW <- GBD_ASW[AGE_5YR]

compute_asr <- function(rates_mat) {
  # rates_mat: rows = age (in AGE_5YR order), cols = period
  # Reweight by GBD ASW
  asr <- colSums(rates_mat * GBD_ASW, na.rm = TRUE)
  asr
}

# ---- MAIN LOOP across sex × measure --------------------------------------
results <- list()
combos <- expand.grid(sex     = c("Both","Male","Female"),
                      measure = c("Deaths","DALYs (Disability-Adjusted Life Years)",
                                  "Incidence"),
                      stringsAsFactors = FALSE)

for (i in seq_len(nrow(combos))) {
  sex_i <- combos$sex[i]
  meas_i <- combos$measure[i]
  lbl <- paste(sex_i, meas_i, sep = " | ")

  mats <- build_matrices(dg, sex_i, meas_i)
  fit  <- run_nordpred(mats, lbl)
  proj <- extract_projection(fit, mats, lbl)

  asr_all <- compute_asr(proj$rates)
  total_counts <- colSums(proj$counts, na.rm = TRUE)

  period_lbl <- c(HIST_PERIOD_LBL, FCST_PERIOD_LBL)
  period_mid <- c(HIST_PERIOD_MID, FCST_PERIOD_MID)

  results[[lbl]] <- data.table(
    sex            = sex_i,
    measure        = meas_i,
    period         = period_lbl,
    period_mid     = period_mid,
    period_type    = ifelse(period_mid <= 2019, "historical", "forecast"),
    asr_per_100k   = as.numeric(asr_all),
    total_count    = as.numeric(total_counts)
  )
}

projections <- rbindlist(results)
fwrite(projections, file.path(DATA_OUT, "nordpred_projections.csv"))
cat("\nProjection table saved to:", file.path(DATA_OUT, "nordpred_projections.csv"), "\n")

# ---- Build the historical "observed" series straight from GBD for plot ----
# Pull GBD's own age-standardized rate from the original `d` frame (before our
# aggregation collapsed the metric_name column).
obs <- d[location_name == "Germany" &
         age_name == "Age-standardized" & metric_name == "Rate" &
         measure_name %in% c("Deaths","DALYs (Disability-Adjusted Life Years)","Incidence") &
         sex_name %in% c("Both","Male","Female"),
         .(sex = sex_name, measure = measure_name, year, asr = val,
           asr_lo = lower, asr_hi = upper)]
fwrite(obs, file.path(DATA_OUT, "germany_asr_observed_1990_2023.csv"))

# ---- Figure 1: Historical + Forecast for Both sexes ----------------------
# Two panels: Deaths and DALYs (the most policy-relevant); skip Incidence here
# (will go to Figure 2 / Appendix).
both_obs  <- obs[sex == "Both" & measure %in% c("Deaths","DALYs (Disability-Adjusted Life Years)")]
both_proj <- projections[sex == "Both" &
                          measure %in% c("Deaths","DALYs (Disability-Adjusted Life Years)") &
                          period_type == "forecast"]

# Pretty labels
nice_measure <- function(x) {
  ifelse(x == "Deaths", "Deaths",
         ifelse(x == "DALYs (Disability-Adjusted Life Years)", "DALYs",
                ifelse(x == "Incidence", "Incident cases", x)))
}
both_obs[,  measure_lbl := nice_measure(measure)]
both_proj[, measure_lbl := nice_measure(measure)]

# DPI: PNG ≤ 1900px longest side. For 18 cm × 12 cm: dpi ≤ 1900*2.54/18 ≈ 268
# We use 250 dpi for PNG; PDF is vector and unrestricted.
PNG_DPI <- 250
W_CM <- 18
H_CM <- 10

p1 <- ggplot() +
  geom_line(data = both_obs,
            aes(x = year, y = asr, colour = "Historical (GBD 2023)"), linewidth = 0.7) +
  geom_ribbon(data = both_obs,
              aes(x = year, ymin = asr_lo, ymax = asr_hi),
              fill = "grey70", alpha = 0.3) +
  geom_line(data = both_proj,
            aes(x = period_mid, y = asr_per_100k,
                colour = "Nordpred projection (this study)"),
            linewidth = 0.7, linetype = "dashed") +
  geom_point(data = both_proj,
             aes(x = period_mid, y = asr_per_100k),
             colour = "#1f77b4", size = 1.5) +
  facet_wrap(~ measure_lbl, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c("Historical (GBD 2023)" = "black",
                                 "Nordpred projection (this study)" = "#1f77b4")) +
  scale_x_continuous(breaks = seq(1990, 2050, 10)) +
  labs(x = "Year", y = "Age-standardized rate per 100,000",
       colour = NULL,
       title  = "Lip and oral cavity cancer burden in Germany",
       subtitle = "Both sexes, age-standardized rate; observed 1990–2023, projected to 2040") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        plot.title.position = "plot",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# Both-sexes figure retained for reference only; the manuscript's Figure 1 is the
# sex-stratified version produced by 02_results_analysis.R (pooled-sex rates are
# uninformative given the sex divergence).
SUP_DIR <- file.path(FIG_DIR, "_superseded"); dir.create(SUP_DIR, showWarnings = FALSE)
ggsave(file.path(SUP_DIR, "bothsexes_asr_trend_forecast.pdf"),
       p1, width = W_CM, height = H_CM, units = "cm", device = "pdf")
ggsave(file.path(SUP_DIR, "bothsexes_asr_trend_forecast.png"),
       p1, width = W_CM, height = H_CM, units = "cm", dpi = PNG_DPI)
cat("\nBoth-sexes reference figure written to:", SUP_DIR, "\n")

# ---- Summary printout -----------------------------------------------------
cat("\n=== Summary: ASR per 100,000 (Both sexes) ===\n")
key_years <- c(1996, 2021, 2026, 2031, 2036, 2041)
key_lbls  <- c("1994-98","2019-23","2024-28","2029-33","2034-38","2039-43")
key <- projections[sex == "Both" &
                    measure %in% c("Deaths","DALYs (Disability-Adjusted Life Years)","Incidence") &
                    period_mid %in% key_years]
key <- key[, .(measure, period, asr_per_100k = round(asr_per_100k, 1),
               total_count = round(total_count))]
print(key)
cat("\nDone.\n")
