# 03_joinpoint_table2.R
# Joinpoint-style change-point detection for trends 1990–2023 (gold-standard
# cancer epidemiology trend analysis) + Table 2 risk-attributable burden.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(tidyr)
  library(scales); library(segmented)
})

DATA_DIR <- "."  # repository root
CAUSE_CSV <- file.path(DATA_DIR, "data/gbd_cause",
                      "IHME-GBD_2023_DATA-684738bc-1.csv")
RISK_CSV  <- file.path(DATA_DIR, "data/gbd_risk",
                      "IHME-GBD_2023_DATA-6a8e77e9-1.csv")
OUT_DIR <- file.path(DATA_DIR, "output")
FIG_DIR  <- file.path(OUT_DIR, "figures")
TBL_DIR  <- file.path(OUT_DIR, "tables")

W_CM <- 18; PNG_DPI <- 250

d  <- fread(CAUSE_CSV)
dr <- fread(RISK_CSV)

# ============================================================================
# Joinpoint regression: detect change-points in log(ASR) ~ year
# ============================================================================
cat("\n[Joinpoint] Detecting trend change-points...\n")

jp_input <- d[location_name == "Germany" &
              age_name == "Age-standardized" & metric_name == "Rate" &
              measure_name %in% c("Deaths","DALYs (Disability-Adjusted Life Years)",
                                  "Incidence") &
              sex_name %in% c("Both","Male","Female"),
              .(measure = measure_name, sex = sex_name, year, asr = val)]

joinpoint_one <- function(df, max_segments = 2) {
  df <- df[asr > 0]
  if (nrow(df) < 8) return(NULL)
  base <- lm(log(asr) ~ year, data = df)
  out <- list()
  # Try 1 and 2 breakpoints; pick best by BIC
  best <- list(model = base, n_breaks = 0, bic = BIC(base))
  for (k in 1:max_segments) {
    fit_k <- tryCatch(
      segmented(base, seg.Z = ~ year, npsi = k),
      error = function(e) NULL
    )
    if (!is.null(fit_k) && BIC(fit_k) < best$bic - 2) {  # require BIC drop > 2
      best <- list(model = fit_k, n_breaks = k, bic = BIC(fit_k))
    }
  }
  m <- best$model
  if (best$n_breaks == 0) {
    b   <- coef(m)["year"]
    se  <- summary(m)$coefficients["year","Std. Error"]
    segs <- data.table(
      seg_start = min(df$year), seg_end = max(df$year),
      apc = (exp(b)-1)*100,
      apc_lo = (exp(b - 1.96*se)-1)*100,
      apc_hi = (exp(b + 1.96*se)-1)*100
    )
  } else {
    psi <- round(m$psi[, "Est."])
    slopes <- slope(m)$year
    boundaries <- c(min(df$year), psi, max(df$year))
    segs <- data.table(
      seg_start = head(boundaries, -1),
      seg_end   = tail(boundaries, -1),
      apc       = (exp(slopes[,"Est."])-1)*100,
      apc_lo    = (exp(slopes[,"CI(95%).l"])-1)*100,
      apc_hi    = (exp(slopes[,"CI(95%).u"])-1)*100
    )
  }
  segs[, n_breaks := best$n_breaks]
  segs[, period   := sprintf("%d-%d", seg_start, seg_end)]
  segs[, apc_fmt  := sprintf("%+.2f (%+.2f to %+.2f)", apc, apc_lo, apc_hi)]
  segs
}

jp_results <- jp_input[, joinpoint_one(.SD), by = .(measure, sex)]
setorder(jp_results, measure, sex, seg_start)

# AAPC: weighted average of segment APCs
aapc_calc <- function(segs) {
  if (nrow(segs) == 1) return(segs$apc)
  weights <- segs$seg_end - segs$seg_start + 1
  # Geometric average rate via log-space
  log_factors <- log(1 + segs$apc / 100) * weights
  total_w     <- sum(weights)
  (exp(sum(log_factors) / total_w) - 1) * 100
}
aapc <- jp_results[, .(aapc = round(aapc_calc(.SD), 2)),
                   by = .(measure, sex)]

fwrite(jp_results, file.path(TBL_DIR, "joinpoint_segments.csv"))
fwrite(aapc,       file.path(TBL_DIR, "joinpoint_aapc.csv"))

cat("Joinpoint segments table (selected, Both sexes):\n")
print(jp_results[sex == "Both", .(measure, period, n_breaks, apc_fmt)])
cat("\nAAPC (average annual % change, all periods 1990-2023):\n")
print(aapc)

# ============================================================================
# Table 2: Risk-attributable burden detail
# ============================================================================
cat("\n[Table 2] Risk-attributable burden detail...\n")

# Pull risk × measure × year × sex, Number metric
risk_n <- dr[location_name == "Germany" &
             rei_name %in% c("Tobacco","High alcohol use") &
             measure_name %in% c("Deaths","DALYs (Disability-Adjusted Life Years)") &
             metric_name == "Number" &
             age_name == "All ages" &
             sex_name %in% c("Both","Male","Female") &
             year %in% c(1990, 2023),
             .(sex = sex_name, measure = measure_name, year, risk = rei_name,
               attributable = val)]

# Pull totals
total_n <- d[location_name == "Germany" &
             measure_name %in% c("Deaths","DALYs (Disability-Adjusted Life Years)") &
             metric_name == "Number" &
             age_name == "All ages" &
             sex_name %in% c("Both","Male","Female") &
             year %in% c(1990, 2023),
             .(sex = sex_name, measure = measure_name, year, total = val)]

t2 <- merge(risk_n, total_n, by = c("sex","measure","year"))
t2[, paf := round(attributable / total * 100, 1)]
t2[, attributable := round(attributable)]
t2[, total := round(total)]
t2[, measure_lbl := ifelse(measure == "Deaths", "Deaths", "DALYs")]
setcolorder(t2, c("measure_lbl","sex","risk","year","attributable","total","paf"))
t2[, c("measure") := NULL]
setorder(t2, measure_lbl, sex, risk, year)
fwrite(t2, file.path(TBL_DIR, "table2_risk_attributable.csv"))
cat("Table 2 risk-attributable burden:\n")
print(t2)

# Wide format for paper-ready
t2_w <- dcast(t2, measure_lbl + sex + risk ~ year,
              value.var = c("attributable","paf"))
setnames(t2_w,
         old = c("attributable_1990","attributable_2023","paf_1990","paf_2023"),
         new = c("Attributable 1990","Attributable 2023","PAF % 1990","PAF % 2023"))
fwrite(t2_w, file.path(TBL_DIR, "table2_risk_attributable_wide.csv"))
cat("\nTable 2 wide:\n")
print(t2_w)

cat("\nAll analyses complete.\n")
