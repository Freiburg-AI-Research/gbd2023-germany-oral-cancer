# 05_extract_headline_numbers.R
# Pull the exact 2023 (and 1990) Germany numbers with 95% UIs for the
# Results section narrative. Output: a compact reference CSV.

suppressPackageStartupMessages({
  library(data.table); library(dplyr)
})

DATA_DIR <- "."  # repository root
CAUSE_CSV <- file.path(DATA_DIR, "data/gbd_cause",
                      "IHME-GBD_2023_DATA-684738bc-1.csv")
RISK_CSV  <- file.path(DATA_DIR, "data/gbd_risk",
                      "IHME-GBD_2023_DATA-6a8e77e9-1.csv")
TBL_DIR  <- file.path(DATA_DIR, "analysis", "output", "tables")

d  <- fread(CAUSE_CSV)
dr <- fread(RISK_CSV)

fmt <- function(v, lo, hi, dp = 2) {
  sprintf(paste0("%.", dp, "f (%.", dp, "f to %.", dp, "f)"), v, lo, hi)
}
fmt_int <- function(v, lo, hi) {
  sprintf("%s (%s to %s)",
          formatC(round(v),  big.mark = ",", format = "d"),
          formatC(round(lo), big.mark = ",", format = "d"),
          formatC(round(hi), big.mark = ",", format = "d"))
}

# ---- 1. Headline 2023 numbers (Germany, all ages, by sex, all 6 measures) ----
headlines_2023 <- d[location_name == "Germany" & year == 2023 &
                    age_name == "All ages" & metric_name == "Number" &
                    sex_name %in% c("Both","Male","Female")]
headlines_2023 <- headlines_2023[, .(sex = sex_name, measure = measure_name,
                                    val, upper, lower)]
headlines_2023[, fmt := mapply(fmt_int, val, lower, upper)]
fwrite(headlines_2023, file.path(TBL_DIR, "headlines_2023_germany.csv"))
cat("2023 Germany headline counts (95% UI):\n")
print(headlines_2023[, .(sex, measure, fmt)])

# ---- 2. Headline 2023 ASRs (age-standardized) ----
asr_2023 <- d[location_name == "Germany" & year == 2023 &
              age_name == "Age-standardized" & metric_name == "Rate" &
              sex_name %in% c("Both","Male","Female")]
asr_2023 <- asr_2023[, .(sex = sex_name, measure = measure_name,
                        val, upper, lower)]
asr_2023[, fmt := mapply(fmt, val, lower, upper, 2)]
fwrite(asr_2023, file.path(TBL_DIR, "asr_2023_germany.csv"))
cat("\n2023 Germany ASR per 100k (95% UI):\n")
print(asr_2023[, .(sex, measure, fmt)])

# ---- 3. Headline 1990 (for comparison) ----
asr_1990 <- d[location_name == "Germany" & year == 1990 &
              age_name == "Age-standardized" & metric_name == "Rate" &
              sex_name %in% c("Both","Male","Female")]
asr_1990 <- asr_1990[, .(sex = sex_name, measure = measure_name,
                        val, upper, lower)]
asr_1990[, fmt := mapply(fmt, val, lower, upper, 2)]
fwrite(asr_1990, file.path(TBL_DIR, "asr_1990_germany.csv"))
cat("\n1990 Germany ASR per 100k (95% UI):\n")
print(asr_1990[, .(sex, measure, fmt)])

# ---- 4. Western Europe + Global 2023 comparators (Both sexes) ----
comp <- d[year == 2023 & age_name == "Age-standardized" & metric_name == "Rate" &
          sex_name == "Both" &
          location_name %in% c("Germany","Western Europe","Global")]
comp <- comp[, .(location = location_name, measure = measure_name,
                val, upper, lower)]
comp[, fmt := mapply(fmt, val, lower, upper, 2)]
fwrite(comp, file.path(TBL_DIR, "comparators_2023.csv"))
cat("\n2023 comparators (Both sexes, ASR per 100k):\n")
print(comp[, .(location, measure, fmt)])

# ---- 5. PAFs in 2023 with UIs (use Percent metric) ----
paf_2023 <- dr[location_name == "Germany" & year == 2023 &
               age_name == "All ages" & metric_name == "Percent" &
               sex_name %in% c("Both","Male","Female") &
               rei_name %in% c("Tobacco","High alcohol use") &
               measure_name %in% c("Deaths","DALYs (Disability-Adjusted Life Years)")]
paf_2023 <- paf_2023[, .(sex = sex_name, measure = measure_name,
                        risk = rei_name, val = val * 100,
                        lower = lower * 100, upper = upper * 100)]
paf_2023[, fmt := mapply(fmt, val, lower, upper, 1)]
fwrite(paf_2023, file.path(TBL_DIR, "paf_2023_with_ui.csv"))
cat("\n2023 PAFs (% with 95% UI):\n")
print(paf_2023[, .(sex, measure, risk, fmt)])

# ---- 6. Compute M:F ratio at 1990 and 2023 (for ASR convergence) ----
mf_ratio <- function(year_target) {
  x <- d[location_name == "Germany" & year == year_target &
         age_name == "Age-standardized" & metric_name == "Rate" &
         sex_name %in% c("Male","Female")]
  m <- x[sex_name == "Male",   .(measure_name, val)]
  f <- x[sex_name == "Female", .(measure_name, val)]
  setnames(m, "val", "male");   setnames(f, "val", "female")
  rr <- merge(m, f, by = "measure_name")
  rr[, ratio := round(male / female, 2)]
  rr[, year := year_target]
  rr
}
mf <- rbind(mf_ratio(1990), mf_ratio(2023))
fwrite(mf, file.path(TBL_DIR, "mf_ratio_1990_2023.csv"))
cat("\nMale:Female ASR ratio:\n")
print(mf[, .(year, measure_name, male, female, ratio)])

cat("\nAll headline number tables written to:", TBL_DIR, "\n")
