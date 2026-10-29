# 14_sex_absolute_relative.R
# Absolute and relative change 1990 -> 2023 by sex, to address the floor-effect
# critique: women started from a much lower baseline, so "men benefited and
# women did not" must be read against both absolute (per 100,000) and relative
# (%) change, and the male-female gap.

suppressPackageStartupMessages(library(data.table))
DATA_DIR <- "."  # repository root
CAUSE_CSV <- file.path(DATA_DIR, "data/gbd_cause", "IHME-GBD_2023_DATA-684738bc-1.csv")
TBL_DIR   <- file.path(DATA_DIR, "output", "tables")

d <- fread(CAUSE_CSV)
asr <- d[location_name == "Germany" & age_name == "Age-standardized" & metric_name == "Rate" &
         sex_name %in% c("Male","Female") & year %in% c(1990, 2023) &
         measure_name %in% c("Incidence","Deaths","DALYs (Disability-Adjusted Life Years)"),
         .(measure = measure_name, sex = sex_name, year, asr = val, lo = lower, hi = upper)]
asr[, measure := fifelse(measure == "DALYs (Disability-Adjusted Life Years)", "DALYs", measure)]
w <- dcast(asr, measure + sex ~ year, value.var = c("asr","lo","hi"))
w[, abs_change := asr_2023 - asr_1990]
w[, rel_change_pct := (asr_2023 / asr_1990 - 1) * 100]
setorder(w, measure, sex)

gap <- dcast(asr[, .(measure, sex, year, asr)], measure + year ~ sex, value.var = "asr")
gap[, `:=`(abs_gap = Male - Female, ratio = Male / Female)]
gapw <- dcast(gap, measure ~ year, value.var = c("abs_gap","ratio"))
gapw[, gap_closed_pct := (1 - abs_gap_2023 / abs_gap_1990) * 100]

# Share of the total (both-sex) absolute ASR decline attributable to men
both <- d[location_name == "Germany" & age_name == "Age-standardized" & metric_name == "Rate" &
          sex_name == "Both" & year %in% c(1990, 2023) &
          measure_name %in% c("Incidence","Deaths","DALYs (Disability-Adjusted Life Years)"),
          .(measure = fifelse(measure_name == "DALYs (Disability-Adjusted Life Years)", "DALYs", measure_name), year, asr = val)]
bothw <- dcast(both, measure ~ year, value.var = "asr")[, both_abs_change := `2023` - `1990`]

out <- merge(w[, .(measure, sex, asr_1990, asr_2023, abs_change, rel_change_pct)],
             gapw[, .(measure, abs_gap_1990, abs_gap_2023, ratio_1990, ratio_2023, gap_closed_pct)], by = "measure")
out <- merge(out, bothw[, .(measure, both_abs_change)], by = "measure")
fwrite(out, file.path(TBL_DIR, "sex_absolute_relative_change.csv"))

cat("=== ASR change 1990 -> 2023 by sex (GBD world standard) ===\n")
print(out[, .(measure, sex, asr_1990 = round(asr_1990, 2), asr_2023 = round(asr_2023, 2),
              abs_change = round(abs_change, 2), rel_change_pct = round(rel_change_pct, 1))])
cat("\n=== Male-female gap ===\n")
print(unique(out[, .(measure, abs_gap_1990 = round(abs_gap_1990, 2), abs_gap_2023 = round(abs_gap_2023, 2),
                     ratio_1990 = round(ratio_1990, 2), ratio_2023 = round(ratio_2023, 2),
                     gap_closed_pct = round(gap_closed_pct, 0), both_sex_abs_change = round(both_abs_change, 2))]))
