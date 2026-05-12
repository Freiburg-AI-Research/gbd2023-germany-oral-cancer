# 06_zfkd_comparison.R
# Cross-validate GBD 2023 Germany ASR estimates against the German national
# cancer registry (ZfKD / RKI). Produces eFigure 3 and a comparison table.
#
# Input:
#   - data/zfkd/zfkd_incidence_c00_c08.csv
#   - data/zfkd/zfkd_mortality_c00_c08.csv
#   - GBD 2023 Germany ASRs (already in cause CSV)
# Output:
#   - eFigure 3: GBD vs ZfKD ASR trajectories (incidence and mortality, M and F)
#   - eTable: comparison numbers + correlation statistics

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(tidyr)
  library(scales)
})

DATA_DIR <- "."  # repository root
GBD_CAUSE <- file.path(DATA_DIR, "data/gbd_cause",
                      "IHME-GBD_2023_DATA-684738bc-1.csv")
ZFKD_INC <- file.path(DATA_DIR, "data/zfkd", "zfkd_incidence_c00_c08.csv")
ZFKD_MOR <- file.path(DATA_DIR, "data/zfkd", "zfkd_mortality_c00_c08.csv")
OUT_DIR <- file.path(DATA_DIR, "output")
FIG_DIR  <- file.path(OUT_DIR, "figures")
TBL_DIR  <- file.path(OUT_DIR, "tables")

W_CM <- 18; PNG_DPI <- 250

# ---- Read ZfKD CSV (semicolon delimited, German decimal comma) -----------
read_zfkd <- function(path, measure_label) {
  # CSV is semicolon-delimited with German decimal commas. Use fread with explicit
  # delimiter and decimal handling.
  dat <- fread(path, sep = ";", dec = ",", encoding = "UTF-8",
               header = TRUE, fill = TRUE)
  # Header row was "[empty];[empty];[empty];Lippe (C00);...;C08;[empty]"
  # First 3 columns have no name -> fread assigns "V1","V2","V3"
  # Last column is the trailing empty -> fread may add a final empty col
  # Rename systematically:
  diag_cols <- grep("\\(C0[0-8]\\)", names(dat), value = TRUE)
  setnames(dat, diag_cols, paste0("C0", seq_along(diag_cols) - 1))
  setnames(dat, names(dat)[1:3], c("standard","sex","year"))
  # Drop any all-NA trailing column
  drop_cols <- names(dat)[sapply(dat, function(x) all(is.na(x)))]
  if (length(drop_cols)) dat[, (drop_cols) := NULL]
  dat <- dat[!is.na(year)]
  num_cols <- paste0("C0", 0:8)
  dat[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
  dat[, year := as.integer(year)]
  dat[, asr_total := rowSums(.SD, na.rm = TRUE), .SDcols = num_cols]
  dat[, measure := measure_label]
  dat[, sex_en  := ifelse(sex == "weiblich", "Female",
                          ifelse(sex == "männlich", "Male", NA_character_))]
  dat[, standard_short := ifelse(grepl("Welt", standard), "Segi",
                                  ifelse(grepl("Europa", standard), "Europe1976", NA_character_))]
  dat[, .(standard_short, sex = sex_en, year, asr_total, measure)]
}

zinc <- read_zfkd(ZFKD_INC, "Incidence")
zmor <- read_zfkd(ZFKD_MOR, "Deaths")
zfkd <- rbind(zinc, zmor)
fwrite(zfkd, file.path(TBL_DIR, "zfkd_long.csv"))
cat("ZfKD long-format data sample:\n")
print(head(zfkd, 6))
cat(sprintf("\nZfKD rows: %d (expected ~ 25 yr * 2 sex * 2 std * 2 measure = 200)\n",
            nrow(zfkd)))

# ---- GBD 2023 Germany ASR (use Welt/Segi-equivalent: GBD World Standard) ----
g <- fread(GBD_CAUSE)
gbd <- g[location_name == "Germany" & age_name == "Age-standardized" &
         metric_name == "Rate" & sex_name %in% c("Male","Female") &
         measure_name %in% c("Deaths","Incidence"),
         .(measure = measure_name, sex = sex_name, year, asr = val,
           asr_lo = lower, asr_hi = upper)]
gbd[, measure := ifelse(measure == "Deaths", "Deaths", "Incidence")]

# ---- Combine for plotting (use Segi for GBD comparability) ----
zfkd_segi <- zfkd[standard_short == "Segi"]
combined <- rbind(
  zfkd_segi[, .(source = "ZfKD (Segi)", measure, sex, year, asr = asr_total,
                asr_lo = NA_real_, asr_hi = NA_real_)],
  gbd[, .(source = "GBD 2023 (World Standard)", measure, sex, year, asr,
          asr_lo, asr_hi)]
)
combined <- combined[year >= 1999]

# ---- Compute comparison statistics ----
wide <- dcast(combined[year >= 1999 & year <= 2022],
              measure + sex + year ~ source, value.var = "asr")
setnames(wide, "GBD 2023 (World Standard)", "GBD")
setnames(wide, "ZfKD (Segi)", "ZfKD")

stats_by_strata <- wide[!is.na(GBD) & !is.na(ZfKD), .(
  n          = .N,
  mean_GBD   = mean(GBD),
  mean_ZfKD  = mean(ZfKD),
  diff_pct   = mean((GBD - ZfKD) / ZfKD * 100),
  abs_diff_pct = mean(abs(GBD - ZfKD) / ZfKD * 100),
  pearson_r  = cor(GBD, ZfKD)
), by = .(measure, sex)]
fwrite(stats_by_strata, file.path(TBL_DIR, "gbd_vs_zfkd_stats.csv"))
cat("\nGBD vs ZfKD comparison statistics (Segi standard, 1999-2022):\n")
print(stats_by_strata)

# ---- eFigure 3: side-by-side trajectories ----
combined[, measure := factor(measure, levels = c("Incidence","Deaths"))]
combined[, sex := factor(sex, levels = c("Male","Female"))]

p3 <- ggplot(combined,
            aes(x = year, y = asr, colour = source, linetype = source)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2, alpha = 0.7) +
  geom_ribbon(data = combined[source == "GBD 2023 (World Standard)"],
              aes(ymin = asr_lo, ymax = asr_hi, fill = source),
              alpha = 0.15, colour = NA) +
  facet_grid(measure ~ sex, scales = "free_y") +
  scale_colour_manual(values = c("GBD 2023 (World Standard)" = "#1f77b4",
                                 "ZfKD (Segi)" = "#d62728")) +
  scale_fill_manual(values  = c("GBD 2023 (World Standard)" = "#1f77b4",
                                 "ZfKD (Segi)" = "#d62728")) +
  scale_linetype_manual(values = c("GBD 2023 (World Standard)" = "solid",
                                   "ZfKD (Segi)" = "dashed")) +
  scale_x_continuous(breaks = seq(2000, 2020, 10)) +
  labs(x = "Year", y = "Age-standardised rate per 100,000",
       colour = NULL, fill = NULL, linetype = NULL,
       title = "Cross-validation: GBD 2023 vs ZfKD (RKI) estimates for Germany",
       subtitle = "Lip and oral cavity cancer (ICD-10 C00–C08), 1999–2023, by sex and measure") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom",
        plot.title.position = "plot",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "figureS3_gbd_vs_zfkd.pdf"),
       p3, width = W_CM, height = 12, units = "cm", device = "pdf")
ggsave(file.path(FIG_DIR, "figureS3_gbd_vs_zfkd.png"),
       p3, width = W_CM, height = 12, units = "cm", dpi = PNG_DPI)
cat("\neFigure 3 written.\n")

# ---- Print summary for the manuscript ----
cat("\n=== Headline comparison numbers (for Methods/Discussion) ===\n\n")
overall <- stats_by_strata[, .(
  mean_pearson_r = mean(pearson_r),
  mean_abs_pct_diff = mean(abs_diff_pct),
  range_pct_diff = sprintf("%.1f to %.1f", min(diff_pct), max(diff_pct))
)]
cat(sprintf("Average Pearson correlation across the four strata: r = %.3f\n",
            overall$mean_pearson_r))
cat(sprintf("Average absolute percentage difference: %.1f%%\n",
            overall$mean_abs_pct_diff))
cat(sprintf("Range of mean % difference (GBD - ZfKD)/ZfKD: %s%%\n",
            overall$range_pct_diff))
