# 12_subsite_analysis.R
# Sub-site decomposition of the GBD "lip and oral cavity cancer" cause (ICD-10
# C00-C08) using ZfKD registry data, which report each 3-digit site separately.
# Addresses the reviewer/co-author critique that C00-C08 pools etiologically
# distinct entities:
#   C00      lip (UV/age-driven, cSCC-like)
#   C01      base of tongue (conventionally oropharyngeal; HPV-related)
#   C02-C06  oral cavity proper (tobacco/alcohol-driven)
#   C07-C08  major salivary glands (distinct etiology)
# Outputs: composition table (5-yr average), sub-site EAPC table, eFigure.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales) })

DATA_DIR <- "."  # repository root
ZFKD_INC <- file.path(DATA_DIR, "data/zfkd", "zfkd_incidence_c00_c08.csv")
ZFKD_MOR <- file.path(DATA_DIR, "data/zfkd", "zfkd_mortality_c00_c08.csv")
OUT_DIR  <- file.path(DATA_DIR, "output")
FIG_DIR  <- file.path(OUT_DIR, "figures"); TBL_DIR <- file.path(OUT_DIR, "tables")
W_CM <- 18; PNG_DPI <- 250

GROUPS <- list(`Lip (C00)` = "C00", `Base of tongue (C01)` = "C01",
               `Oral cavity (C02-C06)` = c("C02","C03","C04","C05","C06"),
               `Major salivary glands (C07-C08)` = c("C07","C08"))
GROUP_LEVELS <- names(GROUPS)

read_zfkd_sites <- function(path, measure_label) {
  dat <- fread(path, sep = ";", dec = ",", encoding = "UTF-8", header = TRUE, fill = TRUE)
  setnames(dat, names(dat)[1:3], c("standard","sex","year"))
  dat <- dat[!is.na(year)]
  site_cols <- grep("C0[0-8]", names(dat), value = TRUE)
  codes <- regmatches(site_cols, regexpr("C0[0-8]", site_cols))
  setnames(dat, site_cols, codes)
  dat[, (codes) := lapply(.SD, as.numeric), .SDcols = codes]
  long <- melt(dat[, c("standard","sex","year", codes), with = FALSE],
               id.vars = c("standard","sex","year"), variable.name = "code", value.name = "asr")
  long[, code := as.character(code)]
  long[, sex := fifelse(sex == "männlich", "Male", "Female")]
  long[, standard_short := fifelse(grepl("Welt", standard), "Segi", "Europe1976")]
  long[, group := NA_character_]
  for (g in names(GROUPS)) long[code %in% GROUPS[[g]], group := g]
  long[, measure := measure_label]
  long[, .(standard_short, sex, year = as.integer(year), code, group, asr, measure)]
}
sites <- rbind(read_zfkd_sites(ZFKD_INC, "Incidence"), read_zfkd_sites(ZFKD_MOR, "Deaths"))
grp <- sites[, .(asr = sum(asr)), by = .(standard_short, sex, year, group, measure)]
grp[, group := factor(group, levels = GROUP_LEVELS)]
grp[, measure := factor(measure, levels = c("Incidence","Deaths"))]
grp[, sex := factor(sex, levels = c("Male","Female"))]
fwrite(grp, file.path(TBL_DIR, "zfkd_subsite_group_asr.csv"))

# ---- Composition: 5-year average 2018-2022 (smooths one-decimal rounding) ---
comp <- grp[standard_short == "Segi" & year %in% 2018:2022,
            .(asr_mean = mean(asr)), by = .(measure, sex, group)]
comp[, total := sum(asr_mean), by = .(measure, sex)]
comp[, share_pct := 100 * asr_mean / total]
setorder(comp, measure, sex, group)
fwrite(comp, file.path(TBL_DIR, "zfkd_subsite_composition_2018_2022.csv"))
cat("=== Sub-site composition of C00-C08, ZfKD, Segi ASR, 2018-2022 mean ===\n")
print(comp[, .(measure, sex, group, asr = round(asr_mean, 2), share = sprintf("%.0f%%", share_pct))])

# ---- Sub-site trends: EAPC 1999-2022 and 2009-2022 ------------------------
eapc <- function(y, x) {
  ok <- y > 0
  if (sum(ok) < 5 || sd(log(y[ok])) == 0) return(c(eapc = NA_real_, lo = NA_real_, hi = NA_real_))
  f <- lm(log(y[ok]) ~ x[ok]); b <- unname(coef(f)[2]); se <- unname(sqrt(vcov(f)[2, 2]))
  c(eapc = (exp(b) - 1) * 100, lo = (exp(b - 1.96 * se) - 1) * 100, hi = (exp(b + 1.96 * se) - 1) * 100)
}
trend <- rbind(
  grp[standard_short == "Segi" & year %in% 1999:2022, c(as.list(eapc(asr, year)), window = "1999-2022", mean_asr = mean(asr)), by = .(measure, sex, group)],
  grp[standard_short == "Segi" & year %in% 2009:2022, c(as.list(eapc(asr, year)), window = "2009-2022", mean_asr = mean(asr)), by = .(measure, sex, group)])
# Rates reported to one decimal: flag strata where mean ASR < 0.3 (rounding-dominated)
trend[, reliable := mean_asr >= 0.3]
setorder(trend, window, measure, sex, group)
fwrite(trend, file.path(TBL_DIR, "zfkd_subsite_eapc.csv"))
cat("\n=== Sub-site EAPC (ZfKD, Segi); strata with mean ASR < 0.3 flagged unreliable ===\n")
print(trend[, .(window, measure, sex, group, mean_asr = round(mean_asr, 2),
                EAPC = ifelse(is.na(eapc), "n/a", sprintf("%+.2f (%+.2f to %+.2f)", eapc, lo, hi)),
                reliable)])

# ---- eFigure: sub-site group trends, 1999-2023, Segi ------------------------
pal <- c(`Lip (C00)` = "#e6ab02", `Base of tongue (C01)` = "#7570b3",
         `Oral cavity (C02-C06)` = "#d95f02", `Major salivary glands (C07-C08)` = "#1b9e77")
p <- ggplot(grp[standard_short == "Segi"], aes(year, asr, colour = group)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1) +
  facet_grid(measure ~ sex, scales = "free_y") +
  scale_colour_manual(values = pal) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .06))) +
  scale_x_continuous(breaks = seq(2000, 2020, 10)) +
  labs(x = "Year", y = "ASR per 100,000 (Segi world standard)", colour = NULL,
       title = "Sub-site composition of ICD-10 C00-C08 in Germany (ZfKD registry data), 1999-2023",
       subtitle = "2023 provisional. Rates are published to one decimal; strata below ~0.3 per 100,000 are rounding-dominated.") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(), plot.title.position = "plot") +
  guides(colour = guide_legend(nrow = 2))
ggsave(file.path(FIG_DIR, "figureS5_zfkd_subsites.pdf"), p, width = W_CM, height = 12, units = "cm", device = "pdf")
ggsave(file.path(FIG_DIR, "figureS5_zfkd_subsites.png"), p, width = W_CM, height = 12, units = "cm", dpi = PNG_DPI)
cat("\neFigure 5 (sub-sites) written.\n")
