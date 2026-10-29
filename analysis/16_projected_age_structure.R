# 16_projected_age_structure.R
# Age distribution of incident cases and deaths, 2019-2023 (observed period)
# vs 2039-2043 (Nordpred projection), under (a) static population and
# (b) UN WPP 2024 medium-variant population. Reports share aged >=70/75/80 and
# mean age, by sex. Service-planning implication: the future patient is older.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(nordpred); library(wpp2024) })
DATA_DIR <- "."  # repository root
CAUSE_CSV <- file.path(DATA_DIR, "data/gbd_cause", "IHME-GBD_2023_DATA-684738bc-1.csv")
OUT_DIR   <- file.path(DATA_DIR, "output"); FIG_DIR <- file.path(OUT_DIR, "figures"); TBL_DIR <- file.path(OUT_DIR, "tables")
W_CM <- 18; PNG_DPI <- 250

AGE_RAW <- c("<5 years","5-9 years","10-14 years","15-19 years","20-24 years","25-29 years","30-34 years","35-39 years",
             "40-44 years","45-49 years","50-54 years","55-59 years","60-64 years","65-69 years","70-74 years",
             "75-79 years","80-84 years","85-89 years","90-94 years","95+ years")
AGE18 <- c(AGE_RAW[1:17], "85+ years"); AGE_MID18 <- c(seq(2.5, 82.5, by = 5), 90)
HIST <- list(c(1994,1998),c(1999,2003),c(2004,2008),c(2009,2013),c(2014,2018),c(2019,2023))
L <- c("1994-1998","1999-2003","2004-2008","2009-2013","2014-2018","2019-2023")
FL <- c("2024-2028","2029-2033","2034-2038","2039-2043"); FMID <- c(2026, 2031, 2036, 2041)

# ---- WPP 2024 person-years per future period (as in 07_population_sensitivity.R)
data(popMprojMed); data(popFprojMed)
wmap <- c("<5 years"="0-4","5-9 years"="5-9","10-14 years"="10-14","15-19 years"="15-19","20-24 years"="20-24","25-29 years"="25-29",
          "30-34 years"="30-34","35-39 years"="35-39","40-44 years"="40-44","45-49 years"="45-49","50-54 years"="50-54","55-59 years"="55-59",
          "60-64 years"="60-64","65-69 years"="65-69","70-74 years"="70-74","75-79 years"="75-79","80-84 years"="80-84")
interp <- function(dt, agew, yr) { row <- dt[age == agew]; yrs <- as.integer(grep("^[0-9]+$", names(dt), value = TRUE))
  lo <- max(yrs[yrs <= yr]); hi <- min(yrs[yrs >= yr]); v1 <- as.numeric(row[[as.character(lo)]]); v2 <- as.numeric(row[[as.character(hi)]])
  if (lo == hi) v1 else v1 + (v2 - v1) * (yr - lo) / (hi - lo) }
wpp_pyr <- function(sex) {
  dt <- as.data.table(if (sex == "Male") popMprojMed else popFprojMed)[name == "Germany"]
  m <- matrix(NA_real_, 18, 4, dimnames = list(AGE18, FL))
  for (a in AGE18) for (i in 1:4) m[a, i] <- if (a == "85+ years") sum(sapply(c("85-89","90-94","95-99","100+"), function(w) interp(dt, w, FMID[i]))) else interp(dt, wmap[[a]], FMID[i])
  m * 1000 * 5 }
WPP <- list(Male = wpp_pyr("Male"), Female = wpp_pyr("Female")); WPP$Both <- WPP$Male + WPP$Female

# ---- GBD cases/pyr matrices
d <- fread(CAUSE_CSV)
mk <- function(sex_i, meas_i) {
  g <- d[location_name == "Germany" & age_name %in% AGE_RAW & sex_name == sex_i & measure_name == meas_i]
  n <- g[metric_name == "Number", .(year, age_name, n = val)]; r <- g[metric_name == "Rate", .(year, age_name, r = val)]
  p <- merge(n, r, by = c("year","age_name"))[, pop := ifelse(r > 0, n / r * 1e5, NA_real_)]
  p[, age := fifelse(age_name %in% c("85-89 years","90-94 years","95+ years"), "85+ years", age_name)]
  a <- p[, .(n = sum(n), pop = sum(pop, na.rm = TRUE)), by = .(year, age)]; a[, period := NA_character_]
  for (i in seq_along(HIST)) a[year >= HIST[[i]][1] & year <= HIST[[i]][2], period := L[i]]
  agg <- a[!is.na(period), .(cases = sum(n), pyr = sum(pop)), by = .(age, period)]
  cm <- dcast(agg, age ~ period, value.var = "cases"); cm <- as.data.frame(cm[match(AGE18, age), L, with = FALSE]); rownames(cm) <- AGE18; cm[is.na(cm)] <- 0; cm <- round(cm)
  pm <- dcast(agg, age ~ period, value.var = "pyr");   pm <- as.data.frame(pm[match(AGE18, age), L, with = FALSE]); rownames(pm) <- AGE18; pm[is.na(pm) | pm <= 0] <- 1
  list(cm = cm, pm = pm) }
run <- function(sex_i, meas_i, future_pyr) {
  m <- mk(sex_i, meas_i); pm <- cbind(m$pm, as.data.frame(future_pyr))
  fit <- nordpred(cases = m$cm, pyr = pm, startestage = 6, startuseage = 6, noperiods = 4, cuttrend = c(0, .25, .5, .75), linkfunc = "power5")
  cnt <- as.matrix(nordpred.getpred(fit, incidence = FALSE, standpop = NULL)); rownames(cnt) <- AGE18   # rows are in age order
  cnt }

rows <- list(); dist_rows <- list()
for (sx in c("Both","Male","Female")) for (ms in c("Incidence","Deaths")) {
  static <- { m <- mk(sx, ms); fut <- matrix(rep(m$pm[, 6], 4), 18, dimnames = list(AGE18, FL)); run(sx, ms, fut) }
  dyn <- run(sx, ms, WPP[[sx]])
  for (scen in c("static","WPP 2024")) {
    cnt <- if (scen == "static") static else dyn
    for (per in c("2019-2023","2039-2043")) {
      v <- cnt[, per]; tot <- sum(v)
      rows[[length(rows) + 1]] <- data.table(sex = sx, measure = ms, scenario = scen, period = per, annual_count = tot / 5,
        share_70plus = sum(v[15:18]) / tot * 100, share_75plus = sum(v[16:18]) / tot * 100, share_80plus = sum(v[17:18]) / tot * 100,
        mean_age = sum(v * AGE_MID18) / tot)
      dist_rows[[length(dist_rows) + 1]] <- data.table(sex = sx, measure = ms, scenario = scen, period = per, age = AGE18, count = v / 5) } } }
tab <- rbindlist(rows); dist <- rbindlist(dist_rows)
fwrite(tab, file.path(TBL_DIR, "projected_age_structure.csv")); fwrite(dist, file.path(TBL_DIR, "projected_age_distribution.csv"))
cat("=== Age structure of annual cases/deaths: 2019-2023 vs 2039-2043 ===\n")
print(tab[order(measure, sex, scenario, period), .(measure, sex, scenario, period, annual = round(annual_count), `>=70` = sprintf("%.0f%%", share_70plus),
      `>=75` = sprintf("%.0f%%", share_75plus), `>=80` = sprintf("%.0f%%", share_80plus), mean_age = round(mean_age, 1))])

# Figure: age distribution of annual incident cases and deaths, 2019-23 vs 2039-43 (WPP), by sex
dd <- dist[scenario == "WPP 2024" & sex != "Both" & age %in% AGE18[6:18]]
dd[, age := factor(age, levels = AGE18)]; dd[, sex := factor(sex, levels = c("Male","Female"))]
dd[, measure := factor(measure, levels = c("Incidence","Deaths"))]
p <- ggplot(dd, aes(age, count, fill = period)) + geom_col(position = position_dodge(width = .8), width = .75) +
  facet_grid(measure ~ sex, scales = "free_y") +
  scale_fill_manual(values = c(`2019-2023` = "grey60", `2039-2043` = "#1f77b4"), name = NULL) +
  scale_x_discrete(labels = function(l) sub(" years", "", l)) +
  scale_y_continuous(expand = expansion(mult = c(0, .06))) +
  labs(x = "Age group (years)", y = "Annual number", title = "Age distribution of lip and oral cavity cancer cases and deaths in Germany",
       subtitle = "Observed period 2019–2023 vs Nordpred projection 2039–2043 with UN WPP 2024 medium-variant population") +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom", strip.text = element_text(face = "bold"), panel.grid.minor = element_blank(),
        plot.title.position = "plot", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIG_DIR, "figure_projected_age_structure.pdf"), p, width = W_CM, height = 11, units = "cm", device = "pdf")
ggsave(file.path(FIG_DIR, "figure_projected_age_structure.png"), p, width = W_CM, height = 11, units = "cm", dpi = PNG_DPI)
ggsave(file.path(FIG_DIR, "figure_projected_age_structure.tif"), p, width = W_CM, height = 11, units = "cm", dpi = 600, device = "tiff", compression = "lzw", type = "cairo")
cat("\nfigure_projected_age_structure written.\n")
