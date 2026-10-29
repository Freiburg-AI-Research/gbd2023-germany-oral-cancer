# Lip and oral cavity cancer in Germany — GBD 2023 national analysis

Analysis code, derived tables, and figures for a national analysis of lip and
oral cavity cancer in Germany, 1990–2023, with projections to 2040, based on
the Global Burden of Disease Study 2023 integrated with national cancer-registry
data (German Centre for Cancer Registry Data, ZfKD).

## Summary of the analysis

- **Population** — Germany, with Western Europe and the global aggregate as comparators.
- **Cause** — Lip and oral cavity cancer (GBD 2023 cause; ICD-10 C00–C08), with
  sub-site decomposition (C00; C01; C02–C06; C07–C08) from registry data.
- **Time horizon** — 1990–2023 observed; 2024–2043 projected.
- **Metrics** — Age-standardised rates (incidence, prevalence, deaths, YLLs, YLDs,
  DALYs), age-specific rates, absolute counts, mean age at diagnosis and death,
  and tobacco-/alcohol-attributable burden, by sex and 5-year age group.
- **Trend analysis** — log-linear regression (EAPC) overall and by age band;
  segmented regression (Joinpoint) for change-points.
- **Age-period-cohort modelling and projection** — Nordpred (power-5 link,
  drift attenuation) with extraction of birth-cohort effects, drift and
  goodness-of-fit diagnostics; projection to 2039–2043 under a static population
  and under the UN World Population Prospects 2024 medium variant, including the
  projected age distribution of cases and deaths.
- **Decomposition** — change in annual cases and deaths 1990→2023 split into
  population growth, population ageing and change in age-specific rates
  (Das Gupta-type, averaged over factor orderings).
- **Agreement with registry data** — GBD age-specific rates re-standardised to
  the Segi and old European standards; ratio, Bland–Altman limits and slope
  comparison against ZfKD, 1999–2022 and 2009–2022.

## Repository layout

```
.
├── README.md
├── LICENSE
├── analysis/                 # R scripts (numbered, run in order)
│   ├── 00_setup.R                     # one-time package install (nordpred, segmented, wpp2024, ...)
│   ├── 01_nordpred_germany.R          # Nordpred projections for deaths, DALYs, incidence
│   ├── 02_results_analysis.R          # EAPC, Figure 1 (sex-stratified trends), age pattern, risk decomposition, Table 1
│   ├── 03_joinpoint_table2.R          # Joinpoint regression and risk-attributable table
│   ├── 04_joinpoint_figure.R          # Joinpoint figure
│   ├── 05_extract_headline_numbers.R  # headline 2023 numbers with 95% UIs
│   ├── 06_zfkd_comparison.R           # ZfKD data import (long format)
│   ├── 07_population_sensitivity.R    # UN WPP 2024 population sensitivity
│   ├── 11_crossvalidation_rebuild.R   # GBD vs ZfKD agreement on common standard populations
│   ├── 12_subsite_analysis.R          # sub-site composition and trends (ZfKD)
│   ├── 13_nordpred_apc_effects.R      # birth-cohort effects, drift, goodness of fit (Figure 3)
│   ├── 14_sex_absolute_relative.R     # absolute and relative change by sex, 1990–2023
│   ├── 15_age_specific_trends.R       # age-band EAPCs, mean age at diagnosis/death (Figure 2)
│   ├── 16_projected_age_structure.R   # age distribution of projected cases and deaths (Figure 4)
│   └── 17_decomposition.R             # growth / ageing / rate decomposition; attributable change
├── data/
│   ├── README.md             # instructions for downloading GBD, ZfKD, and WPP inputs
│   ├── gbd_cause/            # place GBD cause CSV here (see data/README.md)
│   ├── gbd_risk/             # place GBD risk CSV here
│   └── zfkd/                 # place ZfKD CSVs here (sub-site analysis and agreement)
└── output/
    ├── figures/              # main and supplementary figures (PDF and PNG)
    ├── tables/               # tables (CSV)
    └── data/                 # derived datasets (CSV)
```

Main-text figures of the manuscript correspond to `figure1_germany_sex_stratified`,
`figure_age_specific_trends`, `figureS6_nordpred_cohort_effects` and
`figure_projected_age_structure`; the remaining figures are supplementary.

## How to reproduce

1. **Install R** (version 4.5.1 or newer was used here).
2. **Clone this repository** and open R at the repository root:

   ```bash
   git clone https://github.com/Freiburg-AI-Research/gbd2023-germany-oral-cancer.git
   cd gbd2023-germany-oral-cancer
   ```

3. **Download the input data** following `data/README.md`. The two GBD CSV
   files go to `data/gbd_cause/` and `data/gbd_risk/`; the two ZfKD CSV files
   go to `data/zfkd/`. The WPP 2024 input is installed automatically as an R
   package from GitHub (PPgp/wpp2024) by `00_setup.R`.

4. **Run the analysis scripts in order**:

   ```bash
   Rscript analysis/00_setup.R
   Rscript analysis/01_nordpred_germany.R
   Rscript analysis/02_results_analysis.R
   Rscript analysis/03_joinpoint_table2.R
   Rscript analysis/04_joinpoint_figure.R
   Rscript analysis/05_extract_headline_numbers.R
   Rscript analysis/06_zfkd_comparison.R
   Rscript analysis/07_population_sensitivity.R
   Rscript analysis/11_crossvalidation_rebuild.R
   Rscript analysis/12_subsite_analysis.R
   Rscript analysis/13_nordpred_apc_effects.R
   Rscript analysis/14_sex_absolute_relative.R
   Rscript analysis/15_age_specific_trends.R
   Rscript analysis/16_projected_age_structure.R
   Rscript analysis/17_decomposition.R
   ```

   All scripts are designed to be run from the repository root. Outputs are
   written to `output/`.

## R packages

The setup script installs the following CRAN packages: `data.table`, `ggplot2`,
`dplyr`, `tidyr`, `readr`, `scales`, `patchwork`, `remotes`, `segmented`. The
`nordpred` package is installed from a GitHub fork because the original CRAN
release is no longer maintained for current R versions. The `wpp2024` package
is installed from [PPgp/wpp2024](https://github.com/PPgp/wpp2024) on GitHub;
it is not on CRAN because of its size.

## Citation

If you use this code or the derived outputs, please cite the corresponding
manuscript (citation details to be added on publication).

Underlying data should be cited as documented in `data/README.md`:

- GBD 2023: Global Burden of Disease Collaborative Network. Global Burden of
  Disease Study 2023 (GBD 2023) Results. Seattle: IHME, 2024.
- ZfKD: Zentrum für Krebsregisterdaten im Robert Koch-Institut.
- WPP 2024: United Nations, Department of Economic and Social Affairs,
  Population Division. World Population Prospects 2024 (released 11 July 2024).

## License

The code in this repository is released under the MIT License (see `LICENSE`).
The underlying GBD, ZfKD, and WPP data remain subject to their respective
terms of use.

## Contact

Babak Saravi · Department of Oral, Maxillofacial and Facial Plastic Surgery,
University Medical Center Düsseldorf, Heinrich Heine University Düsseldorf,
Germany.
