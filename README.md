# Lip and oral cavity cancer in Germany — GBD 2023 national analysis

Analysis code, derived tables, and figures for a national analysis of the burden,
attributable risk factors, and projections to 2040 of lip and oral cavity cancer
in Germany, based on the Global Burden of Disease Study 2023.

## Summary of the analysis

- **Population** — Germany, with Western Europe and the global aggregate as comparators.
- **Cause** — Lip and oral cavity cancer (GBD 2023 cause; ICD-10 C00–C08).
- **Time horizon** — 1990–2023 observed; 2024–2043 projected.
- **Metrics** — Age-standardised rates (incidence, prevalence, deaths, YLLs, YLDs,
  DALYs) and absolute counts, by sex and 5-year age group.
- **Statistical methods** — log-linear regression for the estimated annual percent
  change (EAPC); segmented regression for joinpoint detection of trend
  change-points; the Nordpred age-period-cohort model with the power-5 link
  function and the standard trend-attenuation sequence for the projection
  (Møller et al., *Stat Med* 2003).
- **Cross-validation** — comparison against the German Centre for Cancer Registry
  Data (ZfKD, Robert Koch-Institute) for matching ICD-10 scope.
- **Sensitivity** — population-projection sensitivity using the UN World
  Population Prospects 2024 medium-variant projection for Germany.

## Repository layout

```
.
├── README.md
├── LICENSE
├── analysis/                 # R scripts (numbered, run in order)
│   ├── 00_setup.R            # one-time package install (nordpred, segmented, wpp2024, ...)
│   ├── 01_nordpred_germany.R # main projection: Nordpred for deaths, DALYs, incidence
│   ├── 02_results_analysis.R # EAPC, Figure 2 (age pattern), Figure 3 (risk decomposition), Table 1, supplements
│   ├── 03_joinpoint_table2.R # Joinpoint regression and Table 2 (risk-attributable burden)
│   ├── 04_joinpoint_figure.R # eFigure 2 (joinpoint visualisation)
│   ├── 05_extract_headline_numbers.R  # headline 2023 numbers with 95% UIs for the manuscript
│   ├── 06_zfkd_comparison.R  # ZfKD vs GBD 2023 cross-validation (eFigure 3)
│   └── 07_population_sensitivity.R    # UN WPP 2024 sensitivity (eFigure 4)
├── data/
│   ├── README.md             # instructions for downloading GBD, ZfKD, and WPP inputs
│   ├── gbd_cause/            # place GBD cause CSV here (see data/README.md)
│   ├── gbd_risk/             # place GBD risk CSV here
│   └── zfkd/                 # place ZfKD CSVs here (optional, for cross-validation)
└── output/
    ├── figures/              # main and supplementary figures (PDF and PNG)
    ├── tables/               # tables (CSV)
    └── data/                 # derived datasets (CSV)
```

## How to reproduce

1. **Install R** (version 4.5.1 or newer was used here).
2. **Clone this repository** and open R at the repository root:

   ```bash
   git clone https://github.com/Freiburg-AI-Research/gbd2023-germany-oral-cancer.git
   cd gbd2023-germany-oral-cancer
   ```

3. **Download the input data** following `data/README.md`. The two GBD CSV
   files go to `data/gbd_cause/` and `data/gbd_risk/`; the two optional ZfKD
   CSV files go to `data/zfkd/`. The WPP 2024 input is installed automatically
   as an R package.

4. **Run the analysis scripts in order**:

   ```bash
   Rscript analysis/00_setup.R
   Rscript analysis/01_nordpred_germany.R
   Rscript analysis/02_results_analysis.R
   Rscript analysis/03_joinpoint_table2.R
   Rscript analysis/04_joinpoint_figure.R
   Rscript analysis/05_extract_headline_numbers.R
   Rscript analysis/06_zfkd_comparison.R   # requires the optional ZfKD CSVs
   Rscript analysis/07_population_sensitivity.R
   ```

   All scripts are designed to be run from the repository root. Outputs are
   written to `output/`.

## R packages

The setup script installs the following CRAN packages: `data.table`, `ggplot2`,
`dplyr`, `tidyr`, `readr`, `scales`, `patchwork`, `remotes`, `segmented`,
`wpp2024`. The `nordpred` package is installed from a GitHub fork because the
original CRAN release is no longer maintained for current R versions.

## Citation

If you use this code or the derived outputs, please cite the corresponding
manuscript (citation details to be added on publication).

Underlying data should be cited as documented in `data/README.md`:

- GBD 2023: Global Burden of Disease Collaborative Network. Global Burden of
  Disease Study 2023 (GBD 2023) Results. Seattle: IHME, 2024.
- ZfKD: Zentrum für Krebsregisterdaten im Robert Koch-Institut.
- WPP 2024: UN Population Division. World Population Prospects 2024.

## License

The code in this repository is released under the MIT License (see `LICENSE`).
The underlying GBD, ZfKD, and WPP data remain subject to their respective
terms of use.

## Contact

Babak Saravi · Department of Oral, Maxillofacial and Facial Plastic Surgery,
University Medical Center Düsseldorf, Heinrich Heine University Düsseldorf,
Germany.
