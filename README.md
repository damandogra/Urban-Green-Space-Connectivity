# Spatial Justice through Urban Green Space Connectivity

Project report of Group E — Applied Spatial Analytics 2026, TU Delft.

**Live report:** https://applied-spatial-analytics.github.io/create-your-report-groupe/

## Prerequisites

- [R](https://cran.r-project.org/) ≥ 4.3
- [RStudio](https://posit.co/download/rstudio-desktop/) or [VS Code](https://code.visualstudio.com/) with the Quarto extension
- [Quarto CLI](https://quarto.org/docs/get-started/) ≥ 1.4

---

## Data

Data files are **not tracked in this repo**. Place them at the paths expected by `R/00_config.R`:

Download the data folder from: [Download link](https://drive.google.com/file/d/1ZwqZ9QyFR51ZMT_FNu4M2gbkrB5Pw6Z9/view?usp=sharing)

```
data/
├── delft/
│   ├── raster/    # delft_worldpop_proj.tif, ndvi_delft_proj.tif, …
│   └── vector/    # delft_boundary_proj.gpkg, delft_osm_green_proj.gpkg, …
└── Yuexiu/
    ├── raster/    # Yuexiu_viirs_proj.tif, …
    └── vector/    # yuexiu_boundary.gpkg, guangzhou_osm_green_proj.gpkg, …
```

> All paths are centralised in `R/00_config.R` — edit that file if your data lives elsewhere.

---

## How to run

### 1. Clone and open

```bash
git clone https://github.com/Applied-Spatial-Analytics/create-your-report-groupe.git
```

Open `asa2025-report.Rproj` in RStudio, or open the folder in VS Code.

### 2. Install required R packages once

```r
install.packages(c("sf", "ggplot2", "dplyr", "ineq", "tiff", "tidyverse", "patchwork", "here"))
```

### 3. Run the pipeline scripts

```r
source("_entire-workflow.R")  # pipeline that runs all the scripts in order
```

This runs `R/01_load_data.R` through `R/06_mcda_nbs.R` in sequence. Intermediate
`.rds` objects and all `fig_*.png` figures are saved automatically to `outputs/`
(path controlled by `OUT_ROOT` in `R/00_config.R`).

### 4. Render the report

In RStudio: open `report.qmd` and click **Render**.
In VS Code: open `report.qmd` and click **Render**, or from the terminal:

```bash
quarto render report.qmd
```

Output: `docs/index.html` (the `output-dir: docs` setting lives in `_quarto.yml`
at the project root). This is also the folder GitHub Pages serves the live
report from — no manual renaming or copying needed.

---

## File structure

```text
APPLIED_SPATIAL_ANALYTICS/
├── data/                     # not tracked — see Data section above
│   ├── delft/
│   │   ├── raster/
│   │   └── vector/
│   └── Yuexiu/
│       ├── raster/
│       └── vector/
│
├── R/
│   ├── _entire-workflow.R    # Sources all scripts in order
│   ├── 00_config.R           # ALL paths, CRS constants, thresholds — loaded by every script
│   ├── 01_load_data.R        # Load + validate all layers, save checked objects to /outputs/
│   ├── 02_accessibility.R    # SQ1 — green space per capita, buffers, nearest distance
│   ├── 03_typology_biodiversity.R  # SQ2 — OSM typology, NDVI zonal stats, GBIF density
│   ├── 04_spatial_justice.R  # SQ3 — Gini, bivariate choropleth, income/VIIRS correlation
│   ├── 05_connectivity.R     # SQ4 — fragmentation metrics, graph connectivity
│   └── 06_mcda_nbs.R         # SQ5 — MCDA scoring, corridor prioritisation
│
├── outputs/                  # all .rds intermediates + fig_*.png, written by the R scripts
│
├── docs/                     # rendered report — index.html + supporting assets
│                              # this is what GitHub Pages serves
│
├── report_files/             # Quarto's own render cache/dependencies only
│                              # (execute-results/, libs/) — not analysis output
│
├── _quarto.yml                # project config: sets output-dir to docs/
├── report.qmd                 # narrative + figures pulled from outputs/
└── README.md
```

---

## Reproducibility notes

- All file paths live in `R/00_config.R` — no hardcoded paths elsewhere.
- Uses only open, globally available datasets (OSM, WorldPop, GBIF, VIIRS, NDVI).
- The workflow is fully reproducible in any city/district by updating `00_config.R`.
- The report is published via GitHub Pages directly from the `docs/` folder on `main` —
  every render-and-push updates the live version automatically.
