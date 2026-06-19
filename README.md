# Running the pipeline
## Prerequisites
 
- [R](https://cran.r-project.org/) ≥ 4.3
- [RStudio](https://posit.co/download/rstudio-desktop/) or [Quarto CLI](https://quarto.org/docs/get-started/)
- [Quarto](https://quarto.org/) ≥ 1.4
Install required R packages once:
 
```r
install.packages(c("sf", "ggplot2", "dplyr", "ineq", "tiff", "tidyverse", "patchwork"))
```
---

## Data
 
Data files are **not tracked in this repo**. Place them at the paths expected by `R/00_config.R`:
Download data folder from: [Download link](https://drive.google.com/file/d/1ZwqZ9QyFR51ZMT_FNu4M2gbkrB5Pw6Z9/view?usp=sharing)
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

## Visualisation settings
The colour palette used in the project is defined in R/00_config.R under COLORS.

---

## Green space filter

Green spaces are filtered to include only publicly accessible and ecologically functional urban green spaces.
Park/Recreation
Forest/Woodland,
Grass/Meadow,
Nature reserve/Scrub

The green areas that are excluded from all analysis are:
Cemetery
Allotment/Agriculture
Sports facility
Other/Unclassified

---

## How to run
 
### 1. Clone and open
 
```bash
git clone https://github.com/Applied-Spatial-Analytics/create-your-report-groupe.git
```
 
Open `asa2025-report.Rproj` in RStudio.
 
### 2. Run the pipeline scripts in order
 
```r
source("R/00_config.R")               # paths & constants
source("R/01_load_data.R")            # load & validate layers
source("R/02_accessibility.R")        # SQ1 — green space access
source("R/03_typology_biodiversity.R")# SQ2 — typology & NDVI
source("R/04_spatial_justice.R")      # SQ3 — Gini, equity
source("R/05_connectivity.R")         # SQ4 — fragmentation & graph
source("R/06_mcda_nbs.R")             # SQ5 — MCDA & corridors
```
 
Figures are saved automatically to `report_files/`.
 
### 3. Render the report
 
In RStudio: open `report.qmd` and click **Render**, or from the terminal:
 
```bash
quarto render report.qmd
```
 
Output: `report.html`
 
---

## File strcture
```text
APPLIED_SPATIAL_ANALYTICS/
├── data/
│   ├── delft/
│   │   ├── raster/          # delft_worldpop_proj.tif, ndvi_delft_proj.tif, etc.
│   │   └── vector/          # delft_boundary_proj.gpkg, delft_wijken_proj.gpkg, etc.
│   └── Yuexiu/
│       ├── raster/          # Yuexiu_viirs_proj.tif, Yuexiu_worldpop_proj.tif, etc.
│       └── vector/          # Yuexiu_boundary_proj.gpkg, Yuexiu_subdistricts_proj.gpkg, etc.
│
├── R/
│   ├── 00_config.R          # ALL paths, CRS constants, thresholds — loaded by every script
│   ├── 01_load_data.R       # Load + validate all layers, save checked objects to /outputs/
│   ├── 02_accessibility.R   # SQ1 — green space per capita, buffers, nearest distance
│   ├── 03_typology_biodiversity.R  # SQ2 — OSM typology, NDVI zonal stats, GBIF density
│   ├── 04_spatial_justice.R # SQ3 — Gini, bivariate choropleth, income/VIIRS correlation
│   ├── 05_connectivity.R    # SQ4 — fragmentation metrics, graph connectivity
│   └── 06_mcda_nbs.R        # SQ5 — MCDA scoring, corridor prioritisation
│
├── report_files/
│   └── [figures auto-saved here by ggsave()]
│
├── report.qmd               # Narrative + renders figures from report_files/
└── README.md
```

## Reproducibility notes
 
- All file paths live in `R/00_config.R` — no hardcoded paths elsewhere.
- Uses only open, globally available datasets (OSM, WorldPop, GBIF, VIIRS, NDVI).
- The workflow is fully reproducible in any city/district by updating `00_config.R`.
