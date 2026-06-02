# Report template for the Applied Spatial Analytics 2026 course

This is a template repository used as a starting point for the group
reports produced in the Applied Spatial Analytics 2025 course at
TU Delft.

By starting this assignment in GitHub Classroom, you created a copy
of this repository that you have write access to. You will continue to
work on your report in that repository throughout the quarter. A great
way to practice and apply what you learned in the "Intro to Git and GitHub"
assignment!

### Getting started with the report

1. In RStudio, create a new project from version control. Use the
   URL of your repository to clone it.
   
2. To start working on the report, open the `report.qmd` file, add the
   names of your group members, and press on the "Render" button. Next,
   stage and commit all changed files and push them to GitHub. You will
   follow this **stage -> commit -> push** workflow every time you make a
   change.

### Feedback

In the **Pull requests** section of your repository, you will find a
**Feedback** pull request. We will use this pull request to provide
feedback on your report throughout the quarter. You can also use this
pull request to ask questions about the feedback.

### Asking for help

If you have questions about the assignment, please ask them in
[Discussions](https://github.com/Applied-Spatial-Analytics/asa2026/discussions).


#### File strcture
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
├── test.R                   # Quick scratch / debugging
└── README.md
```
