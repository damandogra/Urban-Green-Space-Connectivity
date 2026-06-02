# R/00_config.R
library(here)

# ── Root-relative paths ──────────────────────────────────────────────────────
DATA_ROOT  <- here("data")
OUT_ROOT   <- here("report_files")
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

# ── Yuexiu (Guangzhou) paths ─────────────────────────────────────────────────
YX <- list(
  boundary     = file.path(DATA_ROOT, "Yuexiu/vector/Yuexiu_boundary_proj.gpkg"),
  subdistricts = file.path(DATA_ROOT, "Yuexiu/vector/Yuexiu_subdistricts_proj.gpkg"),
  green        = file.path(DATA_ROOT, "Yuexiu/vector/Yuexiu_osm_green_proj.gpkg"),
  roads        = file.path(DATA_ROOT, "Yuexiu/vector/Yuexiu_osm_roads_proj.gpkg"),
  water_poly   = file.path(DATA_ROOT, "Yuexiu/vector/Yuexiu_osm_water_polygon_proj.gpkg"),
  water_line   = file.path(DATA_ROOT, "Yuexiu/vector/Yuexiu_osm_water_proj.gpkg"),
  gbif         = file.path(DATA_ROOT, "Yuexiu/vector/Yuexiu_gbif_proj.gpkg"),
  viirs_pts    = file.path(DATA_ROOT, "Yuexiu/vector/Yuexiu_VIIRS_Points_Joined_proj.gpkg"),
  worldpop     = file.path(DATA_ROOT, "Yuexiu/raster/Yuexiu_worldpop_proj.tif"),
  worldcover   = file.path(DATA_ROOT, "Yuexiu/raster/Yuexiu_worldcover_proj.tif"),
  ndvi         = file.path(DATA_ROOT, "Yuexiu/raster/Yuexiu_ndvi_proj.tif"),
  viirs_rast   = file.path(DATA_ROOT, "Yuexiu/raster/Yuexiu_viirs_proj.tif")
)

# ── Delft paths ──────────────────────────────────────────────────────────────
DL <- list(
  boundary   = file.path(DATA_ROOT, "delft/vector/delft_boundary_proj.gpkg"),
  wijken     = file.path(DATA_ROOT, "delft/vector/delft_wijken_proj.gpkg"),   # NOTE: confirm layer name
  income     = file.path(DATA_ROOT, "delft/vector/delft_income_proj.gpkg"),
  green      = file.path(DATA_ROOT, "delft/vector/delft_osm_green_proj.gpkg"),  # you'll need to add this
  roads      = file.path(DATA_ROOT, "delft/vector/delft_osm_roads_proj.gpkg"),
  gbif       = file.path(DATA_ROOT, "delft/vector/delft_gbif_proj.gpkg"),
  worldpop   = file.path(DATA_ROOT, "delft/raster/delft_worldpop_proj.tif"),
  worldcover = file.path(DATA_ROOT, "delft/raster/delft_worldcover_proj.tif"),
  ndvi       = file.path(DATA_ROOT, "delft/raster/ndvi_delft_proj.tif")
)

# ── CRS ──────────────────────────────────────────────────────────────────────
CRS_DELFT  <- 28992   # Amersfoort RD New
CRS_YX     <- 4547    # UTM Zone 49N (adjust to whatever your data uses)

# ── Analysis thresholds ──────────────────────────────────────────────────────
BUFFER_300M <- 300
BUFFER_500M <- 500
MIN_PATCH_HA <- 0.1        # minimum green patch size to include
DISPERSAL_THRESH_M <- 150  # for connectivity graph edges

# ── MCDA weights (sum to 1) ───────────────────────────────────────────────────
# Adjust these weights after group discussion
MCDA_WEIGHTS <- c(
  accessibility  = 0.30,
  biodiversity   = 0.25,
  connectivity   = 0.25,
  social_equity  = 0.20
)