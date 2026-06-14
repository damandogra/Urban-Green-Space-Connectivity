library(here)

# ── Output directory ──────────────────────────────────────────────────────────
OUT_ROOT <- here("report_files")
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

DATA_ROOT <- here("data")

# ── Yuexiu (Guangzhou) ────────────────────────────────────────────────────────
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

# ── Delft (Netherlands) ───────────────────────────────────────────────────────
DL <- list(
  boundary   = file.path(DATA_ROOT, "delft/vector/delft_boundary_proj.gpkg"),
  wijken     = file.path(DATA_ROOT, "delft/vector/delft_wijken_proj.gpkg"),
  buurten    = file.path(DATA_ROOT, "delft/vector/delft_buurten_proj.gpkg"),
  income     = file.path(DATA_ROOT, "delft/vector/delft_income_proj.gpkg"),
  green      = file.path(DATA_ROOT, "delft/vector/delft_osm_green_proj.gpkg"),
  roads      = file.path(DATA_ROOT, "delft/vector/delft_osm_roads_proj.gpkg"),
  water      = file.path(DATA_ROOT, "delft/vector/delft_osm_water_proj.gpkg"),
  gbif       = file.path(DATA_ROOT, "delft/vector/delft_gbif_proj.gpkg"),
  worldpop   = file.path(DATA_ROOT, "delft/raster/delft_worldpop_proj.tif"),
  worldcover = file.path(DATA_ROOT, "delft/raster/delft_worldcover_proj.tif"),
  ndvi       = file.path(DATA_ROOT, "delft/raster/ndvi_delft_proj.tif")
)

# ── CRS ───────────────────────────────────────────────────────────────────────
CRS_YX   <- 32649   # WGS 84 / UTM Zone 49N
CRS_DELFT <- 28992  # Amersfoort / RD New

# ── Analysis constants ────────────────────────────────────────────────────────
BUFFER_300M        <- 300
BUFFER_500M        <- 500
MIN_PATCH_HA       <- 0.1   # drop green patches smaller than this
DISPERSAL_THRESH_M <- 150   # max gap for connectivity graph edges
WATER_LINE_BUFFER  <- 5     # metres to buffer water lines into polygons

# ── MCDA weights (must sum to 1) ──────────────────────────────────────────────
MCDA_WEIGHTS <- c(
  accessibility = 0.30,
  biodiversity  = 0.25,
  connectivity  = 0.25,
  equity        = 0.20
)

