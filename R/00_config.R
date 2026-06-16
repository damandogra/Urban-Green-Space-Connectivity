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

# ── Network Entry Access Helper Function ──────────────────────────────────────
generate_network_access <- function(green_sf, roads_sf, local_crs) {
  # 1. Properly project inputs to local CRS so metric distance functions work
  green_m <- sf::st_transform(green_sf, local_crs)
  roads_m <- sf::st_transform(roads_sf, local_crs)

  # 2. Extract road intersection points (vertices) to simulate network entry access
  # (Resolves the 'green_m' not found error by establishing variables sequentially)
  road_points <- sf::st_cast(roads_m, "POINT")

  # 3. Create the multi-ring walking catchments (100m, 300m, 500m)
  b100 <- sf::st_union(sf::st_buffer(green_m, 100))
  b300 <- sf::st_union(sf::st_buffer(green_m, 300))
  b500 <- sf::st_union(sf::st_buffer(green_m, 500))

  # 4. Filter for access points that intersect along the 300m walking perimeter perimeter
  intersecting_points <- road_points[sf::st_intersects(road_points, b300, sparse = FALSE), ]

  # Return a structured list matching your main script's expected data bindings
  list(
    roads = roads_m,
    green = green_m,
    b100  = b100,
    b300  = b300,
    b500  = b500,
    pts   = intersecting_points
  )
}

