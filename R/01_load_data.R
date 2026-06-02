source("R/00_config.R")
library(sf); library(terra)

# Load all vector layers
yx_bnd  <- read_sf(YX$boundary)
yx_sub  <- read_sf(YX$subdistricts)   # 18 Jiēdào
yx_grn  <- read_sf(YX$green)
yx_rds  <- read_sf(YX$roads)
yx_gbif <- read_sf(YX$gbif)
yx_viirs_pts <- read_sf(YX$viirs_pts)

dl_bnd  <- read_sf(DL$boundary)
dl_wijk <- read_sf(DL$wijken)          # Wijken (confirm: buurten or wijken?)
dl_inc  <- read_sf(DL$income)
dl_grn  <- read_sf(DL$green)
dl_rds  <- read_sf(DL$roads)
dl_gbif <- read_sf(DL$gbif)

# Load rasters
yx_pop   <- rast(YX$worldpop)
yx_ndvi  <- rast(YX$ndvi)
yx_cover <- rast(YX$worldcover)
dl_pop   <- rast(DL$worldpop)
dl_ndvi  <- rast(DL$ndvi)
dl_cover <- rast(DL$worldcover)

# Validate CRS
stopifnot(st_crs(yx_bnd)$epsg == CRS_YX)
stopifnot(st_crs(dl_bnd)$epsg == CRS_DELFT)

# Filter: drop tiny slivers below minimum patch size
yx_grn <- yx_grn |> mutate(area_ha = as.numeric(st_area(geom)) / 10000) |>
  filter(area_ha >= MIN_PATCH_HA)
dl_grn <- dl_grn |> mutate(area_ha = as.numeric(st_area(geom)) / 10000) |>
  filter(area_ha >= MIN_PATCH_HA)

# Save validated objects for use in downstream scripts
saveRDS(list(yx_bnd=yx_bnd, yx_sub=yx_sub, yx_grn=yx_grn, yx_rds=yx_rds,
             yx_gbif=yx_gbif, yx_viirs_pts=yx_viirs_pts,
             yx_pop=yx_pop, yx_ndvi=yx_ndvi, yx_cover=yx_cover),
        file.path(OUT_ROOT, "yuexiu_data.rds"))

saveRDS(list(dl_bnd=dl_bnd, dl_wijk=dl_wijk, dl_inc=dl_inc, dl_grn=dl_grn,
             dl_rds=dl_rds, dl_gbif=dl_gbif,
             dl_pop=dl_pop, dl_ndvi=dl_ndvi, dl_cover=dl_cover),
        file.path(OUT_ROOT, "delft_data.rds"))

message("Data loaded and validated.")