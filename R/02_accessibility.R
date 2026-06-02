Sys.setenv(PROJ_LIB = "C:/Program Files/R/R-4.6.0/library/sf/proj")

source("R/00_config.R")
library(sf); library(terra); library(exactextractr)
library(dplyr); library(ggplot2)
sf_use_s2(FALSE)

d  <- readRDS(file.path(OUT_ROOT, "yuexiu_data.rds"))
dl <- readRDS(file.path(OUT_ROOT, "delft_data.rds"))

# ── Reload rasters and reproject everything to WGS84 ─────────────────────────
d$yx_pop  <- rast(YX$worldpop)
dl$dl_pop <- rast(DL$worldpop)

d$yx_sub  <- st_transform(st_make_valid(d$yx_sub), 4326)
d$yx_sub  <- d$yx_sub[!is.na(st_dimension(d$yx_sub)), ]
d$yx_sub  <- d$yx_sub[st_dimension(d$yx_sub) == 2, ]
d$yx_sub  <- st_collection_extract(d$yx_sub, "POLYGON")
d$yx_grn  <- st_transform(st_make_valid(d$yx_grn), 4326)
d$yx_bnd  <- st_transform(d$yx_bnd, 4326)
d$yx_pop  <- project(d$yx_pop, "EPSG:4326")

dl$dl_wijk <- st_transform(st_make_valid(dl$dl_wijk), 4326)
dl$dl_grn  <- st_transform(st_make_valid(dl$dl_grn),  4326)
dl$dl_bnd  <- st_transform(dl$dl_bnd,  4326)
dl$dl_pop  <- project(dl$dl_pop, "EPSG:4326")

# ── Helper functions ──────────────────────────────────────────────────────────
pct_pop_within_buffer <- function(green_sf, pop_rast, boundary_sf, dist_m) {
  buf        <- st_union(st_buffer(green_sf, dist_m))
  pop_total  <- exact_extract(pop_rast, boundary_sf, "sum")
  pop_within <- exact_extract(pop_rast, st_intersection(boundary_sf, buf), "sum")
  pop_within / pop_total * 100
}

calc_green_per_capita <- function(admin_sf, green_sf, pop_rast, id_col) {
  green_union <- st_union(green_sf)
  green_area_m2 <- sapply(st_geometry(admin_sf), function(g) {
    g     <- st_set_crs(st_sfc(g), st_crs(admin_sf))
    inter <- suppressWarnings(st_intersection(g, green_union))
    if (length(inter) == 0) return(0)
    as.numeric(st_area(inter))
  })
  admin_sf |>
    mutate(
      green_area_m2 = green_area_m2,
      pop_count     = exact_extract(pop_rast, admin_sf, "sum"),
      green_pc_m2   = green_area_m2 / pop_count
    )
}

nearest_park_dist <- function(admin_sf, green_sf) {
  centroids   <- st_centroid(admin_sf)
  nearest_idx <- st_nearest_feature(centroids, green_sf)
  dists       <- st_distance(centroids, green_sf[nearest_idx, ], by_element = TRUE)
  admin_sf |> mutate(mean_dist_to_green_m = as.numeric(dists))
}

# ── 1A. Green space per capita ────────────────────────────────────────────────
yx_sub_access  <- calc_green_per_capita(d$yx_sub,   d$yx_grn,  d$yx_pop,  "name")
dl_wijk_access <- calc_green_per_capita(dl$dl_wijk, dl$dl_grn, dl$dl_pop, "wijknaam")

# ── 1B. % population within 300m and 500m ────────────────────────────────────
yx_300 <- pct_pop_within_buffer(d$yx_grn,  d$yx_pop,  d$yx_bnd,  300)
yx_500 <- pct_pop_within_buffer(d$yx_grn,  d$yx_pop,  d$yx_bnd,  500)
dl_300 <- pct_pop_within_buffer(dl$dl_grn, dl$dl_pop, dl$dl_bnd, 300)
dl_500 <- pct_pop_within_buffer(dl$dl_grn, dl$dl_pop, dl$dl_bnd, 500)

# ── 1C. Nearest park distance ─────────────────────────────────────────────────
yx_sub_access  <- nearest_park_dist(yx_sub_access,  d$yx_grn)
dl_wijk_access <- nearest_park_dist(dl_wijk_access, dl$dl_grn)

# ── Save outputs ──────────────────────────────────────────────────────────────
saveRDS(yx_sub_access,  file.path(OUT_ROOT, "yx_sub_access.rds"))
saveRDS(dl_wijk_access, file.path(OUT_ROOT, "dl_wijk_access.rds"))

# ── Figures (Data Intact with Visual Cap) ─────────────────────────────────────
p1 <- ggplot(yx_sub_access) +
  geom_sf(aes(fill = green_pc_m2)) +
  scale_fill_viridis_c(
    name = "m² per person\n(Log Scale)", 
    option = "G", 
    trans = "log10",
    labels = function(x) sprintf("%.1f", x)
  ) +
  theme_minimal() +
  labs(title = "Green Space per Capita — Yuexiu Jiēdào")

p2 <- ggplot(dl_wijk_access) +
  geom_sf(aes(fill = green_pc_m2)) +
  scale_fill_viridis_c(
    name = "m² per person", 
    option = "G",
    limits = c(0, 100),            # Caps the color ramp scale at 100 m²
    oob = scales::squish           # Compresses outliers into the top color block
  ) +
  theme_minimal() +
  labs(title = "Green Space per Capita — Delft Wijken")

ggsave(file.path(OUT_ROOT, "fig_access_per_capita.png"),
       p1 + p2, width = 14, height = 6)

message("Accessibility analysis complete.")