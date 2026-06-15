source("R/00_config.R")
library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

sf_use_s2(FALSE)

d  <- readRDS(file.path(OUT_ROOT, "yuexiu_data.rds"))
dl <- readRDS(file.path(OUT_ROOT, "delft_data.rds"))

# Unwrapping raster
d$yx_pop   <- rast(d$yx_pop)
d$yx_ndvi  <- rast(d$yx_ndvi)
d$yx_cover <- rast(d$yx_cover)
d$yx_viirs <- rast(d$yx_viirs)

dl$dl_pop   <- rast(dl$dl_pop)
dl$dl_ndvi  <- rast(dl$dl_ndvi)
dl$dl_cover <- rast(dl$dl_cover)

# ── 1A. Green space per capita per subdistrict ────────────────────────────────
# For each admin unit: intersect green polygons, sum area, divide by population

calc_green_per_capita <- function(admin_sf, green_sf, pop_rast, local_crs) {
  # 1. Ensure geometries are valid and non-empty
  admin_sf <- admin_sf[!st_is_empty(admin_sf), ]
  admin_sf <- st_make_valid(admin_sf)

  admin_m <- st_transform(admin_sf, local_crs)
  green_m <- st_transform(green_sf, local_crs)
  green_union <- st_union(green_m)

  results <- lapply(seq_len(nrow(admin_m)), function(i) {
    unit    <- admin_m[i, ]
    clipped <- tryCatch(
      st_intersection(green_union, st_geometry(unit)),
      error = function(e) st_sfc(crs = local_crs)
    )
    area_m2 <- if (length(clipped) > 0 && !st_is_empty(clipped)) {
      as.numeric(st_area(clipped))
    } else { 0 }
    area_m2
  })

  # 2. Project the admin layer ONCE before extracting to avoid issues inside mutate
  admin_for_extract <- st_transform(admin_sf, st_crs(pop_rast))

  admin_sf |>
    mutate(
      green_area_m2 = unlist(results),
      pop_count     = exact_extract(pop_rast, admin_for_extract, "sum"), # Use the pre-transformed layer
      pop_count     = pmax(pop_count, 1, na.rm = TRUE),
      green_pc_m2   = green_area_m2 / pop_count
    )
}

message("Calculating Yuexiu green per capita...")
yx_sub_access <- calc_green_per_capita(d$yx_sub, d$yx_grn, d$yx_pop, CRS_YX)

message("Calculating Delft green per capita...")
dl_wijk_access <- calc_green_per_capita(dl$dl_wijk, dl$dl_grn, dl$dl_pop, CRS_DELFT)

# ── 1B. % population within 300m and 500m of any green space ─────────────────
pct_pop_within_buffer <- function(green_sf, pop_rast, boundary_sf, dist_m, local_crs) {
  green_m    <- st_transform(green_sf, local_crs)
  bnd_m      <- st_transform(boundary_sf, local_crs)
  bnd_pop    <- st_transform(boundary_sf, st_crs(pop_rast))
  buf        <- st_union(st_buffer(green_m, dist_m))
  buf_pop    <- st_transform(st_sf(geometry = buf, crs = local_crs), st_crs(pop_rast))
  clipped    <- tryCatch(st_intersection(bnd_pop, buf_pop), error = function(e) NULL)

  pop_total  <- exact_extract(pop_rast, bnd_pop,  "sum")
  pop_within <- if (!is.null(clipped) && nrow(clipped) > 0) exact_extract(pop_rast, clipped, "sum") else 0

  round(sum(pop_within, na.rm = TRUE) / sum(pop_total,  na.rm = TRUE) * 100, 1)
}

message("Calculating buffer coverage...")
yx_300 <- pct_pop_within_buffer(d$yx_grn,  d$yx_pop, d$yx_bnd, 300, CRS_YX)
yx_500 <- pct_pop_within_buffer(d$yx_grn,  d$yx_pop, d$yx_bnd, 500, CRS_YX)
dl_300 <- pct_pop_within_buffer(dl$dl_grn, dl$dl_pop, dl$dl_bnd, 300, CRS_DELFT)
dl_500 <- pct_pop_within_buffer(dl$dl_grn, dl$dl_pop, dl$dl_bnd, 500, CRS_DELFT)

message(sprintf("Yuexiu: %.1f%% within 300m, %.1f%% within 500m", yx_300, yx_500))
message(sprintf("Delft:  %.1f%% within 300m, %.1f%% within 500m", dl_300, dl_500))

# ── 1C. Mean nearest green space distance per subdistrict ────────────────────
nearest_park_dist <- function(admin_sf, green_sf, local_crs) {
  admin_m  <- st_transform(admin_sf, local_crs)
  green_m  <- st_transform(green_sf, local_crs)
  cents    <- st_centroid(admin_m)
  idx      <- st_nearest_feature(cents, green_m)
  dists    <- st_distance(cents, green_m[idx, ], by_element = TRUE)
  admin_sf |> mutate(nearest_green_m = as.numeric(dists))
}

yx_sub_access  <- nearest_park_dist(yx_sub_access,  d$yx_grn,  CRS_YX)
dl_wijk_access <- nearest_park_dist(dl_wijk_access, dl$dl_grn, CRS_DELFT)

# ── Save ──────────────────────────────────────────────────────────────────────
buffer_summary <- data.frame(
  city     = c("Yuexiu", "Yuexiu", "Delft", "Delft"),
  buffer_m = c(300, 500, 300, 500),
  pct_pop  = c(yx_300, yx_500, dl_300, dl_500)
)

saveRDS(yx_sub_access,   file.path(OUT_ROOT, "yx_sub_access.rds"))
saveRDS(dl_wijk_access,  file.path(OUT_ROOT, "dl_wijk_access.rds"))
saveRDS(buffer_summary,  file.path(OUT_ROOT, "buffer_summary.rds"))

# ── Figures ───────────────────────────────────────────────────────────────────

# Figure 1A: Green space per capita maps (log scale, same for both)
p1 <- ggplot(yx_sub_access) +
  geom_sf(aes(fill = green_pc_m2)) +
  scale_fill_gradient(
    name = "m² per person\n(pseudo-log)",
    trans = scales::pseudo_log_trans(sigma = 1),
    breaks = c(0, 1, 10, 100, 400),
    low = "#f7fcf5",
    high = "#00441b",
    na.value = "grey80") +
  theme_minimal() +
  labs(title = "Green Space per Capita — Yuexiu Jiēdào",
       subtitle = "Source: OSM green polygons + WorldPop")

p2 <- ggplot(dl_wijk_access) +
  geom_sf(aes(fill = green_pc_m2)) +
  scale_fill_gradient(
    name = "m² per person\n(pseudo-log)",
    trans = scales::pseudo_log_trans(sigma = 1),
    breaks = c(0, 1, 10, 100, 400),
    low = "#f7fcf5",
    high = "#00441b",
    na.value = "grey80") +
  theme_minimal() +
  labs(title = "Green Space per Capita — Delft Wijken",
       subtitle = "Source: OSM green polygons + WorldPop")

ggsave(file.path(OUT_ROOT, "fig_access_per_capita.png"),
       p1 + p2, width = 14, height = 6, dpi = 300)

# Figure 1B: Mean nearest green distance maps
p3 <- ggplot(yx_sub_access) +
  geom_sf(aes(fill = nearest_green_m)) +
  scale_fill_viridis_c(name = "Distance (m)", option = "B", direction = -1) +
  theme_minimal() +
  labs(title = "Mean Distance to Nearest Green Space — Yuexiu",
       subtitle = "From subdistrict centroid")

p4 <- ggplot(dl_wijk_access) +
  geom_sf(aes(fill = nearest_green_m)) +
  scale_fill_viridis_c(name = "Distance (m)", option = "B", direction = -1) +
  theme_minimal() +
  labs(title = "Mean Distance to Nearest Green Space — Delft",
       subtitle = "From wijk centroid")

ggsave(file.path(OUT_ROOT, "fig_nearest_green_distance.png"),
       p3 + p4, width = 14, height = 6, dpi = 300)

# Figure 1C: Buffer coverage bar chart
p5 <- ggplot(buffer_summary, aes(x = factor(buffer_m), y = pct_pop, fill = city)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  scale_fill_manual(values = c("Yuexiu" = "#21918c", "Delft" = "#440154")) +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  theme_minimal(base_size = 12) +
  labs(title = "Population within Walking Distance of Green Space",
       subtitle = "300m and 500m Euclidean buffer catchments",
       x = "Buffer distance (m)", y = "% of population covered",
       fill = "City")

ggsave(file.path(OUT_ROOT, "fig_buffer_coverage.png"), p5, width = 8, height = 5, dpi = 300)

# ── Figure 1D: Network Entry Gates & Multi-Ring Walking Buffers Panel ───────────
library(sf)
library(ggplot2)
library(patchwork)

# 1. Redefine the function completely and correctly
generate_network_access <- function(green_sf, road_sf, local_crs) {
  # Create projected spatial datasets sequentially
  green_m <- st_transform(green_sf, local_crs)
  road_m  <- st_transform(road_sf, local_crs)

  # Build multi-ring buffers
  buf_500  <- st_union(st_buffer(green_m, 500))
  buf_300  <- st_union(st_buffer(green_m, 300))
  buf_100  <- st_union(st_buffer(green_m, 100))

  # Extract point intersections along the 300m threshold boundary line
  buf_300_line <- st_cast(buf_300, "MULTILINESTRING")
  access_pts   <- st_intersection(road_m, buf_300_line)
  access_pts   <- access_pts[st_geometry_type(access_pts) %in% c("POINT", "MULTIPOINT"), ]

  list(b100 = buf_100, b300 = buf_300, b500 = buf_500, roads = road_m, green = green_m, pts = access_pts)
}

# 2. Run calculations using the ACTUAL dataset list names (rds instead of roads)
yx_net <- generate_network_access(d$yx_grn, d$yx_rds, CRS_YX)
dl_net <- generate_network_access(dl$dl_grn, dl$dl_rds, CRS_DELFT)

# 3. Build the Plots
p_yx_network <- ggplot() +
  geom_sf(data = yx_net$roads, color = "grey85", size = 0.3) +
  geom_sf(data = yx_net$b500, fill = "#238b45", alpha = 0.08, color = NA) +
  geom_sf(data = yx_net$b300, fill = "#238b45", alpha = 0.15, color = "#238b45", linetype = "dashed", size = 0.4) +
  geom_sf(data = yx_net$b100, fill = "#238b45", alpha = 0.25, color = NA) +
  geom_sf(data = yx_net$green, fill = "#00441b", color = NA) +
  geom_sf(data = yx_net$pts, color = "#d95f02", size = 1.2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Network Entry Thresholds — Yuexiu",
       subtitle = "Orange points indicate road intersections at the 300m buffer boundary")

p_dl_network <- ggplot() +
  geom_sf(data = dl_net$roads, color = "grey85", size = 0.3) +
  geom_sf(data = dl_net$b500, fill = "#238b45", alpha = 0.08, color = NA) +
  geom_sf(data = dl_net$b300, fill = "#238b45", alpha = 0.15, color = "#238b45", linetype = "dashed", size = 0.4) +
  geom_sf(data = dl_net$b100, fill = "#238b45", alpha = 0.25, color = NA) +
  geom_sf(data = dl_net$green, fill = "#00441b", color = NA) +
  geom_sf(data = dl_net$pts, color = "#d95f02", size = 1.2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Network Entry Thresholds — Delft",
       subtitle = "Orange points indicate road intersections at the 300m buffer boundary")

# 4. Save combined map side-by-side
ggsave(file.path(OUT_ROOT, "fig_network_walk_access.png"),
       p_yx_network + p_dl_network, width = 14, height = 6, dpi = 300)

message("Figure 1D saved successfully!")

