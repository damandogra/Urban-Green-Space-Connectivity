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
  admin_for_extract <- st_transform(admin_sf, local_crs)

  admin_sf |>
    mutate(
      green_area_m2 = unlist(results),
      pop_count     = exact_extract(pop_rast, admin_for_extract, "sum"),
      pop_count     = pmax(pop_count, 1, na.rm = TRUE),
      green_pc_m2   = green_area_m2 / pop_count
    )
}

# Unpack population rasters after loading
d$yx_pop <- terra::rast(d$yx_pop)
dl$dl_pop <- terra::rast(dl$dl_pop)

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

# Figure 0: Population density maps
# Residents per km², calculated from extracted population and polygon area

yx_sub_access <- yx_sub_access %>%
  mutate(
    area_km2 = as.numeric(st_area(.)) / 1e6,
    pop_density_km2 = pop_count / area_km2
  )

dl_wijk_access <- dl_wijk_access %>%
  mutate(
    area_km2 = as.numeric(st_area(.)) / 1e6,
    pop_density_km2 = pop_count / area_km2
  )

density_breaks <- c(0, 1000, 2500, 5000, 10000, 20000, 40000, 80000)

density_labels <- c(
  "0–1k",
  "1k–2.5k",
  "2.5k–5k",
  "5k–10k",
  "10k–20k",
  "20k–40k",
  "40k–80k"
)

density_cols <- c(
  COLORS$beige,
  "#e7d6bc",
  "#d8b48a",
  "#c98f5f",
  COLORS$red_light,
  "#7f3f3e",
  COLORS$red
)

yx_sub_access$pop_density_class <- cut(
  yx_sub_access$pop_density_km2,
  breaks = density_breaks,
  labels = density_labels,
  include.lowest = TRUE,
  right = TRUE
)

dl_wijk_access$pop_density_class <- cut(
  dl_wijk_access$pop_density_km2,
  breaks = density_breaks,
  labels = density_labels,
  include.lowest = TRUE,
  right = TRUE
)

p_den_yx <- ggplot(yx_sub_access) +
  geom_sf(aes(fill = pop_density_class), show.legend = TRUE) +
  scale_fill_manual(
    name = "residents/km²",
    values = setNames(density_cols, density_labels),
    limits = density_labels,
    drop = FALSE,
    na.value = COLORS$grey85,
    guide = guide_legend(reverse = TRUE)
  ) +
  theme_minimal() +
  labs(
    title = "Population Density – Yuexiu Jiēdào",
    subtitle = "Source: WorldPop + administrative polygons"
  )

p_den_dl <- ggplot(dl_wijk_access) +
  geom_sf(aes(fill = pop_density_class), show.legend = TRUE) +
  scale_fill_manual(
    name = "residents/km²",
    values = setNames(density_cols, density_labels),
    limits = density_labels,
    drop = FALSE,
    na.value = COLORS$grey85,
    guide = guide_legend(reverse = TRUE)
  ) +
  theme_minimal() +
  labs(
    title = "Population Density – Delft Wijken",
    subtitle = "Source: WorldPop + wijk polygons"
  )

ggsave(
  file.path(OUT_ROOT, "fig_population_density.png"),
  p_den_yx + p_den_dl,
  width = 14,
  height = 6,
  dpi = 300
)

# Figure 1A: Green space per capita maps (log scale, same for both)
p1 <- ggplot(yx_sub_access) +
  geom_sf(aes(fill = green_pc_m2)) +
  scale_fill_stepsn(
    name = "m² per person",

    colours = c(
      "#f4ecd0",  # 0–0.1
      "#e6dfbb",  # 0.1–0.5
      "#d5d0a5",  # 0.5–1
      "#c3c193",  # 1–2.5
      "#b1b181",  # 2.5–5
      "#92975f",  # 5–10
      "#5c612f",  # 10–25
      "#252f18",  # 25–100
      "#111809"   # 100–500
    ),

    values = scales::rescale(
      c(0, 0.1, 0.5, 1, 2.5, 5, 10, 25, 100, 500)
    ),

    breaks = c(0, 0.1, 0.5, 1, 2.5, 5, 10, 25, 100, 500),
    limits = c(0, 500),
    na.value = "grey85"
  ) +
  theme_minimal() +
  labs(title = "Green Space per Capita — Yuexiu Jiēdào",
       subtitle = "Source: OSM green polygons + WorldPop")

p2 <- ggplot(dl_wijk_access) +
  geom_sf(aes(fill = green_pc_m2)) +
  scale_fill_stepsn(
    name = "m² per person",

    colours = c(
      "#f4ecd0",  # 0–0.1
      "#e6dfbb",  # 0.1–0.5
      "#d5d0a5",  # 0.5–1
      "#c3c193",  # 1–2.5
      "#b1b181",  # 2.5–5
      "#92975f",  # 5–10
      "#5c612f",  # 10–25
      "#252f18",  # 25–100
      "#111809"   # 100–500
    ),

    values = scales::rescale(
      c(0, 0.1, 0.5, 1, 2.5, 5, 10, 25, 100, 500)
    ),

    breaks = c(0, 0.1, 0.5, 1, 2.5, 5, 10, 25, 100, 500),
    limits = c(0, 500),
    na.value = "grey85"
  ) +
  theme_minimal() +
  labs(title = "Green Space per Capita — Delft Wijken",
       subtitle = "Source: OSM green polygons + WorldPop")

ggsave(file.path(OUT_ROOT, "fig_access_per_capita.png"),
       p1 + p2, width = 14, height = 6, dpi = 300)

# Figure 1B: Mean nearest green distance maps
p3 <- ggplot(yx_sub_access) +
  geom_sf(aes(fill = nearest_green_m)) +
  scale_fill_gradient(
    name = "Distance (m)",
    low  = COLORS$pink_light,
    high = COLORS$red
  ) +
  theme_minimal() +
  labs(title = "Mean Distance to Nearest Green Space — Yuexiu",
       subtitle = "From subdistrict centroid")

p4 <- ggplot(dl_wijk_access) +
  geom_sf(aes(fill = nearest_green_m)) +
  scale_fill_gradient(
    name = "Distance (m)",
    low  = COLORS$pink_light,
    high = COLORS$red
  ) +
  theme_minimal() +
  labs(title = "Mean Distance to Nearest Green Space — Delft",
       subtitle = "From wijk centroid")

ggsave(file.path(OUT_ROOT, "fig_nearest_green_distance.png"),
       p3 + p4, width = 14, height = 6, dpi = 300)

# Figure 1C: Buffer coverage bar chart
p5 <- ggplot(buffer_summary, aes(x = factor(buffer_m), y = pct_pop, fill = city)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  scale_fill_manual(values = c("Yuexiu" = COLORS$orange, "Delft" = COLORS$blue)) +
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

  # Entry points: roads near green-space boundaries
  entry_buffer <- st_buffer(st_boundary(green_m), 3)

  road_m_filtered <- road_m[
    road_m$highway %in% c(
      "primary", "secondary", "tertiary",
      "residential", "living_street", "unclassified",
      "service",
      "pedestrian", "cycleway", "footway"
    ),
  ]

  # Remove only deep interior park paths, but keep edge paths
  green_core <- st_buffer(st_union(green_m), -10)

  roads_outside_green <- st_difference(road_m_filtered, green_core)

  # Road segments close to green-space boundary
  entry_segments <- st_intersection(roads_outside_green, entry_buffer)

  # Convert candidate edge road segments to points
  access_pts_raw <- st_centroid(entry_segments)

  # Cluster nearby candidate points into representative entry points
  access_pts <- access_pts_raw |>
    st_buffer(35) |>
    st_union() |>
    st_cast("POLYGON") |>
    st_centroid() |>
    st_as_sf()

  list(
    b100 = buf_100,
    b300 = buf_300,
    b500 = buf_500,
    entry_buffer = entry_buffer,
    entry_segments = entry_segments,
    roads = road_m,
    green = green_m,
    pts = access_pts
  )
}

# 2. Run calculations using the ACTUAL dataset list names (rds instead of roads)
yx_net <- generate_network_access(d$yx_grn, d$yx_rds, CRS_YX)
dl_net <- generate_network_access(dl$dl_grn, dl$dl_rds, CRS_DELFT)

# 3. Build the Plots
p_yx_network <- ggplot() +
  geom_sf(data = yx_net$roads, color = "grey85", size = 0.3) +
  geom_sf(data = yx_net$b500, fill = COLORS$pink_light, alpha = 0.08, color = NA) +
  geom_sf(data = yx_net$entry_buffer, fill = COLORS$blue_light, alpha = 0.25, color = COLORS$beige, linetype = "dashed", size = 0.4) +
  geom_sf(data = yx_net$b100, fill = COLORS$pink, alpha = 0.25, color = NA) +
  geom_sf(data = yx_net$green, fill = COLORS$green_dark, color = NA) +
  geom_sf(data = yx_net$pts, color = COLORS$red, size = 0.4, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Network Entry Thresholds — Yuexiu",
       subtitle = "Red points indicate clustered road-boundary intersections as potential green-space entry points"
  )

p_dl_network <- ggplot() +
  geom_sf(data = dl_net$roads, color = "grey85", size = 0.3) +
  geom_sf(data = dl_net$b500, fill = COLORS$pink_light, alpha = 0.08, color = NA) +
  geom_sf(data = dl_net$entry_buffer, fill = COLORS$blue_light, alpha = 0.25, color = COLORS$beige, linetype = "dashed", size = 0.4) +
  geom_sf(data = dl_net$b100, fill = COLORS$pink, alpha = 0.25, color = NA) +
  geom_sf(data = dl_net$green, fill = COLORS$green_dark, color = NA) +
  geom_sf(data = dl_net$pts, color = COLORS$red, size = 0.4, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Network Entry Thresholds — Delft",
       subtitle = "Red points indicate clustered road-boundary intersections as potential green-space entry points"
  )

# 4. Save combined map side-by-side
ggsave(file.path(OUT_ROOT, "fig_network_walk_access.png"),
       p_yx_network + p_dl_network, width = 14, height = 6, dpi = 300)

message("Figure 1D saved successfully!")

