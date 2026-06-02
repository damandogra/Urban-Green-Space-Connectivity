source("R/00_config.R")
library(sf); library(terra); library(exactextractr)
library(dplyr); library(ggplot2); library(patchwork); library(scales)

# Enforce planar calculations for accurate metric projections
sf_use_s2(FALSE)

d  <- readRDS(file.path(OUT_ROOT, "yuexiu_data.rds"))
dl <- readRDS(file.path(OUT_ROOT, "delft_data.rds"))

# Local metric systems for planar spatial arithmetic
yx_crs <- "EPSG:32649" # UTM Zone 49N (Guangzhou - meters)
dl_crs <- "EPSG:28992" # Amersfoort / RD New (Dutch Grid - meters)

# ── Reload NDVI Rasters fresh from source paths ──────────────────────────────
d_ndvi_rast  <- rast(YX$ndvi)
dl_ndvi_rast <- rast(DL$ndvi)

# ── 2A. Green space typology from OSM tags (Safe from missing columns) ───────
classify_green_type <- function(green_sf) {
  # Dynamically pull columns if they exist, fallback to empty characters if missing
  leisure_vec <- if ("leisure" %in% names(green_sf)) green_sf$leisure else ""
  landuse_vec <- if ("landuse" %in% names(green_sf)) green_sf$landuse else ""
  
  green_sf |>
    mutate(green_type = case_when(
      grepl("park",         leisure_vec, ignore.case = TRUE) ~ "Park",
      grepl("garden",       leisure_vec, ignore.case = TRUE) ~ "Garden",
      grepl("nature_reserve|forest", landuse_vec, ignore.case = TRUE) ~ "Nature/Forest",
      grepl("grass|meadow", landuse_vec, ignore.case = TRUE) ~ "Grass/Meadow",
      grepl("cemetery",     landuse_vec, ignore.case = TRUE) ~ "Cemetery",
      TRUE ~ "Other green"
    ))
}

# ── 2B. NDVI zonal statistics per green patch ────────────────────────────────
d$yx_grn  <- st_transform(st_make_valid(d$yx_grn), st_crs(d_ndvi_rast))
dl$dl_grn <- st_transform(st_make_valid(dl$dl_grn), st_crs(dl_ndvi_rast))

yx_grn_ndvi <- d$yx_grn |>
  mutate(ndvi_mean = exact_extract(d_ndvi_rast, d$yx_grn, "mean", progress = FALSE),
         ndvi_sd   = exact_extract(d_ndvi_rast, d$yx_grn, "stdev", progress = FALSE))

dl_grn_ndvi <- dl$dl_grn |>
  mutate(ndvi_mean = exact_extract(dl_ndvi_rast, dl$dl_grn, "mean", progress = FALSE),
         ndvi_sd   = exact_extract(dl_ndvi_rast, dl$dl_grn, "stdev", progress = FALSE))

# ── 2C. GBIF species density per green patch (species per hectare) ───────────
species_density <- function(green_sf, gbif_sf, local_metric_crs) {
  green_metric <- st_transform(st_make_valid(green_sf), local_metric_crs)
  gbif_metric  <- st_transform(st_make_valid(gbif_sf), local_metric_crs)
  
  counts <- lengths(st_intersects(green_metric, gbif_metric))
  
  green_sf |>
    mutate(
      gbif_count     = counts,
      area_ha        = as.numeric(st_area(st_geometry(green_metric))) / 10000, 
      species_per_ha = ifelse(area_ha > 0, gbif_count / area_ha, 0)
    )
}

yx_grn_bio <- species_density(yx_grn_ndvi, d$yx_gbif, yx_crs)
dl_grn_bio <- species_density(dl_grn_ndvi, dl$dl_gbif, dl_crs)

# ── 2D. Blue-green ratio per subdistrict ─────────────────────────────────────
cover_rast_yx <- rast(YX$worldcover)
cover_rast_dl <- rast(DL$worldcover)

d$yx_sub <- d$yx_sub |> 
  st_make_valid() |> 
  filter(!st_is_empty(st_geometry(d$yx_sub))) |> 
  st_transform(st_crs(cover_rast_yx))

dl$dl_wijk <- dl$dl_wijk |> 
  st_make_valid() |> 
  filter(!st_is_empty(st_geometry(dl$dl_wijk))) |> 
  st_transform(st_crs(cover_rast_dl))

bg_ratio <- function(admin_sf, cover_rast) {
  admin_sf |>
    mutate(
      pct_green = exact_extract(cover_rast, admin_sf,
        fun = function(values, cov) mean(values %in% c(10, 20, 30, 40), na.rm = TRUE) * 100, progress = FALSE),
      pct_water = exact_extract(cover_rast, admin_sf,
        fun = function(values, cov) mean(values == 80, na.rm = TRUE) * 100, progress = FALSE),
      bg_ratio  = pct_green / (pct_water + 0.001)
    )
}

# Process calculations across spatial administrative bounds
yx_sub_bg  <- bg_ratio(d$yx_sub, cover_rast_yx)
dl_wijk_bg <- bg_ratio(dl$dl_wijk, cover_rast_dl)

# Apply qualitative labels safely using our updated classification step
yx_grn_bio <- classify_green_type(yx_grn_bio)
dl_grn_bio <- classify_green_type(dl_grn_bio)

# ── Save Completed R Objects ──────────────────────────────────────────────────
saveRDS(yx_grn_bio, file.path(OUT_ROOT, "yx_grn_bio.rds"))
saveRDS(dl_grn_bio, file.path(OUT_ROOT, "dl_grn_bio.rds"))
saveRDS(yx_sub_bg,  file.path(OUT_ROOT, "yx_sub_bg.rds"))
saveRDS(dl_wijk_bg, file.path(OUT_ROOT, "dl_wijk_bg.rds"))


# ── 2E. GRAPHICS GENERATION ENGINE (APPENDED) ─────────────────────────────────

# 1. Plot NDVI Boxplots (Vegetation Health Across Typologies)
yx_temp <- yx_grn_bio |> st_drop_geometry() |> mutate(City = "Yuexiu (Guangzhou)")
dl_temp <- dl_grn_bio |> st_drop_geometry() |> mutate(City = "Delft (Netherlands)")
combined_patches <- bind_rows(yx_temp, dl_temp)

p_boxplot <- ggplot(combined_patches, aes(x = green_type, y = ndvi_mean, fill = City)) +
  geom_boxplot(outlier.size = 1, alpha = 0.8, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("Yuexiu (Guangzhou)" = "#21918c", "Delft (Netherlands)" = "#440154")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  ) +
  labs(
    title = "Vegetation Quality Index (NDVI) Across Green Space Typologies",
    subtitle = "Comparing canopy vibrancy benchmarks across structural classifications",
    x = "OSM Typology Classification",
    y = "Mean Patch NDVI Profile",
    fill = "Study Region:"
  )

ggsave(file.path(OUT_ROOT, "fig_ndvi_typology_boxplot.png"), p_boxplot, width = 10, height = 6, dpi = 300)

# 2. Plot Biodiversity Scatter Plot (Island Biogeography Evaluation)
p_scatter <- ggplot(combined_patches |> filter(area_ha > 0), aes(x = area_ha, y = gbif_count + 1, color = green_type)) +
  geom_point(aes(shape = City), size = 2.5, alpha = 0.7) +
  geom_smooth(aes(group = City, linetype = City), method = "lm", color = "black", se = FALSE, size = 0.8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_color_viridis_d(option = "D") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Species-Area Relationship Matrix",
    subtitle = "Evaluating biodiversity density logs against patch sizing boundaries",
    x = "Log Patch Area Size (Hectares)",
    y = "Log Total Species Counts (GBIF Record Count + 1 Offset)",
    color = "Green Typology:",
    shape = "Study Region:",
    linetype = "Urban Trend:"
  ) +
  theme(legend.box = "vertical")

ggsave(file.path(OUT_ROOT, "fig_biodiversity_area_scaling.png"), p_scatter, width = 11, height = 6.5, dpi = 300)

# 3. Plot Blue-Green Ratio Spatial Maps
m1 <- ggplot(yx_sub_bg) +
  geom_sf(aes(fill = bg_ratio), color = "white", size = 0.1) +
  scale_fill_viridis_c(
    name = "Blue-Green\nRatio",
    option = "E",
    trans = "log10",
    labels = function(x) sprintf("%.2f", x)
  ) +
  theme_minimal() +
  labs(title = "Yuexiu Matrix Composition")

m2 <- ggplot(dl_wijk_bg) +
  geom_sf(aes(fill = bg_ratio), color = "white", size = 0.1) +
  scale_fill_viridis_c(
    name = "Blue-Green\nRatio",
    option = "E",
    limits = c(0, 150),
    oob = scales::squish
  ) +
  theme_minimal() +
  labs(title = "Delft Matrix Composition")

# Save side-by-side spatial comparative layout
ggsave(file.path(OUT_ROOT, "fig_blue_green_landscape_maps.png"), m1 + m2, width = 14, height = 6, dpi = 300)

message("Script 03 completed successfully. All figures and background models processed.")