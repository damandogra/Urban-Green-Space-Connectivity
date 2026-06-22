# Context maps for presentation and report introduction
# Setup -----------------------------------------------------
source("R/00_config.R")

library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggspatial)

# Study area units
dl_wijk_context <- dl_wijk
yx_sub_context  <- yx_sub

# Guangzhou municipality boundary
yx_bnd <- st_read("data/Yuexiu/vector/Yuexiu_boundary_proj.gpkg")

# Make equal-size map extents
make_square_extent <- function(x, size_m = 9000) {

  bb <- st_bbox(x)

  cx <- (bb$xmin + bb$xmax) / 2
  cy <- (bb$ymin + bb$ymax) / 2

  half <- size_m / 2

  c(
    xmin = cx - half,
    xmax = cx + half,
    ymin = cy - half,
    ymax = cy + half
  )
}

dl_equal <- make_square_extent(dl_bnd, size_m = 9000)

yx_equal <- make_square_extent(yx_bnd, size_m = 9000)

# Shared map theme -----------------------------------------

theme_context <- theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 18, face = "bold", color = "grey20")
  )

## figures
# Delft and Guangzhou context

p_context_gz <- ggplot() +
  geom_sf(
    data = gz_bnd,
    fill = COLORS$orange_light,
    color = NA
  ) +
  geom_sf(
    data = yx_bnd,
    fill = COLORS$orange,
    color = "white",
    linewidth = 0.12
  ) +
  geom_sf(
    data = gz_bnd,
    fill = NA,
    color = COLORS$red_light,
    linewidth = 0.55
  ) +
  theme_context +
  labs(title = "Guangzhou") +

  annotation_scale(
    location = "bl",
    style = "ticks",
    line_width = 0.6,
    text_cex = 0.7,
    width_hint = 0.2
  )

p_context_dl <- ggplot() +
  geom_sf(
    data = dl_bnd,
    fill = COLORS$blue_light,
    color = NA
  ) +
  geom_sf(
    data = dl_wijk_context,
    fill = NA,
    color = "white",
    linewidth = 0.12
  ) +
  geom_sf(
    data = dl_bnd,
    fill = NA,
    color = COLORS$red_light,
    linewidth = 0.55
  ) +
  theme_context +
  labs(title = "Delft") +

  annotation_scale(
    location = "bl",
    style = "ticks",
    unit_category = "metric",
    width_hint = 0.3,
    text_cex = 0.7,
    line_width = 0.6
  )

p_context_dl + p_context_gz

# spatial unit (su) Delft and Yuexiu context

p_su_dl <- ggplot() +
  geom_sf(
    data = dl_bnd,
    fill = COLORS$blue_light,
    color = NA
  ) +
  geom_sf(
    data = dl_wijk_context,
    fill = NA,
    color = "white",
    linewidth = 0.12
  ) +
  geom_sf(
    data = dl_bnd,
    fill = NA,
    color = COLORS$red_light,
    linewidth = 0.55
  ) +
  theme_context +
  labs(title = "Delft") +
  annotation_scale(
    location = "bl",
    style = "ticks",
    unit_category = "metric",
    width_hint = 0.30,
    text_cex = 0.7,
    line_width = 0.6
  ) +
  coord_sf(
    xlim = c(dl_equal["xmin"], dl_equal["xmax"]),
    ylim = c(dl_equal["ymin"], dl_equal["ymax"]),
    expand = FALSE
  )

p_su_yx <- ggplot() +
  geom_sf(
    data = yx_bnd,
    fill = COLORS$orange_light,
    color = NA
  ) +
  geom_sf(
    data = yx_sub_context,
    fill = NA,
    color = "white",
    linewidth = 0.12
  ) +
  geom_sf(
    data = yx_bnd,
    fill = NA,
    color = COLORS$red_light,
    linewidth = 0.55
  ) +
  theme_context +
  labs(title = "Yuexiu") +

  coord_sf(
    xlim = c(dl_equal["xmin"], dl_equal["xmax"]),
    ylim = c(dl_equal["ymin"], dl_equal["ymax"]),
    expand = FALSE
  )

p_su_dl + p_su_yx

# Urban green space ------------------------------------

# Prepare roads and green-space layers
dl_roads <- dl_net$roads
yx_roads <- yx_net$roads

dl_green <- dl_net$green
yx_green <- yx_net$green

# Prepare water layers
dl_water_context <- st_intersection(dl_water, dl_bnd)

yx_water_context <- st_transform(yx_water_poly, st_crs(yx_bnd))
yx_water_context <- st_intersection(yx_water_context, yx_bnd)

# Clip roads and green space to study area boundaries
dl_roads <- st_intersection(dl_roads, dl_bnd)
dl_green <- st_intersection(dl_green, dl_bnd)

yx_roads <- st_transform(yx_roads, st_crs(yx_bnd))
yx_green <- st_transform(yx_green, st_crs(yx_bnd))

yx_roads <- st_intersection(yx_roads, yx_bnd)
yx_green <- st_intersection(yx_green, yx_bnd)

p_ugs_dl <- ggplot() +
  geom_sf(data = dl_bnd,
          fill = COLORS$beige,
          color = NA) +
  geom_sf(
    data = dl_water_context,
    fill = COLORS$blue_light,
    color = NA,
    alpha = 0.6
  ) +
  geom_sf(data = dl_roads,
          color = COLORS$grey85,
          linewidth = 0.10) +
  geom_sf(data = dl_green,
          fill = COLORS$green_dark,
          color = NA,
          alpha = 0.95) +
  geom_sf(data = dl_wijk_context,
          fill = NA,
          color = COLORS$red,
          linewidth = 0.16) +
  geom_sf(data = dl_bnd,
          fill = NA,
          color = COLORS$red_light,
          linewidth = 0.55) +
  theme_context +
  labs(title = "Delft") +
  annotation_scale(
    location = "bl",
    style = "ticks",
    unit_category = "metric",
    width_hint = 0.22,
    text_cex = 0.7,
    line_width = 0.6
  ) +
  coord_sf(
    xlim = c(dl_equal["xmin"], dl_equal["xmax"]),
    ylim = c(dl_equal["ymin"], dl_equal["ymax"]),
    expand = FALSE
  )

p_ugs_yx <- ggplot() +
  geom_sf(data = yx_bnd,
          fill = COLORS$beige,
          color = NA) +
  geom_sf(
    data = yx_water_context,
    fill = COLORS$blue_light,
    color = NA,
    alpha = 0.6
  ) +
  geom_sf(data = yx_roads,
          color = COLORS$grey85,
          linewidth = 0.15) +
  geom_sf(data = yx_green,
          fill = COLORS$green_dark,
          color = NA,
          alpha = 0.9) +
  geom_sf(data = yx_sub_context,
          fill = NA,
          color = COLORS$red,
          linewidth = 0.16) +
  geom_sf(data = yx_bnd,
          fill = NA,
          color = COLORS$red_light,
          linewidth = 0.55) +
  theme_context +
  labs(title = "Yuexiu") +
  coord_sf(
    xlim = c(yx_equal["xmin"], yx_equal["xmax"]),
    ylim = c(yx_equal["ymin"], yx_equal["ymax"]),
    expand = FALSE
  )

p_ugs_dl + p_ugs_yx

ggsave(
  file.path(report_files, "context_urban_green_space.png"),
  p_ugs_dl + p_ugs_yx,
  width = 14,
  height = 7,
  dpi = 300
)
