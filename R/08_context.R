# Context maps for presentation and report introduction
# Setup ------------------------------------------------------
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
# 0A Delft and Guangzhou context

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
  theme_map_clean() +
  labs( title = "Guangzhou", subtitle = "Source: ??")

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
  theme_map_clean() +
  labs( title = "Delft", subtitle = "Source: CBS wijk data")

  annotation_scale(
    location = "bl",
    style = "ticks",
    unit_category = "metric",
    width_hint = 0.3,
    text_cex = 0.7,
    line_width = 0.6
  )

fig_context_gz_dl <- p_context_gz + p_context_dl +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(title = "Context", theme = theme(plot.title = element_text(face = "bold", size = 14)))

fig_context_gz_dl

ggsave( file.path(OUT_ROOT, "fig_context_gz_dl.png"), fig_context_gz_dl, width = 14, height = 6, dpi = 300)

# 0B spatial unit (su) Delft and Yuexiu context
# transform both maps to metre-based CRS
yx_plot <- st_transform(yx_sub_access, CRS_YX)
dl_plot <- st_transform(dl_wijk_access, 28992)  # RD New

# calculate real maps widths
yx_bb <- st_bbox(yx_plot)
dl_bb <- st_bbox(dl_plot)

yx_w <- yx_bb$xmax - yx_bb$xmin
dl_w <- dl_bb$xmax - dl_bb$xmin

p_su_yx <- ggplot() +
  geom_sf(data = yx_bnd, fill = COLORS$orange_light, color = NA) +
  geom_sf(data = yx_sub_context,  fill = NA, color = "white", linewidth = 0.12) +
  geom_sf(data = yx_bnd, fill = NA, color = COLORS$red_light, linewidth = 0.55) +
  annotation_scale( location = "bl", style = "ticks", width_hint = 0.25, text_cex = 0.7, line_width = 0.4 ) +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs( title = "Yuexiu", subtitle = "Source: WorldPop + administrative polygons")


p_su_dl <- ggplot() +
  geom_sf(data = dl_bnd, fill = COLORS$blue_light, color = NA) +
  geom_sf(data = dl_wijk_context, fill = NA, color = "white", linewidth = 0.12 ) +
  geom_sf(data = dl_bnd, fill = NA, color = COLORS$red_light, linewidth = 0.55 ) +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs(title = "Delft", subtitle = "Source: CBS wijk data")


fig_context_yx_dl <- p_su_yx +  p_su_dl +
  plot_layout(widths = c(yx_w, dl_w), guides = "collect")+
  plot_annotation( title = "Context", theme = theme( plot.title = element_text(face = "bold", size = 14))) &
  theme(legend.position = "right")

fig_context_yx_dl

ggsave(file.path(OUT_ROOT, "fig_context_yx_dl.png"), fig_context_yx_dl, width = 14, height = 6, dpi = 300)


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

# transform both maps to metre-based CRS
yx_plot <- st_transform(yx_sub_access, CRS_YX)
dl_plot <- st_transform(dl_wijk_access, 28992)  # RD New

# calculate real maps widths
yx_bb <- st_bbox(yx_plot)
dl_bb <- st_bbox(dl_plot)

yx_w <- yx_bb$xmax - yx_bb$xmin
dl_w <- dl_bb$xmax - dl_bb$xmin

p_ugs_yx <- ggplot() +
  geom_sf(data = yx_bnd, aes(fill = "Urban fabric"), color = NA) +
  geom_sf(data = yx_water_context, aes(fill = "Water"), color = NA, alpha = 0.6) +
  geom_sf(data = yx_roads, aes(color = "Roads"), linewidth = 0.15) +
  geom_sf(data = yx_green, aes(fill = "Urban green space"), color = NA, alpha = 0.9) +
  geom_sf(data = yx_sub_context, aes(color = "Administrative boundaries"), fill = NA, linewidth = 0.16) +
  geom_sf(data = yx_bnd, aes(color = "City boundary"), fill = NA, linewidth = 0.55) +
  annotation_scale( location = "bl", style = "ticks", width_hint = 0.25, text_cex = 0.7, line_width = 0.4 ) +
  coord_sf(expand = FALSE, datum = NA) +
  scale_fill_manual(name = NULL, values = c("Urban fabric" = COLORS$beige, "Water" = COLORS$blue_light, "Urban green space" = COLORS$green_dark)) +
  scale_color_manual( name = NULL, values = c( "Roads" = COLORS$grey85, "Administrative boundaries" = COLORS$red, "City boundary" = COLORS$red_light)) +
  theme_map_clean() +
  labs(title = "Yuexiu - Jiēdào", subtitle = "Source: WorldPop + administrative polygons")

p_ugs_dl <- ggplot() +
  geom_sf(data = dl_bnd, aes(fill = "Urban fabric"), color = NA) +
  geom_sf(data = dl_water_context, aes(fill = "Water"), color = NA, alpha = 0.6) +
  geom_sf(data = dl_roads, aes(color = "Roads"), linewidth = 0.10) +
  geom_sf(data = dl_green, aes(fill = "Urban green space"), color = NA, alpha = 0.95) +
  geom_sf(data = dl_wijk_context, aes(color = "Administrative boundaries"), fill = NA, linewidth = 0.16) +
  geom_sf(data = dl_bnd, aes(color = "City boundary"), fill = NA, linewidth = 0.55) +
  coord_sf(expand = FALSE, datum = NA) +
  scale_fill_manual(name = NULL, values = c( "Urban fabric" = COLORS$beige, "Water" = COLORS$blue_light, "Urban green space" = COLORS$green_dark), guide = "none") +
  scale_color_manual( name = NULL, values = c( "Roads" = COLORS$grey85, "Administrative boundaries" = COLORS$red, "City boundary" = COLORS$red_light), guide = "none") +
  theme_map_clean() +
  labs(title = "Delft - Wijken", subtitle = "Source: CBS wijk data")

fig_ugs_yx_dl <- p_ugs_yx + p_ugs_dl +
  plot_layout(
    widths = c(yx_w, dl_w),
    guides = "collect") +
  plot_annotation(title = "Context", theme = theme( plot.title = element_text(face = "bold", size = 14))) &
  theme(legend.position = "right")

fig_ugs_yx_dl

ggsave(file.path(OUT_ROOT, "fig_ugs_yx_dl.png"), fig_ugs_yx_dl, width = 14, height = 6, dpi = 300)
