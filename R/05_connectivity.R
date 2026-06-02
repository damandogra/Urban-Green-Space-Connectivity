source("R/00_config.R")
library(sf); library(igraph); library(dplyr)
library(landscapemetrics); library(terra); library(ggplot2)

d  <- readRDS(file.path(OUT_ROOT, "yuexiu_data.rds"))
dl <- readRDS(file.path(OUT_ROOT, "delft_data.rds"))

# ── 4A. Landscape fragmentation metrics (landscapemetrics) ───────────────────
# Requires binary raster: green = 1, non-green = 0
# Use the WorldCover raster, reclassify tree/grass/shrub as green
make_green_raster <- function(cover_rast) {
  m <- matrix(c(10, 1,   # tree cover
                20, 1,   # shrubland
                30, 1,   # grassland
                40, 0,   # cropland
                50, 0,   # built-up
                60, 0,   # bare
                80, 0,   # water
                95, 1),  # mangrove
              ncol = 2, byrow = TRUE)
  classify(cover_rast, m, others = 0)
}

yx_green_rast <- make_green_raster(d$yx_cover)
dl_green_rast <- make_green_raster(dl$dl_cover)

# Core metrics: number of patches (NP), mean patch size (MPS),
#               euclidean nearest-neighbour distance (ENN)
yx_frag <- calculate_lsm(yx_green_rast, what = c("lsm_c_np", "lsm_c_area_mn", "lsm_c_enn_mn"))
dl_frag <- calculate_lsm(dl_green_rast, what = c("lsm_c_np", "lsm_c_area_mn", "lsm_c_enn_mn"))

# ── 4B. Functional connectivity graph ────────────────────────────────────────
# Nodes = green patch centroids; edges where centroids are within dispersal threshold
build_connectivity_graph <- function(green_sf, threshold_m) {
  cents <- st_centroid(green_sf)
  dist_mat <- st_distance(cents)          # n × n distance matrix

  edges <- which(dist_mat > units::as_units(0, "m") &
                   dist_mat <= units::as_units(threshold_m, "m"),
                 arr.ind = TRUE)

  g <- graph_from_edgelist(edges, directed = FALSE)
  V(g)$patch_id   <- seq_len(nrow(green_sf))
  V(g)$area_ha    <- green_sf$area_ha
  E(g)$dist_m     <- as.numeric(dist_mat[edges])

  # Betweenness centrality — identifies stepping-stone patches
  V(g)$betweenness <- betweenness(g, normalized = TRUE)
  g
}

yx_graph <- build_connectivity_graph(d$yx_grn,  DISPERSAL_THRESH_M)
dl_graph <- build_connectivity_graph(dl$dl_grn, DISPERSAL_THRESH_M)

# Key metrics
cat("Yuexiu — components:", components(yx_graph)$no,
    "| avg betweenness:", mean(V(yx_graph)$betweenness), "\n")
cat("Delft  — components:", components(dl_graph)$no,
    "| avg betweenness:", mean(V(dl_graph)$betweenness), "\n")

# Add betweenness back to sf for mapping
d$yx_grn <- d$yx_grn |>
  mutate(betweenness = V(yx_graph)$betweenness,
         component   = components(yx_graph)$membership)

saveRDS(list(yx_graph=yx_graph, dl_graph=dl_graph,
             yx_frag=yx_frag, dl_frag=dl_frag,
             yx_grn_conn=d$yx_grn),
        file.path(OUT_ROOT, "connectivity_results.rds"))