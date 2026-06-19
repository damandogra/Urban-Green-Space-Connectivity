# not debugged and ran

source("R/00_config.R")
library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

sf_use_s2(FALSE)

# ── Load all prior outputs ────────────────────────────────────────────────────
d  <- readRDS(file.path(OUT_ROOT, "yuexiu_data.rds"))
dl <- readRDS(file.path(OUT_ROOT, "delft_data.rds"))

yx_sub_access  <- readRDS(file.path(OUT_ROOT, "yx_sub_access.rds"))
dl_wijk_access <- readRDS(file.path(OUT_ROOT, "dl_wijk_access.rds"))

yx_sub_ndvi    <- readRDS(file.path(OUT_ROOT, "yx_sub_ndvi.rds"))
dl_wijk_ndvi   <- readRDS(file.path(OUT_ROOT, "dl_wijk_ndvi.rds"))

yx_bivar       <- readRDS(file.path(OUT_ROOT, "yx_bivar.rds"))
dl_bivar       <- readRDS(file.path(OUT_ROOT, "dl_bivar.rds"))

yx_graph       <- readRDS(file.path(OUT_ROOT, "yx_graph.rds"))
dl_graph       <- readRDS(file.path(OUT_ROOT, "dl_graph.rds"))

# MCDA weights from config
w <- MCDA_WEIGHTS   # accessibility=0.30, biodiversity=0.25, connectivity=0.25, equity=0.20

# ── Helper: min-max normalise to [0, 1] ───────────────────────────────────────
norm_minmax <- function(x, invert = FALSE) {
  mn <- min(x, na.rm = TRUE)
  mx <- max(x, na.rm = TRUE)
  if (mx == mn) return(rep(0.5, length(x)))
  s <- (x - mn) / (mx - mn)
  if (invert) 1 - s else s
}

# ── 5A. Build MCDA score table per admin unit ─────────────────────────────────
# For each city separately (different admin units), assemble four sub-scores
# then compute the weighted MCDA composite.

build_mcda <- function(access_sf, ndvi_sf, bivar_sf, graph_obj,
                        local_crs, city_label) {

  # --- Accessibility sub-score ---
  # Higher green_pc_m2 → better; lower nearest_green_m → better
  acc <- access_sf |>
    st_drop_geometry() |>
    mutate(
      s_green_pc   = norm_minmax(green_pc_m2,   invert = FALSE),
      s_nearest    = norm_minmax(nearest_green_m, invert = TRUE),   # closer = better
      score_access = (s_green_pc + s_nearest) / 2
    ) |>
    select(matches("^(name|naam|BU_|WK_|GEO_|jiedao|subdistrict)", ignore.case = TRUE),
           score_access)

  # Use geometry from access_sf as the spatial backbone
  backbone <- access_sf |>
    mutate(row_id = row_number())

  # --- Biodiversity sub-score ---
  # Use mean NDVI per admin unit as proxy (already computed in 03)
  bio <- ndvi_sf |>
    st_drop_geometry() |>
    mutate(score_bio = norm_minmax(ndvi_mean, invert = FALSE)) |>
    select(score_bio)

  # --- Equity sub-score ---
  # Pop density in a LOW-green tertile → HIGH equity need → invert so high score = needs intervention
  # For MCDA we treat equity score as "intervention urgency from equity lens":
  #   High pop + Low green → high score
  equity <- bivar_sf |>
    st_drop_geometry() |>
    mutate(
      equity_score_raw = case_when(
        bivar_class == "Low-High"  ~ 1.00,
        bivar_class == "Mid-High"  ~ 0.75,
        bivar_class == "Low-Mid"   ~ 0.65,
        bivar_class == "Mid-Mid"   ~ 0.40,
        bivar_class == "Low-Low"   ~ 0.30,
        bivar_class == "High-Mid"  ~ 0.20,
        bivar_class == "High-High" ~ 0.15,
        bivar_class == "Mid-Low"   ~ 0.10,
        bivar_class == "High-Low"  ~ 0.10,
        TRUE ~ 0.10
      ),
      score_equity = equity_score_raw
    ) |>
    select(score_equity)

  # --- Connectivity sub-score ---
  # For each admin unit: mean betweenness of patches whose centroid falls within it
  admin_m  <- st_transform(access_sf, local_crs)
  nodes_m  <- st_transform(graph_obj$nodes, local_crs)
  node_cents <- st_centroid(nodes_m)

  joined_conn <- st_join(admin_m, node_cents["degree"],
                         join = st_contains) |>
    st_drop_geometry() |>
    group_by(row_number()) |>
    summarise(mean_degree = mean(degree, na.rm = TRUE), .groups = "drop")

  conn_scores <- joined_conn |>
    mutate(
      mean_degree = replace_na(mean_degree, 0),
      score_conn = 1 - norm_minmax(mean_degree)
    ) |>
    select(score_conn)
  # Patches with low betweenness in an admin unit → connectivity gap → invert
  conn_scores <- joined_conn |>
    mutate(
      mean_degree = replace_na(mean_degree, 0),
      score_conn = 1 - norm_minmax(mean_degree)   # low connectivity = high need
    ) |>
    select(score_conn)
  message(city_label, " connectivity mean betweenness:")
  print(summary(joined_conn$mean_betweenness))

  message(city_label, " connectivity score:")
  print(summary(conn_scores$score_conn))

  # --- Assemble & compute weighted composite ---
  n <- nrow(access_sf)
  # Ensure all sub-score vectors are length n
  pad <- function(x, len) { if (length(x) < len) c(x, rep(NA_real_, len - length(x))) else x[seq_len(len)] }

  score_access  <- pad(acc$score_access,    n)
  score_bio     <- pad(bio$score_bio,       n)
  score_equity  <- pad(replace_na(equity$score_equity, 0.10), n)
  score_equity  <- replace_na(score_equity, 0.10)
  score_conn    <- pad(conn_scores$score_conn, n)

  access_sf |>
    mutate(
      score_access   = score_access,
      score_bio      = score_bio,
      score_equity   = score_equity,
      score_conn     = score_conn,
      mcda_composite = w["accessibility"] * score_access  +
                       w["biodiversity"]  * score_bio     +
                       w["connectivity"]  * score_conn    +
                       w["equity"]        * score_equity,
      city           = city_label,
      # Priority tier
      priority_tier  = case_when(
        mcda_composite >= quantile(mcda_composite, 0.67, na.rm = TRUE) ~ "High priority",
        mcda_composite >= quantile(mcda_composite, 0.33, na.rm = TRUE) ~ "Medium priority",
        TRUE                                                            ~ "Low priority"
      )
    )
}

message("Building MCDA scores — Yuexiu...")
yx_mcda <- build_mcda(yx_sub_access, yx_sub_ndvi, yx_bivar, yx_graph,
                       CRS_YX, "Yuexiu")

message("Building MCDA scores — Delft...")
dl_mcda <- build_mcda(dl_wijk_access, dl_wijk_ndvi, dl_bivar, dl_graph,
                       CRS_DELFT, "Delft")

# ── 5B. NbS corridor prioritisation ───────────────────────────────────────────
# Identify potential NbS intervention zones:
#   Criteria: HIGH priority MCDA tier + admin units with low connectivity patches
#             → candidate corridors between isolated patches

identify_nbs_corridors <- function(mcda_sf, graph_obj, local_crs, thresh_m = DISPERSAL_THRESH_M) {
  high_priority <- mcda_sf |>
    filter(priority_tier == "High priority") |>
    st_transform(local_crs)

  # Isolated patches: enn_m > dispersal threshold (not currently connected)
  isolated_patches <- graph_obj$nodes |>
    st_transform(local_crs) |>
    filter(enn_m > thresh_m)

  if (nrow(isolated_patches) == 0) {
    message("  No isolated patches at threshold ", thresh_m, " m")
    return(list(high_priority = high_priority, isolated_patches = isolated_patches,
                corridor_lines = NULL))
  }

  # Potential corridors: lines from each isolated patch to its nearest neighbour
  # regardless of threshold — showing the gap that needs bridging
  cents <- st_centroid(isolated_patches)
  all_cents <- st_centroid(st_transform(graph_obj$nodes, local_crs))

  nearest_idx  <- st_nearest_feature(cents, all_cents)
  corridor_lines <- lapply(seq_len(nrow(cents)), function(i) {
    a <- st_coordinates(cents[i, ])
    b <- st_coordinates(all_cents[nearest_idx[i], ])
    st_linestring(rbind(a, b))
  })
  corridors_sf <- st_sf(
    from_patch   = isolated_patches$area_ha,
    gap_m        = isolated_patches$enn_m,
    geometry     = st_sfc(corridor_lines, crs = local_crs)
  )

  list(high_priority   = high_priority,
       isolated_patches = isolated_patches,
       corridor_lines  = corridors_sf)
}

message("Identifying NbS corridors — Yuexiu...")
yx_nbs <- identify_nbs_corridors(yx_mcda, yx_graph, CRS_YX)

message("Identifying NbS corridors — Delft...")
dl_nbs <- identify_nbs_corridors(dl_mcda, dl_graph, CRS_DELFT)

# ── Save outputs ──────────────────────────────────────────────────────────────
saveRDS(yx_mcda, file.path(OUT_ROOT, "yx_mcda.rds"))
saveRDS(dl_mcda, file.path(OUT_ROOT, "dl_mcda.rds"))
saveRDS(yx_nbs,  file.path(OUT_ROOT, "yx_nbs.rds"))
saveRDS(dl_nbs,  file.path(OUT_ROOT, "dl_nbs.rds"))

# MCDA weight table for reporting
weight_df <- data.frame(
  Criterion    = c("Accessibility", "Biodiversity / NDVI", "Connectivity", "Equity"),
  Weight       = unname(w),
  Sub_scores   = c("Green pc + nearest distance",
                   "Mean NDVI per unit",
                   "Inverse mean betweenness",
                   "Pop × green density matrix")
)
write.csv(weight_df, file.path(OUT_ROOT, "mcda_weights.csv"), row.names = FALSE)
message("MCDA weights saved to mcda_weights.csv")

# ── Figures ───────────────────────────────────────────────────────────────────

# Figure 5A: MCDA composite score maps
p_mcda_yx <- ggplot(yx_mcda) +
  geom_sf(aes(fill = mcda_composite), colour = "white", linewidth = 0.3) +
  scale_fill_viridis_c(name = "MCDA score\n(0–1)", option = "A",
                        direction = -1, na.value = "grey80") +
  theme_minimal(base_size = 11) +
  labs(title = "MCDA Composite Score — Yuexiu Subdistricts",
       subtitle = sprintf("Weights: Acc %.0f%% | Bio %.0f%% | Conn %.0f%% | Eq %.0f%%",
                          w["accessibility"]*100, w["biodiversity"]*100,
                          w["connectivity"]*100, w["equity"]*100))

p_mcda_dl <- ggplot(dl_mcda) +
  geom_sf(aes(fill = mcda_composite), colour = "white", linewidth = 0.3) +
  scale_fill_viridis_c(name = "MCDA score\n(0–1)", option = "A",
                        direction = -1, na.value = "grey80") +
  theme_minimal(base_size = 11) +
  labs(title = "MCDA Composite Score — Delft Wijken",
       subtitle = sprintf("Weights: Acc %.0f%% | Bio %.0f%% | Conn %.0f%% | Eq %.0f%%",
                          w["accessibility"]*100, w["biodiversity"]*100,
                          w["connectivity"]*100, w["equity"]*100))

ggsave(file.path(OUT_ROOT, "fig_mcda_maps.png"),
       p_mcda_yx + p_mcda_dl, width = 14, height = 7, dpi = 300)

# Figure 5B: Priority tier maps
tier_colours <- c("High priority"   = "#d73027",
                   "Medium priority" = "#fee090",
                   "Low priority"    = "#EBF6F9")

p_tier_yx <- ggplot(yx_mcda) +
  geom_sf(aes(fill = priority_tier), colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = tier_colours, name = "NbS Priority") +
  theme_minimal(base_size = 11) +
  labs(title = "NbS Intervention Priority — Yuexiu",
       subtitle = "Tertile classification of MCDA composite score")

p_tier_dl <- ggplot(dl_mcda) +
  geom_sf(aes(fill = priority_tier), colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = tier_colours, name = "NbS Priority") +
  theme_minimal(base_size = 11) +
  labs(title = "NbS Intervention Priority — Delft",
       subtitle = "Tertile classification of MCDA composite score")

ggsave(file.path(OUT_ROOT, "fig_priority_tiers.png"),
       p_tier_yx + p_tier_dl, width = 14, height = 7, dpi = 300)

# Figure 5C: NbS corridor maps (isolated patches + proposed links)
plot_nbs_map <- function(nbs_obj, bnd_sf, local_crs, city_label) {
  bnd <- st_transform(bnd_sf, local_crs)

  p <- ggplot() +
    geom_sf(data = bnd, fill = "grey95", colour = "grey60", linewidth = 0.5) +
    geom_sf(data = st_transform(nbs_obj$high_priority, local_crs),
            fill = "#d73027", alpha = 0.25, colour = NA)

  if (!is.null(nbs_obj$corridor_lines) && nrow(nbs_obj$corridor_lines) > 0) {
    p <- p +
      geom_sf(data = st_transform(nbs_obj$corridor_lines, local_crs),
              colour = "#e6550d", linewidth = 0.7, linetype = "dashed", alpha = 0.8)
  }

  p +
    geom_sf(data = nbs_obj$isolated_patches,
            fill = "#74c476", colour = "white", alpha = 0.85, linewidth = 0.2) +
    theme_minimal(base_size = 11) +
    labs(title = paste("NbS Corridor Prioritisation —", city_label),
         subtitle = paste0("Red fill = high-priority admin zones | ",
                           "Green = isolated patches | Dashed = proposed corridors"))
}

p_nbs_yx <- plot_nbs_map(yx_nbs, d$yx_bnd, CRS_YX,    "Yuexiu")
p_nbs_dl <- plot_nbs_map(dl_nbs, dl$dl_bnd, CRS_DELFT, "Delft")

ggsave(file.path(OUT_ROOT, "fig_nbs_corridors.png"),
       p_nbs_yx + p_nbs_dl, width = 14, height = 7, dpi = 300)

# Figure 5D: Sub-score radar / spider chart (city-level means)
score_means <- bind_rows(
  yx_mcda |> st_drop_geometry() |>
    summarise(across(c(score_access, score_bio, score_equity, score_conn), mean, na.rm = TRUE)) |>
    mutate(city = "Yuexiu"),
  dl_mcda |> st_drop_geometry() |>
    summarise(across(c(score_access, score_bio, score_equity, score_conn), mean, na.rm = TRUE)) |>
    mutate(city = "Delft")
) |>
  pivot_longer(-city, names_to = "criterion", values_to = "score") |>
  mutate(criterion = recode(criterion,
    "score_access"  = "Accessibility",
    "score_bio"     = "Biodiversity",
    "score_equity"  = "Equity",
    "score_conn"    = "Connectivity"
  ))

p_radar_bar <- ggplot(score_means, aes(x = criterion, y = score, fill = city)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dotted", colour = "grey70") +
  scale_fill_manual(values = c("Yuexiu" = "#21918c", "Delft" = "#440154")) +
  scale_y_continuous(limits = c(0, 1), labels = label_percent()) +
  coord_polar() +
  theme_minimal(base_size = 12) +
  labs(title = "MCDA Sub-score Profile — City Comparison",
       subtitle = "Mean normalised score (0=worst, 1=best) across admin units",
       x = NULL, y = NULL, fill = "City")

ggsave(file.path(OUT_ROOT, "fig_mcda_radar.png"),
       p_radar_bar, width = 7, height = 7, dpi = 300)

# Figure 5E: MCDA score distribution comparison
score_dist <- bind_rows(
  yx_mcda |> st_drop_geometry() |> transmute(city = "Yuexiu", mcda_composite),
  dl_mcda |> st_drop_geometry() |> transmute(city = "Delft",  mcda_composite)
)

p_dist <- ggplot(score_dist, aes(x = mcda_composite, fill = city)) +
  geom_density(alpha = 0.6) +
  geom_rug(aes(colour = city), alpha = 0.6) +
  scale_fill_manual(values  = c("Yuexiu" = "#21918c", "Delft" = "#440154")) +
  scale_colour_manual(values = c("Yuexiu" = "#21918c", "Delft" = "#440154"),
                      guide = "none") +
  theme_minimal(base_size = 12) +
  labs(title = "Distribution of MCDA Composite Scores",
       subtitle = "Density plot — all admin units",
       x = "MCDA composite score (0–1)", y = "Density", fill = "City")

ggsave(file.path(OUT_ROOT, "fig_mcda_distribution.png"),
       p_dist, width = 8, height = 5, dpi = 300)

message("Script 06 complete — MCDA & NbS figures saved.")
message("")
message("=== All scripts complete. Output files in: ", OUT_ROOT, " ===")
