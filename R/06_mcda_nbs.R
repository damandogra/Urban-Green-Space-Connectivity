source("R/00_config.R")
library(sf); library(dplyr); library(ggplot2); library(viridis)

# Load all prior outputs
access  <- readRDS(file.path(OUT_ROOT, "yx_sub_access.rds"))
bio     <- readRDS(file.path(OUT_ROOT, "yx_grn_bio.rds"))
sj      <- readRDS(file.path(OUT_ROOT, "spatial_justice_results.rds"))
conn    <- readRDS(file.path(OUT_ROOT, "connectivity_results.rds"))

# ── Step 1: Aggregate patch-level indicators to subdistrict level ─────────────
# Join green patch biodiversity and connectivity metrics to Jiēdào
yx_sub <- access  # start from subdistrict access layer

# Mean NDVI per subdistrict (from patch means, area-weighted)
grn_bio <- bio |> st_drop_geometry()
yx_sub_bio <- yx_sub |>
  st_join(bio |> select(ndvi_mean, species_per_ha, betweenness=betweenness),
          join = st_intersects) |>
  st_drop_geometry() |>
  group_by(across(-c(ndvi_mean, species_per_ha, betweenness))) |>
  summarise(ndvi_mean      = mean(ndvi_mean, na.rm = TRUE),
            species_per_ha = mean(species_per_ha, na.rm = TRUE),
            mean_betweenness = mean(betweenness, na.rm = TRUE)) |>
  ungroup()

# Rejoin geometry
yx_mcda <- yx_sub |>
  left_join(st_drop_geometry(yx_sub_bio), by = "geometry")

# ── Step 2: Min-max normalise each indicator to [0, 1] ───────────────────────
normalise <- function(x) (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))

yx_mcda <- yx_mcda |>
  mutate(
    score_access  = normalise(green_pc_m2),
    score_bio     = normalise(ndvi_mean) * 0.5 + normalise(species_per_ha) * 0.5,
    score_conn    = normalise(mean_betweenness),
    # For equity: INVERT green_pc_m2 rank so LOW access = HIGH priority
    score_equity  = 1 - normalise(green_pc_m2)
  )

# ── Step 3: Weighted MCDA score ──────────────────────────────────────────────
yx_mcda <- yx_mcda |>
  mutate(
    mcda_score = MCDA_WEIGHTS["accessibility"]  * score_access  +
                 MCDA_WEIGHTS["biodiversity"]   * score_bio     +
                 MCDA_WEIGHTS["connectivity"]   * score_conn    +
                 MCDA_WEIGHTS["social_equity"]  * score_equity
  )

# ── Step 4: NBS Corridor Prioritisation ──────────────────────────────────────
# Priority = subdistricts with LOW access AND isolated green patches
# Isolation is flagged by patch belonging to a small graph component
conn_data <- conn$yx_grn_conn

# Flag isolated patches (component size < 3 patches)
comp_sizes <- table(conn_data$component)
conn_data <- conn_data |>
  mutate(isolated = comp_sizes[as.character(component)] < 3)

# Spatially join: which Jiēdào contain isolated patches?
isolated_patches <- conn_data |> filter(isolated)
yx_mcda <- yx_mcda |>
  mutate(has_isolated_patch = lengths(st_intersects(geometry, isolated_patches)) > 0)

# Priority tier: high MCDA score + has isolated patches = Tier 1 intervention
yx_mcda <- yx_mcda |>
  mutate(priority_tier = case_when(
    mcda_score > quantile(mcda_score, 0.66) & has_isolated_patch ~ "Tier 1 — Urgent",
    mcda_score > quantile(mcda_score, 0.33)                       ~ "Tier 2 — Moderate",
    TRUE                                                           ~ "Tier 3 — Monitor"
  ))

# ── Step 5: Final map ─────────────────────────────────────────────────────────
p_mcda <- ggplot(yx_mcda) +
  geom_sf(aes(fill = mcda_score)) +
  scale_fill_viridis_c(option = "inferno", name = "MCDA Score\n(higher = more urgent)") +
  theme_minimal() +
  labs(title = "MCDA Priority Score — Yuexiu Jiēdào",
       subtitle = paste0("Weights: Accessibility ",
                         MCDA_WEIGHTS["accessibility"]*100, "% | Biodiversity ",
                         MCDA_WEIGHTS["biodiversity"]*100, "% | Connectivity ",
                         MCDA_WEIGHTS["connectivity"]*100, "% | Equity ",
                         MCDA_WEIGHTS["social_equity"]*100, "%"))

p_priority <- ggplot(yx_mcda) +
  geom_sf(aes(fill = priority_tier)) +
  scale_fill_manual(values = c("Tier 1 — Urgent"   = "#d73027",
                               "Tier 2 — Moderate" = "#fee090",
                               "Tier 3 — Monitor"  = "#91cf60"),
                    name = "NBS Priority") +
  theme_minimal() +
  labs(title = "NBS Corridor Prioritisation — Yuexiu")

ggsave(file.path(OUT_ROOT, "fig_mcda_yuexiu.png"),    p_mcda,     width=10, height=8)
ggsave(file.path(OUT_ROOT, "fig_priority_yuexiu.png"), p_priority, width=10, height=8)

saveRDS(yx_mcda, file.path(OUT_ROOT, "mcda_final.rds"))
message("MCDA complete.")