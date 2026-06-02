source("R/00_config.R")
library(sf); library(ineq); library(dplyr); library(ggplot2); library(biscale); library(cowplot)

# Load underlying data products
yx_sub_access  <- readRDS(file.path(OUT_ROOT, "yx_sub_access.rds"))
dl_wijk_access <- readRDS(file.path(OUT_ROOT, "dl_wijk_access.rds"))
dl             <- readRDS(file.path(OUT_ROOT, "delft_data.rds"))
d              <- readRDS(file.path(OUT_ROOT, "yuexiu_data.rds"))

# Ensure matching coordinate references
yx_sub_access  <- st_transform(yx_sub_access, 4326)
dl_wijk_access <- st_transform(dl_wijk_access, 4326)

# ── 3A. Gini coefficient (Equity Metric Evaluation) ──────────────────────────
yx_gini <- Gini(yx_sub_access$green_pc_m2, na.rm = TRUE)
dl_gini <- Gini(dl_wijk_access$green_pc_m2, na.rm = TRUE)

cat("\n===========================================================\n")
cat("            ENVIRONMENTAL EQUITY INDEX MATRIX              \n")
cat("===========================================================\n")
cat("Yuexiu Green Space Gini Coefficient (Inequality):", round(yx_gini, 3), "\n")
cat("Delft Green Space Gini Coefficient (Inequality):  ", round(dl_gini, 3), "\n")
cat("===========================================================\n\n")


# ── 3B. Bivariate Choropleth Layer Processing (CORRECTED AXIS MAPPING) ───────
# We map x = pop_count (Demand) and y = green_pc (Supply) to correctly match 
# how the bi_legend grid structures and displays the "GrPink" color matrix.
yx_bi <- bi_class(yx_sub_access, 
                  x = pop_count, 
                  y = green_pc_m2, 
                  style = "quantile", dim = 3)

# Process Delft using identical axis logic and cap outliers at 150m² to preserve color contrast
dl_wijk_capped <- dl_wijk_access |>
  mutate(green_pc_capped = ifelse(green_pc_m2 > 150, 150, green_pc_m2))

dl_bi <- bi_class(dl_wijk_capped, 
                  x = pop_count, 
                  y = green_pc_capped, 
                  style = "quantile", dim = 3)


# ── 3C. Income Correlation (Delft) / VIIRS proxy (Yuexiu) ────────────────────
# 1. Delft Socio-Economic Assessment Matrix
dl_income_source <- if (!is.null(dl$dl_inc)) dl$dl_inc else dl$income
dl_income_source <- st_transform(dl_income_source, st_crs(dl_wijk_access))

dl_joined <- dl_wijk_access |>
  st_join(dl_income_source, join = st_intersects, largest = TRUE)

inc_col_name <- grep("inc", names(dl_joined), value = TRUE, ignore.case = TRUE)[1]

cor_delft <- cor(dl_joined$green_pc_m2, dl_joined[[inc_col_name]],
                 use = "complete.obs", method = "pearson")
cat("Delft green-income correlation Pearson r =", round(cor_delft, 3), "\n")

# 2. Yuexiu Nightlight Economics Proxy Matrix
yx_viirs <- if (!is.null(d$yx_viirs_pts)) d$yx_viirs_pts else d$viirs_pts
yx_sub_layer <- if (!is.null(d$yx_sub)) d$yx_sub else d$subdistricts

yx_viirs     <- st_transform(yx_viirs, st_crs(yx_sub_layer))
yx_sub_layer <- st_transform(yx_sub_layer, st_crs(yx_sub_layer))

# Detect and normalize the light intensity column name inside point file
rad_col_name <- grep("rad|val|int|nit|zen|avg", names(yx_viirs), value = TRUE, ignore.case = TRUE)[1]
if (!is.na(rad_col_name) && rad_col_name != "radiance") {
  yx_viirs <- yx_viirs |> rename(radiance = !!rad_col_name)
} else if (is.na(rad_col_name)) {
  first_data_col <- setdiff(names(yx_viirs), attr(yx_viirs, "sf_column"))[1]
  yx_viirs <- yx_viirs |> rename(radiance = !!first_data_col)
}

# Find the zone identifier column name directly on yx_sub_access FIRST to prevent suffix conflicts
yx_sub_key <- grep("name|jiedao|district|id|sub", names(yx_sub_access), value = TRUE, ignore.case = TRUE)[1]
if (is.na(yx_sub_key)) yx_sub_key <- setdiff(names(yx_sub_access), c("green_pc_m2", "pop_count", "geometry"))[1]

# Spatial Join execution
yx_viirs_joined <- st_join(yx_viirs, yx_sub_layer)

# Find where that matching column name ended up in the joined table (handling .x/.y suffixes automatically)
yx_joined_key <- grep(paste0("^", yx_sub_key), names(yx_viirs_joined), value = TRUE)[1]
if (is.na(yx_joined_key)) yx_joined_key <- yx_sub_key

# Aggregate nightlight radiance using the dynamically matched join key
yx_viirs_agg <- yx_viirs_joined |>
  st_drop_geometry() |>
  group_by(across(all_of(yx_joined_key))) |>   
  summarise(viirs_mean = mean(as.numeric(radiance), na.rm = TRUE)) |>
  rename(join_key = !!yx_joined_key)

# Sync the original main dataset using its localized key name to merge smoothly
yx_sub_access_clean <- yx_sub_access |> rename(join_key = !!yx_sub_key)

yx_sub_joined <- yx_sub_access_clean |>
  left_join(yx_viirs_agg, by = "join_key")

cor_yx <- cor(yx_sub_joined$green_pc_m2, yx_sub_joined$viirs_mean,
              use = "complete.obs", method = "pearson")
cat("Yuexiu green-VIIRS Economic Proxy Pearson r =", round(cor_yx, 3), "\n\n")


# ── Save Completed Workspace Results ──────────────────────────────────────────
saveRDS(list(yx_gini = yx_gini, dl_gini = dl_gini,
             cor_delft = cor_delft, cor_yx = cor_yx,
             yx_sub_joined = yx_sub_joined, dl_joined = dl_joined),
        file.path(OUT_ROOT, "spatial_justice_results.rds"))


# ── 3D. GRAPHICS GENERATION ENGINE (CORRECTED MAPS & LEGEND AXES) ─────────────

# 1. Render Bivariate Spatial Map — Yuexiu
bi_map_yx <- ggplot(yx_bi) +
  geom_sf(aes(fill = bi_class), show.legend = FALSE, color = "white", size = 0.2) +
  bi_scale_fill(pal = "GrPink", dim = 3) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Yuexiu Spatial Justice Grid", 
    subtitle = "Population Demand vs. Green Space per Capita"
  )

# 2. Render Bivariate Spatial Map — Delft
bi_map_dl <- ggplot(dl_bi) +
  geom_sf(aes(fill = bi_class), show.legend = FALSE, color = "white", size = 0.2) +
  bi_scale_fill(pal = "GrPink", dim = 3) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Delft Spatial Justice Grid", 
    subtitle = "Population Demand vs. Green Space per Capita"
  )

# 3. Build the Bivariate 3x3 Coordinate Legend Matrix Square
# X-axis represents horizontal change (Demand), Y-axis represents vertical change (Supply)
bi_legend_cube <- bi_legend(
  pal = "GrPink", dim = 3,
  xlab = "Population Demand (High) ", 
  ylab = "Green Supply (High) ",       
  size = 7
)

# 4. Composite layout components into a unified thematic panel frame
combined_bivariate_layout <- plot_grid(
  bi_map_yx, bi_map_dl, bi_legend_cube,
  labels = c("A", "B", ""),
  rel_widths = c(1, 1, 0.45),
  nrow = 1
)

# Export layout to working directory output paths
ggsave(file.path(OUT_ROOT, "fig_bivariate_equity_maps.png"), 
       combined_bivariate_layout, width = 15, height = 6.5, dpi = 300)

message("Spatial Justice script pipeline processed and graphics saved successfully.")