# not debugged and ran

source("R/00_config.R")
library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

sf_use_s2(FALSE)

# ── Load pre-computed accessibility outputs + raw data ────────────────────────
d  <- readRDS(file.path(OUT_ROOT, "yuexiu_data.rds"))
dl <- readRDS(file.path(OUT_ROOT, "delft_data.rds"))

yx_sub_access  <- readRDS(file.path(OUT_ROOT, "yx_sub_access.rds"))
dl_wijk_access <- readRDS(file.path(OUT_ROOT, "dl_wijk_access.rds"))

d$yx_pop  <- rast(d$yx_pop)
d$yx_viirs <- rast(d$yx_viirs)
dl$dl_pop  <- rast(dl$dl_pop)

# ── 3A. Gini coefficient of green space per capita ────────────────────────────
# Lorenz curve: rank neighbourhoods by green_pc_m2; Gini = 1 - 2*AUC(Lorenz)

gini_index <- function(x) {
  x <- x[!is.na(x) & x >= 0]
  if (length(x) == 0) return(NA_real_)
  x <- sort(x)
  n <- length(x)
  # Standard Gini formula
  2 * sum(seq_len(n) * x) / (n * sum(x)) - (n + 1) / n
}

lorenz_df <- function(x, city_label) {
  x <- sort(x[!is.na(x) & x >= 0])
  cumx <- cumsum(x) / sum(x)
  data.frame(
    city    = city_label,
    cum_pop = seq_along(x) / length(x),
    cum_grn = cumx
  )
}

yx_gini <- gini_index(yx_sub_access$green_pc_m2)
dl_gini <- gini_index(dl_wijk_access$green_pc_m2)

message(sprintf("Gini — Yuexiu: %.3f  |  Delft: %.3f", yx_gini, dl_gini))

lorenz_data <- bind_rows(
  lorenz_df(yx_sub_access$green_pc_m2, "Yuexiu"),
  lorenz_df(dl_wijk_access$green_pc_m2, "Delft")
)
gini_labels <- data.frame(
  city  = c("Yuexiu", "Delft"),
  gini  = c(yx_gini, dl_gini),
  label = sprintf("Gini = %.3f", c(yx_gini, dl_gini))
)

# ── 3B. Bivariate choropleth: green density vs. population density ────────────
# Classify each unit into a 3×3 matrix (low/mid/high green × low/mid/high pop)

tertile_label <- function(x, labels = c("Low", "Mid", "High")) {
  cut(x, breaks = quantile(x, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
      labels = labels, include.lowest = TRUE)
}

# Yuexiu: use green_area_m2 and pop_count from accessibility output
yx_bivar <- yx_sub_access |>
  mutate(
    admin_area_m2 = as.numeric(st_area(st_transform(geometry, CRS_YX))),
    green_density = green_area_m2 / admin_area_m2,  # fraction of area
    pop_density   = pop_count / (admin_area_m2 / 10000),  # per ha
    grn_class     = tertile_label(green_density),
    pop_class     = tertile_label(pop_density),
    bivar_class   = paste(grn_class, pop_class, sep = "-")
  )

dl_bivar <- dl_wijk_access |>
  mutate(
    admin_area_m2 = as.numeric(st_area(st_transform(geometry, CRS_DELFT))),
    green_density = green_area_m2 / admin_area_m2,
    pop_density   = pop_count / (admin_area_m2 / 10000),
    grn_class     = tertile_label(green_density),
    pop_class     = tertile_label(pop_density),
    bivar_class   = paste(grn_class, pop_class, sep = "-")
  )

# 3×3 bivariate colour palette (green axis × population axis)
bivar_palette <- c(
  "Low-Low"   = "#e8f4e8",  # low green, low pop  → pale green
  "Low-Mid"   = "#b8d4b8",
  "Low-High"  = "#5aaa5a",  # low green, high pop → equity concern
  "Mid-Low"   = "#d4e8f4",
  "Mid-Mid"   = "#7fb0cc",
  "Mid-High"  = "#2277aa",
  "High-Low"  = "#f4d4e8",
  "High-Mid"  = "#cc7faa",
  "High-High" = "#8833aa"   # high green, high pop → green-rich dense area
)

# ── 3C. Green space – income / nightlight correlation ─────────────────────────
# Delft: Pearson r between green_pc_m2 and mean disposable income
# Yuexiu: Pearson r between green_pc_m2 and mean VIIRS radiance

# --- Delft income ---
dl_inc_sf <- dl$dl_inc
# Extract a numeric income column (look for any column with "income", "inkomen", or "besteed")
inc_col <- names(dl_inc_sf)[grepl("income|inkomen|besteed|gemiddeld",
                                   names(dl_inc_sf), ignore.case = TRUE)][1]
if (is.na(inc_col)) {
  message("Warning: no income column auto-detected in dl_inc. Check column names:")
  print(names(dl_inc_sf))
  inc_col <- names(dl_inc_sf)[2]  # fallback: second column
}
message(sprintf("Using income column: %s", inc_col))

# Spatial join income to wijken
dl_inc_joined <- st_join(
  dl_wijk_access |> select(green_pc_m2, nearest_green_m),
  dl_inc_sf[, inc_col],
  join = st_intersects,
  largest = TRUE
) |>
  rename(income = all_of(inc_col)) |>
  mutate(income = as.numeric(income))

dl_corr <- cor(dl_inc_joined$green_pc_m2, dl_inc_joined$income,
                use = "complete.obs", method = "pearson")
message(sprintf("Delft Pearson r (green pc ~ income): %.3f", dl_corr))

# --- Yuexiu VIIRS ---
yx_viirs_zonal <- yx_sub_access |>
  mutate(
    viirs_mean = exact_extract(d$yx_viirs,
                                st_transform(yx_sub_access, st_crs(d$yx_viirs)),
                                "mean")
  )
yx_corr <- cor(yx_viirs_zonal$green_pc_m2, yx_viirs_zonal$viirs_mean,
                use = "complete.obs", method = "pearson")
message(sprintf("Yuexiu Pearson r (green pc ~ VIIRS): %.3f", yx_corr))

# ── Save outputs ──────────────────────────────────────────────────────────────
saveRDS(lorenz_data,     file.path(OUT_ROOT, "lorenz_data.rds"))
saveRDS(gini_labels,     file.path(OUT_ROOT, "gini_labels.rds"))
saveRDS(yx_bivar,        file.path(OUT_ROOT, "yx_bivar.rds"))
saveRDS(dl_bivar,        file.path(OUT_ROOT, "dl_bivar.rds"))
saveRDS(dl_inc_joined,   file.path(OUT_ROOT, "dl_inc_joined.rds"))
saveRDS(yx_viirs_zonal,  file.path(OUT_ROOT, "yx_viirs_zonal.rds"))

corr_summary <- data.frame(
  city      = c("Delft", "Yuexiu"),
  proxy     = c("Household income (CBS)", "Nighttime light (VIIRS)"),
  pearson_r = c(dl_corr, yx_corr)
)
saveRDS(corr_summary, file.path(OUT_ROOT, "corr_summary.rds"))

# ── Figures ───────────────────────────────────────────────────────────────────

# Figure 3A: Lorenz curves with Gini annotations
p_lorenz <- ggplot(lorenz_data, aes(x = cum_pop, y = cum_grn, colour = city)) +
  geom_line(linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_text(data = gini_labels,
            aes(x = 0.25, y = c(0.82, 0.72), label = label, colour = city),
            size = 4, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = c("Yuexiu" = "#21918c", "Delft" = "#440154")) +
  scale_x_continuous(labels = label_percent()) +
  scale_y_continuous(labels = label_percent()) +
  theme_minimal(base_size = 12) +
  labs(title = "Lorenz Curve — Green Space per Capita Inequality",
       subtitle = "Diagonal = perfect equality; further below = more unequal",
       x = "Cumulative share of neighbourhoods (ranked by green pc)",
       y = "Cumulative share of total green space",
       colour = "City")

ggsave(file.path(OUT_ROOT, "fig_lorenz_gini.png"),
       p_lorenz, width = 8, height = 6, dpi = 300)

# Figure 3B: Bivariate choropleth maps
p_bv_yx <- ggplot(yx_bivar) +
  geom_sf(aes(fill = bivar_class), colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = bivar_palette, na.value = "grey80",
                    name = "Green–Population\nmatrix") +
  theme_minimal() +
  labs(title = "Bivariate: Green Density × Population Density — Yuexiu",
       subtitle = "3×3 tertile classification")

p_bv_dl <- ggplot(dl_bivar) +
  geom_sf(aes(fill = bivar_class), colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = bivar_palette, na.value = "grey80",
                    name = "Green–Population\nmatrix") +
  theme_minimal() +
  labs(title = "Bivariate: Green Density × Population Density — Delft",
       subtitle = "3×3 tertile classification")

ggsave(file.path(OUT_ROOT, "fig_bivariate_choropleth.png"),
       p_bv_yx + p_bv_dl, width = 14, height = 6, dpi = 300)

# 3×3 legend tile
legend_df <- expand.grid(
  grn_class = c("Low", "Mid", "High"),
  pop_class = c("Low", "Mid", "High")
) |>
  mutate(bivar_class = paste(grn_class, pop_class, sep = "-"),
         grn_class = factor(grn_class, c("Low", "Mid", "High")),
         pop_class = factor(pop_class, c("Low", "Mid", "High")))

p_legend <- ggplot(legend_df, aes(x = grn_class, y = pop_class, fill = bivar_class)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_manual(values = bivar_palette, guide = "none") +
  theme_minimal(base_size = 10) +
  labs(x = "← Green density →", y = "← Population density →",
       title = "Bivariate legend")

ggsave(file.path(OUT_ROOT, "fig_bivariate_legend.png"),
       p_legend, width = 3.5, height = 3.5, dpi = 300)

# Figure 3C: Scatter plots — green pc vs. income/VIIRS
p_corr_dl <- ggplot(dl_inc_joined |> st_drop_geometry(),
                     aes(x = income, y = green_pc_m2)) +
  geom_point(colour = "#440154", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "#440154", fill = "#440154", alpha = 0.15) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("r = %.3f", dl_corr),
           hjust = 1.1, vjust = 1.5, size = 5, fontface = "bold") +
  scale_y_continuous(labels = label_comma()) +
  theme_minimal(base_size = 12) +
  labs(title = "Green Space per Capita vs. Household Income — Delft",
       subtitle = "Pearson correlation; wijk level",
       x = "Mean disposable income (€)", y = "Green space (m² per person)")

p_corr_yx <- ggplot(yx_viirs_zonal |> st_drop_geometry(),
                     aes(x = viirs_mean, y = green_pc_m2)) +
  geom_point(colour = "#21918c", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "#21918c", fill = "#21918c", alpha = 0.15) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("r = %.3f", yx_corr),
           hjust = 1.1, vjust = 1.5, size = 5, fontface = "bold") +
  scale_y_continuous(labels = label_comma()) +
  theme_minimal(base_size = 12) +
  labs(title = "Green Space per Capita vs. VIIRS Nighttime Light — Yuexiu",
       subtitle = "Pearson correlation; subdistrict level",
       x = "Mean VIIRS radiance (nW/cm²/sr)", y = "Green space (m² per person)")

ggsave(file.path(OUT_ROOT, "fig_equity_correlations.png"),
       p_corr_dl + p_corr_yx, width = 14, height = 6, dpi = 300)

message("Script 04 complete — spatial justice figures saved.")