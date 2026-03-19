############################################################
## LOAD LIBRARIES
############################################################
library(terra)
library(ggplot2)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(cmocean)
library(patchwork)

############################################################
## 1. WORLD SHAPEFILE
############################################################
world_sf <- ne_countries(scale = "medium", returnclass = "sf")
world_vect <- vect(world_sf)

world_bold <- world_sf 
  
############################################################
## 2. LOAD FUNCTION
############################################################
load_and_prepare <- function(path) {
  r <- rast(path)
  r <- mask(crop(r, ext(world_vect)), world_vect)
  as.data.frame(r, xy = TRUE, na.rm = TRUE) %>%
    rename(vc = 3)
}

############################################################
## 3. LOAD ALL RASTERS (UPDATED PATHS)
############################################################

# -------- FL --------
R0_FL_CF <- load_and_prepare('outputs/R0_FL_counterfactual.tif')
R0_FL_F  <- load_and_prepare('outputs/R0_FL_factual.tif')
R0_FL_D  <- load_and_prepare('outputs/R0_FL_delta.tif')

# -------- NC --------
R0_NC_CF <- load_and_prepare('outputs/R0_NC_counterfactual.tif')
R0_NC_F  <- load_and_prepare('outputs/R0_NC_factual.tif')
R0_NC_D  <- load_and_prepare('outputs/R0_NC_delta.tif')

# -------- CA --------
R0_CA_CF <- load_and_prepare('outputs/R0_CA_counterfactual.tif')
R0_CA_F  <- load_and_prepare('outputs/R0_CA_factual.tif')
R0_CA_D  <- load_and_prepare('outputs/R0_CA_delta.tif')


############################################################
## DEFINE COLOR SCALE + PLOT FUNCTION
############################################################

vc_breaks <- c(-Inf, 1, 5, 10, 20, 40, Inf)

vc_labels <- c("<1","1–5","5–10","10–20","20–40",">40")

vc_colors <- c(
  "#fff5eb","#fdd0a2","#fdae6b",
  "#fd8d3c","#e6550d","#a63603"
)

world_bold <- world_sf

make_vc_plot <- function(df, title, legend_title = "R0 class") {
  
  df$vc_class <- cut(
    df$vc,
    breaks = vc_breaks,
    labels = vc_labels,
    include.lowest = TRUE,
    ordered_result = TRUE
  )
  
  ggplot() +
    geom_tile(data = df, aes(x = x, y = y, fill = vc_class)) +
    
    geom_sf(
      data = world_bold,
      fill = NA,
      color = "black",
      linewidth = 0.4,
      inherit.aes = FALSE
    ) +
    
    coord_sf(xlim = c(-180,180), ylim = c(-60,60), expand = FALSE) +
    
    scale_fill_manual(
      values = vc_colors,
      drop = FALSE,
      name = legend_title
    ) +
    
    labs(title = title) +
    
    theme_void() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = c(0.8, 0.2),
      plot.margin = ggplot2::margin(2,2,2,2)
    )
}

## FIGURE 1: FACTUAL + ΔR0 (MAIN)
############################################################
############################################################

# -------- FL --------
p1 <- make_vc_plot(R0_FL_F,  "a) FL – Factual")
p2 <- make_vc_plot(R0_FL_D,  "b) FL – ΔR0")

# -------- NC --------
p3 <- make_vc_plot(R0_NC_F,  "c) NC – Factual")
p4 <- make_vc_plot(R0_NC_D,  "d) NC – ΔR0")

# -------- CA --------
p5 <- make_vc_plot(R0_CA_F,  "e) CA – Factual")
p6 <- make_vc_plot(R0_CA_D,  "f) CA – ΔR0")

final_plot <- (p1 | p2) /
  (p3 | p4) /
  (p5 | p6) +
  plot_layout(guides = "collect")

print(final_plot)

ggsave(
  "figure_R0_main.png",
  final_plot,
  width = 14,
  height = 10,
  dpi = 300,
  bg = "white"
)

############################################################
############################################################
## FIGURE 2: AGREEMENT MAP (NO RESAMPLING)
############################################################
############################################################

# -------- DELTA (USE ORIGINAL RASTERS) --------
dCA <- rast('outputs/R0_CA_delta.tif')
dFL <- rast('outputs/R0_FL_delta.tif')
dNC <- rast('outputs/R0_NC_delta.tif')

# -------- MASK TO LAND --------
mask_land <- function(r){
  mask(crop(r, ext(world_vect)), world_vect)
}

dCA <- mask_land(dCA)
dFL <- mask_land(dFL)
dNC <- mask_land(dNC)

# -------- SIGN --------
sign_CA <- sign(dCA)
sign_FL <- sign(dFL)
sign_NC <- sign(dNC)

# -------- AGREEMENT SCORE --------
agreement <- (sign_CA == sign_FL) + 
  (sign_CA == sign_NC) + 
  (sign_FL == sign_NC)

names(agreement) <- "agreement"

# -------- GLOBAL AGREEMENT METRIC (FIXED) --------
df_vals <- data.frame(
  v_CA = as.vector(values(dCA)),
  v_FL = as.vector(values(dFL)),
  v_NC = as.vector(values(dNC))
)

# remove NA
df_vals <- na.omit(df_vals)

# sign agreement (all three same sign)
agree_all <- mean(
  sign(df_vals$v_CA) == sign(df_vals$v_FL) &
    sign(df_vals$v_CA) == sign(df_vals$v_NC)
) * 100
# format label
agree_label <- paste0("Global agreement = ", round(agree_all,1), "%")

# -------- DATAFRAME --------
df_agree <- as.data.frame(agreement, xy = TRUE, na.rm = TRUE)

p_agree <- ggplot(df_agree, aes(x = x, y = y, fill = factor(agreement))) +
  geom_tile() +
  
  geom_sf(data = world_bold, fill = NA, color = "black",
          linewidth = 0.4, inherit.aes = FALSE) +
  
  coord_sf(xlim = c(-180,180), ylim = c(-60,60), expand = FALSE) +
  
  scale_fill_manual(
    values = c(
      "0" = "#f7f7f7",
      "1" = "#fee8c8",
      "2" = "#fdbb84",
      "3" = "#e34a33"
    ),
    name = "Agreement",
    labels = c("None","Low","Medium","High")
  ) +
  
  # 🔥 ADD THIS
  annotate("text",
           x = -170, y = 55,
           label = agree_label,
           hjust = 0,
           size = 4,
           fontface = "bold") +
  
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold"),
    plot.margin = ggplot2::margin(2,2,2,2)
  ) +
  
  labs(title = "Agreement in ΔR0 across parameterizations")

ggsave(
  "figure_R0_agreement.png",
  p_agree,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)
############################################################
############################################################
## FIGURE 3: ΔR0 COMPARISON (UPDATED & MEANINGFUL)
############################################################
############################################################

############################################################
## 7. LOAD DELTA RASTERS (NO RESAMPLING)
############################################################

dCA <- rast('outputs/R0_CA_delta.tif')
dFL <- rast('outputs/R0_FL_delta.tif')
dNC <- rast('outputs/R0_NC_delta.tif')

# mask to land
mask_land <- function(r){
  mask(crop(r, ext(world_vect)), world_vect)
}

dCA <- mask_land(dCA)
dFL <- mask_land(dFL)
dNC <- mask_land(dNC)

names(dCA) <- "CA"
names(dFL) <- "FL"
names(dNC) <- "NC"

stack_all <- c(dCA, dFL, dNC)

############################################################
## 8. DATAFRAME + CONTINENTS (ROBUST FIX)
############################################################

df <- as.data.frame(stack_all, xy = TRUE, na.rm = TRUE)

points_sf <- st_as_sf(df, coords = c("x", "y"), crs = st_crs(world_sf))

world_cont <- ne_countries(scale = "medium", returnclass = "sf")[, "continent"]

points_sf <- suppressWarnings(st_join(points_sf, world_cont))

df <- points_sf %>%
  st_drop_geometry()

df <- df %>%
  dplyr::filter(continent %in% c(
    "Africa","Asia","South America",
    "North America","Europe","Oceania"
  ))

############################################################
## 9. WIDE FORMAT
############################################################

df_wide <- df %>%
  dplyr::select(CA, FL, NC, continent) %>%
  tidyr::drop_na()

############################################################
## METRICS FUNCTION
############################################################

get_metrics <- function(x, y){
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  
  r   <- cor(x, y, method = "pearson")
  rho <- cor(x, y, method = "spearman")
  agree <- mean(sign(x) == sign(y)) * 100
  
  paste0(
    "r = ", round(r,2),
    "\nρ = ", round(rho,2),
    "\nAgree = ", round(agree,1), "%"
  )
}

############################################################
## 10. PAIRWISE COMPARISON PLOTS
############################################################

# FL vs CA
p1 <- ggplot(df_wide, aes(x = CA, y = FL, color = continent)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate("text",
           x = min(df_wide$CA, na.rm=TRUE),
           y = max(df_wide$FL, na.rm=TRUE),
           label = get_metrics(df_wide$CA, df_wide$FL),
           hjust = 0, vjust = 1,
           size = 3.5, fontface = "bold") +
  labs(x = "ΔR0 — CA", y = "ΔR0 — FL", title = "a) FL vs CA") +
  scale_color_brewer(palette = "Set2") +
  theme_classic()

# FL vs NC
p2 <- ggplot(df_wide, aes(x = NC, y = FL, color = continent)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate("text",
           x = min(df_wide$NC, na.rm=TRUE),
           y = max(df_wide$FL, na.rm=TRUE),
           label = get_metrics(df_wide$NC, df_wide$FL),
           hjust = 0, vjust = 1,
           size = 3.5, fontface = "bold") +
  labs(x = "ΔR0 — NC", y = "ΔR0 — FL", title = "b) FL vs NC") +
  scale_color_brewer(palette = "Set2") +
  theme_classic()

# CA vs NC
p3 <- ggplot(df_wide, aes(x = NC, y = CA, color = continent)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate("text",
           x = min(df_wide$NC, na.rm=TRUE),
           y = max(df_wide$CA, na.rm=TRUE),
           label = get_metrics(df_wide$NC, df_wide$CA),
           hjust = 0, vjust = 1,
           size = 3.5, fontface = "bold") +
  labs(x = "ΔR0 — NC", y = "ΔR0 — CA", title = "c) CA vs NC") +
  scale_color_brewer(palette = "Set2") +
  theme_classic()

############################################################
## 11. COMBINE PANEL
############################################################

final_scatter <- p1 | p2 | p3

print(final_scatter)

ggsave(
  "figure_delta_comparison.png",
  final_scatter,
  width = 12,
  height = 4.5,
  dpi = 300,
  bg = "white"
)
