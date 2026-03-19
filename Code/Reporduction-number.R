############################################################
## 1. Clear workspace
############################################################
rm(list = ls())
cat("\014")

############################################################
## 2. Libraries
############################################################
library(dplyr)
library(tidyr)
library(ggplot2)

############################################################
## 3. Development functions
## Temperature-only Brière functions fitted from digitized data
## RH is kept in the function signature for compatibility
############################################################

briere_rate <- function(T, a, Tmin, Tmax, m){
  val <- a * T * (T - Tmin) * (Tmax - T)^(1 / m)
  val[T <= Tmin | T >= Tmax] <- 0
  pmax(0, val)
}

# -------- FL --------
egg_rate_FL <- function(T, RH = NULL){
  briere_rate(T, a = 8.17e-5, Tmin = 11.8, Tmax = 35.8, m = 10)
}
larva_rate_FL <- function(T, RH = NULL){
  briere_rate(T, a = 3.54e-7, Tmin = 17.0, Tmax = 57.2, m = 0.5)
}
nymph_rate_FL <- function(T, RH = NULL){
  briere_rate(T, a = 2.38e-4, Tmin = 17.0, Tmax = 35.5, m = 10)
}

# -------- NC --------
egg_rate_NC <- function(T, RH = NULL){
  briere_rate(T, a = 2.02e-7, Tmin = 20.0, Tmax = 60.0, m = 0.5)
}
larva_rate_NC <- function(T, RH = NULL){
  briere_rate(T, a = 1.88e-4, Tmin = 12.8, Tmax = 35.0, m = 10)
}
nymph_rate_NC <- function(T, RH = NULL){
  briere_rate(T, a = 4.60e-5, Tmin = 10.2, Tmax = 43.1, m = 2.58)
}

# -------- CA --------
egg_rate_CA <- function(T, RH = NULL){
  briere_rate(T, a = 7.61e-5, Tmin = 14.0, Tmax = 36.8, m = 3.81)
}
larva_rate_CA <- function(T, RH = NULL){
  briere_rate(T, a = 2.31e-4, Tmin = 15.3, Tmax = 35.0, m = 10)
}
nymph_rate_CA <- function(T, RH = NULL){
  briere_rate(T, a = 9.85e-5, Tmin = 9.79, Tmax = 39.1, m = 10)
}

############################################################
## 4. Mortality functions
## Table 2 gives weekly mortality probability:
## MP = b1 + b2*T + b3*RH + b4*(T*RH)
## Convert to daily rate: mu = -log(1 - MP)/7
############################################################

mp_to_rate <- function(MP){
  MP <- pmax(0.001, pmin(1 - 1e-6, MP))  # 🔥 minimum mortality
  -log(1 - MP) / 7
}

# -------- FL --------
muL_FL <- function(T, RH){
  mp_to_rate(2.42 + 0.03*T - 0.03*RH - 4.58e-5*(T*RH))
}
muP_FL <- function(T, RH){
  mp_to_rate(-2.61 + 0.12*T + 0.03*RH - 1.45e-3*(T*RH))
}
muA_FL <- function(T, RH){
  mp_to_rate(-0.11 + 0.01*T + 9.01e-4*RH - 4.93e-5*(T*RH))
}

# -------- NC --------
muL_NC <- function(T, RH){
  mp_to_rate(-0.94 + 0.12*T + 0.01*RH - 1.27e-3*(T*RH))
}
muP_NC <- function(T, RH){
  mp_to_rate(-0.57 + 0.05*T + 0.01*RH - 5.26e-4*(T*RH))
}
muA_NC <- function(T, RH){
  mp_to_rate(-0.05 + 0.01*T + 6.00e-4*RH - 5.03e-5*(T*RH))
}

# -------- CA --------
muL_CA <- function(T, RH){
  mp_to_rate(2.88 + 0.01*T - 0.04*RH + 2.94e-4*(T*RH))
}
muP_CA <- function(T, RH){
  mp_to_rate(-1.58 + 0.08*T + 0.01*RH - 7.22e-4*(T*RH))
}
muA_CA <- function(T, RH){
  mp_to_rate(-0.31 + 0.02*T + 3.05e-3*RH - 1.43e-4*(T*RH))
}

############################################################
## 5. Fecundity
############################################################

oviposition_fun <- function(T){
  ifelse(T < 15, 0,
         ifelse(T <= 35, 1, 
                ifelse(T <= 40, 0.5, 0)))
}

F_fun <- function(T, F0 = 805.4, p_female = 0.5){
  
  ovip_days <- 14.3   # fixed biological constant
  
  (F0 / ovip_days) * oviposition_fun(T) * p_female
}

############################################################
## 6. R0 function
## Simplified form with M_E = 0:
## R0 = (F * G2 * G3) / ((G2 + M_L) * (G3 + M_P) * M_A)
############################################################

R0_strain <- function(T, RH,
                      dev_fun_egg, dev_fun_larva, dev_fun_nymph,
                      mort_fun_L, mort_fun_P, mort_fun_A,
                      F0 = 805.4,
                      p_female = 0.5) {
  
  # Development rates
  G1 <- dev_fun_egg(T)    # kept for completeness
  G2 <- dev_fun_larva(T)
  G3 <- dev_fun_nymph(T)
  
  # Mortality rates
  ML <- mort_fun_L(T, RH)
  MP <- mort_fun_P(T, RH)
  MA <- mort_fun_A(T, RH)
  
  # Effective reproduction
  F <- F_fun(T, F0 = F0, p_female = p_female)
  
  # Simplified R0
  R0 <- (F * G2 * G3) / ((G2 + ML) * (G3 + MP) * MA)
  
  # Numerical safety
  bad <- is.nan(R0) | is.infinite(R0) | R0 < 0
  R0[bad] <- 0
  
  return(R0)
}

############################################################
## 7. Grid
############################################################

T_vec  <- seq(-5, 40, by = 0.5)
RH_vec <- seq(17, 92, by = 5)

grid <- expand.grid(
  T  = T_vec,
  RH = RH_vec
)

############################################################
## 8. Compute R0
############################################################

grid$R0_FL <- R0_strain(
  T = grid$T, RH = grid$RH,
  dev_fun_egg = egg_rate_FL,
  dev_fun_larva = larva_rate_FL,
  dev_fun_nymph = nymph_rate_FL,
  mort_fun_L = muL_FL,
  mort_fun_P = muP_FL,
  mort_fun_A = muA_FL
)

grid$R0_NC <- R0_strain(
  T = grid$T, RH = grid$RH,
  dev_fun_egg = egg_rate_NC,
  dev_fun_larva = larva_rate_NC,
  dev_fun_nymph = nymph_rate_NC,
  mort_fun_L = muL_NC,
  mort_fun_P = muP_NC,
  mort_fun_A = muA_NC
)

grid$R0_CA <- R0_strain(
  T = grid$T, RH = grid$RH,
  dev_fun_egg = egg_rate_CA,
  dev_fun_larva = larva_rate_CA,
  dev_fun_nymph = nymph_rate_CA,
  mort_fun_L = muL_CA,
  mort_fun_P = muP_CA,
  mort_fun_A = muA_CA
)

############################################################
## 9. Prepare plotting data
############################################################

grid_long <- grid %>%
  pivot_longer(
    cols = starts_with("R0_"),
    names_to = "Scenario",
    values_to = "R0"
  ) %>%
  mutate(
    Scenario = recode(
      Scenario,
      R0_FL = "FL",
      R0_NC = "NC",
      R0_CA = "CA"
    ),
    R0_class = cut(
      R0,
      breaks = c(-Inf, 1, 5, 10, 20, 40, Inf),
      labels = c("<1", "1–5", "5–10", "10–20", "20–40", ">40"),
      right = FALSE
    )
  )

############################################################
## 10. Plot
############################################################

cols <- c(
  "#fff5eb", "#fdd0a2", "#fdae6b",
  "#fd8d3c", "#e6550d", "#a63603"
)

p <- ggplot(grid_long, aes(x = T, y = RH, fill = R0_class)) +
  geom_tile() +
  scale_fill_manual(values = cols, name = expression(R[0])) +
  facet_wrap(~Scenario, ncol = 3) +
  labs(
    x = "Temperature (°C)",
    y = "Relative Humidity (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90", colour = NA),
    legend.position = "right",
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

print(p)

############################################################
## 11. Save outputs
############################################################

ggsave(
  filename = "R0_heatmap_simple.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "R0_heatmap_simple.tiff",
  plot = p,
  width = 10,
  height = 6,
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = "R0_heatmap_simple.pdf",
  plot = p,
  width = 10,
  height = 6,
  bg = "white"
)

############################################################
## 12. Optional: print fitted development equations
############################################################

cat("\nDevelopment equations used:\n")
cat("FL Egg   : a=8.17e-5, Tmin=11.8, Tmax=35.8, m=10\n")
cat("FL Larva : a=3.54e-7, Tmin=17.0, Tmax=57.2, m=0.5\n")
cat("FL Nymph : a=2.38e-4, Tmin=17.0, Tmax=35.5, m=10\n\n")

cat("NC Egg   : a=2.02e-7, Tmin=20.0, Tmax=60.0, m=0.5\n")
cat("NC Larva : a=1.88e-4, Tmin=12.8, Tmax=35.0, m=10\n")
cat("NC Nymph : a=4.60e-5, Tmin=10.2, Tmax=43.1, m=2.58\n\n")

cat("CA Egg   : a=7.61e-5, Tmin=14.0, Tmax=36.8, m=3.81\n")
cat("CA Larva : a=2.31e-4, Tmin=15.3, Tmax=35.0, m=10\n")
cat("CA Nymph : a=9.85e-5, Tmin=9.79, Tmax=39.1, m=10\n")
