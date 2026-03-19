

############################################################
## 0. SETUP
############################################################

setwd('/Users/kagboka/Desktop/Anaplasma/')

library(terra)
library(ggplot2)
library(tidyverse)
library(sf)
library(rnaturalearth)

############################################################
## 1. LOAD RASTERS
############################################################

load_mean_raster <- function(folder){
  files <- list.files(folder, pattern="\\.tif$", full.names=TRUE)
  mean(rast(files), na.rm=TRUE)
}

# factual
r_T  <- load_mean_raster('/Users/kagboka/Desktop/Ritter_work/LTM_tif/factual/tas/')
r_RH <- load_mean_raster('/Users/kagboka/Desktop/Ritter_work/LTM_tif/factual/hurs/')

# counterfactual
r_T_cf  <- load_mean_raster('/Users/kagboka/Desktop/Ritter_work/LTM_tif/counterfactual/tas/')
r_RH_cf <- load_mean_raster('/Users/kagboka/Desktop/Ritter_work/LTM_tif/counterfactual/hurs/')

############################################################
## 2. HARMONISE
############################################################

align_rasters <- function(rT, rRH){
  rRH <- project(rRH, rT, method="bilinear")
  ext_common <- terra::ext(
    max(xmin(rT), xmin(rRH)),
    min(xmax(rT), xmax(rRH)),
    max(ymin(rT), ymin(rRH)),
    min(ymax(rT), ymax(rRH))
  )
  list(
    T  = crop(rT, ext_common),
    RH = crop(rRH, ext_common)
  )
}

f  <- align_rasters(r_T, r_RH)
cf <- align_rasters(r_T_cf, r_RH_cf)

r_in    <- c(f$T,  f$RH)
r_in_cf <- c(cf$T, cf$RH)

names(r_in)    <- c("T","RH")
names(r_in_cf) <- c("T","RH")

############################################################
## 3. R0 FUNCTIONS
############################################################

R0_fun_factory <- function(devE, devL, devN, muL, muP, muA){
  function(T, RH){
    R0_strain(
      T, RH,
      devE, devL, devN,
      muL, muP, muA
    )
  }
}

R0_FL_fun <- R0_fun_factory(egg_rate_FL, larva_rate_FL, nymph_rate_FL, muL_FL, muP_FL, muA_FL)
R0_NC_fun <- R0_fun_factory(egg_rate_NC, larva_rate_NC, nymph_rate_NC, muL_NC, muP_NC, muA_NC)
R0_CA_fun <- R0_fun_factory(egg_rate_CA, larva_rate_CA, nymph_rate_CA, muL_CA, muP_CA, muA_CA)

############################################################
## 4. APPLY MODEL (ALL STRAINS)
############################################################

run_model <- function(fun){
  list(
    factual = lapp(r_in,    fun=function(T,RH) fun(T,RH)),
    cf      = lapp(r_in_cf, fun=function(T,RH) fun(T,RH))
  )
}

FL <- run_model(R0_FL_fun)
NC <- run_model(R0_NC_fun)
CA <- run_model(R0_CA_fun)

# delta
FL$delta <- FL$factual - FL$cf
NC$delta <- NC$factual - NC$cf
CA$delta <- CA$factual - CA$cf

############################################################
## 5. MASK LAND
############################################################

world_vect <- vect(ne_countries(scale="medium", returnclass="sf"))

mask_all <- function(r){
  mask(crop(r, ext(world_vect)), world_vect)
}

FL <- lapply(FL, mask_all)
NC <- lapply(NC, mask_all)
CA <- lapply(CA, mask_all)


############################################################
## 6. SAVE RASTERS (CORRECT)
############################################################

out_dir <- "outputs"
dir.create(out_dir, showWarnings=FALSE)

save_all <- function(obj, prefix){
  writeRaster(obj$factual, file.path(out_dir, paste0(prefix,"_factual.tif")), overwrite=TRUE)
  writeRaster(obj$cf,      file.path(out_dir, paste0(prefix,"_counterfactual.tif")), overwrite=TRUE)
  writeRaster(obj$delta,   file.path(out_dir, paste0(prefix,"_delta.tif")), overwrite=TRUE)
}

save_all(FL, "R0_FL")
save_all(NC, "R0_NC")
save_all(CA, "R0_CA")
############################################################
## 7. PLOT (FL ONLY EXAMPLE)
############################################################

to_df <- function(r, label){
  df <- as.data.frame(r, xy=TRUE, na.rm=TRUE)
  names(df)[3] <- "value"
  df$scenario <- label
  df
}

df_all <- bind_rows(
  to_df(FL$factual, "Factual"),
  to_df(FL$cf,      "Counterfactual"),
  to_df(FL$delta,   "Delta")
)

p <- ggplot(df_all, aes(x=x, y=y, fill=value)) +
  geom_tile() +   # FIXED WARNING
  facet_wrap(~scenario) +
  coord_equal() +
  theme_bw()

print(p)

############################################################
## 13. SAVE
############################################################

ggsave("R0_FL_climate_comparison.png", p, width = 12, height = 5, dpi = 300)
ggsave("R0_FL_climate_comparison.tiff", p, width = 12, height = 5, dpi = 300, compression = "lzw")

