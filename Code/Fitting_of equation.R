############################################################
## 1. LIBRARIES
############################################################
library(dplyr)
library(purrr)
library(minpack.lm)
library(Metrics)
library(pracma)

############################################################
## 2. DATA (DIGITIZED)
############################################################

dev_data <- data.frame(
  Strain = rep(c("FL","NC","CA"), each = 15),
  Stage  = rep(rep(c("Egg","Larva","Nymph"), each = 5), times = 3),
  T = rep(c(20,23,27,30,35), times = 9),
  rate = c(
    # FL
    0.02,0.025,0.04,0.055,0.065,
    0.03,0.055,0.09,0.10,0.11,
    0.02,0.04,0.08,0.11,0.14,
    # NC
    0.02,0.035,0.045,0.055,0.06,
    0.04,0.055,0.08,0.12,0.10,
    0.03,0.045,0.06,0.075,0.09,
    # CA
    0.02,0.03,0.05,0.06,0.065,
    0.03,0.05,0.09,0.12,0.10,
    0.03,0.035,0.06,0.075,0.10
  )
)

############################################################
## 3. MODEL FUNCTIONS
############################################################

briere_fun <- function(T, a, Tmin, Tmax, m){
  a * T * (T - Tmin) * (Tmax - T)^(1/m)
}

performance_fun <- function(T, b, Tmin, c, Tmax){
  b*(T - Tmin)*(1 - exp(c*(T - Tmax)))
}

############################################################
## 4. SAFE FITTING FUNCTIONS
############################################################

fit_briere_safe <- function(T, y){
  try(nlsLM(
    y ~ briere_fun(T, a, Tmin, Tmax, m),
    start = list(a=1e-4, Tmin=min(T)-5, Tmax=max(T)+5, m=2),
    lower = c(a=1e-8, Tmin=0, Tmax=max(T), m=0.5),
    upper = c(a=1, Tmin=min(T), Tmax=60, m=10),
    control = nls.lm.control(maxiter=500)
  ), silent=TRUE)
}

fit_perf_safe <- function(T, y){
  try(nlsLM(
    y ~ performance_fun(T, b, Tmin, c, Tmax),
    start = list(b=0.01, Tmin=min(T)-5, c=-0.1, Tmax=max(T)+5),
    lower = c(b=1e-6, Tmin=0, c=-10, Tmax=max(T)),
    upper = c(b=1, Tmin=min(T), c=0, Tmax=60),
    control = nls.lm.control(maxiter=500)
  ), silent=TRUE)
}

############################################################
## 5. FIT MODELS + STORE PARAMETERS
############################################################

fit_models <- function(df){
  
  T <- df$T
  y <- df$rate
  
  fit_briere <- fit_briere_safe(T, y)
  fit_perf   <- fit_perf_safe(T, y)
  
  res_list <- list()
  
  # -------- BRIERE --------
  if(!inherits(fit_briere, "try-error")){
    
    pred <- predict(fit_briere)
    coefs <- coef(fit_briere)
    
    res_list[[1]] <- data.frame(
      Model="Briere",
      RMSE=rmse(y, pred),
      R2=cor(y, pred)^2,
      AUC=trapz(T, pred),
      Converged=TRUE,
      a=coefs["a"],
      Tmin=coefs["Tmin"],
      Tmax=coefs["Tmax"],
      m=coefs["m"],
      b=NA, c=NA
    )
    
  } else {
    
    res_list[[1]] <- data.frame(
      Model="Briere",
      RMSE=NA, R2=NA, AUC=NA, Converged=FALSE,
      a=NA, Tmin=NA, Tmax=NA, m=NA,
      b=NA, c=NA
    )
  }
  
  # -------- PERFORMANCE --------
  if(!inherits(fit_perf, "try-error")){
    
    pred <- predict(fit_perf)
    coefs <- coef(fit_perf)
    
    res_list[[2]] <- data.frame(
      Model="Performance",
      RMSE=rmse(y, pred),
      R2=cor(y, pred)^2,
      AUC=trapz(T, pred),
      Converged=TRUE,
      a=NA, m=NA,
      Tmin=coefs["Tmin"],
      Tmax=coefs["Tmax"],
      b=coefs["b"],
      c=coefs["c"]
    )
    
  } else {
    
    res_list[[2]] <- data.frame(
      Model="Performance",
      RMSE=NA, R2=NA, AUC=NA, Converged=FALSE,
      a=NA, Tmin=NA, Tmax=NA, m=NA,
      b=NA, c=NA
    )
  }
  
  bind_rows(res_list) %>% arrange(RMSE)
}

############################################################
## 6. APPLY TO ALL STRAIN × STAGE
############################################################

results <- dev_data %>%
  group_by(Strain, Stage) %>%
  group_modify(~fit_models(.x)) %>%
  ungroup()

############################################################
## 7. RANK MODELS
############################################################

results_ranked <- results %>%
  group_by(Strain, Stage) %>%
  arrange(RMSE, .by_group = TRUE) %>%
  mutate(Rank = row_number()) %>%
  ungroup()

############################################################
## 8. BEST MODELS (WITH PARAMETERS)
############################################################

best_models <- results_ranked %>%
  filter(Rank == 1)

############################################################
## 9. PRINT RESULTS
############################################################

print(results_ranked)

cat("\n================ BEST MODELS (WITH PARAMETERS) ================\n")
print(best_models)
