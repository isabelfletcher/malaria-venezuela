#######################################################################

######## Spatiotemporal model of malaria incidence in Bolivar ########

#######################################################################

## Load libraries
library(raster)
library(INLA)
library(sf)
library(sp)
library(dplyr)
library(rgdal)
library(MASS)
library(Metrics)
library(INLAutils)

### Define neighbourhood matrix
venezuela <- rgdal::readOGR("data/ven_shp", "Parroquias_Venezuela")
bolivar <- subset(venezuela, venezuela$ESTADO == "Bolívar")

nb.map <- poly2nb(bolivar)
nb2INLA("map.graph", nb.map)

## Read in data
data <- read.csv("data/inla_input/data.csv")

# Fixed effects
n    <- length(data$PF)

# Add API so modify offset to include it 
e  <- (data$Population/12)/1000

# Land use variables
mines           <- scale(data$total_mines)[,1]
forest_decrease <- scale(data$forest_decrease_cum)[,1]
urbanization    <- scale(data$urban_increase_cum)[,1]

# Socioeconomic data
health_travel <- scale(data$health_travel_time)[,1]

# Spatial random effects
s1 <- rep(1:46, 252)
Parroquia <- data$Parroquia
Municipio <- data$Municipio

# Monthly/seasonal random effects
t1 <- rep(rep(1:12, each = 46), 21)

# Yearly random effects
t2 <- rep(1996:2016, each = 552)

############################################################################
## Function to run spatiotemporal models of malaria incidence

run_models <- function(parasite, final_model){
  
  ## Condition on parasite
  if (parasite == "P. falciparum") {
    
    # Cases
    y <- data$PF
    
    # Climate variables - choosing most appropriate monthly time lags
    temp <- scale(data$temp_lag2)[,1]
    prcp <- scale(data$prcp_lag0)[,1]
    
    # Model filename
    mod_filename <- "pf"
    
    } else { # P. vivax
    
    # Cases
    y <- data$PV
    
    # Climate variables
    temp <- scale(data$temp_lag2)[,1]
    prcp <- scale(data$prcp_lag0)[,1]

    
    # Model filename
    mod_filename <- "pv"
    
  }
  
  df_inla <- data.frame(y, e,
                        Parroquia, Municipio,
                        s1, t1, t2,
                        temp, prcp, 
                        mines, forest_decrease,
                        urbanization, 
                        health_travel)
  
  ## Model formulas to be tested
  formulas <- c(# Spatial and monthly random effects
    y ~ 1 + f(s1, model = "bym", graph = "map.graph") +
      f(t1, model = "rw1"),
    # Spatial, monthly and annual random effects
    y ~ 1 + f(s1, model = "bym", graph = "map.graph") +
      f(t1, model = "rw1") +
      f(t2, model = "iid"),
    # ... + temperature
    y ~ 1 + f(s1, model = "bym", graph = "map.graph") +
      f(t1, model = "rw1") +
      f(t2, model = "iid") +
      f(inla.group(temp), model = "rw1"),
    # ... + precipitation
    y ~ 1 + f(s1, model = "bym", graph = "map.graph") +
      f(t1, model = "rw1") +
      f(t2, model = "iid") +
      f(inla.group(temp), model = "rw1") +
      f(inla.group(prcp), model = "rw1"),
    # ... + mines
    y ~ 1 + f(s1, model = "bym", graph = "map.graph") +
      f(t1, model = "rw1") +
      f(t2, model = "iid") +
      f(inla.group(temp), model = "rw1") +
      f(inla.group(prcp), model = "rw1") +
      mines,
    # ... + forest decrease
    y ~ 1 + f(s1, model = "bym", graph = "map.graph") +
      f(t1, model = "rw1") +
      f(t2, model = "iid") +
      f(inla.group(temp), model = "rw1") +
      f(inla.group(prcp), model = "rw1") +
      mines + 
      forest_decrease,
    # ... + urbanization
    y ~ 1 + f(s1, model = "bym", graph = "map.graph") +
      f(t1, model = "rw1") +
      f(t2, model = "iid") +
      f(inla.group(temp), model = "rw1") +
      f(inla.group(prcp), model = "rw1") +
      mines + 
      forest_decrease + 
      urbanization,
    # ... + health travel time 
    y ~ 1 + f(s1, model = "bym", graph = "map.graph") +
      f(t1, model = "rw1") +
      f(t2, model = "iid") +
      f(inla.group(temp), model = "rw1") +
      f(inla.group(prcp), model = "rw1") +
      mines + 
      forest_decrease + 
      urbanization +
      health_travel)
  
    ## Run models
    model_dic <- NULL
  
    for (i in 1:length(formulas)){
    
    mod <- inla(formulas[[i]], data = df_inla, family = "zeroinflatednbinomial0", 
                offset = log(e), verbose = TRUE,
                control.inla = list(strategy = 'adaptive'),
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                       config = TRUE, 
                                       return.marginals = TRUE), 
                control.predictor = list(link = 1, compute = TRUE), 
                control.family = list(link = "log"))
    
    # Save model
    file_path <- paste0("models/mod", i, "_", mod_filename, ".RData")
    save(mod, file = file_path)
    
    model_formula <- as.character(row.names(mod$summary.fixed)[2:length(row.names(mod$summary.fixed))])
    model_formula <- paste(model_formula, collapse = ", ")
    model_formula_nl <- paste(as.character(names(mod$summary.random)[5:6]), collapse = ", ")
    model_formula <- paste(model_formula_nl, model_formula, sep = ",")
    
    ## DIC values for all models
    df_dic <- data.frame(Model     = i,
                         Formula   = model_formula,
                         Parasite  = parasite, 
                         DIC       = mod$dic$dic,
                         WAIC      = mod$waic$waic,
                         log_score = format(-mean(log(mod$cpo$cpo), na.rm = T),scientific=F))

    ## Model RMSE
    df_rmse <- data.frame(cases    = y,
                          fit      = mod$summary.fitted.values$`0.5quant`) %>%
               dplyr::mutate(e = (fit - cases)^2) %>%
            dplyr::summarise(mae  = mae(cases, fit),
                             mse  = mse(cases, fit),
                             rmse = sqrt(mean(e, na.rm = T)))
    
    # Add to dic df
    df_dic <- df_dic %>% mutate(mae  = df_rmse$mae,
                                mse  = df_rmse$mse,
                                rmse = df_rmse$rmse)
    model_dic <- rbind(model_dic, df_dic)

    }
    
    # Get estimates, fitted values and predictions for the final model
    files <- list.files("models", ".RData", full.names = TRUE)
    model <- files[grepl(paste0(final_model, "_", mod_filename), files)]
    load(model)
    
    ## Covariate estimates for selected model
    model_estimates <- mod$summary.fixed %>%
      dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
      mutate(Model    = final_model,
             Parasite = parasite,
             Variable = rownames(mod$summary.fixed)) %>%
      dplyr::slice(-1) %>% 
      dplyr::rename(lci = `0.025quant`,
                    uci = `0.975quant`)
    
    ## Non-linear climate relationships for selected model
    model_nl <- rbind(mod$summary.random$`inla.group(temp)` %>% 
                     mutate(Variable = "Temperature",
                            Parasite = parasite) %>%
                     dplyr::rename(lci = `0.025quant`,
                                   uci = `0.975quant`,
                                   value = ID) %>%
                     dplyr::select(Parasite, Variable, value, mean, uci, lci),
                   mod$summary.random$`inla.group(prcp)` %>% 
                     mutate(Variable = "Precipitation",
                            Parasite = parasite) %>%
                     dplyr::rename(lci = `0.025quant`,
                                   uci = `0.975quant`,
                                   value = ID) %>%
                     dplyr::select(Parasite, Variable, value, mean, uci, lci))
    
    ## Monthly random effects for selected model
    model_t1 <- mod$summary.random$t1 %>%
                mutate(Model     = final_model,
                       Formula   = model_formula,
                       Parasite  = parasite) %>%
                dplyr::rename(Month = ID,
                              lci   = `0.025quant`,
                              uci   = `0.975quant`) %>%
                dplyr::select(Month, mean, lci, uci, Parasite)
    
    ## Interannual random effects for selected model
    model_t2 <- mod$summary.random$t2 %>%
      mutate(Model     = final_model,
             Formula   = model_formula,
             Parasite  = parasite) %>%
      dplyr::rename(Year = ID,
                    lci   = `0.025quant`,
                    uci   = `0.975quant`) %>%
      dplyr::select(Year, mean, lci, uci, Parasite)
    
    ## Spatial random effects for selected model
    model_s <- mod$summary.random$s1[1:46,] %>%
      mutate(Model     = final_model,
             Formula   = model_formula,
             Parasite  = parasite,
             Municipio = bolivar$MUNICIPIO,
             Parroquia = bolivar$PARROQUIA) %>%
      dplyr::rename(lci   = `0.025quant`,
                    uci   = `0.975quant`,) %>%
      dplyr::select(Parroquia, Municipio, mean, lci, uci, Parasite)
    
    ## Fitted values for selected model
    model_fitted <- mod$summary.fitted.values %>%
      dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
      mutate(Model     = final_model,
             Formula   = model_formula,
             Parasite  = parasite,
             Year      = t2,
             Month     = t1,
             Parroquia = Parroquia,
             Municipio = Municipio,
             Population = (e*12)*1000, 
             cases      = y,
             cpo        = mod$cpo$cpo) %>% 
      dplyr::rename(lci = `0.025quant`,
                    uci = `0.975quant`)
    
    
    ## Posterior predictive distribution for selected model
    s = 1000
    
    samples = inla.posterior.sample(s, mod)
    samples.m <- sapply(samples, function(x) x$latent[1:n])
    samples.m <- exp(samples.m)
    kappa <- as.vector(sapply(samples, function(x) x$hyperpar[1]))
    pred.samples=matrix(NA,nrow=dim(samples.m)[1],ncol=s)
    
    for (l in 1:dim(samples.m)[1])
    {
      pred.samples[l,]=rnegbin(s,samples.m[l,],kappa)
    }
    
    model_predictions <- data.frame(Model      = final_model,
                                 Formula    = model_formula, 
                                 Parasite   = parasite,
                                 Year       = t2,
                                 Month      = t1,
                                 Parroquia  = Parroquia,
                                 Municipio  = Municipio,
                                 Population = (e*12)*1000, 
                                 cases      = y,
                                 pred       = apply(pred.samples,1,mean),
                                 lci        = apply(pred.samples,1,quantile,probs=c(0.025),na.rm=T),
                                 uci        = apply(pred.samples,1,quantile,probs=c(0.975),na.rm=T))
  
    model_results <- list("dic" = model_dic, 
                          "estimates" = model_estimates, 
                          "non-linear climate" = model_nl,
                          "fitted" = model_fitted, 
                          "post_predictions" = model_predictions,
                          "monthly_random_effects" = model_t1,
                          "interannual_random_effects" = model_t2,
                          "spatial_random_effects" = model_s)
  
    file_path <- paste0("models/model_results_", mod_filename, ".RData")
    save(model_results, file = file_path)
   
}

############################################################################
## Run models for each parasite

run_models("P. falciparum", 8)

run_models("P. vivax", 8)
