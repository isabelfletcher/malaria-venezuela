## Load libraries
library(raster)
library(INLA)
library(sf)
library(sp)
library(dplyr)
library(rgdal)
library(MASS)
library(Metrics)
library(hydroGOF)
library(INLAutils)
library(spdep)

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

# Climate variables
temp <- scale(data$temp_avg)[,1]
prcp <- scale(data$prcp_avg)[,1]
nino <- scale(data$nino34_lag8)[,1]

# Land use variables
mines           <- scale(data$total_mines)[,1]
forest_decrease <- scale(data$deforestation_lag0)[,1]
urbanization    <- scale(data$urban_increase_cum)[,1]

# Spatial random effects
s1 <- rep(1:46, 252)
Parroquia <- data$Parroquia
Municipio <- data$Municipio

# Monthly/seasonal random effects
t1 <- rep(rep(1:12, each = 46), 21)

# Yearly random effects
t2 <- rep(1996:2016, each = 552)

###############################################################
# Non-linear models

# Cut data with median mines 
data$mining_cat[data$total_mines <= 2] <- 1 # low mining
data$mining_cat[data$total_mines > 2]  <- 2 # high mining

mining_cat <- data$mining_cat

#### Loop through parasites
parasites <- c("P. falciparum", "P. vivax")

for (j in parasites){
  
  # Condition on parasite
  if (j == "P. falciparum"){
    y <- data$PF
    mod_filename <- "_pf"
  } else {
    y <- data$PV
    mod_filename <- "_pv"
  }
  
  df_inla <- data.frame(y, e,
                        Parroquia, Municipio,
                        s1, t1, t2,
                        temp, prcp, nino,
                        mines, forest_decrease,
                        mining_cat,
                        urbanization)
  
  formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
    f(t1, model = "rw1") +
    f(t2, model = "iid") +
    f(inla.group(temp), model = "rw1", replicate = mining_cat) +
    f(inla.group(prcp), model = "rw1", replicate = mining_cat) +
    nino +
    forest_decrease + 
    mines +
    urbanization 
  
  # Run model
  mod <- inla(formula, data = df_inla, family = "zeroinflatednbinomial0",
              offset = log(e), verbose = TRUE,
              control.inla = list(strategy = 'adaptive'),
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
                                     config = TRUE,
                                     return.marginals = TRUE),
              control.predictor = list(link = 1, compute = TRUE),
              control.family = list(link = "log"))
  
  # Save model
  save(mod, file = paste0("models/models_nonlinear/mod", mod_filename, ".RData"))
  
  # Model without mining
  formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
    f(t1, model = "rw1") +
    f(t2, model = "iid") +
    f(inla.group(temp), model = "rw1", replicate = mining_cat) +
    f(inla.group(prcp), model = "rw1", replicate = mining_cat) +
    nino +
    forest_decrease + 
    #mines +
    urbanization
  
  # Run model
  mod <- inla(formula, data = df_inla, family = "zeroinflatednbinomial0",
              offset = log(e), verbose = TRUE,
              control.inla = list(strategy = 'adaptive'),
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
                                     config = TRUE,
                                     return.marginals = TRUE),
              control.predictor = list(link = 1, compute = TRUE),
              control.family = list(link = "log"))
  
  # Save model
  save(mod, file = paste0("models/models_nonlinear/mod_w_mining", mod_filename, ".RData"))
  
  # Model without interaction
  formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
    f(t1, model = "rw1") +
    f(t2, model = "iid") +
    f(inla.group(temp), model = "rw1") +
    f(inla.group(prcp), model = "rw1") +
    nino +
    forest_decrease + 
    mines +
    urbanization 
  
  # Run model
  mod <- inla(formula, data = df_inla, family = "zeroinflatednbinomial0",
              offset = log(e), verbose = TRUE,
              control.inla = list(strategy = 'adaptive'),
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
                                     config = TRUE,
                                     return.marginals = TRUE),
              control.predictor = list(link = 1, compute = TRUE),
              control.family = list(link = "log"))
  
  # Save model
  save(mod, file = paste0("models/models_nonlinear/mod_w_interaction", mod_filename, ".RData"))
  
  # Model without interaction or mining
  formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
    f(t1, model = "rw1") +
    f(t2, model = "iid") +
    f(inla.group(temp), model = "rw1") +
    f(inla.group(prcp), model = "rw1") +
    nino +
    forest_decrease + 
    # mines +
    urbanization 
  
  # Run model
  mod <- inla(formula, data = df_inla, family = "zeroinflatednbinomial0",
              offset = log(e), verbose = TRUE,
              control.inla = list(strategy = 'adaptive'),
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
                                     config = TRUE,
                                     return.marginals = TRUE),
              control.predictor = list(link = 1, compute = TRUE),
              control.family = list(link = "log"))
  
  # Save model
  save(mod, file = paste0("models/models_nonlinear/mod_w_interaction_w_mining", mod_filename, ".RData"))
  
  
}

###############################################################
# Linear models
#### Loop through parasites
parasites <- c("P. falciparum", "P. vivax")

for (j in parasites){
  
  # Condition on parasite
  if (j == "P. falciparum"){
    y <- data$PF
    mod_filename <- "_pf"
  } else {
    y <- data$PV
    mod_filename <- "_pv"
  }
  
  df_inla <- data.frame(y, e,
                        Parroquia, Municipio,
                        s1, t1, t2,
                        temp, prcp, nino,
                        mines, 
                        forest_decrease,
                        urbanization)
  
  formula <- y ~ 1 + f(s1, model = "bym2", graph = "map.graph") +
    f(t1, model = "rw1") +
    f(t2, model = "iid") +
    temp +
    prcp +
    nino +
    mines +
    forest_decrease + 
    urbanization
  
  # Run model
  mod <- inla(formula, data = df_inla, family = "zeroinflatednbinomial0",
              offset = log(e), verbose = TRUE,
              control.inla = list(strategy = 'adaptive'),
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
                                     config = TRUE,
                                     return.marginals = TRUE),
              control.predictor = list(link = 1, compute = TRUE),
              control.family = list(link = "log"))
  
  # Save model
  save(mod, file = paste0("models/models_linear/mod", mod_filename, ".RData"))
  
}
