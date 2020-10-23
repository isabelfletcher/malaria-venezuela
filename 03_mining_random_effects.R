## Load libraries
pacman::p_load("raster", "INLA", "sf", "sp", "spdep", "dplyr", "ggplot2",
               "scales", "knitr", "kableExtra", "doBy", "ggpubr",
               "RColorBrewer", "Hmisc", "cowplot", "bookdown", "raster",
               "lattice", "rasterVis", "brinla", "miceadds", "gridExtra",
               "ggpubr", "gganimate", "gifski", "transformr", "tidyr",
               "ggregplot", "grid", "MASS", "INLAutils", "reshape2", "zoo",
               "viridis", "ggcorrplot", "Metrics", "lubridate")

### Define neighbourhood matrix
setwd("/Volumes/Isabel Backup/PHID Dropbox/Isabel Fletcher/Isabel_PhD/projects/venezuela/spatio_temporal_model_era5_land")
venezuela <- rgdal::readOGR("data/ven_shp", "Parroquias_Venezuela")
bolivar <- subset(venezuela, venezuela$ESTADO == "Bolívar")

nb.map <- poly2nb(bolivar)
nb2INLA("map.graph", nb.map)

##### Read in data
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

################################################
# Empty df to fill data
model_dic       <- NULL
model_estimates <- NULL
model_fit       <- NULL
model_sp         <- NULL

# Loop through parasites
parasites <- c("P. falciparum", "P. vivax")

for (i in parasites){

# Condition on parasite
if (j == "P. falciparum"){
  y <- data$PF
  temp <- scale(data$temp_lag2)[,1]
  prcp <- scale(data$prcp_lag0)[,1]
  mod_filename <- "_pf"
} else {
  y <- data$PV
  mod_filename <- "_pv"
  temp <- scale(data$temp_lag2)[,1]
  prcp <- scale(data$prcp_lag0)[,1]
}

df_inla <- data.frame(y, e,
                      Parroquia, Municipio,
                      s1, t1, t2,
                      temp, prcp,
                      forest_decrease,
                      mines, urbanization, 
                      health_travel)

# Run model
formula <- y ~ 1 + f(s1, model = "bym", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(temp), model = "rw1") +
  f(inla.group(prcp), model = "rw1") +
  mines + 
  forest_decrease + 
  health_travel + 
  urbanization 

mod <- inla(formula, data = df_inla, family = "zeroinflatednbinomial0", 
            offset = log(e), verbose = TRUE,
            control.inla = list(strategy = 'adaptive'),
            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                   config = FALSE, 
                                   return.marginals = FALSE), 
            control.predictor = list(link = 1, compute = TRUE), 
            control.family = list(link = "log"))

save(mod, file = paste0("models/mod_mining", mod_filename, ".RData"))

# Run model
formula <- y ~ 1 + f(s1, model = "bym", graph = "map.graph") +
  f(t1, model = "rw1") +
  f(t2, model = "iid") +
  f(inla.group(temp), model = "rw1") +
  f(inla.group(prcp), model = "rw1") +
  #mines +
  forest_decrease + 
  health_travel + 
  urbanization 

mod_w_mining <- inla(formula, data = df_inla, family = "zeroinflatednbinomial0", 
                     offset = log(e), verbose = TRUE,
                     control.inla = list(strategy = 'adaptive'),
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, 
                                            config = FALSE, 
                                            return.marginals = FALSE), 
                     control.predictor = list(link = 1, compute = TRUE), 
                     control.family = list(link = "log"))

save(mod_w_mining, file = paste0("models/mod_w_mining", mod_filename, ".RData"))

}
