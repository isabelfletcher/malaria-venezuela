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


j <- "P. falciparum"

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

save(mod, file = "variable_testing/mod_mining.RData")

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

save(mod_w_mining, file = "variable_testing/mod_w_mining.RData")


sp_df <- rbind(mod$summary.random$s1[1:46,] %>%
        dplyr::select(mean, `0.025quant`, `0.5quant`, `0.975quant`) %>%
        mutate(Parasite = j,
               Model      = "With mining",
               Parroquia  = bolivar$PARROQUIA,
               Municipio  = bolivar$MUNICIPIO) %>%
        dplyr::rename(lci = `0.025quant`,
                      uci = `0.975quant`,
                      mci = `0.5quant`) %>%
        tibble::remove_rownames(),
      mod_w_mining$summary.random$s1[1:46,] %>%
        dplyr::select(mean, `0.025quant`, `0.5quant`, `0.975quant`) %>%
        mutate(Parasite = j,
               Model      = "Without mining",
               Parroquia  = bolivar$PARROQUIA,
               Municipio  = bolivar$MUNICIPIO) %>%
        dplyr::rename(lci = `0.025quant`,
                      uci = `0.975quant`,
                      mci = `0.5quant`) %>%
        tibble::remove_rownames())

# Rename Dalla Costa parroquias (two municipalities)
sp_df$Parroquia[sp_df$Municipio == "Sifontes" & sp_df$Parroquia == "Dalla Costa"] <- "Dalla Costa*"
sp_df$Parroquia[sp_df$Municipio == "Caroní" & sp_df$Parroquia == "Dalla Costa"] <- "Dalla Costa†"

ggplot(sp_df, aes(Parroquia)) + 
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_linerange(aes(ymin = lci, ymax = uci, group = Model, colour = Model), size = 0.6, position = position_dodge(width = 0.5)) +
  geom_point(aes(Parroquia, mean, group = Model, colour = Model), position = position_dodge(width = 0.5)) +
  theme_classic() + 
  ylab("Relative risk (log API)") + xlab("Parroquia") +
  facet_wrap(~Parasite, scales = "free") +
  theme(strip.text = element_text(face = "italic"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(colour = "black", fill = NA, size=0.5)) +
  scale_x_discrete(labels = wrap_format(10)) +
  xlim(unique(sp_df$Parroquia)) +
  scale_fill_manual(values = c("#6872A1", "#AFBCE8")) +
  scale_colour_manual(values = c("#6872A1", "#AFBCE8")) +
  coord_flip()
