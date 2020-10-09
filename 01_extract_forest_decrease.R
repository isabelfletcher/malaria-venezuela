library(raster)
library(dplyr)
library(rgdal)

venezuela <- rgdal::readOGR("ven_shp", "Parroquias_Venezuela")
bolivar <- subset(venezuela, venezuela$ESTADO == "Bolívar")

files <- list.files(path = "input/bolivar", pattern = ".tif", full.names = TRUE)

polygon <- subset(bolivar, bolivar$PARROQUIA == polygon_name & bolivar$MUNICIPIO == municipio_name)

######################################################
# 1. Create new rasters of forest change 1996-2016
revalue_ras <- function(file_1, file_2){
  
  # Read in file1
  ras_1 <- file_1 %>% raster() %>% crop(extent(bolivar)) %>% raster::mask(bolivar)
  
  # Extract all values labelled as forest, based on legend from ESA CCI
  ras_1 <- ras_1 == 40 | ras_1 == 50 |
    ras_1 == 60 | ras_1 == 70 | 
    ras_1 == 80 | ras_1 == 90 | 
    ras_1 == 100 | ras_1 == 110
  
  # Read in file2
  ras_2 <- file_2 %>% raster() %>% crop(extent(bolivar)) %>% raster::mask(bolivar)
  
  # Extract all values labeled as forest, based on legend from ESA CCI
  ras_2 <- ras_2 == 40 | ras_2 == 50 |
    ras_2 == 60 | ras_2 == 70 | 
    ras_2 == 80 | ras_2 == 90 | 
    ras_2 == 100 | ras_2 == 110
  
  # Loop through all cells in polygon, relabelling values 
  for (i in 1:ncell(ras_1)){
    
    # Only use values not NA
    if(is.na(ras_2[i]) == TRUE){
    } else if (ras_2[i] == ras_1[i]){
      # If value is the same, then revalue to 0 (no change)
      ras_2[i] <- 0 
      # If value is 1 then revalue to 1 (increase)
    } else if (ras_2[i] == 1) {
      ras_2[i] <- 1
    } else {
      # If else revalue to -1 (decrease)
      ras_2[i] <- -1
    }
  }
  # Save new raster to file
  writeRaster(ras_2, gsub("input/bolivar/", paste0("output/"), file_2), overwrite = TRUE)
}

# Loop through all files
nfiles <- length(files)-1
for (i in 1:nfiles){
    revalue_ras(files[i], files[i+1])
}

######################################################
# 2. Extract area where forest cover has decreased

files     <- list.files(path = "output", pattern = ".tif", full.names = TRUE)
ras_stack <- stack(files)

# Extract cells where forest has decreased (-1)
ras_stack <- ras_stack == -1

# Function to calculate annual total area where forest has decreased, for each polygon

# List polygons 
parroquias <- bolivar$PARROQUIA
municipios <- bolivar$MUNICIPIO

`%notin%` <- Negate(`%in%`)

# Create empty df to fill 
ras_df <- NULL

calc_area <- function(parroquia, municipio){
  
  polygon     <- subset(bolivar, bolivar$PARROQUIA == parroquia & bolivar$MUNICIPIO == municipio)
  ras_polygon <- crop(ras_stack, extent(polygon))
  
  # Loop through all files
  for (i in 1:nlayers(ras_polygon)) {
    
    ras <- ras_polygon[[i]]
    a   <- area(ras)
    
    # Ignore NA values
    if (all(is.na(c(values(ras)))) == TRUE){
      a_df <- data.frame(forest_decrease = 0,
                         Parroquia       = parroquia,
                         Municipio       = municipio,
                         Year            = gsub("X", "", names(ras)[1]),
                         zone            = 1)
    } else {
      # Calculate total area where forest has decreased
      a_df <- as.data.frame(zonal(a, ras, "sum")) %>%
        subset(zone == 1) %>%
        dplyr::rename(forest_decrease = sum) %>% 
        mutate(Parroquia = parroquia,
               Municipio = municipio,
               Year      = gsub("X", "", names(ras)[1]))
      # Include years with no forest cover decrease
      if (1 %notin% a_df$zone){
        a_df <- data.frame(forest_decrease = 0,
                           Parroquia       = parroquia,
                           Municipio       = municipio,
                           Year            = gsub("X", "", names(ras)[1]),
                           zone            = 1)
      } 
    }
    ras_df <- rbind(ras_df, a_df)
  }
  ras_df <- ras_df %>% dplyr::select(-zone)
  return(ras_df)
}

# Loop through all polygons and combine data
data <- NULL
for (i in 1:length(parroquias)){
  
  df <- calc_area(parroquias[i], municipios[i])
  
  data <- rbind(data, df)
}

######################################################
# 3. Calculate cumulative forest decrease and save data

data %>% mutate(forest_decrease_cum = cumsum(forest_decrease)) %>% write.csv("input/data.csv")
