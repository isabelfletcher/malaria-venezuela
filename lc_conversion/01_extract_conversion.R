library(raster)
library(rgdal)
library(dplyr)

venezuela <- rgdal::readOGR("ven_shp", "Parroquias_Venezuela")
bolivar <- subset(venezuela, venezuela$ESTADO == "Bolívar")
files <- list.files(path = "input/bolivar", pattern = ".tif", full.names = TRUE)

# Function to extract land cover conversions for a polygon
extract_gradient <- function(parroquia, municipio) {
  
  polygon <- subset(bolivar, bolivar$PARROQUIA == parroquia & bolivar$MUNICIPIO == municipio)
  
  # Create empty df
  data <- NULL
  
  # Loop through each file, revalue and extract cell values
  for (i in files){
    
    # Revalue raster by creating classification matrix
    reclass_mat <- matrix(c(0, 10, NA,
                            10, 40, 1,   # Cropland
                            40, 120, 2,  # Tree cover
                            120, 160, 3, # Grass-/Shrubland
                            160, 190, 4, # Flooded vegetation
                            190, 200, 7, # Urban
                            200, 210, 6, # Bare 
                            210, 220, 5, # Water bodies
                            220, Inf, NA), ncol = 3, byrow = TRUE)
    
    # Read in raster, crop and revalue
    ras <- raster(i) %>% crop(extent(polygon)) %>% raster::mask(polygon) %>% raster::reclassify(reclass_mat)
    
    # Loop through all cells in polygon and get cell values
    # Condition on NA, ignoring NA values
    for (j in 1:ncell(ras)){
      
      if(!is.na(ras[j]) == TRUE){
        
        # Extract cell coordinates and matching polygon
        pts <- as.data.frame(xyFromCell(ras, j))
        pts <- SpatialPoints(pts, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
        
        # Extract cell value and combine into df
        df <- data.frame(Year = gsub("[^0-9.-]", "", i),
                         cell = j,
                         value = ras[j],
                         x = xFromCell(ras, j),
                         y = yFromCell(ras, j)) %>%
          mutate(Year = gsub("[.]", "", Year)) %>%
          # Extract over polygon
          cbind(over(pts, bolivar)) %>%
          dplyr::select(-ESTADO) %>%
          dplyr::rename(Parroquia = PARROQUIA,
                        Municipio = MUNICIPIO)
        data <- rbind(data, df)
        
      } else {
        # do nothing
      }
    } # end loop of cells
  } # end loop of files
  
  # Relabel values
  data$value <- as.factor(data$value)
  levels(data$value)[levels(data$value)=="1"] <- "Cropland"
  levels(data$value)[levels(data$value)=="2"] <- "Tree cover"
  levels(data$value)[levels(data$value)=="3"] <- "Grass-/Shrubland"
  levels(data$value)[levels(data$value)=="4"] <- "Flooded vegetation"
  levels(data$value)[levels(data$value)=="5"] <- "Water bodies"
  levels(data$value)[levels(data$value)=="7"] <- "Urban"
  levels(data$value)[levels(data$value)=="6"] <- "Bare"
  
  # Create empty df
  change_df <- NULL
  
  # Loop through all cells in df and identify conversion
  for (k in unique(data$cell)){
    
    df <- data %>% subset(cell == k)
    
    # If land cover type is different to previous year, create a variable of this change
    if(length(unique(df$value)) > 1){
      lc_df <- df %>% mutate(Conversion = paste(subset(df, df$Year == 1996)$value, "to", subset(df, df$Year == 2016)$value),
                             Change     = "Change")
    } else {
      # If land cover doesn't change label as 'no change'
      lc_df <- df %>% mutate(Conversion = paste(unique(df$value), "to", unique(df$value)),
                             Change = "No change")
    }
    change_df <- rbind(change_df, lc_df)
  }
  
  # Save data
  save(change_df, file = gsub(" ", "", paste0("output/lc_gradient_", tolower(gsub(" ", "", parroquia)), "_", tolower(gsub(" ", "", municipio)), ".RData")))
  
  # end function
}

#############################
# List polygons to extract for 
parroquias <- bolivar$PARROQUIA
municipios <- bolivar$MUNICIPIO

# Loop through all polygons
for (p in 1:length(parroquias)){
  
  extract_gradient(parroquias[p], municipios[p])
  
}

#############################
# Combine data into single file
files <- list.files("output", full.names = TRUE, pattern = ".RData")

load(files[1])
data <- change_df

for (i in 2:length(files)){
  load(files[i])
  data <- rbind(data, change_df)
}

data <- subset(data, data$Year == 2016) %>% dplyr::select(-Year)
save(data, file = "output/data.RData")
