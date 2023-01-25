## ----setup, include=FALSE-----------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

library(tidyverse)
library(here)
library(terra)
library(sf)
library(tmap)
library(geos)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyterra)

#File path to the HWC_data folder

HWC_data <- "/Volumes/GoogleDrive/.shortcut-targets-by-id/1YB-Hz3L-kWyiZMg2UM89GQkvqXyZUW1H/HWC_data"


## ----eval = FALSE-------------------------------------------------------------------------------------------------------
## #year <- 2030
## 
## #ssp <- 1


## -----------------------------------------------------------------------------------------------------------------------
popfolder <- paste0("FPOP_SSP", ssp)

popfile <- paste0("FPOP_SSP", ssp, "_" , year, ".tif")

pop <- rast(here(HWC_data, "/Geospatial Data/Pop_dens/fpop_data/", popfolder, popfile))


## -----------------------------------------------------------------------------------------------------------------------
#turn off spherical geometry to simplify joins, etc.
sf_use_s2(FALSE)

#read in asia lulc for reprojection
asia <- rast(here(HWC_data, "/Geospatial Data/Diminin_Replication_Data/cropped_lulc/lulc_asia.tif")) 

#reproject the designated population raster
popdens <- project(pop, asia)
  
#mask pop density to asia extent

popdens_asia <- mask(popdens, asia)
  
#aggregate pop density
asia_popdens_agg <- terra::aggregate(popdens_asia, fact = 10, fun = sum, na.rm = TRUE)
  
#cutoff value for highest decile (from baseline)

topdec <- 48851.86875
  
#reclassify so that the upper decile is 1 and all other values are NA
  
top10_pop <- asia_popdens_agg
  
top10_pop[top10_pop < topdec] <- NA
  
top10_pop[top10_pop >= topdec] <- 1


## ----eval = FALSE-------------------------------------------------------------------------------------------------------
## #year <- 2030
## 
## #ssp <- 1
## 
## #rcp <- 26 #2.6 - no decimals in the file path


## -----------------------------------------------------------------------------------------------------------------------

cropfolder <- paste0("SSP", ssp, "_", "RCP", rcp)

cropfile <- paste0("global_", "SSP", ssp, "_", "RCP", rcp, "_", year, ".tif")

crop <- rast(here(HWC_data, "/Geospatial Data/Chen_LULC_data", cropfolder, cropfile)) %>% 
  project(asia)


## -----------------------------------------------------------------------------------------------------------------------
#turn off spherical geometry to simplify joins, etc.
sf_use_s2(FALSE)

#mask lulc to asia extent

lulc_asia <- mask(crop, asia)

#filter just for cropland by reclassifying (cropland is 1, everything else is 0)

asia_cropland <- lulc_asia

asia_cropland[asia_cropland != 5] <- 0

asia_cropland[asia_cropland == 5] <- 1

#aggregate cropland

asia_cropland_agg <- terra::aggregate(asia_cropland, fact = 10, fun = mean, na.rm = TRUE)

#reclassify aggregated agriculture into 3 categories

# keep only the top decile (as in Di Minin et al.) - cutoff from baseline is 1 for 90th percentile

top10_crop <- asia_cropland_agg

top10_crop[top10_crop < 1] <- NA

top10_crop[top10_crop >= 1] <- 1


## -----------------------------------------------------------------------------------------------------------------------
#read in the data, set NA values to 0

human_dens <- top10_pop

human_dens[is.na(human_dens)] <- 0 

crop_dens <- top10_crop

crop_dens[is.na(crop_dens)] <- 0 

#raster math to create ranked conflict layer

ranked_conflict <- (human_dens + crop_dens) %>% 
  project(asia)

#re-mask to Asia

ranked_conflict_masked <- mask(ranked_conflict, asia)


## -----------------------------------------------------------------------------------------------------------------------
crs <- crs(asia)

#read in the extended range map and reproject to the correct crs

ext_range <- read_sf(here(HWC_data, "/Geospatial Data/Diminin_Replication_Data/extended_range/extended_range_asia.shp")) %>% 
  st_transform(crs = crs)

#buffer the extended range map by 10 km

ext_range_buffered <- st_buffer(ext_range, 10000)

#intersect the pressure layer with extended range map

#convert range map to a spatvector and then a spatraster

range_buff_rast <- vect(ext_range_buffered) %>% 
  rasterize(asia)

#cut out internal polygons of protected areas so that only the buffer remains

ext_range_rast <- vect(ext_range) %>% 
  rasterize(asia)

buffers_only <- mask(range_buff_rast, ext_range_rast, inverse = TRUE)

#read in the ranked pressures map

ranked_press <- ranked_conflict_masked %>% 
  project(asia, method = "near")

conflict_int <- mask(ranked_press, buffers_only)

#create a basemap

basemap_asia <- ne_countries(
  continent = "Asia",
  scale = "medium", 
  returnclass = "sf") %>% 
  st_transform(crs = crs(ranked_press))  # make sure crs is same as raster

basemap_asia_vect <- basemap_asia %>% 
  vect()

#plot conflict boundaries

pal <- c("darkgreen", "goldenrod", "orangered")

plot(basemap_asia_vect, col = "grey98")
plot(conflict_int, add = TRUE, pal = pal)


## -----------------------------------------------------------------------------------------------------------------------

name <- paste0("SSP", ssp, "_" , "RCP", rcp, "_", year, ".tif")

writeRaster(conflict_int, filename = here(HWC_data, "/Geospatial Data/Diminin_Replication_Data/conflict_boundaries_rasters_asia", name), overwrite = TRUE)

