# For this to run, you need the 2016 Cropland Data Layer (CDL), a vector shapefile of US state boundaries, and a pesticide reclass table
# CDL available here: https://www.nass.usda.gov/Research_and_Science/Cropland/Release/index.php
# shapefile available here: https://www.igismap.com/united-states-shapefile-download-free-map-boundary-states-and-county/
# reclass table available by contacting Dr. Maggie Douglas ()

library(sf)
library(raster)
library(logger)
library(tidyr)

# bring in cropland data layer (CDL)

cdl_2016 <- raster("E:\\FireflyAnalysis_October2020\\CDL_rasters\\2016_30m_cdls\\2016_30m_cdls.img") # read in CDL
plot(cdl_2016, main="CDL 2016") # map CDL

# bring in shapefile of states

us_states1 <- st_read("E:\\FireflyAnalysis_October2020\\states\\us_states.shp") #  read shapefile of US states
us_states1 <- sf::st_transform(us_states1, crs = crs(cdl_2016)) # reproject to match CDL

# overlay maps

plot(us_states1$geometry, add = TRUE) # they look good!

# extract North Carolina polygon

NorthCarolinaBoundary <- subset(us_states1, name == "North Carolina") # extract state boundary polygon
plot(NorthCarolinaBoundary$geometry) # nice!

# use NC polygon to extract CDL for NC

CroppedCDL <- crop(cdl_2016, NorthCarolinaBoundary) # crop CDL using NC boundary polygon
MaskedCDL <- raster::mask(CroppedCDL, NorthCarolinaBoundary) # mask cropped CDL using NC bounds
plot(MaskedCDL, main="CDL, North Carolina 2016") # this takes a minute or two
writeRaster(MaskedCDL, "E:\\2021_NC_Bees\\NorthCarolinaBees2021\\NC_CDL2016_Raster.tif") # takes a minute or two

# prepare pesticide reclass table for NC 2016

ReclassTable1 <- read.csv("E:\\FireflyAnalysis_October2020\\insecticide_reclass_table\\beetox_I_cdl_reclass_20200717.csv") # read in reclass table
Subset_ReclassTable <- subset(ReclassTable1, state_alpha == "NC" & year == "2016") # subset table
write.csv(Subset_ReclassTable, "E:\\2021_NC_Bees\\NorthCarolinaBees2021\\NC2016_reclassTable.csv") # export (necessary for SpeedyBeeModel)

# run SpeedyBeeModel on NC raster

SpeedyBeeModel::insecticide_index(
  output_dir = "E:\\2021_NC_Bees\\NorthCarolinaBees2021\\NC2016_PestMap", # output loc
  pesticide_path = "E:\\2021_NC_Bees\\NorthCarolinaBees2021\\NC2016_reclassTable.csv", # path to FILTERED reclass table
  landcover_path = "E:\\2021_NC_Bees\\NorthCarolinaBees2021\\NC_CDL2016_Raster.tif", # path to NC CDL map
  forage_range = 2000, # range in meters
  guild_table = NA, # No define guild table - not running multi-species stuff
  ins_method = "mean", # mean is default; could do oral/contact or whatev
  agg_factor = NA,
  normalize = F,
  check_pesttable = F # leave useW, and rastertag as default
)

# Did it work??

NC2016map <- raster("E:\\2021_NC_Bees\\NorthCarolinaBees2021\\NC2016_PestMap\\NC_CDL2016_Raster_insecticide.tif")
plot(NC2016map) # hell yea it did!