# For this to run, you need the 2016 Cropland Data Layer (CDL), a vector shapefile of US state boundaries, and a pesticide reclass table
# CDL available here: https://www.nass.usda.gov/Research_and_Science/Cropland/Release/index.php
# shapefile available here: https://www.igismap.com/united-states-shapefile-download-free-map-boundary-states-and-county/
# reclass table available by contacting Dr. Maggie Douglas ()

# install.packages("devtools")
# install.packages("raster")
library(devtools)
devtools::install_github("land-4-bees/SpeedyBeeModel")

library(sf); library(raster); library(logger); library(tidyr); library(SpeedyBeeModel)
rm(list=ls())
# bring in cropland data layer (CDL)

cdl_2020 <- raster("E:\\2021_NC_Bees\\2020_30m_cdls\\2020_30m_cdls.img") # read in CDL
#cdl_2020 <- raster::raster('Y:/SpatialData/NASS_CDL/CDL2020/2020_30m_cdls.img') #file path on Melanie's computer
plot(cdl_2020, main="CDL 2020") # map CDL

# bring in shapefile of states (use an R package to grab this shapefile to avoid having to change the file path)
us_states1 <- tigris::states() %>% 
  sf::st_transform(us_states1, crs = crs(cdl_2020))

# us_states1 <- st_read("E:\\FireflyAnalysis_October2020\\states\\us_states.shp") #  read shapefile of US states
# us_states1 <- sf::st_transform(us_states1, crs = crs(cdl_2020)) # reproject to match CDL

# overlay maps

plot(us_states1$geometry, add = TRUE) # they look good!

# extract Pennsylvania polygon

PennsylvaniaBoundary <- subset(us_states1, NAME == "Pennsylvania") # extract state boundary polygon
plot(PennsylvaniaBoundary$geometry) # nice!

# use polygon to extract CDL for PA

CroppedCDL <- raster::crop(cdl_2020, PennsylvaniaBoundary) # crop CDL using PA boundary polygon
MaskedCDL <- raster::mask(CroppedCDL, PennsylvaniaBoundary) # mask cropped CDL using PA bounds
plot(MaskedCDL, main="CDL, North Carolina 2016") # this takes a minute or two
writeRaster(MaskedCDL, "./PA_CDL2020_Raster.tif", overwrite = TRUE) # takes a minute or two
# unique(MaskedCDL)
# dev.off()
# hist(MaskedCDL)

# prepare pesticide reclass table for PA 2020

ReclassTable1 <- read.csv("E:\\FireflyAnalysis_October2020\\insecticide_reclass_table\\beetox_I_cdl_reclass_20200717.csv") # read in reclass table
Subset_ReclassTable <- subset(ReclassTable1, state_alpha == "PA" & year == "2014") # subset table
write.csv(Subset_ReclassTable, "E:\\2021_NC_Bees\\NorthCarolinaBees2021\\PA2020_reclassTable.csv") # export (necessary for SpeedyBeeModel)

##################################### [PESTICIDES]

# run SpeedyBeeModel on NC raster

SpeedyBeeModel::insecticide_index(
  output_dir = "E:\\2021_NC_Bees\\NorthCarolinaBees2021\\PA2020_PestMap", # output loc
  pesticide_path = "E:\\2021_NC_Bees\\NorthCarolinaBees2021\\PA2020_reclassTable.csv", # path to FILTERED reclass table
  landcover_path = "E:\\2021_NC_Bees\\NorthCarolinaBees2021\\PA_CDL2020_Raster.tif", # path to NC CDL map
  forage_range = 2000, # range in meters
  guild_table = NA, # No define guild table - not running multi-species stuff
  ins_method = "mean", # mean is default; could do oral/contact or whatev
  agg_factor = NA,
  normalize = F,
  check_pesttable = F # leave useW, and rastertag as default
)

# Did it work??

PA2020map <- raster("E:\\2021_NC_Bees\\NorthCarolinaBees2021\\PA2020_PestMap\\PA_CDL2020_Raster_insecticide.tif")
plot(PA2020map) # hell yea it did!

##################################### [LONSDORF INDICES]

ForageTable1 <- read.csv("./cdl_reclass_forage_quality_updated.csv")
# View(ForageTable1)

SpeedyBeeModel::forage_index( 
    output_dir = "./PA2020_LonsdorfIndices", # output loc
    landcover_path  = "./PA_CDL2020_Raster.tif", # path to PA CDL map
    forage_table = ForageTable1,
    seasons = c("Floral_Spring", "Nest_Cavity"),
    forage_range = 2000,
    guild_table = NA,
    agg_factor = NA,
    normalize = T,
    compress_rasters = T
)

PAfloral <- raster("./PA2020_LonsdorfIndices/PA_CDL2020_Raster_Floral_Spring.tif")
plot(PAfloral)

PAnest <- raster("./PA2020_LonsdorfIndices/PA_CDL2020_Raster_Nest_Cavity.tif")
plot(PAnest)

PAfloral2 <- raster("./PA2020_LonsdorfIndices/PA_CDL2020_Raster_Floral_Summer.tif")
plot(PAfloral2)

# getting an error about missing raster classes
# this is text from mel to figure out what's missing

forage_table <- ForageTable1
hab.r <- raster::raster("E:\\2021_NC_Bees\\NorthCarolinaBees2021\\NC_CDL2016_Raster.tif")
same <- unique(raster::values(hab.r)) %in% forage_table$LULC
missing <- unique(raster::values(hab.r))[!same]
missing <- missing[!is.na(missing)]
missing

