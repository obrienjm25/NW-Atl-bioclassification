# Project Description and Code Overview ====
#
# Project:      Bioclassification
# Description:  Identify and delineate predominant assemblages of demersal fish
#               and benthic invertebrates in 4 marine bioregions in Atlantic 
#               Canada: 1) Northern Gulf of St. Lawrence, 2) Southern Gulf of St.
#               Lawrence, 3) Newfoundland & Labrador, 4) Maritimes

# Contact:      e-mail: John.Obrien@dfo-mpo.gc.ca | tel: +1.782.640.1522
# Publication:  O'Brien JM, Stanley RRE, Jeffery NW, Heaslip SG, DiBacco C, Wang Z
#               (2021) Modelling demersal fish and benthic invertebrate assemblages 
#               in support of marine conservation planning. Ecol Appl
#
# Code Overview:
# 1. Generate polygons for Study Area Boundaries, Fishnet Grids, Coastline
# 2. Clean up and subsetting of bottom trawl survey data
# 3. Cast to wide format
# 4. Aggregate taxa observations to 4-km grid cells
# 5. Remove rare taxa and sites with only a single taxon recorded
#    Convert taxon abundance to PA
#    Construct Site X Species Matrix
# 6. Figure Prep
#
# Requirements:
# R version 3.6.3
# Position and catch composition data from Department of Fisheries and Oceans Canada 
# annual multispecies bottom trawl surveys from 2007 to 2017

# House Keeping ====

# Load required packages

library(dplyr)
library(tidyr)
library(sf)
library(purrr)
library(data.table)
library(raster)
library(maptools)
library(ggplot2)
library(wesanderson)

# Source functions

source("Code/BioclassificationFunctions.R")

# Set controls

# String identifying the regional study area (i.e. NGSL, SGSL, NL, or MAR)
region <- "MAR"

# Habitat associations of taxa to be included
BottomHabitats <- c("bathydemersal", "benthic", "benthopelagic", 
                    "demersal", "sessile", "reef-associated", 
                    "bathypelagic")

# Combination of variables in survey datasets that uniquely define a trawl set
IDvars <- c("year", "month", "day", "set")

# Coordinate reference system of desired projection (integer EPSG code)
crs_proj <- ifelse(region == "NL", 26921, 26920)

# 1. Prepare polygons for study area boundaries and land borders ====

# Get polygons for countries, provinces & states within/bordering Atlantic Canada
# from database of global administrative boundaries

CAN <- raster::getData("GADM", country = "CAN", level = 1) # Canada
USA <- raster::getData("GADM", country = "USA", level = 1) # United States
SPM <- raster::getData("GADM", country = "SPM", level = 1) # St. Pierre & Miquelon

# Subset and bind polygons

AtlCan <- bind(CAN, USA) %>% # join Canada and USA polygons
  st_as_sf() %>% # convert to sf object
  # subset provinces and states in Atlantic region
  filter(NAME_1 %in% c("Maine", "Newfoundland and Labrador", "Prince Edward Island", 
                       "Qu\u{e9}bec", "Nova Scotia", "New Brunswick")) %>% 
  bind_rows(., st_as_sf(SPM)) %>% # join with St. Pierre Miquelon
  st_transform(4326) # transform to WGS84

# In each region: a) Aggregate survey strata to create boundaries of study area
#                 b) Crop land polygons to buffered extent around region
#                 c) Create fishnet grid covering the study area polygon

## NGSL ----

load("Data/Shapefiles/NGSL strata borders.rda")
NGSL_RV <- PolySet2SpatialPolygons(strat.N) %>% 
  st_as_sf() %>% 
  st_transform(26920) %>% 
  mutate(ID = 1:nrow(.)) %>% 
  filter(!(ID %in% c(1:5))) %>% 
  st_union() 

st_write(obj = NGSL_RV,
         dsn = "Data/Shapefiles",
         layer = "StudyArea_NGSL",
         driver = "ESRI Shapefile")

# b) Crop land polygons

NGSL_Land <- AtlCan %>% 
  st_transform(26920) %>% # transform to UTM20
  st_crop(., st_buffer(st_as_sfc(st_bbox(NGSL_RV)), 50000))

st_write(obj = NGSL_Land,
         dsn = "Data/Shapefiles",
         layer = "LandBorders_NGSL", 
         driver = "ESRI Shapefile")

# c) Create fishnet grid

SGSL_Grid <- GridFilter(regionshape = NGSL_RV, 
                        resol = 4000, 
                        landshape = NGSL_Land)

st_write(obj = NGSL_Grid,
         dsn = "Data/Shapefiles",
         layer = "Grid4km_NGSL",
         driver = "ESRI Shapefile")

## SGSL ----

# a) Aggregate survey strata

load("Data/Shapefiles/SGSL strata borders.rda")
SGSL_RV <- PolySet2SpatialPolygons(strat.S) %>% 
  st_as_sf() %>% 
  st_transform(26920) %>% 
  mutate(ID = 1:nrow(.)) %>% 
  st_union() 

st_write(obj = SGSL_RV,
         dsn = "Data/Shapefiles",
         layer = "StudyArea_SGSL",
         driver = "ESRI Shapefile")

# b) Crop land polygons

SGSL_Land <- AtlCan %>% 
  st_transform(26920) %>% # transform to UTM20
  st_crop(., st_buffer(st_as_sfc(st_bbox(SGSL_RV)), 50000))

st_write(obj = SGSL_Land,
         dsn = "Data/Shapefiles",
         layer = "LandBorders_SGSL", 
         driver = "ESRI Shapefile")

# c) Create fishnet grid

SGSL_Grid <- GridFilter(regionshape = SGSL_RV, 
                      resol = 4000, 
                      landshape = SGSL_Land)

st_write(obj = SGSL_Grid,
         dsn = "Data/Shapefiles",
         layer = "Grid4km_SGSL",
         driver = "ESRI Shapefile")

## NL ----

# a) Aggregate survey strata

NAFO <- st_read("Data/Shapefiles/NAFO_Divisions.shp") %>% # NAFO Divisions
  st_set_crs(., 4326) %>%  # set crs to WGS84
  st_transform(26921) %>% # transform to UTM21
  filter(., !(ZONE %in% c("2G","3M","0B"))) # exclude divisions unsampled by RV survey

NL_RV <- st_read("Data/Shapefiles/NF_SamplingStrata_20140514.shp") %>% # NL Survey strata
  st_transform(26921) %>% # transform to UTM21
  st_intersection(., NAFO) %>% # intersect with NAFO divisions
  st_union() # aggregate

st_write(obj = NL_RV,
         dsn = "Data/Shapefiles",
         layer = "StudyArea_NL",
         driver = "ESRI Shapefile")

# b) Crop land polygons

NL_Land <- AtlCan %>% 
  st_transform(26921) %>% # transform to UTM21
  st_crop(., st_buffer(st_as_sfc(st_bbox(NL_RV)), 50000))

st_write(obj = NL_Land,
         dsn = "Data/Shapefiles",
         layer = "LandBorders_NL", 
         driver = "ESRI Shapefile")


# c) Create fishnet grid

NL_Grid <- GridFilter(regionshape = NL_RV, 
                      resol = 4000, 
                      landshape = st_buffer(NL_Land, 5000))

st_write(obj = NL_Grid,
         dsn = "Data/Shapefiles",
         layer = "Grid4km_NL",
         driver = "ESRI Shapefile")

## MAR ----

# a) Aggregate survey strata

MAR_Reg <- st_read("Data/Shapefiles/MaritimesPlanningArea.shp") %>% # Maritimes Planning Region
  st_transform(26920) %>%  # transform to UTM20
  as(., "Spatial")

MAR_RV <- st_read("Data/Shapefiles/MaritimesRegionStrataBoundaries.shp") %>% # Maritimes Survey Strata
  st_transform(26920) %>% # transform to UTM20
  as(., "Spatial") %>% 
  cover(MAR_Reg, .) %>% 
  aggregate() %>% 
  st_as_sf() 
 
st_write(obj = MAR_RV,
         dsn = "Data/Shapefiles",
         layer = "StudyArea_MAR",
         driver = "ESRI Shapefile")

# b) Crop land polygons

MAR_Land <- st_crop(AtlCan, c(xmin = -71, ymin = 39, xmax = -45, ymax = 59)) %>% 
  st_transform(26920) # transform to UTM20

st_write(obj = MAR_Land,
         dsn = "Data/Shapefiles",
         layer = "LandBorders_MAR", 
         driver = "ESRI Shapefile")

# c) Create fishnet grid

MAR_Grid <- GridFilter(regionshape = MAR_RV, resol = 4000, landshape = MAR_Land)

st_write(obj = MAR_Grid,
         dsn = "Data/Shapefiles",
         layer = "Grid4km_MAR",
         driver = "ESRI Shapefile")

# 2. Select and Subset Trawl Data Sources ====

trawl_data <- readRDS(file = paste0("Data/TrawlData_", region, ".rds")) %>% 
  filter(vertical.position %in% BottomHabitats) %>%
  filter(year >= 2007) %>% 
  filter(!grepl("atlantic herring", common_name, ignore.case = TRUE)) %>%
  rowwise() %>% 
  mutate(SetID = paste(c_across(all_of(IDvars)), collapse = "_")) %>% 
  ungroup()

cat("Starting with", nrow(trawl_data), "observations of",
    length(unique(trawl_data$common_name)), "taxa from",
    length(unique(trawl_data$SetID)), "unique trawl sets in", region, "region", sep = " ")

# 3. Cast to wide format and convert to sf object ====

# reshape and convert to simple features (point geometry)

trawl_wide <- pivot_wider(trawl_data, 
                          id_cols = SetID, 
                          names_from = species,
                          values_from = totwgt,
                          values_fn = sum,
                          values_fill = 0) %>% 
  inner_join(., trawl_data[!duplicated(trawl_data$SetID),], by = "SetID") %>% 
  dplyr::select(SetID, longitude, latitude, all_of(IDvars), unique(trawl_data$species)) %>% 
  st_as_sf(., coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(., crs = crs_proj)

# Calculate average distance to nearest neighbour

distMat <- st_distance(trawl_wide)
diag(distMat) <- NA
NearNeigh <- apply(distMat, 2, min, na.rm = TRUE)

cat("For trawl survey points in", region, "region, the average distance to nearest
    neighbout is ~", round(mean(NearNeigh)/1000, digits = 2), "km", sep = " ")

# write to shapefile

st_write(obj = trawl_wide,
         dsn = "Data/Shapefiles",
         layer = paste0("TrawlWide_pt_", region),
         driver = "ESRI Shapefile")

# 4. Aggregate Trawl Observations to 4-km grid cells ====

# spatial join of grid with trawl points

grid <- st_read(paste0("Data/Shapefiles/Grid4km_", region, ".shp")) %>% 
  st_set_crs(crs_proj)

joinedgrid <- st_join(grid, trawl_wide) %>%
  filter(!is.na(SetID))

# List of populated grid cells

populated <- unique(joinedgrid$GridID)

cat("There are", length(populated), "4-km grid cells in", region,
    "region populated with trawl survey data from", nrow(joinedgrid),
    "unique sets", sep = " ")

# Write shapefile with trawl set frequency in populated grid cells

TrawlFreq <- joinedgrid %>% 
  group_by(GridID) %>% 
  summarise(Freq = n())

st_write(obj = TrawlFreq,
         dsn = "Data/Shapefiles",
         layer = paste0("TrawlFreq_Subgrid_", region),
         driver = "ESRI Shapefile")

# Examine distribution of sampling effort among grid cells

# Frequency table of trawl sets by grid cell

TrawlFreq <- table(joinedgrid$GridID) %>%
  data.frame() %>% 
  `colnames<-`(c("Grid","Frequency"))

# Histogram with frequency distribution

hist(TrawlFreq$Frequency, 
     main = "Allocation of sets among grid cells", 
     xlab = "Sets per grid", 
     col = "grey90")

cat("In", region, "region:", "\n\n", 
    round((nrow(TrawlFreq[TrawlFreq$Frequency > 1, ])/nrow(TrawlFreq)) * 100),
    "% of grid cells with trawl data contained more than one set",
    "\n\n  mode =", raster::modal(TrawlFreq$Frequency),
    "\n\n  95th percentile =", quantile(TrawlFreq$Frequency, 0.95),
    "\n\n  max =", max(TrawlFreq$Frequency))

# 5. Construct Site X Species Matrix of taxa P-A ====

# List of taxa

taxa <- st_drop_geometry(joinedgrid) %>% 
  names() %>% 
  setdiff(., c("GridID", "SetID", IDvars))

# Site X Species matrix (0 = absent, 1 = present)

SiteXSpecies <- joinedgrid %>%
  # remove geometry column
  st_drop_geometry() %>% 
  # select GridID and taxa columns
  dplyr::select(GridID, all_of(taxa)) %>% 
  group_by(GridID) %>%
  # sum catch data where 2+ sets fall in same grid cell
  summarise(across(all_of(taxa), sum)) %>% 
  ungroup() %>% 
  # convert to P-A
  mutate(across(all_of(taxa), ~if_else(.x > 0, 1, 0)))

# Identify barren sites (i.e. those with 1 or fewer taxa recorded)

BarrenSites <- SiteXSpecies %>% 
  rowwise() %>% 
  # calculate taxa richness at each site
  mutate(TaxaCount = sum(c_across(all_of(taxa)))) %>% 
  ungroup() %>% 
  filter(TaxaCount < 2) %>% 
  pull(GridID)

# Identify rare taxa (i.e. those observed at fewer than 1% of sites)

dropThreshold <- (ceiling(nrow(SiteXSpecies)*0.01))

cat("Identifying rare taxa observed at fewer than 1 % of sites",
    "\n\nDrop threshold is", dropThreshold, "out of",
    nrow(SiteXSpecies), "sites")

rareTaxa <- SiteXSpecies %>% 
  dplyr::select(all_of(taxa)) %>% 
  colSums() %>% 
  Filter(function(x) {x < dropThreshold}, .) %>% 
  names()

# Trim the data

cat("Dropping", length(BarrenSites), "barren sites and", 
    length(rareTaxa), "rare taxa from Site x Species matrix")

SiteXSpecies_trim <- SiteXSpecies %>% 
  filter(!(GridID %in% BarrenSites)) %>% 
  dplyr::select(-all_of(rareTaxa)) %>% 
  data.frame() %>% 
  `rownames<-`(.$GridID) %>% 
  dplyr::select(-GridID)

cat("Trimmed Site X Species matrix for", region, "region with",
    ncol(SiteXSpecies_trim), "taxa from", nrow(SiteXSpecies_trim), "sites")

# Save the data

write.csv(SiteXSpecies_trim, 
          file = paste0("Data/ClusterData_", region, ".csv"), 
          row.names = TRUE)

# 6. Figure Prep ====

## Fig. 1 ----

# Study Domain Map

# Coastline polygon

coast_atl <- AtlCan %>% 
  st_crop(., c(xmin = -69.5, ymin = 39, xmax = -45, ymax = 59)) %>% 
  st_simplify(., dTolerance = 0.025)


# read in shapefiles with study domains

names <- c("MAR","NGSL","NL","SGSL")
StudyArea <- list.files(pattern = "^StudyArea.+.shp", 
                        path = "Data/Shapefiles", 
                        full.names = TRUE) %>% 
  map(., st_read) %>%
  map(., ~ st_transform(., 4326)) %>%
  set_names(., names) %>% 
  map2(., names, ~ mutate(.x, region = .y)) %>% 
  map(., ~dplyr::select(., region, geometry)) %>% 
  map(., ~st_cast(., "MULTIPOLYGON")) %>% 
  rbindlist() %>% 
  st_as_sf() %>% 
  arrange((region))

# points for trawl sets

RVdata <- readRDS("Data/RVdata_4regions.rds") %>% 
  st_as_sf(., coords = c("longitude","latitude"), crs = 4326) %>% 
  filter(., year > 2007) %>%  
  distinct(., region, year, geometry, .keep_all = TRUE) %>% 
  filter(., st_intersects(., st_union(StudyArea), sparse = FALSE)) %>% 
  mutate(., region2 = case_when(
    region == "GULF" ~ "SGSL",
    region == "MARITIME" ~ "MAR",
    region == "QUEBEC" ~ "NGSL",
    region == "NEWFOUNDLAND" ~ "NL"
  ))

# map study area and trawl locations

pal <- wes_palette("Darjeeling1", n = 4, type = "discrete") #colour palette
names(pal) = c("NGSL","NL","SGSL", "MAR") # name palette elements

study.dom.gg <- ggplot() +
  geom_sf(data = coast_atl, fill = "grey30", col = "grey20",lwd = 0.25) +
  geom_sf(data = StudyArea, aes(fill = region), col = NA, alpha = 0.3) +
  geom_sf(data = RVdata, aes(col = region2), size = 0.25, shape = 16) + 
  scale_color_manual(values = pal, 
                     breaks = c("NGSL","SGSL","NL","MAR"),
                     aesthetics = c("fill","colour"), 
                     guide = guide_legend("Region", 
                                          override.aes = list(alpha = 0.3),
                                          title.theme = element_text(face = "bold", size = 7))) +
  coord_sf(xlim = c(-69.5,-45.5), 
           ylim = c(40.5,58),
           expand = FALSE) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        plot.margin = unit(c(0,0.5,0,0), "lines"),
        panel.grid.major = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(7, "points"),
        legend.position = c(0.85, 0.85),
        legend.background = element_blank()); study.dom.gg

ggsave(plot = study.dom.gg, 
       "Output/StudyDomain_600dpi.tiff", 
       width = 3, 
       unit = "in", 
       compression = "lzw",
       dpi = 600)

cat("Proceed with 2ClusterAnalysis.R")
