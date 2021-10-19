# Project and Code Overview ====
#
# Project:      Bioclassification
# Contact:      e-mail: John.Obrien@dfo-mpo.gc.ca | tel: +1.782.640.1522
# Publication:  O'Brien JM, Stanley RRE, Jeffery NW, Heaslip SG, DiBacco C, Wang Z
#               (2021) Modelling demersal fish and benthic invertebrate assemblages 
#               in support of marine conservation planning. Ecol Appl
#
# Overview:
# 
# 1. Evaluate support for combining spring and fall survey data
#    in NL by comparing community structure between seasons
# 2. Evaluate potential biases introduced by trawl data aggregation
#    and use of higher-level taxa (> genus)
# 3. Examine distribution of grid cells with higher effort
# 4. Explore correlation between indicator taxa richness and species
#    richness
# 5. Examine spatial distribution of minor assemblages
# 6. Figure Prep

# Requirements:
# R version 3.6.3
# Position and catch composition data from Department of Fisheries and Oceans Canada 
# annual multispecies bottom trawl surveys from 2007 to 2017
# Polygons for NAFO divisions
# Saved classification of grid cells from 2ClusterAnalysis.R
# Shapefiles for subset of populated grid cells from 1DataPrep.R giving trawl frequency
# Rasters with predicted assemblage distributions from 4RandomForestClassification.R
# csv files with grid cluster assignments from 2ClusterAnalysis.R
# txt files with indicator taxa from 3IndicatorAnalysis.R
# rds files with sf objects containing classified grid cells from 2ClusterAnalysis.R

# House Keeping ====

# Load required packages

library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(rlang)
library(sf)
library(raster)
library(rasterVis)
library(vegan)
library(simba)
library(ggplot2)
library(tmap)
library(cowplot)
library(wesanderson)

# Source functions

source("Code/BioclassificationFunctions.R")

# 1. Evaluate support for combining spring and fall survey data in NL ====

# Set controls

# String identifying the regional study area (i.e. NGSL, SGSL, NL, or MAR)
region <- "NL"

# Habitat associations of taxa to be included
BottomHabitats <- c("bathydemersal", "benthic", "benthopelagic", 
                    "demersal", "sessile", "reef-associated", 
                    "bathypelagic")

# Combination of variables in survey datasets that uniquely define a trawl set
IDvars <- c("year", "month", "day", "longitude", "latitude")

# Coordinate reference system of desired projection (integer EPSG code)
crs_proj <- ifelse(region == "NL", 26921, 26920)

# Select and Subset Trawl Data

trawl_data <- readRDS(file = paste0("Data/TrawlData_", region, ".rds")) %>% 
  filter(vertical.position %in% BottomHabitats) %>%
  filter(year >= 2007) %>% 
  filter(!grepl("atlantic herring", common_name, ignore.case = TRUE)) %>%
  rowwise() %>% 
  mutate(SetID = paste(c_across(all_of(IDvars)), collapse = "_")) %>% 
  ungroup()

# Split observations into Spring and Fall surveys

trawl_split <- split(trawl_data, trawl_data$season)

# reshape and convert to simple features (point geometry)

trawl_wide <- trawl_split %>% 
  map(~pivot_wider(.x, 
                   id_cols = SetID, 
                   names_from = species,
                   values_from = totwgt,
                   values_fn = sum,
                   values_fill = 0)) %>% 
  map(~inner_join(.x, trawl_data[!duplicated(trawl_data$SetID),], by = "SetID")) %>% 
  map(~dplyr::select(.x, SetID, longitude, latitude, all_of(IDvars), 
                     any_of(unique(trawl_data$species)))) %>% 
  map(~st_as_sf(.x, coords = c("longitude", "latitude"), crs = 4326)) %>% 
  map(~st_transform(.x, crs = crs_proj))

# Spatial subset to NAFO divisions surveyed in both seasons

NAFO <- st_read("Data/Shapefiles/NAFO_Divisions.shp") %>% 
  st_set_crs(4326) %>% 
  st_transform(crs_proj) %>% 
  filter(ZONE %in% c("3L", "3O", "3N"))

trawl_seasons <- trawl_wide %>% 
  map(., ~st_filter(.x, NAFO)) %>% 
  rbindlist(idcol = "Season", fill = TRUE) %>% 
  #replace NA values (species not caught in both seasons) with 0 
  replace_na(list(Cryptopsaras.couesii = 0, 
                  Petromyzon.marinus = 0, 
                  Ulvaria.subbifurcata = 0))

# Convert to Site X Species matrix with P-A

SiteXSpecies <- trawl_seasons %>%
  dplyr::select(any_of(unique(trawl_data$species))) %>% 
  mutate(across(everything(), ~if_else(.x > 0, 1, 0)))

# Visualize multivariate differences in composition between seasons with nMDS

# Initial start
mds_seasons <- metaMDS(SiteXSpecies, 
                       k = 2, 
                       distfun = sim, 
                       distance = "simpson", 
                       zerodist = "ignore")
# Start search from previous solution
mds_seasons2 <- metaMDS(SiteXSpecies, 
                        previous.best = mds_seasons, 
                        trymax = 100)
# Start search from previous solution
mds_seasons3 <- metaMDS(SiteXSpecies, 
                        previous.best = mds_seasons2,
                        trymax = 100,
                        sfgrmin = 1e-8)

# Save results
saveRDS(mds_seasons3, 'Output/NL_mds_k2.rds')

# Convert results to dataframe
df_mds <- data.frame(mds_seasons3[["points"]], 
                     season = trawl_seasons$Season) %>% 
  dplyr::rename(NMDS1 = MDS1, NMDS2 = MDS2)

# 2. Evaluate potential biases introduced by trawl data aggregation ====

# Proceed with 1DataPrep.R

# To remove higher taxa groups, modify L244-250 as follows:

# trawl_data <- readRDS(file = paste0("Data/TrawlData_", region, ".rds")) %>%
#   filter(vertical.position %in% BottomHabitats) %>%
#   filter(year >= 2007) %>%
#   filter(!grepl("atlantic herring", common_name, ignore.case = TRUE)) %>%
#   rowwise() %>%
#   mutate(SetID = paste(c_across(all_of(IDvars)), collapse = "_")) %>%
#   ungroup() %>% 
#   filter(!is.na(genus))
 
# To select most recent set from grid cells with multiple trawl observations,
# modify L294-295 as follows:

# joinedgrid <- st_join(grid, trawl_wide) %>%
#   filter(!is.na(SetID)) %>%
#   group_by(GridID) %>%
#   arrange(desc(year)) %>%
#   summarise(across(everything(), first)) %>% 
#   ungroup()

# Proceed with 2ClusterAnalysis.R to L266

# To save new classification with non-aggregated data, modify L266 as follows:

# saveRDS(colour_grid, file = paste0("Data/Grid_ClusterAssignment_NoAgg_", region, ".rds"))

# Check for correlation among original & new classification

# Read in classified grid cells
class_orig <- readRDS(paste0("Data/Grid_ClusterAssignment_", region, ".rds"))
class_new <- readRDS(paste0("Data/Grid_ClusterAssignment_NoAgg_", region, ".rds"))

# Ensure same grid cells are present in both classifications
all.equal(class_new$GridID, class_orig$GridID) # TRUE

# Create frequency table with counts assigned to each cluster

tbl_orig <- table(class_orig$name)

if(region == "NGSL"){
  tbl_new <- c(table(class_new$name), Unclassified = 0)
} else{
  tbl_new <- table(class_new$name)
} 

tbl <- matrix(data = c(tbl_orig, tbl_new), 
              nrow = 2, 
              ncol = length(unique(class_orig$name)), 
              byrow = TRUE)
dimnames(tbl) <- list(classification = c("original","new"),
                      assemblage = sort(unique(class_orig$name)))

# Chisq test for independence of assemblage and classification dataset
chisq <- chisq.test(tbl)
print(chisq)
cat("chisq =", round(chisq$statistic, 3), 
    ", p <", signif(chisq$p.value, 3), 
    ", df =", chisq$parameter)

# Cramer's V correlation
CramersV <- sqrt(chisq$statistic / sum(tbl))
cat("Cramer's V =", round(CramersV, 3))

# 3. Examine distribution of grid cells with higher effort ====

# Load subgrids for each region with populated grid cells

files <- list.files(path = "Data/Shapefiles", 
                    pattern = "^TrawlFreq_Subgrid_.+shp$", 
                    full.names = TRUE)

# Coordinate reference systems for subgrids
crs <- c(26920, 26920, 26921, 26920) 

# Read in shapefiles and filter cells with 2+ sets only
sub_grids <- map2(files, crs, ~st_read(.x, crs = .y)) %>% 
  map(., ~filter(.x, Freq > 1))

# Examine spatial distribution of these cells

# Switch to interactive view
tmap_mode("view") 

# Plot grid cells with red fill
map(sub_grids, ~tm_shape(.) + tm_polygons(col = "red")) # No particular spatial pattern 

# Load raster layers with predicted assemblage type

class_tif <- list.files(path = "Output", 
                        pattern = "PredClust_map\\.tif$", 
                        full.names = TRUE) %>%  
  map(., raster)

# Calculate area of each assemblage type

class_area <- map(class_tif, freq) %>%
  map(., data.frame) %>% 
  map(., ~filter(., !is.na(value))) %>% 
  map(., ~mutate(., area = count*16))

# Determine which assemblage type that cells with 2+ sets fall within

which_class <- map2(class_tif, sub_grids, ~raster::extract(.x, .y)) %>% 
  map(., unlist) %>% 
  map(., table) %>% 
  map(., data.frame) %>% 
  map(., ~mutate(., value = as.numeric(Var1)))

# Calculate number of cells with extra effort per cluster area

effort_area <- map2(class_area, which_class, ~left_join(.x, .y)) %>% 
  map(., ~mutate(., per_area = Freq/area)) %>% 
  map(., ~dplyr::select(.,-Var1)) %>% 
  data.table::rbindlist()

effort_area$Region <- c(rep("MAR",6), rep("NGSL",3), rep("NL",5), rep("SGSL",4))

# 4. Explore correlation between indicator taxa richness & species richness ====

# List of files with cluster assignments

files <- list.files(path = "Data/",
                    pattern = "^Colourgrid4km.+csv$",
                    full.names = TRUE)

# Calculate richness in each major assemblage type

richness <- map(files, fread) %>% 
  map(., ~ dplyr::select(., -V1, -GridID, -cl, -assigned)) %>% 
  map(., ~filter(., name != 'Unclassified')) %>% 
  map(., ~group_by(., name)) %>% 
  map(., ~summarise_all(., .funs = sum)) %>%
  map(., ~mutate_if(., is.integer, function(x){if_else(x > 0, 1, 0)})) %>% 
  map(., rowwise) %>%
  map(., ~group_by(., name)) %>% 
  map(., ~mutate(., richness = sum(c_across(everything())))) %>% 
  map(., ~dplyr::select(., name, richness)) %>% 
  set_names(., nm = c("MAR", "NGSL", "NL", "SGSL"))

# Calculate number of indicator taxa

ind_taxa <- list.files(path = "Output/", 
                       pattern = "IndicatorSpecies.+txt$",
                       full.names = TRUE) %>% 
  set_names(., nm = c("MAR", "NGSL", "NL", "SGSL")) %>% 
  map(., fread)

# Calculate indicator taxa richness by cluster

ind_rich <- map(ind_taxa, ~group_by(.x, Assemblage)) %>% 
  map(., ~summarise(.x, ind_rich = n()))

# join ind_richness and taxa richness

ind_join <- map2(richness, ind_rich, left_join)

# 5. Examine spatial distribution of minor assemblages ====

# Read in grid cells with cluster assignments

saveRDS(colour_grid, file = paste0("Data/Grid_ClusterAssignment_", region, ".rds"))

files <- list.files("Data/", 
                    pattern = "Grid_ClusterAssignment_.+rds", 
                    full.names = TRUE)

cl_grid <- map(files, readRDS)

# filter observations to only unclassified grid cells (i.e. not assigned to a major cluster)

minor_cl <- map(cl_grid, ~filter(., name == "Unclassified"))
map(minor_cl, nrow) # 5 - 10 cells unclassified in each region

# Look at number in individual minor cluster

N_minor_cl <- list.files(path = "Data/",
                         pattern = "^Colourgrid4km.+csv$",
                         full.names = TRUE) %>% 
  map(., fread) %>% 
  set_names(., nm = c("MAR", "NGSL", "NL", "SGSL")) %>%
  map(., ~group_by(., cl)) %>% 
  map(., ~summarise(., N = n()))
# mostly 1-2 in each minor cluster, max (4-8)

# 6. Figure Prep ====

## Fig. S1 ----

# nMDS comparing community structure between Spring and Fall survey trawl sets in NL

# Display with ggplot
p_mds <- ggplot(df_mds, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(col = season)) +
  guides(color = guide_legend("Season"))

# Save plot
ggsave(plot = p_mds, 
       filename = "Output/NL_mds.tiff", 
       width = 5, 
       height = 4, 
       units = "in", 
       dpi = 300)

## Fig. S10 ----

# Display number of cells with extra effort per cluster area

pal <- wes_palette("Darjeeling1", n = 4, type = "discrete") # colour palette
names(pal) = c("NGSL", "NL", "SGSL", "MAR") # name palette elements

effort_plot <- ggplot(effort_area, aes(x = area, y = Freq)) +
  geom_point(aes(col = Region)) +
  scale_color_manual(values = pal, 
                     breaks = c("NGSL", "SGSL", "NL", "MAR"),
                     aesthetics = c("colour"), 
                     guide = guide_legend("Region'", 
                                          title.theme = element_text(face = "bold"))) +
  geom_smooth(method = "gam") +
  labs(x = expression(Area~(km^2)), y = "Cells with >1 sets (N)") +
  theme_grey()

ggsave(plot = effort_plot, 
       filename = "Output/effort_vs_area.tiff", 
       dpi = 300)

## Fig. S11 ----

# Plot indicator richess ~ taxa richness relationship using all 4 regions

ind_all <- rbindlist(ind_join, idcol = "Region")

# Fit linear model to indicator richess ~ taxa richness

lm_sm <- broom::glance(lm(ind_all$ind_rich ~ ind_all$richness))

p1ot_rich <- ggplot(ind_all, aes(x = richness, y = ind_rich)) +
  geom_point(aes(col = Region)) + 
  scale_color_manual(values = pal, 
                     breaks = c("NGSL", "SGSL", "NL", "MAR"),
                     aesthetics = c("colour"), 
                     guide = guide_legend("Region", 
                                          title.theme = element_text(face = "bold"))) +
  geom_smooth(method = "lm") +
  labs(x = "Taxa richness", y = "Indicator taxa richness") +
  annotate("text", x = 63, y = 30, 
           label = paste(paste0('Rsq = ', round(lm_sm$r.squared, 3)), 
                         paste0('p = ', round(lm_sm$p.value, 3)), 
                         sep = "\n")) +
  theme_grey()

ggsave(plot = plot_rich, 
       filename = "Output/IndRich.tiff", 
       dpi = 300)

## Fig. S12 ----

# Plot distribution of minor clusters with respect to major cluster

# Plot tags
labels = c("MAR", "NGSL", "NL", "SGSL")

# Dislay with ggplot
plots_minor <- pmap(list(class_tif, minor_cl, labels), 
                    ~gplot(..1) + 
                      geom_tile(aes(fill = value),
                                show.legend = FALSE) + 
                      geom_sf(data = ..2, 
                              inherit.aes = FALSE, 
                              fill = "red", 
                              col = "red") +
                      labs(x = NULL, y = NULL) +
                      annotate("text", 
                               label = paste0("bold(",..3,")"),
                               parse = T,
                               size = 6,
                               x = ((extent(..1)[2]-extent(..1)[1])*0.85)+extent(..1)[1],
                               y = ((extent(..1)[4]-extent(..1)[3])*0.9)+extent(..1)[3])) %>% 
  set_names(labels)

# Save individual plots
map2(plots_minor, labels, 
     ~ggsave(plot = .x, 
             filename = paste0("Output/minor_cluster_distN_", .y ,".tiff")))
