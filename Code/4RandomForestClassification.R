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
# 1. Remove highly correlated environmental predictors
# 2. Fit Random Forest Classification to predict assemblage membership
# 3. Evaluation statistics with 10-fold cross-validation
# 4. Predict distribution of major assemblage types
# 5. Highlight areas of greater uncertainty
# 6. Relative influence of environmental predictors (variable importance)
# 7. Figure Prep
#
# Requirements:
# R version 3.6.3
# Land and study area polygons from 1DataPrep.R
# Save sf objects with grid cells assigned cluster membership from
# hierarchical classification in 2ClusterAnalysis.R
# Environmental predictor layers resampled to resolution
# as fishnet grid (sources described in text)
#
# House Keeping ====

# Load required packages

library(raster)
library(rasterVis)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(forcats)
library(corrplot)
library(tcltk2)
library(randomForest)
library(pROC)
library(ggplot2)
library(tmap)
library(cowplot)

# Source functions

source("Code/BioclassificationFunctions.R")

# Set controls

# String identifying the regional study area (i.e. NGSL, SGSL, NL, or MAR)
region <- "MAR"

# Spatial layers

# Polygon for coastline
Land <- st_read(paste0("Data/Shapefiles/LandBorders_", region, ".shp"))
# 5-km land buffer
LandBuffer <- Land %>% 
  st_buffer(5000) 
# Study area boundaries
StudyArea <- st_read(paste0("Data/Shapefiles/StudyArea_", region, ".shp"))
# Stack of environmental raster layers buffered 5-km from shore
env_pred <- list.files(path = "Data/Rasters/",
                       pattern = paste0(region, ".+tif$"),
                       full.names = TRUE,
                       ignore.case = TRUE) %>% 
  stack() %>% 
  mask(., LandBuffer, inverse = TRUE)

# 1. Remove highly correlated environmental predictors ====

# Examine correlation between predictors

# Pearson correlation matrix
vals <-  c() #create empty list for values
for (i in 1:nlayers(env_pred)){ #for each raster layer
  v <-  getValues(raster(env_pred, i)) #get 1D array of values
  vals <-  cbind(vals, v) #add 1D array to col of vals
}
colnames(vals) <- gsub(paste0(region, "_"),"", names(env_pred)) #rename columns

corr_matrix <- cor(vals, use="pairwise.complete.obs", method = "pearson")
diag(corr_matrix) <- NA #set same-variable correlations to NA
corr_abs <- abs(corr_matrix) # absolute values

# graphical display of correlation matrix
corr_plot <- corrplot(corr_matrix) 

# Run GUI to select variable pair with highest correlation and eliminate one
# Continue until no variables remain with correlation higher than threshold (0.7)

threshold <- 0.7

# Set up filenames
tmpfilename <- paste0("Output/", region, "_temporary.txt")
removefilename <- paste0("Output/", region, "_eliminate_list.txt")
retainfilename <- paste0("Output/", region, "_retain_list.txt")

eliminated <- c("") #empty list for eliminated variables
maxcor <- max(corr_abs, na.rm = TRUE) # find highest correlation b/w predictors
maxrow <- which(corr_abs == maxcor, arr.ind = TRUE)[1,1] # row index of maxcor
maxcol <- which(corr_abs == maxcor, arr.ind = TRUE)[1,2] # col index of maxcor

# Write-functions for 2 buttons (variable names of most highly correlated variables)

write1 <- function(){ #if button 1 is pushed
  tmpfile <- file(tmpfilename) #create temporary file
  #write name of maxrow (variable to be eliminated) to tmpfile
  writeLines(colnames(corr_abs)[maxrow], tmpfile) 
  close(tmpfile)
}
write2 <- function(){ #if button 2 is pushed
  tmpfile <- file(tmpfilename)
  # write name of maxcol (variable to be elimnated) to tmpfile
  writeLines(colnames(corr_abs)[maxcol], tmpfile) 
  close(tmpfile)
}

# Create GUI widget
root <- tktoplevel()
btn1 <- tk2button(root, text = paste(colnames(corr_abs)[maxrow]), command = write1)
tkpack(btn1)
btn2 <- tk2button(root, text = paste(colnames(corr_abs)[maxcol]), command = write2)
tkpack(btn2)

# Create 'while' loop that stops running when maxcor falls below 0.7
while(maxcor > threshold){
  tkconfigure(btn1, text = paste(colnames(corr_abs)[maxrow]), command = write1)
  tkconfigure(btn2, text = paste(colnames(corr_abs)[maxcol]), command = write2)
  # Wait until button is pushed/tmpfile is written
  while(!file.exists(tmpfilename)){a = 1} 
  tmpfile <- file(tmpfilename, "rw") # open tmpfile
  eliminate <- readLines(tmpfile, 1) # read one line
  close(tmpfile) # close tmpfile
  file.remove(tmpfilename) # delete tmpfile
  eliminated <- c(eliminated, eliminate) # eliminate variable from correlation matrix
  corr_abs <- corr_abs[-which(colnames(corr_abs) == eliminate), -which(colnames(corr_abs) == eliminate)]
  if(dim(corr_abs)[1]>1){ # if at least 2 variables remain, recalculate maxcor/row/col
    maxcor = max(corr_abs, na.rm = TRUE)
    maxrow <- which(corr_abs == maxcor, arr.ind = TRUE)[1,1] 
    maxcol <- which(corr_abs == maxcor, arr.ind = TRUE)[1,2]
  } else maxcor = 0 #otherwise end the loop
}

tkdestroy(root) # remove widget when finished
cornames <- colnames(corr_abs) # names of retained predictors
# write retained predictor names to file
write.table(cornames, file = retainfilename, row.names = FALSE) 
rmlist <- file(removefilename) # filename of elimiated predictors
writeLines(eliminated, rmlist) # write elimanted predictor names to file
close(rmlist)

# 2. Fit Random Forest Classification to predict assemblage membership ====

# Subset environmental predictors retained for modelling from above

# varlist
retained_vars <- readLines(file(paste0("Output/", region, "_retain_list.txt"))) %>% 
  gsub(pattern = "\\\"",
       replacement = "", 
       x = .) %>% 
  .[-1] %>% 
  gsub(pattern = paste0(region, "_"),
       replacement = "",
       x = .,
       ignore.case = TRUE)
  

env_pred <- env_pred %>% 
  `names<-`(., gsub(pattern = paste0(region, "_"), 
                    replacement = "" , 
                    x = names(.),
                    ignore.case = TRUE)) %>% # rename layers
  .[[retained_vars]] %>% # subset raster stack to include only retained variables
  mask(., StudyArea) # mask raster cells of predictors outside study area boundaries

# Import populated grid cells (i.e. sites) with cluster assignments (i.e. assemblage types)

# Read in shapefiles
Grid_populated <- readRDS(paste0("Data/Grid_ClusterAssignment_", region, ".rds"))

# Names of major assemblages (i.e. with more than 20 sites)
top_n <- table(Grid_populated$name) %>% 
  sort(decreasing = TRUE) %>% 
  Filter(f = function(x) x >= 20) %>% 
  names()

# Exclude smaller clusters
Grid_populated <- Grid_populated %>% 
  filter(name %in% top_n) 

# Extract raster values at each populated grid cell

RF_data <- raster::extract(x = env_pred,
                           y = Grid_populated, 
                           factors = TRUE, 
                           nl = nlayers(env_pred), 
                           df = TRUE) %>% # extract values into df
  dplyr::select(-ID) %>% 
  bind_cols(Grid_populated, .) %>% 
  drop_na() %>% 
  mutate(name = as.factor(name)) %>% 
  droplevels()

saveRDS(RF_data, file = paste0("Data/RF_data_", region, ".rds"))

# Specify model formula

formula <-  as.formula(paste("name","~",paste(retained_vars, collapse = "+"))) 

# Fit model

RF_mod <- randomForest(formula, 
                       data = RF_data, 
                       ntree = 10000, 
                       importance = TRUE)
print(RF_mod)

# Save output

RF_output <- saveRDS(RF_mod, file = paste0("Output/", region, "_RF_model.rds"))

# 3. Evaluation statistics with 10-fold cross-validation ====

nfold <-  10 # number of splits of the data
n <- nrow(RF_data) # number of observations in data
groups <- sample(rep(1:nfold, length = n), n) # assign observation to 
                                              # 1 of 10 splits of the data (groups)
CalVal.list <- 1:nfold %>% # indices for 10 calibration and validation datasets
  map(~ list(cal = which(groups != .), 
             val = which(groups == .))) 

# Compute multi-class AUC using each calibration/validation pair 

cv10 <- CalVal.list %>% 
  map(~ list(model = randomForest(formula, 
                                  data = RF_data[.$cal,], 
                                  ntree = 10000), 
             val = RF_data[.$val,])) %>% 
  map(~ list(obs = .$val$name, 
             predictions = predict(.$model, newdata = .$val, type = "prob"))) %>%  
  map(~ multiclass.roc(.$obs, .$predictions)$auc) 

# Average AUC for each 90:10 split of data

AUC <- mean(unlist(cv10)) 

cat("10-fold cross-validated Area Under the Receiver 
    Operating Characteristic Curve is", AUC) 

# 4. Predict distribution of major assemblage types ====

# Predict assemblage membership over entire raster surfaces

predict_map <- raster::predict(env_pred, RF_mod) 
plot(predict_map)

# Write raster layer to file

writeRaster(predict_map, 
            filename = paste0("Output/", region, "_PredClust_map.tif"), 
            overwrite = TRUE)

# Aggregate raster cells and vectorize boundary for figure display

raster_boundary <- rasterToPolygons(predict_map) %>% 
  st_as_sfc() %>% 
  st_union()

st_write(obj = raster_boundary, 
         dsn = "Data/Shapefiles", 
         layer = paste0(region, "_RasterBoundary"), 
         driver = "ESRI Shapefile")

# 5. Highlight areas of greater uncertainty ====

# Extract environmental raster cell values into dataframe

env_pred_df <- data.frame(getValues(env_pred)) 

# Get random forest predictions for environmental layers cell values as probability

RF_pred <- predict(object = RF_mod, 
                   newdata = env_pred_df, 
                   type = "prob")

# For each raster cell, extract proportion of vote counts 
# for the assemblage type predicted by model (i.e. extract maximum probability)

VoteCounts <- data.frame(RF_pred) %>%
  dplyr::rowwise() %>% 
  mutate(maxVC = max(c_across(everything()))) %>% 
  dplyr::select(maxVC) %>%
  data.frame()

# Create raster layer with maximum proportion of vote counts

VoteCountsRaster <- raster(env_pred[[1]]) 
VoteCountsRaster <- setValues(VoteCountsRaster, as.vector(VoteCounts$maxVC))

# Identify cells with assignment probability less than 0.70 and vectorize

VoteCounts_0.7 <- rasterToPolygons(VoteCountsRaster, 
                                   fun = function(x){x < 0.7}, 
                                   na.rm = TRUE) %>% 
  st_as_sfc()

# Write to shapefile

st_write(obj = VoteCounts_0.7, 
         dsn = "Output", 
         layer = paste0(region, "_RF_uncertainty_0.7"), 
         driver = "ESRI Shapefile")

# 6. Relative influence of environmental predictors (variable importance) ====

# Get permutation variable importance from RF model object

# Variable names

if(region == "NGSL") {
  var_names <- c("Aspect","Avg max PP (spring/summer)", "Avg mean chl (winter)","BPI (1 km)", "BPI (20 km)",
                 "DO","Avg max temperature", "Avg max SST", "Avg min temperature", 
                 "Avg min SST","Avg mean bottom stress", "Avg mean MLD","Avg mean current velocity (EW)", 
                 "Avg mean current velocity (NS)", "pH", "Slope")
} else{
  if(region == "SGSL"){
    var_names <- c("Aspect","Avg max PP (spring/summer)", "BPI (1 km)", "BPI (20 km)",
                   "Avg max temperature", "Avg max current velocity (NS)",
                   "Avg min temperature","Avg min SST","Avg mean temperature",
                   "Avg mean bottom stress", "Avg mean MLD", "Avg mean current velocity (EW)", 
                   "Avg mean current velocity (NS)", "Avg mean MLD (summer)",
                   "Avg mean SST (winter)", "Slope")
  } else {
    if(region == "NL"){
      var_names <- c("Aspect","Avg max PP (spring/summer)","Avg min chl", "Avg mean PP",
                     "Avg mean chl (fall)", "Avg mean chl (spring)","Avg mean chl (summer)",
                     "Avg mean chl (winter)", "Depth", "BPI (0.5 km)", "BPI (20 km)","DO",
                     "Avg max temperature", "Avg max current velocity (EW)", "Avg max current velocity (NS)",
                     "Avg min temperature","Avg mean bottom stress", "Average mean MLD", 
                     "Avg mean current velocity (EW)", "pH", "Avg range salinity", "Slope")
    } else{
      if(region == "MAR"){
        var_names <- c("Aspect","Avg max PP (spring/summer)","Depth","DO","BPI (0.5 km)", "BPI (10 km)",
                       "Avg max salinity","Avg max temperature","Avg min temperature",
                       "Avg mean bottom stress","Slope")
      }
    }
  }
}

# Column names

if(region == "NGSL") {
  col_names <- c("Channel Heads &\nSlopes","Deep Channels","Shallow Banks &\nStraits")
} else{
  if(region == "SGSL"){
    col_names <- c("Inshore/\nMagdalen Is.","Laurentian Channel","Magdalen Shallows", 
                   "Northumberland Strait/\nSt. George's Bay")
  } else {
    if(region == "NL"){
      col_names <- c("Grand Banks","Inner Shelf","Laurentian Channel/\nShelf Break",  
                     "Outer Shelf","Slope")
    } else{
      if(region == "MAR"){
        col_names <- c( "ESS","ESS: Banks","Laurentian Channel/\nShelf Break",
                        "Slope","WSS/Outer BoF", "WSS: Banks/\nInner BoF")
      }
    }
  }
}

# Create dataframe with variable importance for whole model and by class

varImp <- data.frame(importance(RF_mod, scale = FALSE)) %>% 
  bind_cols(., predictor = as.factor(var_names)) %>%
  mutate(predictor = fct_reorder(predictor, MeanDecreaseAccuracy, .desc = FALSE)) %>% 
  arrange(desc(MeanDecreaseAccuracy)) %>%
  rename("Whole Model" = "MeanDecreaseAccuracy") %>% 
  rename(!!!set_names(names(.)[1:(length(RF_mod$classes))], col_names)) %>%
  dplyr::select(-MeanDecreaseGini) %>% 
  pivot_longer(cols = !predictor,
               names_to = "class",
               values_to = "MeanDecreaseAccuracy") %>% 
  mutate(class = relevel(as_factor(class), "Whole Model", 1))

# 8. Figure Prep ====

## Fig. 3 ----

# Colour palette

pal <- eval(parse(text = paste0(region, ".palette"))) %>% 
  arrange(name) %>% 
  dplyr::select(-cl) %>% 
  filter(name != "Unclassified")

# Map with predicted distribution of major clusters

# Raster with predicted distribution of major clusters from random forest
RF_pred <- raster(paste0("Output/", region, "_PredClust_map.tif")) %>% 
  ratify()

# Add other categories to Raster Attribute Table

rat <- levels(RF_pred) %>% 
  bind_cols(., pal)

levels(RF_pred) <- rat

# Convert to sf
RF_pred_sf <- as(RF_pred, "SpatialPolygonsDataFrame") %>%
  st_as_sf() %>% 
  `names<-`(., c("ID", "geometry")) %>% 
  left_join(., levels(RF_pred)[[1]])

# Polygon of raster boundaries
raster_boundary <- st_read(paste0("Data/Shapefiles/", region, "_RasterBoundary.shp")) %>% 
  st_set_crs(crs_proj)

# Polygons highlighting grid cells with assignment probability < 0.70
uncertainty <- st_read(paste0("Output/", region, "_RF_uncertainty_0.7.shp")) %>% 
  st_set_crs(crs_proj)

if(region == "NGSL") {
  labelpos <-  "NW"
  xlim <- c(30174, 1026175)
  ylim <- c(5224423, 5820424)
  xbreak <- c(-70,-66,-62,-58)
  ybreak <- c(48,50,52)
  regiontag <- region
  tagpos <- c(900000, 5800000)
  ypad <- 0.5
} else{
  if(region == "SGSL") {
    labelpos = "SW"
    xlim <- c(280798, 720798)
    ylim <- c(5063640, 5447640)
    xbreak <- c(-65,-63,-61)
    ybreak <- c(49,48,47,46)
    regiontag <- region
    tagpos <- c(685000, 5420000)
    ypad <- 1.2
  } else{
    if(region == "NL") {
      labelpos = "NE"
      xlim <- c(270169.8, 1314169.8)
      ylim <- c(4742510.1, 6398510.1)
      xbreak <- c(-48,-54,-60)
      ybreak <- c(44,48,52,56)
      regiontag <- region
      tagpo <- c(1200000, 6300000)
      ypad <- 0.25
    } else{
      if(region == "MAR") {
        labelpos = "SE"
        xlim <- c(0, 1040000)
        ylim <- c(4530000, 5300000)
        xbreak <- c(-67,-63,-59)
        ybreak <- c(47,45,43,41)
        regiontag <- region
        tagpos <- c(920000, 5220000)
        ypad <- 0.5
      }
    }
  }
}

if(region == "MAR") {
  expand <- FALSE
  expansion <- expr(waiver())
} else{
  expand <- TRUE
  expansion <- c(0.02,0)
}

if(region %in% c("NSGL", "SGSL")) {
  yrot <- 90
} else{
  yrot <- 270
}

if(region %in% c("SGSL", "NL")) {
  scalepos <- "bl"
} else{
  scalepos <- "tl"
}

if(region == "SGSL"){
  scalewid <- 0.2
} else{
  scalewid <- 0.25
}

map_distrib <- plot_colourmap(col_map = RF_pred_sf,
                              boundaries_poly = raster_boundary,
                              land_poly = Land,
                              uncertainty = uncertainty,
                              label_position = labelpos,
                              expand_scale = expand,
                              xlim = xlim,
                              ylim = ylim,
                              x_breaks = xbreak,
                              y_breaks = ybreak,
                              expansion_const = expansion,
                              region_tag = regiontag,
                              tag_position = tagpos,
                              ytext_angle = yrot,
                              scalebar_position = scalepos,
                              scalebar_width = scalewid,
                              scalebar_ypad = ypad)

# Save map
saveRDS(map_distrib, 
        file = paste0("Output/", region, "_RF_predictions_ggmap_uncertainty_0.7.rds"))

if(region == "NGSL") {
  ploth <- 5.4
  plotw <- 9
} else{
  if(region == "SGSL") {
    ploth <- 5.8
    plotw <- 6.7
  } else{
    if(region == "NL") {
      ploth <- 8.8
      plotw <- 5.6
    } else{
      if(region == "MAR") {
        ploth <- 5.8
        plotw <- 8
      }
    }
  }
}

ggsave(plot = map_distrib, 
       filename = paste0("Output/", region, "_RF_predictions_ggmap_uncertainty_0.7.tiff"),
       compression = "lzw", 
       dpi = 300,
       width = plotw, 
       height = ploth, 
       unit = "in")

# Create legend for higher uncertainty areas

legend <- tm_shape(uncertainty) +
  tm_borders(col = "black") +
  tm_add_legend(type = "fill",
                col = "white",
                labels = "Prob. of Assignment < 0.70") +
  tm_layout(legend.text.size = 3,
            legend.title.size = 4,
            legend.only = TRUE,
            legend.width = 10,
            legend.height = 1,
            legend.title.fontface = "bold"); legend

saveRDS(legend, "Output/tmap_legend_RF_prediction.rds")

tmap_save(legend, 
          filename = "Output/tmap_legend_RF_prediction.tiff", 
          width = 10, 
          height = 1, 
          unit = "in", 
          compression = "lzw")

## Fig. 4 ----

# Plot variable importance - whole model only

varImp_plot <- filter(varImp, class == "Whole Model") %>% 
  ggplot(., aes(x = predictor, y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity") +
  facet_grid(~class, labeller = as_labeller(c(`Whole Model` = region))) +
  theme_bw() +
  coord_flip() +
  labs(x = NULL, y = "Mean Decrease in Accuracy") + 
  theme(axis.text.x = element_text(angle = 90), 
        text = element_text(size = 10)); varImp_plot

# Save individual results for each region
saveRDS(varImp_plot, file = paste0("Output/", region,"_VarImpPlot_WM.rds"))

# Combine plots for 4 regions

# Import plots from all 4 regions
varImp.wm <- list.files("Output/", 
                        pattern = "WM.rds$", 
                        full.names = TRUE) %>%
  map(readRDS) %>% # read in ggplots
  map(., ~ . + theme(axis.title.x = element_blank(),# remove x-axis label
                     strip.text = element_text(size = 14, face = 'bold'), # adjust font of strip text
                     axis.text = element_text(size = 12)) + # increase axis text size 
        scale_y_continuous(limits = c(-0.05,0.25))) %>% # apply consistent scale across regions
  set_names(c("MAR", "NGSL", "NL", "SGSL")) %>%  # name list elements by region
  .[c("NGSL", "NL", "SGSL", "MAR")] # re-arrange list elements

# Plot on 2 X 2 grid

varImp_4regions <- plot_grid(plotlist = varImp.wm, 
                             nrow = 2, ncol = 2, 
                             align = "hv", 
                             axis = "lb", 
                             labels = "Mean Decrease in Accuracy",
                             label_x = c(0.7),
                             label_y = c(-1),
                             hjust = -0.5) + 
  theme(plot.margin = ggplot2::margin(t = 0, r = 0,b = 40, l = 0))

# Save output

ggsave(plot = varImp_4regions, 
       filename = "Output/VarImp_4regions_WholeModel.tiff", 
       height = 9, 
       width = 11, 
       unit = "in", 
       compression = "lzw", 
       dpi = 300)

## Fig. S5 ----

# Plot variable importance - faceted by cluster

varImpPlot_grid <- ggplot(varImp, aes(x = predictor, 
                                      y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity") +
  facet_grid(~class) +
  theme_bw() +
  coord_flip() +
  labs(x = NULL, y = "Mean Decrease in Accuracy") + 
  theme(axis.text.x = element_text(angle = 90), 
        text = element_text(size = 10)); varImpPlot_grid

# Save individual results for each region
ggsave(plot = varImpPlot_grid,
       filename = paste0("Output/", region, "_varImpPlot_grid.tiff"),
       width = 9,
       height = 3.5,
       units = "in",
       dpi = 300,
       compression = "lzw")

saveRDS(varImpPlot_grid, 
        file = paste0("Output/", region, "_varImpPlot_grid.rds"))

# Import variable importance plots from all 4 regions

varImp.ls <- list.files("Output/", 
                        pattern = 'ImpPlot_grid.rds$', 
                        full.names = TRUE) %>% 
  map(readRDS) %>% # read in ggplots
  map(., ~ . + theme(axis.title.x = element_blank(), # remove x-axis label
                     axis.text = element_text(size = 12), # increase axis text size
                     strip.text = element_text(size = 12)) + # increase strip text size
        scale_y_continuous(limits = c(-0.05,0.55))) %>% # apply consistent scale across regions
  set_names(c("MAR", "NGSL", "NL", "SGSL")) # name list elements by region

# Plot variable importance plots in 4 rows; 1 per region
varImp_facet_4regions <- ggdraw() +
  draw_plot(varImp.ls$NGSL, x = 0, y = 0.75, width = 0.625, height = 0.25) +
  draw_plot(varImp.ls$SGSL, x = 0, y = 0.5, width = 0.75, height = 0.25) +
  draw_plot(varImp.ls$NL, x = 0, y = 0.25, width = 0.875, height = 0.25) +
  draw_plot(varImp.ls$MAR, x = 0.01, y = 0, width = 0.99, height = 0.25) +
  draw_plot_label(label = c('NGSL','SGSL','NL','MAR'), size = 15, hjust = 1,
                  x = c(0.12,0.12,0.12,0.12), y = c(1,0.75,0.5,0.25))
# Save output

ggsave(plot = varImp_facet_4regions, 
       filename = "Output/VarImp_4regions_ByClass.tiff", 
       height = 17, 
       width = 16, compression = "lzw")

## Fig. S6 ----

# Examine variation in important predictors across clusters - NGSL

# boxplot for min annual bottom temperature
box_minT <- boxplot_colour(.data = RF_data,
                            x = factor(cl, levels = c("1","6","2")),
                            y = min_ann_BT,
                            colour_var = assigned,
                            ytitle = "Avg min temperature (\u00B0C)",
                            xlabs = c("Deep Channels", "Shallow Banks & Straits",
                                      "Channel Heads & Slopes"))

# boxplot for max annual bottom temperature 
box_maxT <- boxplot_colour(.data = RF_data,
                           x = factor(cl, levels = c("1","6","2")),
                           y = max_ann_BT,
                           colour_var = assigned,
                           ytitle = "Avg max temperature (\u00B0C)",
                           xlabs = c("Deep Channels", "Shallow Banks & Straits",
                                      "Channel Heads & Slopes"))

# boxplot for dissolved oxygen
box_DO <- boxplot_colour(.data = RF_data,
                         x = factor(cl, levels = c("1","6","2")),
                         y = dissox,
                         colour_var = assigned,
                         ytitle = expression(DO~'('*mg~l^{-1}*')'),
                         xlabs = c("Deep Channels", "Shallow Banks & Straits",
                                   "Channel Heads & Slopes"))

# boxplot for BPI-20km
box_bpi20 <- boxplot_colour(.data = RF_data,
                            x = factor(cl, levels = c("1","6","2")),
                            y = bpi20,
                            colour_var = assigned,
                            ytitle = "BPI (20 km)",
                            xlabs = c("Deep Channels", "Shallow Banks & Straits",
                                      "Channel Heads & Slopes"))

# combine boxplots in 4 panel figure

boxplot_4panel(tl = box_minT,
               tr = box_DO,
               bl = box_maxT,
               br = box_bpi20,
               bioregion = "NGSL")

## Fig. S7 ----

# Examine variation in important predictors across clusters - SGSL

# boxplot for max annual bottom temperature
box_maxT <- boxplot_colour(.data = RF_data,
                           x = factor(cl, levels = c("6","7","10","1")),
                           y = max_ann_BT,
                           colour_var = assigned,
                           ylim = c(0,25),
                           ytitle = "   Avg max\ntemperature (\u00B0C)",
                           xlabs = c("Magdalen Shallows", "Inshore/Magdalen Is.",
                                     "Laurentian Channel","N. Strait &\n St. George's Bay"))

# boxplot for mean annual bottom temperature
box_meanT <- boxplot_colour(.data = RF_data,
                           x = factor(cl, levels = c("6","7","10","1")),
                           y = mn_ann_BT,
                           colour_var = assigned,
                           ylim = c(0,10),
                           ytitle = "   Avg mean\ntemperature (\u00B0C)",
                           xlabs = c("Magdalen Shallows", "Inshore/Magdalen Is.",
                                     "Laurentian Channel","N. Strait &\n St. George's Bay"))

# boxplot for min annual bottom temperature
box_minT <- boxplot_colour(.data = RF_data,
                           x = factor(cl, levels = c("6","7","10","1")),
                           y = min_ann_BT,
                           colour_var = assigned,
                           ytitle = "   Avg min\ntemperature (\u00B0C)",
                           xlabs = c("Magdalen Shallows", "Inshore/Magdalen Is.",
                                     "Laurentian Channel","N. Strait &\n St. George's Bay"))

# boxplot for Avg max Primary Production (spring/summer)  
box_PP <- boxplot_colour(.data = RF_data,
                           x = factor(cl, levels = c("6","7","10","1")),
                           y = avg_max_sprsum_PP,
                           colour_var = assigned,
                           ytitle = expression(Avg~max~PP~'('*mg~C~m^{-2}~day^{-1}*')'),
                           xlabs = c("Magdalen Shallows", "Inshore/Magdalen Is.",
                                     "Laurentian Channel","N. Strait &\n St. George's Bay"))

# combine boxplots in 4 panel figure

boxplot_4panel(tl = box_meanT,
               tr = box_maxT,
               bl = box_minT,
               br = box_PP,
               bioregion = "SGSL")

## Fig. S8 ----

# Examine variation in important predictors across clusters - NL

# boxplot for depth 
box_depth <- boxplot_colour(.data = RF_data,
                            x = factor(cl, levels = c("7","8","6","1","4")),
                            y = bathy,
                            colour_var = assigned,
                            ytitle = "Depth (m)",
                            xlabs = c("Inner Shelf", "Outer Shelf", "Grand Banks", 
                                      "Slope", "LC/Shelf Break"))

# boxplot for min annual bottom temperature
box_minT <- boxplot_colour(.data = RF_data,
                            x = factor(cl, levels = c("7","8","6","1","4")),
                            y = min_ann_BT,
                            colour_var = assigned,
                            ytitle = "   Average min\ntemperature (\u00B0C)",
                            xlabs = c("Inner Shelf", "Outer Shelf", "Grand Banks", 
                                      "Slope", "LC/Shelf Break"))

# boxplot for max annual bottom temperature
box_maxT <- boxplot_colour(.data = RF_data,
                           x = factor(cl, levels = c("7","8","6","1","4")),
                           y = max_ann_BT,
                           colour_var = assigned,
                           ylim = c(0,11),
                           ytitle = "   Average max\ntemperature (\u00B0C)",
                           xlabs = c("Inner Shelf", "Outer Shelf", "Grand Banks", 
                                     "Slope", "LC/Shelf Break"))

# boxplot for range annual bottom salinity
box_ranSal <- boxplot_colour(.data = RF_data,
                           x = factor(cl, levels = c("7","8","6","1","4")),
                           y = range_ann_BSal,
                           colour_var = assigned,
                           ytitle = "Average range salinity (\u2030)",
                           xlabs = c("Inner Shelf", "Outer Shelf", "Grand Banks", 
                                     "Slope", "LC/Shelf Break"))

# combine boxplots in 4 panel figure

boxplot_4panel(tl = box_depth,
               tr = box_minT,
               bl = box_maxT,
               br = box_ranSal,
               bioregion = "NL")

## Fig. S9 ----

# Examine variation in important predictors across clusters - MAR

# boxplot for depth
box_depth <- boxplot_colour(.data = RF_data,
                            x = factor(cl, levels = c("1","5","4","2","9","8")),
                            y = bathy,
                            colour_var = assigned,
                            ytitle = "Depth (m)",
                            xlabs = c("Slope","LC/Shelf Break","ESS","ESS: Banks",
                                      "WSS/Outer BoF", "WSS: Banks/Inner BoF"))
  
# boxplot for min annual bottom temperature 
box_minT <- boxplot_colour(.data = RF_data,
                           x = factor(cl, levels = c("1","5","4","2","9","8")),
                           y = min_ann_BT,
                           colour_var = assigned,
                           ytitle = "   Average min\ntemperature (\u00B0C)",
                           xlabs = c("Slope","LC/Shelf Break","ESS","ESS: Banks",
                                     "WSS/Outer BoF", "WSS: Banks/Inner BoF"))

# boxplot for max annual bottom temperature
box_maxT <- boxplot_colour(.data = RF_data,
                           x = factor(cl, levels = c("1","5","4","2","9","8")),
                           y = max_ann_BT,
                           colour_var = assigned,
                           ylim = c(0,16),
                           ytitle = "   Average max\ntemperature (\u00B0C)",
                           xlabs = c("Slope","LC/Shelf Break","ESS","ESS: Banks",
                                     "WSS/Outer BoF", "WSS: Banks/Inner BoF"))
  
# boxplot for max annual bottom salinity
box_maxSal <- boxplot_colour(.data = RF_data,
                             x = factor(cl, levels = c("1","5","4","2","9","8")),
                             y = max_ann_BSal,
                             colour_var = assigned,
                             ytitle = "Average max salinity (\u2030)",
                             xlabs = c("Slope","LC/Shelf Break","ESS","ESS: Banks",
                                       "WSS/Outer BoF", "WSS: Banks/Inner BoF"))

# combine boxplots in 4 panel figure

boxplot_4panel(tl = box_depth,
               tr = box_minT,
               bl = box_maxSal,
               br = box_maxT,
               bioregion = "MAR")

cat("Proceed to 5ClimateVulnerability.R")
