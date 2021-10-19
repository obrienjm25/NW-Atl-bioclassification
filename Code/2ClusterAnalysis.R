# Project and Code Overview ====
#
# Project:      Bioclassification
# Contact:      e-mail: John.Obrien@dfo-mpo.gc.ca | tel: +1.782.640.1522
# Publication:  O'Brien JM, Stanley RRE, Jeffery NW, Heaslip SG, DiBacco C, Wang Z
#               (2021) Modelling demersal fish and benthic invertebrate assemblages 
#               in support of marine conservation planning. Ecol Appl
#
# Overview:
# 1. Choose dissimilarity measure
# 2. Choose clustering algorithm
# 3. Construct dendrograms from Site X Species Matrices with chosen dissimilarity
#    measure and clustering algorithm
# 4. Choose dissimilarity threshold by which to cut dendrograms
# 5. Identify major clusters of demersal fish and benthic invertebrates
# 6. Colour-code dendrograms with cluster assignments
# 7. Match cluster assignments to grid locations
# 8. Figure Prep
#
# Requirements:
# R version 3.6.3
# Site X Species matrix derived from 1DataPrep.R
# 4-km fishnet grid spanning study area
#

# House Keeping ====

# Load required packages

library(simba)
library(dendroextras)
library(dendextend)
library(NbClust)
library(dplyr)
library(tidyr)
library(purrr)
library(sf)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggdendro)
library(ggspatial)
library(tmap)

if(!require(dendextendRcpp)) {
  if (!require('devtools')) install.packages('devtools'); 
  devtools::install_github('talgalili/dendextendRcpp')
}

library(dendextendRcpp)

# Source functions

source("Code/BioclassificationFunctions.R")

# Set controls

# String identifying the regional study area (i.e. NGSL, SGSL, NL, or MAR)
region <- "MAR"

# Coordinate reference system of desired projection (integer EPSG code)
crs_proj <- ifelse(region == "NL", 26921, 26920)

# Spatial layers

# Land polygons
Land <- st_read(paste0("Data/Shapefiles/LandBorders_", region, ".shp"))

# Study area boundaries
StudyArea <- st_read(paste0("Data/Shapefiles/StudyArea_", region, ".shp"))
                      
# Load fishnet grid
grid <- st_read(paste0("Data/Shapefiles/Grid4km_", region, ".shp")) %>% 
  st_set_crs(crs_proj)

# Load Site X Species matrix
SiteXSpecies <- read.csv(file = paste0("Data/ClusterData_", region, ".csv"),
                         stringsAsFactors = FALSE,
                         row.names = 1)

# Get colour palettes
source("Code/ColourPalettes.R")

# 1. Choose dissimilarity measure ====

# Choices for dissimilarity measure

diss_measure <-c('simpson','soerensen','jaccard','gower')

# Compute distance matrix using each choice

dist_matrix <- map(diss_measure, ~sim(SiteXSpecies, method = .x)) %>% 
  set_names(nm = diss_measure)

# Create dendrogram object from each distance matrix using
# hierachical clustering

dendro_obj <- map(dist_matrix, ~hclust(., method = 'average')) 

# Compare correlation between resulting dendrograms and original
# distance matrix with cophenetic correlation coefficient

copCor_diss <- map(dendro_obj, cophenetic) %>% 
  map2_dfr(dist_matrix, . , 
       ~ data.frame(diss_measure = attributes(.x)$method, 
                    correlation = cor(.x, .y))) %>% 
  arrange(desc(correlation)) %>% 
  print()

best_diss <- as.character(copCor_diss$diss_measure[[1]])

cat("Best choice of dissimilarity measure is:", "\n\n", 
    best_diss, "=", max(copCor_diss$correlation))

# 2. Choose clustering algorithm ====

# Choice of 8 common clustering algorithm

algo <- c('ward.D', 'ward.D2','sin','com', 'ave','mcq','med','cen') 

# Compare correlation between resulting dendrograms and original
# distance matrix with cophenetic correlation coefficient

copCor_algo <- map(algo, # create dist. matrix and dendrogram with each algorithm
                   ~ list(dm = dist_matrix$simpson, 
                          do = hclust(dist_matrix$simpson, method = .x))) %>% 
  #calculate and print cophenetic correlation for each combination
  map_dfr( ~ data.frame(algorithm = .x$do$method, 
                        correlation = cor(.x$dm, cophenetic(.x$do)))) %>% 
  arrange(desc(correlation)) %>% 
  print()

best_algo <- as.character(copCor_algo$algorithm[1])

cat("Best choice of clustering algorithm is:", "\n\n", 
    best_algo, "=", max(copCor_algo$correlation))

# 3. Get dendrogram from best dissimilarity measure/clustering algorithm ====

best_dendro <- hclust(d = dist_matrix[[best_diss]],
                      method = best_algo)

saveRDS(best_dendro, file = paste0("Data/BestDendrogram_", region, ".rds"))

# 4. Choose dissimilarity threshold to cut dendrogram ====

# Calculate various internal cluster validity indices 
# as a function of number of clusters (k = 2 - 20)

cvi <- c('ch', 'cindex', 'ptbiserial','db','silhouette')
nclust <- c(2:20)
vals <- list()

for (i in cvi) {
  
  vals[[i]] <- NbClust(data = SiteXSpecies, 
                       diss = dist_matrix[[best_diss]], 
                       distance = NULL, 
                       method = 'average',
                       index = i, 
                       min.nc = 2, 
                       max.nc = 20)$All.index
  
}

# Combine results for all indices into one dataframe

cvi_df <- data.frame(NC = rep(nclust, length(cvi)), 
                     Value = unlist(vals))
cvi_df$Index <- gsub('[.][0-9]+','',row.names(cvi_df))

# Plot results

plot_cvi <- ggplot(cvi_df, aes(x = NC, y = Value)) +
  geom_line(size = 1, col = 'blue3') +
  geom_point() +
  labs(x = 'Number of clusters', y = 'Index value') +
  facet_wrap(~ Index, ncol = 1, scales = 'free_y')+
  theme_bw(); plot_cvi

# Examine plots and choose consensus k among indices; if multiple local optima,
# favour solutions with greater number of clusters

# Save plot

saveRDS(plot_cvi, file = paste0("Output/CVI_", region, ".rds"))

# Determine which height will result in which k for dendrogram

heights <- dendextendRcpp_heights_per_k.dendrogram(as.dendrogram(best_dendro)) 
heights <- heights[as.character(c(2:20))] %>% 
  round(., digits = 3) %>% 
  print()

# Create table with cutoff threshold for each region

cutoff_tbl <- tibble(region = c("NGSL", "SGSL", "NL", "MAR"),
                     cutoff = c(0.494, 0.494, 0.512, 0.587))

# 5. Identify major clusters of demersal fish and benthic invertebrates ====

# Cut tree at chosen threshold to generate clusters

thresh <- cutoff_tbl$cutoff[cutoff_tbl$region == region]
cl <- dendroextras::slice(best_dendro, 
                          h = thresh) 

# Get frequency table of cluster memberships

cl_tbl <-as.data.frame(table(cl)) %>%  
  arrange(desc(Freq)) %>% # order by frequency
  print() #Number of sites in each cluster

# Identify number of major clusters (i.e. with 20+ sites)

top_n <- nrow(filter(cl_tbl, Freq >= 20))

# Count the number of sites in the top clusters

top_sites <- sum(cl_tbl$Freq[cl_tbl$Freq >= 20]) 

cat("Using a threshold of", thresh, "there are", top_n, 
    "major clusters (i.e. with 20+ sites).", "\n\n", top_sites,
    "sites were assigned to a major cluster.")

# 6. Colour-code dendrograms with cluster assignments ====

# Assign colours to clusters

pal <- eval(parse(text = paste0(region, ".palette")))

colour_tbl <- left_join(cl_tbl, pal) %>% 
  replace_na(list(assigned = "#e0e0e0", name = "Unclassified")) %>% #set all clusters that occur < cut-off to grey
  arrange(as.numeric(cl)) %>% #order by cluster ID
  print() #Each cluster is assigned a color

# Plot color-coded dendrogram

# Assign colours to branches

colourtree <- colour_branches(best_dendro, 
                             h = thresh,
                             col = colour_tbl$assigned) %>%   
  as.ggdend() %>% 
  segment() %>% 
  replace_na(list(col = '#000000')) %>% 
  left_join(., colour_tbl, by = c('col' = 'assigned')) %>% 
  mutate(lwd = if_else(col == '#000000', 0.75, 0.5))

# 7. Match cluster assignments to grid locations ====

# Join cluster and colour assignments with grids

colour_grid <- tibble(GridID = names(cl), cl) %>%
  mutate(GridID = as.integer(GridID),
         cl = as.factor(cl)) %>% 
  left_join(., grid) %>% 
  left_join(., colour_tbl) %>% 
  filter(!is.na(cl)) %>% 
  arrange(GridID) %>% 
  st_sf(crs = crs_proj) %>% 
  st_cast("MULTIPOLYGON") %>% 
  filter(lengths(st_within(., StudyArea)) > 0) %>% 
  filter(lengths(st_disjoint(., st_buffer(Land, 5000))) > 0)

saveRDS(colour_grid, file = paste0("Data/Grid_ClusterAssignment_", region, ".rds"))

# Join Site X Species matrix with cluster assignment and to write csv

SiteXSpecies_ClustAssigned <- SiteXSpecies %>% 
  mutate(GridID = as.integer(rownames(.))) %>% 
  left_join(., colour_grid, by= "GridID") %>%  		
  dplyr::select(-any_of(c("layer", "Freq", "SP_ID", "geometry")))

write.csv(SiteXSpecies_ClustAssigned, 
          file =  paste0("Data/Colourgrid4km_cluster_assignment_", region, ".csv"))

# 8. Figure Prep ====

## Fig. 2 ----

# Dendrograms - 4 regions

# Display with colour-coded dengrogram with ggplot

if(region %in% c("NSGL", "NL")) {
  tagpos <- c(0.8, 0.5)
  legendpos <- "top"
  legendjust <- c(0,0)
} else{
  tagpos <- c(0.8, 0.92)
  legendpos <- "bottom"
  legendjust <- c(0,1)
}

if(region %in% c("NSGL", "SGSL")) {
  yticks <- element_text(size = 12)
  ytitle <- element_text(size = 14)
} else{
  yticks <- ytitle <- element_blank()
}

gg.colourtree <- plot_ggdendro(dendro_data = colourtree,
                               palette = pal,
                               region_tag = region,
                               tag_position = tagpos,
                               ytick_args = yticks,
                               ylab_args = ytitle,
                               legend_position = legendpos,
                               legend_justification = legendjust)

# Save colour-coded dendrogram
saveRDS(gg.colourtree, file = paste0("Output/ggDendrogram_", region, ".rds"))

# Read in dendrograms for 4 regions
dendro.ls <- list.files(path = "Output/", 
                        pattern = "ggDendro", 
                        full.names = TRUE) %>% # file list
  map(readRDS) %>% # Read files in list
  set_names(c("MAR","NGSL","NL","SGSL")) %>% # name elements in list
  .[c("NGSL","NL","SGSL","MAR")] # re-order list

# convert gg objects to grobs
dendro.ls <- map(dendro.ls, ggplotGrob) 

# Create first 2-panel row
g.r1 <- cbind(dendro.ls$NGSL, dendro.ls$NL, size = "last") # bind panels

# Create 2nd 2-panel row
g.r2 <- cbind(dendro.ls$SGSL, dendro.ls$MAR, size = "last")

# Bind columns
dendro.4panel <- rbind(g.r1, g.r2, size = "first")

# Save 4 panel figure as tiff with lzw compression
ggsave(dendro.4panel, 
       filename = "Output/Dendrograms_4regions.tiff", 
       width = 7, 
       height = 10, 
       compression = "lzw", 
       dpi = 600)

## Fig. S2 ----

# Cluster Validity Indices - 4 regions

# Read in files

saveRDS(plot_cvi, file = paste0("Output/CVI_", region, ".rds"))

p.cvi <- list.files("Output/", 
                    pattern = "CVI.+rds", 
                    full.names = TRUE) %>% # list of rds objects
  map(readRDS) %>% # read rds objects
  map(., ~ . + theme(axis.title = element_blank())) %>% # remove axis labels
  set_names(c("MAR", "NGSL", "NL", "SGSL")) %>%  # name list elements by region
  .[c("NGSL", "SGSL", "NL", "MAR")] # re-arrange list elements'

# Plot on 4 column grid

p <- plot_grid(plotlist = p.cvi, 
               nrow = 1, ncol = 4, 
               align = "hv", 
               axis = "lb", 
               labels = c("Number of Clusters"),
               label_x = c(0.5),
               label_y = c(0),
               hjust = -1) + 
  theme(plot.margin = ggplot2::margin(40,0,40,40)) +
  draw_label("Index value", 
             x = 0, 
             y = 0.5,
             angle = 90, 
             fontface = "bold", 
             size = 14, 
             vjust = -1) + 
  annotate(geom = "text", 
           x = c(0.125,0.375,0.635,0.875), 
           y = c(1,1,1,1), 
           label = c("NGSL","SGSL","NL","MAR"),
           vjust = -1,
           hjust = 0,
           fontface = "bold")

# Save output

ggsave(plot = p, 
       "Output/CVI_4regions.tiff", 
       height = 7.5, 
       width = 7, 
       unit = "in", 
       compression = "lzw", 
       dpi = 300)

## Fig. S3 ----

# Plot distribution of major clusters (i.e. predominant assemblage types)

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

map_clusters <- plot_colourgrid(col_grid = colour_grid,
                                boundaries_poly = StudyArea,
                                land_poly = Land,
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
saveRDS(map_clusters, file = paste0("Output/ColourGrid_ggmap_", region, ".rds"))

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

ggsave(plot = map_clusters, 
       filename = paste0("Output/ColourGrid_ggmap_", region, ".tiff"),
       compression = "lzw", dpi = 300,
       width = plotw, height = ploth, unit = "in")

# Create and save separate map legend

legend_map <- tm_shape(colour_grid) +
  tm_polygons() +
  tm_add_legend(type = 'fill',
                labels = pal$name,
                col = pal$assigned,
                title = region) +
  tm_layout(legend.text.size = 3,
            legend.title.size = 4,
            legend.only = T,
            legend.width = 10,
            legend.height = 5,
            legend.title.fontface = 'bold'); legend_map

saveRDS(legend_map, file = paste0("Output/legend_map_", region, ".rds"))

if(region == "SGSL") {
  legendh <- 5
  legendw <- 12
} else{
  legendh <- 5
  legendw <- 10
}

tmap_save(legend_map, filename = paste0("Output/tmap_legend_", region, ".tiff"), 
          width = legendw, 
          height = legendh, 
          unit = "in", 
          compression = "lzw")

cat("Proceed to 3IndicatorAnalysis.R")