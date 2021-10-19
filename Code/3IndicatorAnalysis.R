# Project and Code Overview ====
#
# Project:      Bioclassification
# Contact:      e-mail: John.Obrien@dfo-mpo.gc.ca | tel: +1.782.640.1522
# Publication:  O'Brien JM, Stanley RRE, Jeffery NW, Heaslip SG, DiBacco C, Wang Z
#               (2021) Modelling demersal fish and benthic invertebrate assemblages 
#               in support of marine conservation planning. Ecol Appl
#
# Overview:
# 1. Some data wrangling with Site X Species matrix with cluster assignments
# 2. Identify taxa emblematic of the predominant assemblage types in each region
# 3. Calculate frequency of occurrence of indicator taxa in associated cluster
# 4. Create table of strong indicators (IndVal > 0.25, p < 0.5)
# 5. Figure Prep
#
# Requirements:
# R version 3.6.3
# Cluster assignments from 2ClusterAnalysis.R
# Taxa list derived from trawl data used in 1DataPrep.R
#

# House Keeping ====

# Load required packages

library(labdsv)
library(dplyr)
library(purrr)
library(data.table)
library(ggplot2)

# Set controls

# String identifying the regional study area (i.e. NGSL, SGSL, NL, or MAR)
region <- "MAR"

# Get colour palettes
source("Code/ColourPalettes.R")

# 1. Data wrangling ====

# Import Site by Species matrix with cluster assignments

colorGridData <- fread(paste0("Data/Colourgrid4km_cluster_assignment_", region, ".csv")) %>% 
  dplyr::select(-V1) #remove superfluous index variable

# List of top clusters

top_n <- table(colorGridData$name) %>% 
  sort(decreasing = TRUE) %>% 
  Filter(f = function(x) x >= 20) %>% 
  names()

# Subset Site X Species matrix to top clusters

colorGridData_top <- colorGridData %>% 
  filter(name %in% top_n) 

PA_data <- colorGridData_top %>% 
  dplyr::select(-c("GridID", "cl", "assigned", "name"))

# 2. Indicator Analysis ====

# Indicator analysis as proposed by Dufrene & Legendre 1997

indicators <- indval(PA_data, colorGridData_top$name) 

# Get cluster with max indval for that taxon
maxcls <- indicators$maxcls 
# Get p-value for each taxon
pval <- indicators$pval 
# Get indicator values for each in its maximum class
indmax <- indicators$indcls 

# 3. Get frequency of occurence of indicators in associated assemblage ====

# Get relative frequency for each taxon in each cluster

relfreq <- indicators$relfrq %>% 
  `names<-`(sort(unique(colorGridData_top$name)))

# Get relative frequency for taxon in its max class

maxfreq <- apply(relfreq, 1, max) 

# 4. Create table of strong indicator taxa in each assemblage type ====

cl_index <- sort(unique(maxcls))
cl_name <- names(relfreq)
cases <- lapply(seq_along(cl_index), function(i){expr(maxcls == cl_index[!!i] ~ cl_name[!!i])})

# Create table with IndVal, p-value, and frequency 
# for strong indicators (IndVal > 0.25, p < 0.05)

IndicatorTable <- tibble(maxcls = indicators$maxcls, pval = indicators$pval, 
                         indmax = indicators$indcls, maxfreq, taxon = names(maxfreq)) %>% 
  filter(pval <= 0.05) %>% 
  filter(indmax >= 0.25) %>% 
  mutate(taxon = gsub("[.]", " ", taxon),
         maxcls = case_when(!!!cases)) %>% 
  arrange(maxcls, desc(indmax)) 

# Match common names with taxa latin names

trawl_data <- readRDS(file = paste0("Data/TrawlData_", region, ".rds"))
TaxaIndex <- trawl_data %>% 
  dplyr::select(species, common_name) %>% 
  distinct()

for (i in 1:length(IndicatorTable$taxon)){
  IndicatorTable$Common_name[i] <- TaxaIndex[agrep(IndicatorTable$taxon[i], TaxaIndex$species),2]
  IndicatorTable$taxon[i] <- TaxaIndex[agrep(IndicatorTable$taxon[i], TaxaIndex$species),1]}

#  Reformat table

IndicatorTable <- IndicatorTable %>% 
  mutate(indmax = round(indmax, 3),
         maxfreq = round(maxfreq*100), 1) %>%
  dplyr::select(maxcls, taxon, Common_name, maxfreq, indmax) %>% 
  `colnames<-`(., c("Assemblage", "Taxa", "Common name", "Freq (%)", "IndVal"))

# Write table to text file 

write.table(IndicatorTable, file = paste0("Output/IndicatorSpecies_", region, ".txt"),
            row.names = F, sep = ",", qmethod = "escape")

# 5. Figure Prep ====

# Fig. S4 ----

# Matrix plot for Indicator Taxa

pal <- bind_rows(GULF.palette, MAR.palette, NL.palette, QC.palette, .id = "Region") %>% 
  mutate(Region = case_when(
    Region == 1 ~ "SGSL",
    Region == 2 ~ "MAR",
    Region == 3 ~ "NL",
    Region == 4 ~ "NGSL"
  ))

# Read in tables with IndVal and Freq

indval.tbl <- list.files("Output/", pattern = "Indi.+txt$", full.names = TRUE) %>%
  map(fread) %>% # read in tables 
  set_names(c("MAR", "NGSL", "NL", "SGSL")) %>% # name list elements
  rbindlist(., idcol = "Region") %>% # row bind & add id column for region
  dplyr::select(-starts_with("Freq")) %>% # remove Frequency
  inner_join(., pal, 
             by = c("Region","Assemblage" = "name")) %>% # join colour palettes
  # fix inconsistent taxa spellings
  mutate(Taxa = gsub("^Gorgonocephalus sp $", "Gorgonocephalus spp.", Taxa),
         Taxa = gsub("sp $", "spp.", Taxa),
         Taxa = gsub("c $", "", Taxa),
         Taxa = gsub("o $", "", Taxa),
         Taxa = gsub("o\\.", "", Taxa),
         Taxa = gsub("f\\.", "", Taxa),
         Taxa = gsub(" f $", "", Taxa),
         Taxa = gsub("c\\.", "", Taxa),
         Taxa = gsub("sp\\.", "spp.", Taxa)) %>%
  mutate(Dupl = duplicated(Taxa)) %>% # Identify shared indicator taxa
  mutate(Region = factor(Region, # re-order factor levels for Region
                         levels = c("NGSL","SGSL","NL","MAR"))) %>% 
  arrange(Taxa) # arrange by species alphabetically # 117 taxa across all regions

# Create list of shared indicator taxa
dupl <- filter(indval.tbl, Dupl == TRUE) %>% 
  pull(Taxa) %>% 
  unique() # 48 unique taxa shared by 2+ regions

# Matrix with indicator taxa shared by at least 2 regions

mat.ind <- ggplot(filter(indval.tbl, Taxa %in% dupl),
                  aes(x = Region, 
                      y = reorder(Taxa, desc(Taxa)),
                      fill = assigned)) +
  geom_tile(width = 1, 
            height = 1, 
            colour = "grey50", 
            size = 0.5) +
  scale_fill_identity(guide = "none") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  labs(y = "Indicator taxa") +
  theme(panel.background = element_rect(fill = "grey85", colour = "grey50"),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "italic")); mat.ind

ggsave(plot = mat.ind, 
       filename = "Output/IndVal_MatrixPlot_SharedTaxa.tiff", 
       width = 5.5, 
       height = 12, 
       compression = "lzw")

cat("Proceed to 4RandomForestClassification.R")