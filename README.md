# Northwest Atlantic Biological Classification

---

## Project Overview

The aim of the __NW-Atl-bioclassification__ project is to develop a __biologically informed__ and __predictive__ classification scheme that delineates the distribution of the predominant assemblages of demersal fish and benthic invertebrates in the Northwest Atlantic. The subdivisions of the classification scheme should define spatially contiguous units of relatively homogeneous community composition (assemblages) along important axes of environmental variation and identify emblematic taxa of each assemblage type. This project will support __planning__ and __monitoring__ of bioregional __conservation networks__. 

Ecological representativity and connectivity are cornerstones of effective MPA network design. The National Framework for Canada's Network of MPAs outlines these design objectives as “facilitating the protection of ecological processes essential for ecosystem functioning”. Ecological processes in this sense reflect the summation and interconnectivity of species and their function within the ecosystem. An ecological classification system that partitions the seascape into homogenous spatial units, based on environmental and biological variables, can be used as a basis to evaluate representativity targets and select model species. 

MPA network design is challenged by the very species diversity it aims to protect. Indeed, the magnitude of biological complexity characterizing marine systems renders incorporating biological attributes of all species intractable. By categorizing this biological diversity hierarchically, and delineating it spatially, we will provide a novel, ecosystem-based tool to achieve representativity targets when designing MPAs. Moreover, we will provide pragmatic approach to identify key species on which to base design and monitoring (i.e., connectivity, size and spacing). 

By integrating taxa-specific observational information and environmental correlates of taxa presence we will delineate taxonomic assemblages and use environmental assignment models (random forest) to project groupings spatially as contiguous ecological classifications. From these classifications, we will identify emblematic species for more detailed connectivity analyses and future monitoring. Using climatological forecasts, we will compare current and future spatial classifications to characterize the risk of climate associated change. 


---

## Code Structure

The R code is intended for classification of assemblage structure from annual multispecies bottom trawl surveys in 4 regions of the Northwest Atlantic:

1. Northern Gulf of St. Lawrence (NGSL)
2. Southern Gulf of St. Lawrence (SGSL)
3. Newfoundland and Labrador (NL)
4. Maritimes (MAR)

The analysis is repeated in each region and the workflow proceeds as follows:

- _1DataPrep.R_: 
	- Defines study area boundaries
	- Creates fishnet grid for aggregation of trawl data
	- Filters and reshapes trawl data
	- Converts abundance to presence-absence
	- Exports Site X Species Matrix
- _2ClusterAnalysis.R_:
       - Select dissimilarity measure and clustering algorithm
       - Hierarchical clustering of Site X Species Matrix
       - Select threshold to slice dendrogram
       - Identify and map distribution of major cluster (assemblages)
- _3IndicatorAnalysis.R_:
	- Identifies emblematic taxa of major assemblage types
- _4RandomForestClassification.R_:
	- Trains and evaluates Random Forest Classifier
	- Characterizes environmental correlates of assemblage distribution
	- Predictively map assemblage distributions
	- Identifies areas of higher uncertainty
- _5ClimateVulnerability.R_:
- Hindcasts and forecasts climate-associated changes in distribution
- Identifies areas and assemblages vulnerable to warming bottom temperatures
- _6AncillaryAnalyses.R_:
	- Other supporting analyses

_BioclassificationFunctions.R_ and _ColourPalettes.R_ are called throughout the preceding 6 code blocks for analyses and figure preparation.

---

## Citation

The code in this repository was used to produce the results and figures in the associated publication:

O'Brien JM, Stanley RRE, Jeffery NW, Heaslip SG, DiBacco C, Wang Z 	(accepted) Modelling demersal fish and benthic invertebrate 	assemblages in support of marine conservation planning. _Ecol Appl_

A Zenodo DOI is also available for the most recent release of __ NW-Atl-bioclassification__:

---

## Contact:

Corresponding Developer: 

John O’Brien 
<https://github.com/obrienjm25>
<John.Obrien@dfo-mpo.gc.ca>

---



