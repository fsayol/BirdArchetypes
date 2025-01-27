## Set up script ## Loading of packages, functions and working directories.

## TO BE SET BY THE USER:
# Define main working directory (Repository folder):
RepoFolder <-"/Users/fsayol/BirdArchetypes/"
# TODO: Delete:
RepoFolder <-"/Users/fsayol/My Drive/2_ResearchProjects/A3_PigotLab/P1_Archetypes/Repository_DataCode/BirdArchetypes_v2/"

# Load packages:
library(archetypes)
library(ape)
library(betareg)
library(circular)
library(cxhull)
library(classInt)
library(cowplot)
library(DescTools)
library(dplyr)
library(entropy)
library(gdata)
library(geiger)
library(geometry)
library(ggplot2)
library(grDevices)
library(hexbin)
library(ks)
library(MASS) 
library(matrixStats)
library(mgcv)
library(MCMCglmm)
library(phytools)
library(plot3D)
library(pracma)
library(proxy)
library(quantreg)
library(raster)
library(RColorBrewer)
library(readxl)
library(REdaS)
library(reshape2)
library(rgl)
library(scales)
library(SimDesign)
library(sf)
library(sp)
library(stringr)
library(vegan)

# Define directories within wd (Relative to main WD, do not require to change)
DataFolder <- paste0(RepoFolder,"data/")
DataTmpFolder <- paste0(RepoFolder,"data/data_tmp/")
OutputFolder <- paste0(RepoFolder,"output/")
FigureFolder <- paste0(RepoFolder,"figures/")

# Create data_tmp, output and figure directories if they do not exist:
if(!file.exists(DataTmpFolder)) dir.create(DataTmpFolder)
if(!file.exists(OutputFolder)) dir.create(OutputFolder)
if(!file.exists(FigureFolder)) dir.create(FigureFolder)

# Load functions:
source(paste0(RepoFolder,"scripts/ArchetypeFunctions.R"))

## End Set up script ##