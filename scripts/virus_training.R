# my firsat R scrtipt

# library(ggplot2) <=> require
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("Bioconductor")
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
library(GenomicFeatures)
library(AnnotationDbi) ## load library
library(tidyverse)
library(AMR)
library(dplyr)

########### Exercise 05 #######
###############################

