# C Library for rocker: apt-get -y install  libz-dev
# apt-get install libglpk-dev

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("DESeq2")
BiocManager::install("AnnotationForge")
BiocManager::install("DEGreport")
BiocManager::install("tximport")
BiocManager::install("rhdf5")




install.packages("ggplot2")
install.packages("reshape2")
install.packages("tidyr")