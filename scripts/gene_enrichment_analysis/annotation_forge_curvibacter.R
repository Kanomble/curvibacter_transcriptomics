#This script reads genome and annotation information from the output of the translation_table.ipynb notebook code.
#It creates an OrgPackage for Curvibacter sp. AEP1-3, that can be used for the enrichGO function of the clusterProfiler
#R package.
library(AnnotationForge)

fSym <- read.csv("data/curvibacter_annotation_files/curvibacter_info_annotation_forge.csv", sep="\t")
fGO <- read.csv("data/curvibacter_annotation_files/curvibacter_go_annotation_forge.csv", sep="\t")
fChr <- read.csv("data/curvibacter_annotation_files/curvibacter_chromosome_annotation_forge.csv", sep="\t")

fSym <- subset(fSym, GID != "none")
fChr <- subset(fChr, GID != "none")


makeOrgPackage(gene_info=fSym, chromosome=fChr, go=fGO,
               version="0.1",
               maintainer="Lukas Becker <lukas.becker@hhu.de>",
               author="Lukas Becker <lukas.becker@hhu.de>",
               outputDir = ".",
               tax_id="1844971",
               genus="Curvibacter",
               species="AEP13",
               goTable="go")

install.packages("org.CAEP13.eg.db", repos=NULL, type="source", unlink=TRUE)
