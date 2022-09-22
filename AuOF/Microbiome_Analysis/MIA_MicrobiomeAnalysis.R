##############################################################################################################################
# A script written as part of the Antibiotics Under Our Feet project
# The script carries out some microbiome analyses using standard R packages
# developed for metagenomic work. 
# Author: Damilola R Oresegun
#
##############################################################################################################################
#                                                        PREAMBLE                                                            #
##############################################################################################################################
# ensure clear workspace
rm(list = ls())
# set the arguments
args = commandArgs(trailingOnly = TRUE)
# ensure at least one argument
#if (length(args) == 0) {
#    stop("At least three arguments are needed. The biom file, the output directory and sample metadata")
#} else if (length(args) < 3) {
#    # use the working directory as the output folder
#    stop("At least three arguments are needed.")
#}
# install and set libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# install the necessary packages
#BiocManager::install(c("microbiome/mia", "devtools","tidyr", "remotes",
#                        "tidyverse","dplyr","cowplot","knitr","phyloseq",
#                        "vegan","ape", "ggplot2","ggplotify", "magrittr",
#                        "data.table"))
# load the libraries
library(mia)
library(dplyr)
library(knitr)
library(phyloseq)
library(vegan)
library(tidyr)
library(magrittr)
library(ape)
library(tidyverse)
library(data.table)
library(miaViz)
library(cowplot)
library(scater)
##############################################################################################################################
#                                               SETUP AND PREPARE DATA                                                       #
##############################################################################################################################
# biom file for testing. REMOVE ONCE COMPLETE
args= c("C:/Users/dro/OneDrive - University of St Andrews/MacKenzie_Institute/BiomFiles/CombinedIsolate.biom", 
        "C:/Users/dro/Dropbox/Work/MacLaptop/CMC_SoilProject/AuOF_microbiome/AuOF/Microbiome_Analysis", 
        "C:/Users/dro/Dropbox/Work/MacLaptop/CMC_SoilProject/AuOF_microbiome/AuOF/Microbiome_Analysis/sample_metadata.csv")
#
# set the argument variables
biomData = args[1]
workDir = args[2]
setwd(workDir)
# load the biom data
CMCData <- loadFromBiom(biomData)
# change the column data and remove unnecessary preamble characters
names(rowData(CMCData)) <- c("Kingdom", "Phylum", "Class", "Order", 
                                  "Family", "Genus", "Species")
modifiedRowData <- BiocParallel::bplapply(rowData(CMCData),
                                      FUN = stringr::str_remove,
                                      pattern = '.*[kpcofgs]__')
# remove any \ in the data
modifiedRowData <- BiocParallel::bplapply(modifiedRowData, 
                                      FUN = stringr::str_remove, 
                                      pattern = '\"')
modifiedRowData <- BiocParallel::bplapply(modifiedRowData, 
                                          FUN = stringr::str_remove, 
                                          pattern = '"""')
# make the list into a dataframe
modifiedRowData <- DataFrame(modifiedRowData)
# see the output 
head(modifiedRowData)
# add back to the data object
rowData(CMCData) <- modifiedRowData
# see the new row data
head(rowData(CMCData))
# add sample metadata to the data object
sampleMeta <- DataFrame(read.table(args[3], sep = ",", header = TRUE))
cat("Sample metadata has been added to the data object \n")
# add the sample names as the row names
rownames(sampleMeta) <- sampleMeta[,1]
# add to the data object
colData(CMCData) <- sampleMeta
# add a taxonomic tree -- we do not have a phylogenetic tree
CMCData <- addTaxonomyTree(CMCData)
cat("Taxonomic tree has been generated from the taxonomic features\n")

# define a function to select an outgroup to root the tree
pick_new_outgroup <- function(tree.unrooted){
  # tablify parts of tree that is needed
  treeDT <- cbind(data.table(tree.unrooted$edge),
                  data.table(length = tree.unrooted$edge.length)
  )[1:Ntip(tree.unrooted)] %>% 
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)
}
# choose the root
outgroup <- pick_new_outgroup(rowTree(CMCData))
cat("The randomly chosen outgroup is ", outgroup, "\n")
# root the tree
rowTree(CMCData) <- ape::root(rowTree(CMCData), 
                              outgroup = outgroup, resolve.root = TRUE)
# calculate relative abundance and add to the data
CMCData <- transformSamples(x = CMCData, abund_values = "counts",
                                 method = "relabundance", name = "Abundance")
# save the abundances and taxa to file
write.csv(meltAssay(CMCData, add_row_data = TRUE,
                    add_col_data = TRUE,
                    abund_values = "Abundance"), "All_Classification.csv",
          row.names = FALSE)
# list the number of reads for each sample
sampleReads = as.data.frame(colSums(assay(CMCData, "counts")))
colnames(sampleReads) = c("Read Count")
write.csv(sampleReads, "Reads_per_sample_per_taxa.csv", row.names = TRUE)
# plot the abundance of the whole dataset
AllAbund1 <- plotAbundanceDensity(CMCData, layout = "jitter", abund_values = "Abundance",
                                  n = as.integer(length(CMCData)), point_size=1, 
                                  point_shape=19, point_alpha=0.7, colour_by="Sample") +
  theme(legend.position = "bottom") +
  scale_x_log10(label=scales::percent)
AllAbund2 <- plotAbundanceDensity(CMCData, layout = "jitter", abund_values = "Abundance",
                                  n = as.integer(length(CMCData)), point_size=1, 
                                  point_shape=19, point_alpha=0.7, colour_by="Location") +
  theme(legend.position = "bottom") +
  scale_x_log10(label=scales::percent)
tiff("Abundance_perSample_Location.tif",width = 4500, height = 3200, res=300, units="px")
png("Abundance_perSample_Location.png",width = 5400, height = 3200, res=300, units="px")
plot_grid(AllAbund1, AllAbund2, labels = c("A", "B"))
dev.off()
# plot abundance density of top 40 species
tiff("Top40_Abundance_perSample_Location.tif",width = 4500, height = 3200, res=300, units="px")
png("Top40_Abundance_perSample_Location.png",width = 4500, height = 3200, res=300, units="px")
plotAbundanceDensity(CMCData, layout = "point", abund_values="Abundance",
                     n = 40, colour_by="Sample", point_alpha=0.9) + 
  scale_x_log10(label=scales::percent) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()
# Plot the abundance of phyla
tiff("Abundance_perPhylum.tif",width = 4500, height = 3200, res=300, units="px")
png("AllData_perPhylum.png",width = 4500, height = 3200, res=300, units="px")
plotAbundance(CMCData, abund_values="Abundance", rank = "Phylum",
              add_x_text = TRUE) +
  theme(legend.key.height = unit(0.2, "cm")) +
  theme(legend.title = element_text(size=8)) +
  theme(legend.text = element_text(size=8)) +
  scale_y_continuous(label = scales::percent) +
  theme(axis.text.x = element_text(angle = 30, hjust=1, size=7)) +
  theme(legend.position = "bottom")
dev.off()
### Calculate prevalence
## prevalence is "How many/how often does this taxa show in the samples"
## detection is the threshold for absence/presence 
# get the top taxa prevalent for each subset
top40Alltaxa <- getTopTaxa(CMCData, rank="Phylum", method = c("median"), 
                           top=40, sort = FALSE)
top40Alltaxa <- rowData(CMCData)[top40Alltaxa, taxonomyRanks(CMCData)]
write.csv(top40Alltaxa, "Top40_Prevalent_Species.csv",
          row.names = FALSE)
top40Alltaxa <- as.data.frame(top40Alltaxa)
View(top40Alltaxa)

cllo <- CMCData[rowData(CMCData)$Kingdom %in% "Eukaryota"& !is.na(rowData(CMCData)$Kingdom), ]
cllo

unique(rowData(CMCData)$Kingdom)
