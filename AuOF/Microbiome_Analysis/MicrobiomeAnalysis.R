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
##############################################################################################################################
#                                               SETUP AND PREPARE DATA                                                       #
##############################################################################################################################
# biom file for testing. REMOVE ONE COMPLETE
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
# quickview of the data
head(rowData(CMCData))
# change the column data and remove unnecessary preamble characters
names(rowData(CMCData)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
modifiedRowData <- BiocParallel::bplapply(rowData(CMCData), 
                                        FUN = stringr::str_remove,
                                        pattern = '.*[kpcofgs]__')
# remove any \ in the data
modifiedRowData <- BiocParallel::bplapply(modifiedRowData, FUN = stringr::str_remove, pattern = '\"')
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
# convert data object to a phyloseq object
CMCData <- makePhyloseqFromTreeSummarizedExperiment(CMCData)
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
outgroup <- pick_new_outgroup(phy_tree(CMCData))
cat("The randomly chosen outgroup is ", outgroup, "\n")
# root the tree
phy_tree(CMCData) <- ape::root(phy_tree(CMCData), outgroup = outgroup, resolve.root = TRUE)
cat("Taxonomic tree has been successfully rooted using the outgroup:", outgroup, "\n")
# check if there are empty features
summary(CMCData@tax_table@.Data=="")
# remove empty, uncharacterised or N.A entries in the taxa features
CMCData <- subset_taxa(CMCData, !is.na(Phylum) & !Phylum %in% c("", "uncharacterised", "uncharacterized"))
CMCData <- subset_taxa(CMCData, !is.na(Class) & !Class %in% c("", "uncharacterised", "uncharacterized"))
CMCData <- subset_taxa(CMCData, !is.na(Family) & !Family %in% c("", "uncharacterised", "uncharacterized"))
CMCData <- subset_taxa(CMCData, !is.na(Order) & !Order %in% c("", "uncharacterised", "uncharacterized"))
CMCData <- subset_taxa(CMCData, !is.na(Species) & !Species %in% c("", "uncharacterised", "uncharacterized"))
# check that they are all false
summary(CMCData@tax_table@.Data=="")
##############################################################################################################################
#                                                      SIMPLE METRICS                                                        #            
##############################################################################################################################
# fix the names on the rownames on the taxa and otu table
taxx <- paste(CMCData@tax_table@.Data[,6], CMCData@tax_table@.Data[,7], sep=" ")
# Get the number of reads for each sample
sample_sums(CMCData)
# add the number of reads back to  the sample metadata in the data object
sample_data(CMCData)$Read_Count <- sample_sums(CMCData)
# Get the total number of reads for each taxa feature across all samples
taxa_sums(CMCData)
# save it to file
taxSums <- data.frame(taxa_sums(CMCData))
rownames(taxSums) <- taxx
colnames(taxSums) <- "Total Number of reads for each species across all samples"
write.csv(taxSums, "Reads_per_taxa.csv", row.names=TRUE)
# get the number of reads per sample per taxa feature
otuTab <- as.data.frame(otu_table(CMCData), row.names = taxx)
write.csv(otuTab, "Reads_per_sample_per_taxa.csv", row.names = TRUE)
# get the number of identified features
table(tax_table(CMCData)[,"Kingdom"])
table(tax_table(CMCData)[,"Phylum"])
# save to file # the frequency is the number of taxa in that kingdom/phylum
write.csv(table(data.frame(tax_table(CMCData)[,"Kingdom"])), "KingdomFeatureCount.csv", row.names = FALSE) 
write.csv(table(data.frame(tax_table(CMCData)[,"Phylum"])), "PhylumFeatureCount.csv", row.names = FALSE)
# Plot absolute abundance of taxa
(plot1 <- plot_bar(CMCData, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Absolute Abundance\n (Reads)\n") +
  facet_wrap(~Location, scales = "free", nrow = 1) +
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(angle = 0, hjust=0.5)))
# save to file
tiff("Taxonomy_AbsoluteAbundance.tif", width=3000, height=2000, units = "px", res=300)
print(plot1)
dev.off()
cat("The total number of reads for each identified taxa across all samples are saved in Reads_per_taxa.csv \n")
cat("The total number of reads for each identified taxa in each sample are saved in Reads_per_sample_per_taxa.csv \n")
cat("The total number of taxa/species in the Kingdom taxonomic classification strata are saved in KingdomFeatureCount.csv \n")
cat("The total number of taxa/species in the Phylum taxonomic classification strata are saved in PhylumFeatureCount.csv \n")
cat("The absolute abundance of the samples using read counts is plotted in Taxonomy_AbsoluteAbundance.tif \n")
##############################################################################################################################
#                                    COMPUTE PREVALENCE AND FILTER IF POSSIBLE                                               #
##############################################################################################################################
# Compute the prevalence of each feature across the samples
CMCDataPrevalence <- apply(X = otu_table(CMCData), 
                            MARGIN = ifelse(taxa_are_rows(CMCData), yes = 1, no = 2),
                            FUN = function(x){sum(x > 0)})
# convert to dataframe
CMCDataPrevalence <- data.frame(Prevalence = CMCDataPrevalence,
                                  TotalAbundance = taxa_sums(CMCData),
                                  tax_table(CMCData))
# bind the taxonomy to the prevalence
CMCDataPrevalence <- cbind(Taxa = taxx, CMCDataPrevalence)
#rownames(CMCDataPrevalence) <- 1:nrow(CMCDataPrevalence)
# save the prevalence information
write.csv(CMCDataPrevalence, "Prevalence_of_all_Taxa.csv", row.names = FALSE)
# check for phyla of low prevalence features by computing the average % prevalence
# and total prevalence of the phyla across all samples
CMCDataPrevalence_Avg <- plyr::ddply(CMCDataPrevalence, "Phylum", function(df){cbind(mean(df$Prevalence), sum(df$Prevalence))})
colnames(CMCDataPrevalence_Avg) <- c("Phylum", "Avg.Prevalence", "Total.Prevalence")
# write to file. This is showing the average number of samples each phylum appears in (Avg.Prevalence) and the total number of
# features associated with that phylum for all samples it appears in (Total Prevalence)
write.csv(CMCDataPrevalence_Avg, "Average_Prevalence.csv", row.names = FALSE)
# plot prevalence and total abundance i.e. this plot show the number of reads and the number of samples found that support each
# identified phylum
(plot2 <- ggplot(CMCDataPrevalence, aes(TotalAbundance, Prevalence/nsamples(CMCData), color = Phylum)) + 
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.3) +
  scale_x_log10() + xlab("Total Abundance/Read count") + ylab("Proportional Prevalence") +
  facet_wrap(~Phylum) + theme(legend.position = "none") +
  labs(caption = "Total number of reads (Total Abundance) plotted against the proportion of samples that support each identified phylum. 
  As such, if a phylum is supported in all samples, the prevalence will be 1.
  Proportional prevalence is calculated by dividing the prevalence of each phylum by the total number of samples.",
       subtitle="Proportional prevalence vs Read count plotted for identified Phyla"))
tiff("Prevalence_TotalAbundance.tif", height = 3200, width = 4500, res = 300, units = "px")
print(plot2)
dev.off()
# filter out low prevalent taxonomic features
keepTaxa = rownames(CMCDataPrevalence)[(CMCDataPrevalence$Prevalence >= (0.05*nsamples(CMCData)))]
if (length(keepTaxa) == length(rownames(CMCDataPrevalence))) {
  cat("No taxonomic features were filtered because all features are supported at least", 0.05*nsamples(CMCData), "of samples.\n")
  cat("This means that all features pass the (0.05 * sample_count) threshold \n")
  CMCData_Filtered = CMCData
} else if (length(keepTaxa) < length(rownames(CMCDataPrevalence))) {
  # filter out the low prevalent taxa
  CMCData_Filtered = prune_taxa(keepTaxa, CMCData)
  # compute prevalence of the filtered taxonomies
  CMCDataPrevalence_Filt = apply(X = otu_table(CMCData_Filtered),
                                  MARGIN = ifelse(taxa_are_rows(CMCData_Filtered), yes = 1, no =2),
                                  FUN = function(x){sum(x > 0)})
  # convert to dataframe
  CMCDataPrevalence_Filt <- data.frame(Prevalence = CMCDataPrevalence_Filt,
                                    TotalAbundance = taxa_sums(CMCData_Filtered),
                                    tax_table(CMCData_Filtered))
  # bind the taxonomy to the prevalence
  filt_taxx <- paste(CMCData_Filtered@tax_table@.Data[,6], CMCData_Filtered@tax_table@.Data[,7], sep=" ")
  CMCDataPrevalence_Filt <- cbind(Taxa = filt_taxx, CMCDataPrevalence_Filt)
  #rownames(CMCDataPrevalence_Filt) <- 1:nrow(CMCDataPrevalence_Filt)
  # save the prevalence information
  write.csv(CMCDataPrevalence_Filt, "Filtered_Prevalence_of_all_Taxa.csv", row.names = FALSE)
  # plot the abundance of the filtered taxa
  (plot3 <- ggplot(CMCDataPrevalence_Filt, aes(TotalAbundance, Prevalence/nsamples(CMCData_Filtered), color = Phylum)) + 
          geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.3) +
          scale_x_log10() + xlab("Total Abundance/Read count") + ylab("Proportional Prevalence") +
          facet_wrap(~Phylum) + theme(legend.position = "none") +
          labs(caption = "Total number of reads (Total Abundance) plotted against the proportion of samples that support each identified phylum after filtering for prevalence. 
    Prevalence filtering threshold was calculated by excluding any taxonomic features supported by <5% of samples. 
    After filtering, if a phylum is supported in all remaining samples, the prevalence will be 1.
    Proportional prevalence is calculated by dividing the prevalence of each phylum by the total number of samples.",
               subtitle="Filtered prevalence vs Read count plotted for identified Phyla"))
  tiff("Filtered_Prevalence_TotalAbundance.tif", height = 3200, width = 4500, res = 300, units = "px")
  print(plot3)
  dev.off()
}
##############################################################################################################################
#                                               COMPUTE RAW RELATIVE ABUNDANCE                                               #
##############################################################################################################################
# calculate and add the relative abundance to the object
CMCData_RelAbund <- transform_sample_counts(CMCData_Filtered, function(x) x/sum(x))
# write to file
write.csv(table(CMCData_RelAbund@otu_table), "Relative_Abundance.csv", row.names = FALSE)
# plot the dataset using phyla to colour the stacked bar chart
(plot4 <- plot_bar(CMCData_RelAbund, fill = "Phylum") + 
    geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") + 
    labs(x = "", y = "Relative Abundance\n") +
    facet_wrap(~Location, scale = "free") +
    theme(legend.position = "bottom") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)))
tiff("Raw_Taxonomy_RelativeAbundance.tif", width = 4500, height = 3200, units = "px", res = 300)
print(plot4)
dev.off()
##############################################################################################################################
#                                               COMPUTE ALPHA DIVERSITY                                                      #
##############################################################################################################################
# Calculate alpha diversity using Shannon and Simpson indices
# normalise reads to an even depth
CMCData_RelAbund_Normal <- rarefy_even_depth(CMCData, rngseed = 12355, replace = FALSE)
# calculate alpha diversity
(alpha_div <- estimate_richness(CMCData_RelAbund_Normal))
write.csv(alpha_div, "AlphaDiversity.csv", row.names = TRUE)
(plot5 <- plot_richness(CMCData_RelAbund_Normal, measures = c("Observed", "Shannon", "Chao1", "Simpson"), color = "Sample") +
    geom_point(aes(color = Sample)) +
    labs(color = "Sample"))
tiff("AlphaDiversity.tif", width = 4500, height = 3200, units = "px", res = 300)
print(plot5)
dev.off()
