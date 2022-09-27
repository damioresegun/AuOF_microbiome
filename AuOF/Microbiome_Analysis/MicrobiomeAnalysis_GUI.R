##############################################################################################################################
# A script written as part of the Antibiotics Under Our Feet project
# The script carries out some microbiome analyses using standard R packages
# developed for metagenomic work. 
# NOTE: This is a script is designed to be run in a GUI not via the command line!
# Author: Damilola R Oresegun
# Usage: In the args variable below, enter the full path to the biom file, the output folder, the sample metadata and the file
#         format to save the output images in. Options for output format are png or tiff
# Example: args=c("my/path/to/CombinedIsolate.biom", 
#                  "my/path/to/the/output/folder",
#                   "my/path/to/my/sample_metadata.csv",
#                    "png")
# Outputs: Several plots and tables holding the composition,abundance, alpha and beta diversity of the input microbiome 
#           community.
##############################################################################################################################
#                                                        PREAMBLE                                                            #
##############################################################################################################################
# ensure clear workspace
rm(list = ls())

# set the arguments
args= c("path/to/CombinedIsolate.biom", 
        "path/to/output/folder", 
        "path/to/sample_metadata.csv",
        "tiff") # Please change this!


# install and set libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# check the r version
BiocManager::install(c("microbiome/mia", "devtools","tidyr", "remotes",
                      "tidyverse","dplyr","cowplot","knitr","phyloseq",
                      "vegan","ape", "ggplot2","ggplotify", "magrittr",
                      "data.table","scater","plyr"))
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
library(scater)
library(cowplot)
library(plyr)
##############################################################################################################################
#                                               SETUP AND PREPARE DATA                                                       #
##############################################################################################################################
# set the argument variables
biomData = args[1]
workDir = args[2]
outImg = args[4]
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
outgroupID <- pick_new_outgroup(phy_tree(CMCData))
taxtablee <- as.data.frame(tax_table(CMCData))
taxtablee <- tibble::rownames_to_column(taxtablee, "ID")
outgroup = paste((taxtablee[grep(outgroupID, taxtablee$ID),][,7]),(taxtablee[grep(outgroupID, taxtablee$ID),])[,8], sep = " ")
cat("The randomly chosen outgroup is ", outgroup, "\n")
# root the tree
phy_tree(CMCData) <- ape::root(phy_tree(CMCData), outgroup = outgroupID, resolve.root = TRUE)
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
pol = sample_sums(CMCData)
cat("The number of reads for", sample_names(CMCData), "are:", sample_sums(CMCData), "respectively\n")
# add the number of reads back to  the sample metadata in the data object
sample_data(CMCData)$Read_Count <- sample_sums(CMCData)
cat("The number of reads for each sample has been added to the data object.\n")
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
cat("The number of features in each kingdom is:\n")
table(tax_table(CMCData)[,"Kingdom"])
cat("The number of features in each phylum is:\n")
table(tax_table(CMCData)[,"Phylum"], exclude = NULL)
# save to file # the frequency is the number of taxa in that kingdom/phylum
write.csv(table(data.frame(tax_table(CMCData)[,"Kingdom"])), "KingdomFeatureCount.csv", row.names = FALSE) 
write.csv(table(data.frame(tax_table(CMCData)[,"Phylum"])), "PhylumFeatureCount.csv", row.names = FALSE)
# Plot absolute abundance of taxa
(plot1 <- plot_bar(CMCData, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Absolute Abundance\n (Reads)\n") +
  facet_wrap(~Location, scales = "free", nrow = 1) +
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(angle = 0, hjust=0.5))+
  theme(legend.text = element_text(size=6))+
  theme(legend.title = element_text(size=7)))
# save to file
if (outImg == "tiff"){
  tiff("Taxonomy_AbsoluteAbundance.tif", width=3000, height=2000, units = "px", res=300)
  print(plot1)
  dev.off()
} else if (outImg == "png") {
  png("Taxonomy_AbsoluteAbundance.png", width=3000, height=2000, units = "px", res=300)
  print(plot1)
  dev.off()
}
cat("The total number of reads for each identified taxa across all samples are saved in Reads_per_taxa.csv \n")
cat("The total number of reads for each identified taxa in each sample are saved in Reads_per_sample_per_taxa.csv \n")
cat("The total number of taxa/species in the Kingdom taxonomic classification strata are saved in KingdomFeatureCount.csv \n")
cat("The total number of taxa/species in the Phylum taxonomic classification strata are saved in PhylumFeatureCount.csv \n")
cat("The absolute abundance of the samples using read counts is plotted in Taxonomy_AbsoluteAbundance.png \n")
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
#rownames(CMCDataPrevalence) <- 1:nrow(CMCDataPrevalence) # works but is problematic downstream. So deactivated for now

# save the prevalence information
write.csv(CMCDataPrevalence, "Prevalence_of_all_Taxa.csv", row.names = FALSE)

cat("The prevalence of each identified taxa has been saved in 'Prevalence_of_all_Taxa.csv'.\n")
cat("Prevalence refers to the number of samples each taxa was identified in.\n")

# check for phyla of low prevalence features by computing the average % prevalence
# and total prevalence of the phyla across all samples
CMCDataPrevalence_Avg <- plyr::ddply(CMCDataPrevalence, "Phylum", function(df){cbind(mean(df$Prevalence), sum(df$Prevalence))})
colnames(CMCDataPrevalence_Avg) <- c("Phylum", "Avg.Prevalence", "Total.Prevalence")

# write to file. This is showing the average number of samples each phylum appears in (Avg.Prevalence) and the total number of
# features associated with that phylum for all samples it appears in (Total Prevalence)
write.csv(CMCDataPrevalence_Avg, "Average_Prevalence_perPhylum.csv", row.names = FALSE)

cat("The average prevalence per phylum has been saved to file.\n")

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
if (outImg == "tiff"){
  tiff("Prevalence_TotalAbundance.tif", width=4000, height=4000, units = "px", res=300)
  print(plot2)
  dev.off()
} else if (outImg == "png") {
  png("Prevalence_TotalAbundance.png", width=4000, height=4000, units = "px", res=300)
  print(plot2)
  dev.off()
}

# filter out low prevalent taxonomic features. Prevalence threshold is 10% of all samples
keepTaxa = rownames(CMCDataPrevalence)[(CMCDataPrevalence$Prevalence >= (0.1 * nsamples(CMCData)))]

if (length(keepTaxa) == length(rownames(CMCDataPrevalence))) {
  
  cat("No taxonomic features were filtered because all features are supported at least", 0.1*nsamples(CMCData), "of samples.\n")
  cat("This means that all features pass the (0.1 * sample_count) threshold \n")
  
  CMCData_Filtered = CMCData
  
} else if (length(keepTaxa) < length(rownames(CMCDataPrevalence))) {
  
  cat("Some taxonomic features are not supported by the threshold number of samples. \n These will be filtered out for low prevalence")
  
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
  
  cat("After filtering, the features that passed the threshold have been saved in Filtered_Prevalence_of_all_Taxa.csv.\n")
  cat("This means the features to be removed were supported by less than", 0.1*nsamples(CMCData), "of samples.\n")
  cat("The filtered dataset will be taken forward.\n")
  
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
  if (outImg == "tiff"){
    tiff("Filtered_Prevalence_TotalAbundance.tif", width=3000, height=2000, units = "px", res=300)
    print(plot3)
    dev.off()
  } else if (outImg == "png") {
    png("Filtered_Prevalence_TotalAbundance.png", width=3000, height=2000, units = "px", res=300)
    print(plot3)
    dev.off()
  }
}
##############################################################################################################################
#                                                   AGGLOMERATE TO GENUS LEVEL                                               #
##############################################################################################################################
# check how many taxa would remain after agglomeration filtering
cat("After agglomertation filtering, the identified taxa will be reduced to", length(get_taxa_unique(CMCData_Filtered, taxonomic.rank = "Genus")), "genus genera\n")
cat("After agglomertation filtering, the identified taxa will be reduced to", length(get_taxa_unique(CMCData_Filtered, taxonomic.rank = "Phylum")), "phylum genera\n")
# agglomerate data to the genus level
CMCData_Genus <- tax_glom(CMCData_Filtered, "Genus", NArm = TRUE)
CMCData_Phylum <- tax_glom(CMCData_Filtered, "Phylum", NArm = TRUE)
cat("The data has been collapsed to just the identified genuses. This is done as a means of reducing the amount of data we are processing\n")
cat("Note that the original identified species still remain, they have just been combined to be their parent genus\n")
# Plot the tree of the pre and post agglomeration data object
(plot4 = plot_tree(CMCData_Filtered, ladderize = "left", color = "Kingdom",
                   base.spacing=0.03, plot.margin=0.2,nodelabf=nodeplotblank) +
  theme(plot.title = element_text(size = 15))+
    labs(subtitle = "Pre-Agglomeration"))

(plot5 = plot_tree(CMCData_Genus, ladderize = "left", color = "Kingdom",
                   base.spacing=0.03, plot.margin=0.2,nodelabf=nodeplotblank) +
    theme(plot.title = element_text(size = 15))+
    labs(subtitle = "Post-Agglomeration by Genus"))

(plot6 = plot_tree(CMCData_Phylum, ladderize = "left", color ="samples", 
                   label.tips = "Phylum", shape = "Kingdom", size = "abundance",base.spacing=0.04, plot.margin=0.4,nodelabf=nodeplotblank) +
    theme(plot.title = element_text(size = 15))+
    labs(caption = "Visual representation of taxonomic trees of identified before and after agglomeration. [A] Pre agglomeration, 4062 species
                    were identified while after agglomerating to the Genus level [B], the species were collapsed to their constituent 1137 genuses.
                    A further agglomeration to the Phylum level [C] results in 33 phyla identified, with the absolute abundance (or read counts)
                    also shown in the size of the coloured labels. All trees are rooted using a randomly picked species in the identified taxa -- ",
         subtitle = "Post agglomeration by Phylum"))
if (outImg == "tiff"){
  tiff("Pre_and_Post_Agglomeration_ByGenus_byPhylum.tif", width=5000, height=4000, units = "px", res=300)
  combp <- plot_grid(plot4, plot5,plot6, labels = c("A", "B", "C"), label_y = 0.9)
  print(combp)
  dev.off()
} else if (outImg == "png") {
  png("Pre_and_Post_Agglomeration_ByGenus_byPhylum.png", width=5000, height=4000, units = "px", res=300)
  combp <- plot_grid(plot4, plot5,plot6, labels = c("A", "B", "C"), label_y = 0.9)
  print(combp)
  dev.off()
}

cat("The data has been agglomerated by Genus and this will be taken forward\n")

##############################################################################################################################
#                                               COMPUTE RAW RELATIVE ABUNDANCE                                               #
##############################################################################################################################
# define function to plot abundances
plot_abundance = function(physeq,title = "",
                          Facet = "Phylum", Color = "Kingdom"){
  # Arbitrary subset, based on Phylum, for plotting
  mphyseq = psmelt(physeq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Sample",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}

# calculate and add the relative abundance to the object
CMCData_RelAbund <- transform_sample_counts(CMCData_Genus, function(x) x/sum(x))

# plot the relative abundance of all phyla after agglomeration
(plot7 = plot_abundance(CMCData_Genus, ""))
(plot8 = plot_abundance(CMCData_RelAbund, ""))
if (outImg == "tiff"){
  tiff("Pre_and_Post_AbundanceValueTransformation.tif", height = 9000, width = 7500, res = 300, units = "px")
  combop <- plot_grid(nrow=2, plot7, plot8, labels = c("A", "B"))
  print(combop)
  dev.off()
  
  tiff("PerPhylum_RelativeAbundance.tif", height = 5000, width = 7500, res = 300, units = "px")
  print(plot8)
  dev.off()
} else if (outImg == "png") {
  png("Pre_and_Post_AbundanceValueTransformation.png", height = 9000, width = 7500, res = 300, units = "px")
  combop <- plot_grid(nrow=2, plot7, plot8, labels = c("A", "B"))
  print(combop)
  dev.off()
  
  png("PerPhylum_RelativeAbundance.png", height = 5000, width = 7500, res = 300, units = "px")
  print(plot8)
  dev.off()
}

# write to relative abundances to file
write.csv(table(CMCData_RelAbund@otu_table), "Relative_Abundance.csv", row.names = FALSE)
# plot the dataset using phyla to colour the stacked bar chart
(plot9 <- plot_bar(CMCData_RelAbund, fill = "Phylum") + 
    geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") + 
    labs(x = "", y = "Relative Abundance\n") +
    facet_wrap(~Location, scale = "free") +
    theme(legend.position = "bottom") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
    labs(caption = "Stacked plot of the calculated relative abundances of all samples after agglomeration to
         the genus level. Plot is further separated based on the location where the samples were extraced.
         Relative abundance was calculated using the equation x/sum(x) where x is the read count for each sample."))
if (outImg == "tiff"){
  tiff("Combined_PerPhylum_RelativeAbundance.tif", width = 4500, height = 3200, units = "px", res = 300)
  print(plot9)
  dev.off()
} else if (outImg == "png") {
  png("Combined_PerPhylum_RelativeAbundance.png", width = 4500, height = 3200, units = "px", res = 300)
  print(plot9)
  dev.off()
}
##############################################################################################################################
#                                               COMPUTE ALPHA DIVERSITY                                                      #
##############################################################################################################################
# Calculate alpha diversity using Shannon and Simpson indices
# normalise reads to an even depth
CMCData_RelAbund_Normal <- rarefy_even_depth(CMCData_Genus, rngseed = 12355, replace = FALSE)
# calculate alpha diversity
(alpha_div <- estimate_richness(CMCData_RelAbund_Normal))
write.csv(alpha_div, "AlphaDiversity.csv", row.names = TRUE)
(plot10 <- plot_richness(CMCData_RelAbund_Normal, 
                         measures = c("Observed", "Shannon", "Chao1", "Simpson"), 
                         color = "Location", shape = "Isolate") +
    geom_point(aes(color = Location)) +
    labs(color = "Location"))
if (outImg == "tiff"){
  tiff("Combined_PerPhylum_RelativeAbundance.tif", width = 4500, height = 3200, units = "px", res = 300)
  print(plot9)
  dev.off()
} else if (outImg == "png") {
  png("Combined_PerPhylum_RelativeAbundance.png", width = 4500, height = 3200, units = "px", res = 300)
  print(plot9)
  dev.off()
}
if (outImg == "tiff"){
  tiff("AlphaDiversity.tif", width = 4500, height = 3200, units = "px", res = 300)
  print(plot10)
  dev.off()
} else if (outImg == "png") {
  png("AlphaDiversity.png", width = 4500, height = 3200, units = "px", res = 300)
  print(plot10)
  dev.off()
}

##############################################################################################################################
#                                               COMPUTE BETA DIVERSITY/ORDINATION                                            #
##############################################################################################################################
# calculate ordination using different methods
dist_methods <- unlist(distanceMethodList)
# no phylogenetic tree is in the data so methods needing phylogenetic trees are removed
dist_methods <- dist_methods[-(1:2)]
# remove user defined method
dist_methods <- dist_methods[-which(dist_methods=="ANY")]
# loop through the distance method and save each plot to a list for later use
plotList <- vector("list", length(dist_methods))
names(plotList) = dist_methods
for (mthd in dist_methods){
  # calculate the distance matrix
  mDist <- phyloseq::distance(CMCData_RelAbund, method = mthd)
  # calculate ordination using MDS/PCoA
  mPcOA <- ordinate(CMCData_RelAbund, "MDS", distance = mDist)
  # make the plot
  dplot <- NULL # make sure the previous plot is cleared
  (dplot <- plot_ordination(CMCData_RelAbund, mPcOA, color = "Location", shape = "Isolate") +
    coord_fixed(sqrt(mPcOA$values$Eigenvalues[2]/mPcOA$values$Eigenvalues[1])) +
    ggtitle(paste("MDS/PCoA using distance method ", mthd, sep="")))
  # save the graphic to file
  plotList[[mthd]] = dplot
}
# show all plots and save to file
plotDF = ldply(plotList, function(x) x$data)
names(plotDF)[1] <- "distance"
(plot11 = ggplot(plotDF, aes(Axis.1, Axis.2, color=Location, shape=Isolate)) +
  geom_point(size=3, alpha=0.5) +
  facet_wrap(~distance, scales="free") +
  ggtitle("MDS/PCoA on various distance metrics"))
# print to file
if (outImg == "tiff"){
  tiff("BetaDiversity_multiMethod.tif", width = 4500, height = 3200, units = "px", res = 300)
  print(plot11)
  dev.off()
} else if (outImg == "png") {
  png("BetaDiversity_multiMethod.png", width = 4500, height = 3200, units = "px", res = 300)
  print(plot11)
  dev.off()
}
