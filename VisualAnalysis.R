###############################################################################
# Script to plot the diversity of microbiome data extracted from macaque
# faeces as part of the Sulawesi project
# Author: Dami Oresegun
#
###############################################################################
# ensure clear workspace
rm(list = ls())
# install and set libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("phyloseq", "ggplot2", "cowplot", "dplyr", "ggplotify"))
library("phyloseq")
library("ggplot2")
library("cowplot")
library("dplyr")
library("ggplotify")
# get arguments
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0) {
  stop("At least one argument is needed to begin")
}
args = c(getwd(), "Combined_BrackenReports.biom", "Species", "knowlesi")
workD = args[1] # the working directory
inpBiom = args[2] # the input biom file
desiredTaxa = args[3] # the taxonomic level of the desired organism
desiredOrg = args[4] # the desired family, organism, family or species to search for
# if species, only the species name is needed. NOT the genus. E.g if you want Plasmodium knowlesi, simply enter knowlesi
# set working directory
setwd(workD)
cat("The working directory is ", getwd())
# read in the data
microData = import_biom(inpBiom)
# change the columns of the taxonomy table
colnames(microData@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# remove the useless characters in the taxa
microData@tax_table@.Data <- substring(microData@tax_table@.Data, 4)
# remove the empty or N.A phyla that does not help us
microData <- subset_taxa(microData, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
microData <- subset_taxa(microData, !is.na(Class) & !Class %in% c("", "uncharacterized"))
microData <- subset_taxa(microData, !is.na(Family) & !Family %in% c("", "uncharacterized"))
microData <- subset_taxa(microData, !is.na(Order) & !Order %in% c("", "uncharacterized"))
microData <- subset_taxa(microData, !is.na(Species) & !Species %in% c("", "uncharacterized"))
summary(microData@tax_table@.Data=="")
###########  Create table, number of features for each phyla ########### 
noFeatures <- as.data.frame(table(tax_table(microData)[, "Phylum"], exclude = NULL))
colnames(noFeatures) <- c("Phylum", "Number of taxa")
write.csv(noFeatures, "PhylumTaxa.csv",row.names = FALSE)
########### Compute prevalence of each feature, store as data.frame ########### 
prevDF = apply(X = otu_table(microData),
               MARGIN = ifelse(taxa_are_rows(microData), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevDF = data.frame(Prevalence = prevDF,
                    TotalAbundance = taxa_sums(microData),
                    tax_table(microData))
prevDF <- cbind(Taxa = rownames(prevDF), prevDF)
rownames(prevDF) <- 1:nrow(prevDF)
write.csv(prevDF, "Prevalence_of_all_Taxa.csv", row.names = FALSE)
###########  Check for phyla of low prevalence features by computing the average % prevalence and total prevalence
          # of the phylum across all samples ########### 
abunDF <- plyr::ddply(prevDF, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
colnames(abunDF) <- c("Phylum", "Avg.Prevalence", "Total.Prevalence")
write.csv(abunDF, "Abundance_Prevalence.csv", row.names = FALSE)
###########  Prevalence filtering ########### 
# plot the abundance
#tiff("Prevalence_Abundance_plot.tif",width = 4500, height = 3200, res=300, units="px",)
png("Prevalence_Abundance_plot.png", height = 800, width=1600)
ggplot(prevDF, aes(TotalAbundance, Prevalence / nsamples(microData),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.3) +
  scale_x_log10() +  xlab("Total Abundance/Reads") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
dev.off()
# Define prevalence threshold as 5% of total samples
prevThres = 0.05 * nsamples(microData)
# use the prevalence threshold to filter taxa
keepTaxa = rownames(prevDF)[(prevDF$Prevalence >= prevThres)] # note here we can also only keep the low abundance/rare sequences
filteredData = prune_taxa(keepTaxa,microData) # this needs more work
########### Compute prevalence of filtered taxonomies ########### 
NewDF = apply(X = otu_table(filteredData),
               MARGIN = ifelse(taxa_are_rows(filteredData), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
NewDF = data.frame(Prevalence = NewDF,
                    TotalAbundance = taxa_sums(filteredData),
                    tax_table(filteredData))
NewDF <- cbind(Taxa = rownames(NewDF), NewDF)
rownames(NewDF) <- 1:nrow(NewDF)
write.csv(NewDF, "Filtered_Prevalence_of_Taxa.csv", row.names = FALSE)
# plot the abundance of filtered
#tiff("Filtered_Prevalence_Abundance_plot.tif",width = 4500, height = 3200, res=300, units="px",)
png("Filtered_Prevalence_Abundance_plot.png", height = 800, width=1600)
ggplot(NewDF, aes(TotalAbundance, Prevalence / nsamples(filteredData),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.3) +
  scale_x_log10() +  xlab("Total Abundance/Reads") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
dev.off()
########### Agglomerate taxa ########### 
# groups all OTUs of the same taxonomy at the species rank
length(get_taxa_unique(filteredData, taxonomic.rank = "Phylum"))
length(get_taxa_unique(microData, taxonomic.rank = "Species"))
aglomeraData = tax_glom(filteredData, "Phylum", NArm = TRUE)
# agglomeration of raw data is needed later
rawAglomeraData = tax_glom(microData,"Phylum", NArm = TRUE)
######################################################################################
# At this point the data has been significantly trimmed and is easier to analyse
# The trimmed and untrimmed data will be constantly compared going forward
######################################################################################
###### Calculate some descriptive statistics ###### 
# number of reads
rawReads <- sample_sums(microData)
filterReads <- sample_sums(aglomeraData)
readList <- list(rawReads, filterReads)
readList = as.data.frame(do.call(cbind,readList))
colnames(readList) <- c("Raw Read Count", "Filtered Read Count")
# ranges and mean
rawSum <- summary(microData@otu_table@.Data)
filtSum <- summary(aglomeraData@otu_table@.Data)
write.table(rawSum, "Raw_DescriptiveStats.tsv", row.names = FALSE, sep = "\t")
write.table(filtSum, "Filtered_DescriptiveStats.tsv", row.names = FALSE, sep = "\t")

###### Calculate alpha diversity for each kingdom ###### 
# Bacteria
wholeBac <- subset_taxa(microData, Kingdom == "Bacteria") # all bacteria for unfiltered data
filtBac <- subset_taxa(aglomeraData, Kingdom == "Bacteria") # all bacteria for filtered data
b1 <- estimate_richness(wholeBac, split = TRUE, measures = c("Shannon"))# alpha diversity for unfiltered data
b2 <- estimate_richness(filtBac, split = TRUE, measures = c("Shannon")) # alpha diversity for filtered dat
blist <- list(b1, b2) # combine into one list
bComb = as.data.frame(do.call(cbind, blist)) # make list into dataframe
colnames(bComb) <- c("Whole Bacteria Dataset", "Filtered Bacteria Dataset") # change column names
bComb <- cbind(Isolate = rownames(bComb), bComb) # remove index column
rownames(bComb) <- 1:nrow(bComb)
#Archaea
wholeArch <- subset_taxa(microData, Kingdom == "Archaea")
filtArch <- subset_taxa(aglomeraData, Kingdom == "Archaea")
a1 <- estimate_richness(wholeArch, split = TRUE, measures = c("Shannon"))
a2 <- estimate_richness(filtArch, split = TRUE, measures = c("Shannon"))
alist <- list(a1, a2)
aComb = as.data.frame(do.call(cbind, alist))
colnames(aComb) <- c("Whole Archaea Dataset", "Filtered Archaea Dataset")
aComb <- cbind(Isolate = rownames(aComb), aComb)
rownames(aComb) <- 1:nrow(aComb)
# Eukaryotes
wholeEuk <- subset_taxa(microData, Kingdom == "Eukaryota")
filtEuk <- subset_taxa(aglomeraData, Kingdom == "Eukaryota")
e1 <- estimate_richness(wholeEuk, split = TRUE, measures = c("Shannon"))
e2 <- estimate_richness(filtEuk, split = TRUE, measures = c("Shannon"))
elist <- list(e1, e2)
eComb = as.data.frame(do.call(cbind, elist))
colnames(eComb) <- c("Whole Eukaryote Dataset", "Filtered Eukaryote Dataset")
eComb <- cbind(Isolate = rownames(eComb), eComb)
rownames(eComb) <- 1:nrow(eComb)
# combine dataframe tables into one and output to csv
AllComb <- list(bComb, aComb, eComb)
AllComb <- Reduce(function(x, y) merge(x, y, all=TRUE), AllComb)
write.csv(AllComb, "Alpha_Diversity.csv", row.names = FALSE)
# plot reads
RawBacReads <- data.frame(Samples = sample_names(wholeBac), Reads = sample_sums(wholeBac))
FiltBacReads <- data.frame(Samples = sample_names(filtBac), Reads = sample_sums(filtBac))
RawArchReads <- data.frame(Samples = sample_names(wholeArch), Reads = sample_sums(wholeArch))
FiltArchReads <- data.frame(Samples = sample_names(filtArch), Reads = sample_sums(filtArch))
RawEukReads <- data.frame(Samples = sample_names(wholeEuk), Reads = sample_sums(wholeEuk))
FiltEukReads <- data.frame(Samples = sample_names(filtEuk), Reads = sample_sums(filtEuk))

rawBacCount <- ggplot(data = RawBacReads, mapping = aes(x = Samples, y= Reads)) + 
  geom_col(aes(color=Samples, fill=Samples)) + theme(legend.position = "none")
filtBacCount <- ggplot(data = FiltBacReads, mapping = aes(x = Samples, y= Reads)) + 
  geom_col(aes(color=Samples, fill=Samples)) + theme(legend.position = "none")
rawArchCount <- ggplot(data = RawArchReads, mapping = aes(x = Samples, y= Reads)) + 
  geom_col(aes(color=Samples, fill=Samples)) + theme(legend.position = "none")
filtArchCount <- ggplot(data = FiltArchReads, mapping = aes(x = Samples, y= Reads)) + 
  geom_col(aes(color=Samples, fill=Samples)) + theme(legend.position = "none")
rawEukCount <- ggplot(data = RawEukReads, mapping = aes(x = Samples, y= Reads)) + 
  geom_col(aes(color=Samples, fill=Samples)) + theme(legend.position = "none")
filtEukCount <- ggplot(data = FiltEukReads, mapping = aes(x = Samples, y= Reads)) + 
  geom_col(aes(color=Samples, fill=Samples)) + theme(legend.position = "none")
legend <- get_legend(filtEukCount + theme(legend.position = "top"))
plot_legend <- plot_grid(legend)
plotdata <- plot_grid(rawBacCount, filtBacCount, rawArchCount, filtArchCount, 
          rawEukCount, filtEukCount,ncol = 2, 
          labels = c("Raw Bacteria data", "Filtered Bacteria data", 
                     "Raw Archaea data", "Filtered Archaea data", 
                     "Raw Eukaryote data", "Filtered Eukaryote data"))
#tiff("ReadCount.tif",width = 4800, height = 3700, res=300, units="px",)
png("ReadCount.png", height = 1200, width=1600)
plot_grid(plotdata, plot_legend, ncol = 2, rel_widths = c(3,1))
dev.off()
###### Plot alpha diversity for each kingdom ######
plotB1 <- plot_richness(physeq = wholeBac, x="samples", title = "Raw Bacteria data",
                        measures = "Shannon", sortby = "Shannon")
plotB1 <- plotB1 + geom_point(size=5, alpha=0.5)
plotB2 <- plot_richness(physeq = filtBac, x="samples", title = "Filtered Bacteria data",
                        measures = "Shannon", sortby = "Shannon")
plotB2 <- plotB2 + geom_point(size=5, alpha=0.5)
plotA1 <- plot_richness(physeq = wholeArch, x="samples", title = "Raw Archaea data",
                        measures = "Shannon", sortby = "Shannon")
plotA1 <- plotA1 + geom_point(size=5, alpha=0.5)
plotA2 <- plot_richness(physeq = filtArch, x="samples", title = "Filtered Archaea data",
                        measures = "Shannon", sortby = "Shannon")
plotA2 <- plotA2 + geom_point(size=5, alpha=0.5)
plotE1 <- plot_richness(physeq = wholeEuk, x="samples", title = "Raw Eukaryote data",
                        measures = "Shannon", sortby = "Shannon")
plotE1 <- plotE1 + geom_point(size=5, alpha=0.5)
plotE2 <- plot_richness(physeq = filtEuk, x="samples", title = "Filtered Eukaryote data",
                        measures = "Shannon", sortby = "Shannon")
plotE2 <- plotE2 + geom_point(size=5, alpha=0.5)
#tiff("Alpha_Diversity.tif",width = 4500, height = 3200, res=300, units="px",)
png("Alpha_Diversity.png", height = 1200, width=1600)
plot_grid(plotB1,plotB2,plotA1,plotA2,plotE1,plotE2, ncol = 2, 
                      labels = c("A", "B", "C", "D", "E", "F"))
dev.off()
#### make plots
# transform the data into a dataframe for manipulation
filteredAgloData <- psmelt(aglomeraData)
RawAgloData <- psmelt(rawAglomeraData)
#tiff("Phyla.tif",width = 5000, height = 3800, res=300, units="px",)
png("Phyla.png", height = 1200, width=1600)
rawPlot <- ggplot(data=RawAgloData, aes(x=Sample, y=Abundance, fill=Phylum)) +
  geom_bar(aes(), stat="identity", position="stack") +
  labs(x = "", y = "Relative Abundance\n") +
  theme(panel.background = element_blank(),legend.position = "bottom")
filtPlot <- ggplot(data=filteredAgloData, aes(x=Sample, y=Abundance, fill=Phylum)) +
  geom_bar(aes(), stat="identity", position="stack") +
  labs(x = "", y = "Relative Abundance\n") +
  theme(panel.background = element_blank(),legend.position = "bottom")
plot_grid(rawPlot, filtPlot, labels = c("Raw Agglomerated Phyla", "Filtered Agglomerated Phyla"))
dev.off()

################
# to come back to
Fyes <- transform_sample_counts(microData, function(x) x/sum(x))
head(otu_table(Fyes))
ps_rare <- rarefy_even_depth(microData, sample.size = 120000, rngseed = 123, replace = FALSE)
sample_sums(ps_rare)
estimate_richness(ps_rare, split = TRUE, measures = "Shannon")
################




# further filtering to change ID of taxa less than 0.5% abundance
#filteredAgloData$Phylum <- as.character(filteredAgloData$Phylum)
#filteredAgloData$Phylum[filteredAgloData$Abundance < 0.5] <- "Phyla < 0.5% abundance"
#unique(filteredAgloData$Phylum)



############################## Comment yard ############################## 
########### Plot alpha diversity with colour ########### 
#tiff("Alpha_Diversity.tif",width = 4500, height = 3200, res=300, units="px",)
#png("Alpha_Diversity.png", height = 1200, width=1600)
# plotB1 <- plot_richness(physeq = wholeBac, x="samples", title = "Raw Bacteria data",
#                         measures = "Shannon", sortby = "Shannon")
# plotB1 <- plotB1 + geom_point(size=5, alpha=0.5,aes(color=samples)) + theme(legend.position = "none")
# plotB2 <- plot_richness(physeq = filtBac, x="samples", title = "Filtered Bacteria data",
#                         measures = "Shannon", sortby = "Shannon")
# plotB2 <- plotB2 + geom_point(size=5, alpha=0.5,aes(color=samples)) + theme(legend.position = "none")
# plotA1 <- plot_richness(physeq = wholeArch, x="samples", title = "Raw Archaea data",
#                         measures = "Shannon", sortby = "Shannon")
# plotA1 <- plotA1 + geom_point(size=5, alpha=0.5,aes(color=samples)) + theme(legend.position = "none")
# plotA2 <- plot_richness(physeq = filtArch, x="samples", title = "Filtered Archaea data",
#                         measures = "Shannon", sortby = "Shannon")
# plotA2 <- plotA2 + geom_point(size=5, alpha=0.5,aes(color=samples)) + theme(legend.position = "none")
# plotE1 <- plot_richness(physeq = wholeEuk, x="samples", title = "Raw Eukaryote data",
#                         measures = "Shannon", sortby = "Shannon")
# plotE1 <- plotE1 + geom_point(size=5, alpha=0.5,aes(color=samples))+ theme(legend.position = "none")
# plotE2 <- plot_richness(physeq = filtEuk, x="samples", title = "Filtered Eukaryote data",
#                         measures = "Shannon", sortby = "Shannon")
# plotE2 <- plotE2 + geom_point(size=5, alpha=0.5,aes(color=samples)) + theme(legend.position = "none")
# legend <- get_legend(plotA2 + theme(legend.position = "bottom"))
# plot_legend <- plot_grid(legend)
# plotdata <- plot_grid(plotB1,plotB2,plotA1,plotA2,plotE1,plotE2, ncol = 2, 
#           labels = c("A", "B", "C", "D", "E", "F"))
# png("Alpha_Diversity.png", height = 1200, width=1600)
# plot_grid(plotdata, plot_legend, ncol = 2, rel_widths = c(3,1))
# dev.off()
############################################################ 