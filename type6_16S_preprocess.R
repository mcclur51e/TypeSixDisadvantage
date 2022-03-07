#!/usr/bin/env Rscript
#####################################################################
########## To run this code in scrum ################################
#####################################################################
# module load R
# Rscript code.R ~/INPUTDIR/

##################################################################################
##################################################################################
### This code takes as input:
### 1) forward and reverse raw fastq files grouped in a directory labeled " raw/ "
### 2) a mapping file named " table_map.csv " containing 1 row for each sample. 
###    a) Where the first column lists the same sample IDs used to label the raw 
###       .fastq files. 
###    b) Must include a column labeled " Sample_or_Control " that labels negative 
###       controls as " Control " and labels samples as " Sample ". 
###    c) Must include a column labeled " InitialConc " that includes a numerical 
###       value for the concentration of DNA measured after initial extraction 
### 3) a file " ~/Masters/silva_nr_v132_train_set_RossMod.fa " as a reference file
###    containing a curated list of 16S sequences for comparison
### 4) a file " ~/Masters/taxa_summary.R " containing code for the fast_melt function
###
### This code will produce many files including:
### 1) A directory labeled " filtered/ " that will include a subset of .fastq files 
###    that contain at least 1 read
### 2) " output/table_noReads.csv " contains a list of samples with no reads after 
###    contaminants and taxa occurring in only one sample have been removed. This 
###    table includes the original number of reads present in the raw data for each 
###    sample listed.
### 3) " output/physeq_start.RData " contains the final phyloseq object that 
###    should be loaded for use in future data analysis
### 4) " output/table_contaminants.csv " contains a list of taxa identified as 
###    contaminants and the statistical results that determines them as such
##################################################################################
##################################################################################

#####################################################################
########## Processing Set −up #######################################
#####################################################################
#Rpath <- "~/R/x86_64-redhat-linux-gnu-library/3.6/"
########## Call libraries for use ##########
# Packages from CRAN
library("ggplot2") # version 3.2.1
library("reshape2") # version 1.4.3
library("data.table") # version 1.12.6
library("RColorBrewer") # design new color palette for data # version 1.1-2
library("scales") # for scientific notation in plots # version 1.0.0
library("tidyr") # version 1.0.0
library("dplyr") # version 0.8.3
library("vegan") # version 2.5-6

# Packages from Bioconductor
library("phyloseq") # version 1.28.0
library("decontam") # identify contaminant ASVs # version 1.4.0
library("DESeq2") # The DESeq2 model internally corrects for library size, so transformed or normalized values such as counts scaled by library size should not be used as input. # version 1.24.0
library("dada2") # version 1.12.1
library("Biostrings") # 2.52.0
library("microbiome") # version 1.6.0

### ggpubr, pheatmap, and pairwiseAdonis currently not installed on discovery

input <- commandArgs(trailingOnly = TRUE) # Read arguments when beginning script in scrum
path <- as.character(input[1])
setwd(path) # set working directory
dir.create("output") # create directory for output files to go

#####################################################################
##### ASV assignment with DADA2 #####################################
#####################################################################
fnFs <- sort(list.files(paste0(path,"raw/"), pattern="R1_001", full.names = TRUE))
fnRs <- sort(list.files(paste0(path,"raw/"), pattern="R2_001", full.names = TRUE))
if(length(fnFs) != length(fnRs)) stop("At least one sample is unpaired. Please check forward and reverse reads are present for all samples")
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:20]) # plot example quality profile(s) of forward sequences
plotQualityProfile(fnRs[1:20]) # plot example quality profile(s) of reverse sequences
save(fnFs,file=("output/output_fnFs.RData")) # Save .RData file 
save(fnRs,file=("output/output_fnRs.RData")) # Save .RData file 

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
save(filtFs,file=("output/output_filtFs.RData")) # Save .RData file
save(filtRs,file=("output/output_filtRs.RData")) # Save .RData file

#####################################################################
##### Starting midway, removing Bimeras #############################
#####################################################################
seqtab <- readRDS("output/output_seqtab.rds")
load("output/output_mergers.RData")
seqtab.all <- seqtab
#seqtab.all <- mergeSequenceTables(seqtab)
seqtab.noBim <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE)
saveRDS(seqtab.noBim, "output/output_seqtabNoBim.rds")

tax <- assignTaxonomy(seqtab.noBim, "~/Masters/silva_nr_v132_train_set_RossMod.fa", multithread=TRUE)
saveRDS(tax, "output/output_taxFinal.rds")

#####################################################################
########## Transfer data into Phyloseq ##############################
#####################################################################
OTU = otu_table(seqtab.noBim, taxa_are_rows=FALSE) # assigns ASV table from DADA output
TAX = tax_table(tax) # reads data.frame into matrix without converting to ASVs # reads csv file into data.frame with row names in column 1 ### .csv file prepped from "otu_table.txt" with taxa data separated by columns and headers added for Kingdom, Phylum, Class, etc. First column = ASV IDs
MAP = sample_data(data.frame(read.csv("table_map.csv", header = TRUE, row.names = 1, check.names=FALSE))) # reads csv file into data.frame with row names in column 1
physeq = phyloseq(OTU, TAX, MAP)
### create and assign ASV numbers to consensus sequences for easier reference ###
ASV <- paste0("ASV", seq(ntaxa(physeq))) 
Sequence <- row.names(tax_table(physeq)) 
bind.asv <- cbind(tax_table(physeq),ASV)
bind.seq <- cbind(bind.asv,Sequence) 
TAX2 = tax_table(as.matrix(bind.seq)) # define new taxonomy table
physeqR = phyloseq(OTU, TAX2, MAP)
save(physeqR,file=("output/physeq_initial.RData")) # Save the phyloseq data object in a .RData file 
write.csv(tax_table(physeqR),"output/table_tax.csv") # Save taxonomy table as .csv
write.csv(otu_table(physeqR),"output/table_otu.csv") # Save ASV table as .csv

############################################
### Identify contaminants using Decontam ###
############################################
physeqSR <- subset_samples(physeqR, Project%in%c("Controls","TypeVImouse") & !Control%in%c("extract"))
sample_data(physeqSR)$InitialConc <- as.numeric(as.character(sample_data(physeqSR)$DNAconc)) + 0.001 # math to make decontam work
sample_data(physeqSR)$is.neg <- sample_data(physeqSR)$Control %in%c("negExtract","negPCR")
df.map <- as.data.frame(sample_data(physeqSR))
df.map$LibrarySize <- sample_sums(physeqSR)
df.map <- df.map[order(df.map$LibrarySize),]
df.map$Index <- seq(nrow(df.map))
ggplot(data=df.map, aes(x=Index, y=LibrarySize, color=Control)) + 
  geom_point() +
  scale_y_log10()

### ID contaminants by either prevalence or frequency
df.contam <- isContaminant(physeqSR, method="either", neg="is.neg", conc="InitialConc",threshold=0.2)
phyN <- prune_taxa(!df.contam$contaminant, physeqSR)
write.csv(subset(df.contam,contaminant=="TRUE"),"output/table_contaminants.csv")
save(phyN,file=("output/physeq_decontam.RData")) # Save the phyloseq data object in a .RData file 

#####################################################################
########## Remove ASVs present in only 1 sample #####################
#####################################################################
source("~/Masters/taxa_summary.R",local=TRUE) # load fast_melt function

dt.phyN = fast_melt(phyN) # make data table from phyloseq object (data table. physeq Rarefied)
prev.phyN = dt.phyN[, list(Prevalence = sum(count > 1),TotalPer = sum(count),
                           MinCount = min(count), MaxCount = max(count)),
                    by = ASV] # make simple table listing 'ASV, Prevalence, TotalPer, MinCount, and MaxCount' (prevalence . physeq Rarefied)
ls.Pres = prev.phyN[(Prevalence > 1), ASV] # Make list of ASVs presMdent in dataset at least once (list . PresMdent)
phyN.pruT <- subset_taxa(phyN,ASV%in%c(ls.Pres)) # remove taxa not in ls.PresMd from phyloseq object (physeq Raw . pruned 2)
phyN.pruS <- subset_samples(phyN.pruT, sample_sums(phyN.pruT)>0)

physeq0 <- subset_samples(physeqSR,!sample_names(physeqSR)%in%sample_names(phyN.pruS))
write.csv(sample_sums(physeq0),"output/table_noReads.csv")
save(phyN.pruS,file=("output/physeq_start.RData")) # Save the phyloseq data object in a .RData file 

### Table outlining how many sequences remain after each pre-processing step
getN <- function(x) sum(getUniques(x))
#filtered <- load("output/output_filtered.RData")
track <- cbind(sapply(mergers, getN), 
               rowSums(seqtab.noBim), 
               rowSums(otu_table(prune_taxa(!df.contam$contaminant, physeqR))), 
               rowSums(otu_table(subset_taxa(physeqR,ASV%in%c(ls.Pres)))))
colnames(track) <- c("merged", "nonchim", "Decontam", "noSingletons")

rownames(track) <- sample.names
#head(track)
write.csv(track,"output/table_trackReads.csv")

phyG.pruS <- tax_glom(phyN.pruS, "Species")
df.otu <- t(as.data.frame(otu_table(phyN.pruS)))
df.otuXtend <- cbind(tax_table(phyN.pruS),df.otu)
rownames(df.otuXtend) <- df.otuXtend[,"ASV"]
write.csv(df.otuXtend,"output/table_otuXtend.csv")

#####################################################################
##################### END OF 16S PRE-PROCESSING #####################
#####################################################################
