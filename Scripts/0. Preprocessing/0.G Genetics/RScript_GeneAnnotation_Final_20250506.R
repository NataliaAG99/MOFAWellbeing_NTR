########
#
# 22/01/2025
#
# project: MOFA Wellbeing
# data: Gene annotation for mofa
# script: Gene annotation for mofa
#
# R version 4.1.1 (2021-08-10)
#
########

rm(list=ls(all=TRUE)) #Delete every previous object in the environment
setwd("/data/natalia/Multi-omics/MOFA")


library(readr)
library(snpsettest)

library(readr)
library(dplyr)
library(data.table)


SNPsfile <- read_csv("GeneticInput/SNPsFilteredReordered_file.csv", 
                                       col_names = FALSE)
info_snp <- data.frame(id = SNPsfile$X2, chr = SNPsfile$X1, pos = SNPsfile$X4)


Genefile <- read_tsv("GeneticInput/annot.bed", col_names = FALSE)

info_gene <- data.frame(gene.id = Genefile$X4, chr = as.numeric(sub("^chr", "", Genefile$X1)), start = Genefile$X2, end = Genefile$X3)
#Chromosomes X, Y and mitochondrial are introduced as NA, since we don't have that data, I can remove those rows

info_gene <- info_gene[!is.na(info_gene$chr), ]

# > sum(is.na(info_gene$gene.id))
# [1] 359388


###################
#Annotate the SNPs
###################

SNPsAnnotated <- map_snp_to_gene(
  info_snp,
  info_gene,
  extend_start = 1,
  extend_end = 0.5,
  only_sets = FALSE
)



# info_gene$start[which(info_gene$gene.id=="MTND1P23")]
# info_gene$end[which(info_gene$gene.id=="MTND1P23")]
# 
# 
# info_gene$start[which(info_gene$gene.id=="MTND2P28")]
# info_gene$end[which(info_gene$gene.id=="MTND2P28")]
# 
# 
# 
# 
# 
# ####Comparing it with MAGMA output
# 
# annot_data <- read_tsv("/data/natalia/Multi-omics/MOFA/GeneticInput/AnnotationMAGMA_with_genenames.txt", col_names = FALSE)



# PGS_df <- fread("/data/natalia/Multi-omics/MOFA/GeneticInput/final_variable_snps.raw") #A subset sample to not overload the environment

PGS_df <- fread("/data/natalia/Multi-omics/MOFA/GeneticInput/final_variable_snps_complete.raw") #A subset sample to not overload the environment
############
#
# Prepare data
# Genomics
#
############
head(PGS_df$IID, 2)
# fisnumbersPGS_num <- as.numeric(PGS_df$IID)
# Save the rownames
fis_numbers_PGS <- data.frame(fisnumber = as.numeric(PGS_df$IID))
# Select the columns from 7 onward
PGS_df <- as.data.frame(PGS_df[,7:ncol(PGS_df)])
# Restore the rownames
rownames(PGS_df) <- fis_numbers_PGS$fisnumber

#Remove the _ATCG from all column names
colnames(PGS_df) <- sub("_.*", "", colnames(PGS_df))


# Create a new data frame to hold the summed values
new_df <- data.frame(row.names = rownames(PGS_df))

for (gene in names(SNPsAnnotated$sets)) {
  rsids <- SNPsAnnotated$sets[[gene]]
  
  # Sum the values of the columns corresponding to the rsIDs for each gene
  new_df[[gene]] <- rowMeans(PGS_df[, rsids, drop = FALSE], na.rm = TRUE)
}

save(new_df, file = "/data/natalia/Multi-omics/MOFA/GeneticInput/GeneData_filtered_20250506.RData")
