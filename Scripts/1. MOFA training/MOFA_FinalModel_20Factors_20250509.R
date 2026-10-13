########
#
# 09/05/2025
#
# project: MOFA Wellbeing
# data: Overlapping sample raw data
# script: MOFA 20 Factors
#
# R version 4.5 (2021-05-09)
#
########

rm(list=ls(all=TRUE)) #Delete every previous object in the environment
setwd("/data/natalia/Multi-omics/MOFA")

#r-packages
# Install MOFA2 from Bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# # BiocManager::install(version = '3.21')
# BiocManager::install("MOFA2")
# BiocManager::install("basilisk")

library(MOFA2)
library(readr)
library(dplyr)
library(data.table)
library(writexl)
library(ggplot2)

############
#
# Read data
#
############

#load data per training data
load("/data/natalia/Multi-omics/MOFA/GeneticInput/GeneData_filtered_20250506.RData") #A subset sample to not overload the environment
PGS_df_1 <- new_df
rm(new_df)
load("/data/natalia/Multi-omics/MOFA/MethylationData_pruned_20241212.RData") #Pruned Cpgs to not overload the env
MS_df <- as.data.frame(beta2_pruned)
rm(beta2_pruned)
load("/data2/PUBLIC/GeneExpression/GODOT_V2_raw_NTR.RData")
trans_df <- genexprdata
rm(genexprdata)
rm(annotation)

############
#
# Prepare data
# Genomics
#
############
head(PGS_df_1[1:5, 1:10])
# fisnumbersPGS_num <- as.numeric(PGS_df$IID)
# Save the rownames
fis_numbers_PGS <- data.frame(fisnumber = rownames(PGS_df_1))
# fis_numbers_PGS <- data.frame(fisnumber = as.numeric(PGS_df$IID))
# Select the columns from 7 onward
# PGS_df <- as.data.frame(PGS_df[,7:ncol(PGS_df)])
# PGS_df <- apply(PGS_df, 2, as.numeric)
# # Restore the rownames
# rownames(PGS_df) <- fis_numbers_PGS$fisnumber


# hist(PGS_df$`TTLL10-AS1`)
# 
# 
# test <- 14 - PGS_df$`TTLL10-AS1` 
# 
# hist(test)
# 
# hist(log(PGS_df$`TTLL10-AS1`))
# hist(sqrt(test))

PGS_df_log <- PGS_df_1
PGS_df_log[sapply(PGS_df_1, is.numeric)] <- lapply(PGS_df_1[sapply(PGS_df_1, is.numeric)], function(x) log(x + 1))

PGS_df <- PGS_df_log

rm(PGS_df_1)
rm(PGS_df_log)

############
#
# Prepare data
# Epigenomics
#
############
load("/data/natalia/Methylation/Blood_450karray_covariates/Covariates_freeze2024.RData")
identical(rownames(MS_df), rownames(covariates)) #TRUE
MS_df <- cbind(MS_df, fisnumber = covariates$fisnumber)
rm(covariates)
MS_df <- MS_df %>% distinct(fisnumber, .keep_all = TRUE) #3049
rownames(MS_df) <- MS_df$fisnumber

#Extract the fisnumbers
fis_numbers_meth <- data.frame(fisnumber = rownames(MS_df))
# # Save the extracted column as a CSV file
#write.table(fis_numbers_meth, "fisnumbersMeth_20240916.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)


############
#
# Prepare data
# Transcriptomics
#
############

#Extract the fisnumbers
fis_numbers_trans <- data.frame(fisnumber = rownames(trans_df))
# # Save the extracted column as a CSV file
# write.csv(fis_numbers, "fisnumbers_Natalia_20240624.csv", row.names = FALSE)

trans_df <- apply(trans_df, 2, as.numeric)
rownames(trans_df) <- fis_numbers_trans$fisnumber


################
#
# Getting the overlapping sample
#
################
length(intersect(fis_numbers_meth$fisnumber, fis_numbers_trans$fisnumber)) #2353
length(Reduce(intersect, list(
  fis_numbers_PGS$fisnumber, 
  fis_numbers_meth$fisnumber, 
  fis_numbers_trans$fisnumber
))) #2320 in common between all the omic layers

total_overlap <- as.data.frame(Reduce(intersect, list(
  fis_numbers_PGS$fisnumber, 
  fis_numbers_meth$fisnumber, 
  fis_numbers_trans$fisnumber
)))

colnames(total_overlap)[1] <- "fisnumber"

# write.csv(as.data.frame(total_overlap), file = "fisnumbers_overlapping_MOFAraw_20241212.csv", row.names = FALSE)

PGS_df_ov <- PGS_df[rownames(PGS_df) %in% total_overlap$fisnumber, ]
MS_df_ov <- MS_df[rownames(MS_df) %in% total_overlap$fisnumber, ]
trans_df_ov <- trans_df[rownames(trans_df) %in% total_overlap$fisnumber, ]


############
#
# Prepare MOFA data
# Genomics
#
############
#PGS_df_ov <- PGS_df_ov[, !colnames(PGS_df_ov) %in% "IID"] #I already deleted it
colnames(PGS_df_ov) <- paste0("GENE_", colnames(PGS_df_ov))
PGS_df_ov <- scale(PGS_df_ov)
head(PGS_df_ov[1:2,1:2])

#GENE_GSTT4 (column 19717) has all NAs (probably not variance at all), so I remove it
PGS_df_ov <- PGS_df_ov[,-19717]

PGS <- t(PGS_df_ov)  # Transpose the matrix if necessary
#counts_df <- as.data.frame(counts)

PGS <- as.matrix(PGS)
# counts <- apply(counts, 2, as.numeric)

# Check the structure of your counts matrix
str(PGS)


############
#
# Prepare data
# Epigenomics
#
############
MS_df_ov <- MS_df_ov[, !colnames(MS_df_ov) %in% "fisnumber"]
colnames(MS_df_ov) <- paste0("CpG_", colnames(MS_df_ov))
MS_df_ov <- scale(MS_df_ov)
head(MS_df_ov[1:2,1:2])

MS <- t(MS_df_ov)  # Transpose the matrix if necessary
#counts_df <- as.data.frame(counts)

MS <- as.matrix(MS)
# counts <- apply(counts, 2, as.numeric)

# Check the structure of your counts matrix
str(MS)


############
#
# Prepare data
# Transcriptomics
#
############
trans_df_ov <- trans_df_ov[, !colnames(trans_df_ov) %in% "fisnumber"]
colnames(trans_df_ov) <- paste0("trans_", colnames(trans_df_ov))
trans_df_ov <- scale(trans_df_ov)
head(trans_df_ov[1:2,1:2])

trans <- t(trans_df_ov)  # Transpose the matrix if necessary
#counts_df <- as.data.frame(counts)

trans <- as.matrix(trans)
# counts <- apply(counts, 2, as.numeric)

# Check the structure of your counts matrix
str(trans)

############
#
# Removing used data
#
############

rm(fis_numbers_PGS)
rm(fis_numbers_meth)
rm(fis_numbers_trans)

rm(PGS_df)
rm(MS_df)
rm(trans_df)

rm(PGS_df_ov)
rm(MS_df_ov)
rm(trans_df_ov)

rm(total_overlap)

############
#
# Create MOFA object
#
############

# Check if the sample names are the same
all.equal(colnames(PGS), colnames(MS))  # Should return TRUE if they are the same
all.equal(colnames(PGS), colnames(trans))  # Should return TRUE if they are the same

# # Rename the sample names if necessary
# colnames(PGS) <- colnames(MS) <- colnames(metab) <- paste0("sample_", 1:ncol(PGS))

# Reorder columns in df2 and df3 to match df1
MS <- MS[, colnames(PGS)]
trans <- trans[, colnames(PGS)]

# Combine the dataframes into a list of views
data_list <- list(View1 = PGS, View2 = MS, View3 = trans)
# data_list <- list(View0 = PGS, View1 = MS, View2 = trans)

# Create a MOFA object
MOFAobject <- create_mofa(data_list)


plot_data_overview(MOFAobject)


############
#
# Define options
#
############

# Retrieve the default data options
data_opts <- get_default_data_options(MOFAobject)
# Set scale_views to TRUE (different ranges among views (omics))
data_opts$scale_views <- TRUE

#Retrieve the model options
#Only change this if you are familiar with the underlying mathematical model
model_opts <- get_default_model_options(MOFAobject)
# model_opts$num_factors <- 10

#Retrieve the train options
train_opts <- get_default_training_options(MOFAobject)
# # Set convergence_mode to slow for definitive results
train_opts$convergence_mode <- "slow"

############
#
# Build and train MOFA object
#
############

#Prepare the MOFA object
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

#Train the MOFA model
outfile = file.path(getwd(),"FinalModel_20Factors_20250509.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)

