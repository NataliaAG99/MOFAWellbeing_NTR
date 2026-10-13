#############################
#### Downstream analyses ####
#############################

#05/08/2025

library(DT)
library(Hmisc)       # For rcorr() function
library(ggcorrplot)  # For visualization
library(ggplot2)     # For plotting
library(MOFA2)
library(gee)
library(readr)
#library(pcpr2)
library(reshape2)
library(lavaan)
library(corrplot)
library(patchwork)

########
# Comparison models
########

MOFA_10F <- load_model("/data/natalia/Multi-omics/MOFA/FinalModel_10Factors_20250509.hdf5")

MOFA_20F <- load_model("/data/natalia/Multi-omics/MOFA/FinalModel_20Factors_20250509.hdf5")


# Extract the factors from the trained MOFA model
factors_20F <- as.data.frame(get_factors(MOFA_20F, factors = "all"))
colnames(factors_20F) <- c("20F_Factor1", "20F_Factor2", "20F_Factor3", "20F_Factor4", "20F_Factor5", "20F_Factor6", "20F_Factor7", "20F_Factor8", "20F_Factor9", "20F_Factor10", "20F_Factor11", "20F_Factor12", "20F_Factor13", "20F_Factor14", "20F_Factor15", "20F_Factor16", "20F_Factor17", "20F_Factor18", "20F_Factor19", "20F_Factor20")

factors_10F <- as.data.frame(get_factors(MOFA_10F, factors = "all"))
colnames(factors_10F) <- c("10F_Factor1", "10F_Factor2", "10F_Factor3", "10F_Factor4", "10F_Factor5", "10F_Factor6", "10F_Factor7", "10F_Factor8", "10F_Factor9", "10F_Factor10")

# Combine factors and covariates into one dataframe
combined_data <- cbind(factors_20F, factors_10F)
# combined_data <- combined_data[, c(1:20, 22:35)]  # Selecting specific columns if needed

# Compute the Spearman correlation matrix
cor_matrix <- cor(combined_data, use = "complete.obs", method = "pearson")

# Plot the correlation matrix using ggcorrplot
ggcorrplot(cor_matrix, 
           hc.order = FALSE,  # Hierarchical clustering
           type = "lower",   # Lower triangle plot
           lab = TRUE,
           lab_size = 2)       # Display correlation values


# Variance decomposition
####################

# Total variance explained per view (10F)
head(get_variance_explained(MOFA_10F)$r2_total[[1]])
# View1      View2      View3 
# 0.1674891 33.2295179 17.3449397  

# Variance explained for every factor in per view (10F)
get_variance_explained(MOFA_10F)$r2_per_factor[[1]]

plot_variance_explained(MOFA_10F, x="view", y="factor", plot_total = T)
plot_variance_explained(MOFA_10F, x="view", y="factor")

# Total variance explained per view (20F)
head(get_variance_explained(MOFA_20F)$r2_total[[1]])
# View1      View2      View3 
# 0.3502786 36.8225455 25.8414328  

# Variance explained for every factor in per view (20F)
get_variance_explained(MOFA_20F)$r2_per_factor[[1]]

plot_variance_explained(MOFA_20F, x="view", y="factor", plot_total = T)
plot_variance_explained(MOFA_20F, x="view", y="factor")

##################
#Remove data and load only model 20F
###################

rm(list = ls())

MOFAobject.trained <- load_model("/data/natalia/Multi-omics/MOFA/FinalModel_20Factors_20250509.hdf5")

############################
# Other functionalities
############################
#Renaming dimensions
views_names(MOFAobject.trained) <- c("Genome", "Epigenome", "Transcriptome")
factors_names(MOFAobject.trained) <- paste("Factor", 1:get_dimensions(MOFAobject.trained)$K, sep=" ")

#Not run the commented part unless necessary
# #Extracting data for downstream analysis
# # "factors" is a list of matrices, one matrix per group with dimensions (nsamples, nfactors)
# factors <- get_factors(model, factors = "all")
# lapply(factors,dim)
# # "weights" is a list of matrices, one matrix per view with dimensions (nfeatures, nfactors)
# weights <- get_weights(model, views = "all", factors = "all")
# lapply(weights,dim)
# # "data" is a nested list of matrices, one matrix per view and group with dimensions (nfeatures, nsamples)
# data <- get_data(model)
# lapply(data, function(x) lapply(x, dim))[[1]]
# 
# ###data in long data.frame format:
# 
# factors <- get_factors(model, as.data.frame = TRUE)
# head(factors, n=3)
# 
# weights <- get_weights(model, as.data.frame = TRUE)
# head(weights, n=3)
# 
# data <- get_data(model, as.data.frame = TRUE)
# head(data, n=3)


########
# Meta data
########

WB_df <- read_csv("/data/natalia/Multi-omics/MOFA/WBData_20250416.csv")

mofa_samples <- samples_names(MOFAobject.trained)[[1]]
# mofa_samples <- samples_names(model)[[1]]

# sample_metadata <- WB_df[WB_df$FISNumbers %in% mofa_samples, ]
# sample_metadata <- as.data.frame(sample_metadata[!duplicated(sample_metadata$FISNumbers), ])
WB_df$Sex_epi <- ifelse(WB_df$Sex_epi == "Male", 1, 0)
colnames(WB_df)[colnames(WB_df) == "FISNumbers"] <- "sample"
# rownames(sample_metadata) <- sample_metadata$FISNumbers
# sample_metadata <- sample_metadata[,c(2:6)]

WB_df <- WB_df[,c("sample", "WBscore", "Age_epi", "timedif", "FamilyNumber", "BMI_epi", "Neut_Perc_epi", "Lymph_Perc_epi", "Mono_Perc_epi", "Eos_Perc_epi", "Baso_Perc_epi", "D", "wbc_trans", "Array_rownum_epi", "Sample_Plate_epi", "smokingnum_epi", "Sex_epi",
                  "PC1_1KG", "PC2_1KG", "PC3_1KG", "PC4_1KG", "PC5_1KG", "PC6_1KG", "PC7_1KG", "PC8_1KG", "PC9_1KG", "PC10_1KG", "PC11_1KG", "PC12_1KG", "PC13_1KG", "PC14_1KG", "PC15_1KG", "PC16_1KG", "PC17_1KG", "PC18_1KG", "PC19_1KG", "PC20_1KG")]
colnames(WB_df) <- c("sample", "WBscore", "Age", "TimeDif", "FamilyNumber", "BMI", "Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils", "Basophils", "D", "WhiteBloodCount", "ArrayRownum", "SamplePlate", "Smoking", "Sex",
                     "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")

metadata_samples <- as.character(WB_df$sample)
# Check which samples are missing
missing_in_metadata <- setdiff(mofa_samples, metadata_samples) #432
missing_in_mofa <- setdiff(metadata_samples, mofa_samples)

sample_metadata_aligned <- as.data.frame(mofa_samples)
colnames(sample_metadata_aligned)[colnames(sample_metadata_aligned) == "mofa_samples"] <- "sample"
sample_metadata_aligned <- merge(sample_metadata_aligned, WB_df[, c("sample", colnames(WB_df)[2:37])], by = "sample", all.x = TRUE)

####
#Add the meta data to the model
####
samples_metadata(MOFAobject.trained) <- sample_metadata_aligned
# samples_metadata(model) <- sample_metadata_aligned
head(samples_metadata(MOFAobject.trained), n=3)


#Table of missingness

# Create a data frame with the missing counts as a single row
#missing_table <- as.data.frame(t(colSums(is.na(WB_df))))
missing_table <- as.data.frame(colSums(is.na(WB_df)))

# Use DT::datatable to create an interactive table
datatable(missing_table, options = list(pageLength = 38, autoWidth = TRUE), 
          caption = 'Missing Values Count by Variable (N=1888)')


#Histogram of WB
hist(WB_df$WBscore)


#######################
# Correlation matrix
########################

# Extract the factors from the trained MOFA model
factors <- as.data.frame(get_factors(MOFAobject.trained, factors = "all"))
colnames(factors) <- c("Factor1", "Factor2", "Factor3", "Factor4", "Factor5", "Factor6", "Factor7", "Factor8", "Factor9", "Factor10", "Factor11", "Factor12", "Factor13", "Factor14", "Factor15", "Factor16", "Factor17", "Factor18", "Factor19", "Factor20")

# Extract the metadata (covariates) from the trained MOFA model
metadata <- samples_metadata(MOFAobject.trained)

# Combine factors and covariates into one dataframe
combined_data_1 <- cbind(factors, metadata)
combined_data <- combined_data_1[, c(1:20, 22:37)]  # Selecting specific columns if needed

combined_data_withPCs <- combined_data_1[, c(1:20, 22:57)]  # Selecting specific columns if needed
colnames(combined_data_withPCs)[colnames(combined_data_withPCs) == "WBscore"] <- "Wellbeing Score"

#####################
# Correlation lavaan
#####################
# Get all numeric variable names
vars <- names(combined_data)
n_vars <- length(vars)

# Initialize matrices to store correlation estimates and p-values
cor_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars, dimnames = list(vars, vars))
p_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars, dimnames = list(vars, vars))


# Loop through variable pairs
for (i in 1:(n_vars - 1)) {
  for (j in (i + 1):n_vars) {
    var1 <- vars[i]
    var2 <- vars[j]
    
    if (var1 != "Sex" && var1 != "Smoking" && var2 != "Sex" && var2 != "Smoking") {
      # Define SEM model for the pair
      model <- paste0(var1, " ~~ ", var2)
    } else if (var1 != "Sex" && (var2 == "Sex" || var2 == "Smoking")) {
      # Define SEM model for the pair
      model <- paste0(var1, " ~ ", var2)
    } else {
      next
    }

    # Fit the SEM model without removing NAs
    fit <- tryCatch(
      sem(model, data = combined_data, cluster = "FamilyNumber", missing = "fiml"), # Full Information Maximum Likelihood (FIML)
      error = function(e) return(NULL)
    )
    
    if (!is.null(fit)) {
      # Extract summary results
      summary_fit <- summary(fit, standardized = TRUE)
      
      # Extract estimate and p-value
      est <- summary_fit$pe$std.all[summary_fit$pe$lhs == var1 & summary_fit$pe$rhs == var2]
      p_val <- summary_fit$pe$pvalue[summary_fit$pe$lhs == var1 & summary_fit$pe$rhs == var2]
      
      # Store results in matrices
      cor_matrix[var1, var2] <- est
      cor_matrix[var2, var1] <- est # Symmetric matrix
      p_matrix[var1, var2] <- p_val
      p_matrix[var2, var1] <- p_val # Symmetric matrix
    }
  }
}


#Remove Family Number, D, Aray Row number and Sample Plate
cor_matrix <- cor_matrix[-c(24,31,33,34), -c(24,31,33,34)]  # Selecting specific columns if needed
p_matrix <- p_matrix[-c(24,31,33,34), -c(24,31,33,34)]  # Selecting specific columns if needed

vars2 <- c("Factor1", "Factor2", "Factor3", "Factor4", "Factor5", "Factor6", "Factor7", "Factor8", "Factor9", "Factor10", "Factor11", "Factor12", "Factor13", "Factor14", "Factor15", "Factor16", "Factor17", "Factor18", "Factor19", "Factor20", "WBscore", "Age", "TimeDif", "BMI", 
           "Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils", "Basophils", "WhiteBloodCount", "Smoking", "Sex")  
# Apply FDR correction
p_values_adj <- matrix(p.adjust(as.vector(p_matrix), method = "fdr"), 
                       nrow = nrow(p_matrix), 
                       ncol = ncol(p_matrix), 
                       dimnames = list(vars2, vars2))

cor_matrix[is.na(cor_matrix)] <- 1
p_values_adj[is.na(p_values_adj)] <- 99

#Rename the variables with the spaces (Wellbeing score and Time difference)
pretty_names <- c(
  "WBscore" = "Wellbeing Score",
  "TimeDif" = "Time Difference"
)

# Change row/column names of correlation and p-value matrices
rownames(cor_matrix) <- ifelse(
  rownames(cor_matrix) %in% names(pretty_names),
  pretty_names[rownames(cor_matrix)],
  rownames(cor_matrix)
)

colnames(cor_matrix) <- ifelse(
  colnames(cor_matrix) %in% names(pretty_names),
  pretty_names[colnames(cor_matrix)],
  colnames(cor_matrix)
)

rownames(p_values_adj) <- ifelse(
  rownames(p_values_adj) %in% names(pretty_names),
  pretty_names[rownames(p_values_adj)],
  rownames(p_values_adj)
)

colnames(p_values_adj) <- ifelse(
  colnames(p_values_adj) %in% names(pretty_names),
  pretty_names[colnames(p_values_adj)],
  colnames(p_values_adj)
)

# Define significance threshold
sig_threshold <- 0.05

# Define the desired variable order
desired_order <- c("Factor1", "Factor2", "Factor3", "Factor4", "Factor5", "Factor6", "Factor7", "Factor8", "Factor9", "Factor10", 
                   "Factor11", "Factor12", "Factor13", "Factor14", "Factor15", "Factor16", "Factor17", "Factor18", "Factor19", 
                   "Factor20", "Wellbeing Score", "Age", "Sex", "Time Difference", "BMI", "Smoking", "Neutrophils", 
                   "Lymphocytes", "Monocytes", "Eosinophils", "Basophils", "WhiteBloodCount")

# Ensure the desired order exists in the dataset
valid_order <- intersect(desired_order, colnames(cor_matrix))

# Reorder matrices
cor_matrix_ordered <- cor_matrix[valid_order, valid_order]
p_values_adj_ordered <- p_values_adj[valid_order, valid_order]

# Plot significant correlations
ggcorrplot(cor_matrix_ordered, 
           hc.order = FALSE,  
           type = "lower",   
           lab = TRUE,
           lab_size = 3,
           p.mat = p_values_adj_ordered,  # Use FDR-corrected p-values
           sig.level = sig_threshold,  # Set significance threshold
           insig = "blank",
           colors = c("darkgreen", "white", "orange"))  # Hide non-significant correlations

# save the image #To guess the inches, divide the px shown in the url of the external viewer by 96
#width 1142 / 96 ≈ 11.9 inches; height =  975 / 96 ≈ 10.2 inches
ggsave("CorrelationMatrix.tiff",last_plot(), path = "/data/natalia/Multi-omics/MOFA/Plots", width = 12, height = 11, dpi=300)

# Plot significant correlations (pvalues)
ggcorrplot(p_values_adj_ordered, 
           hc.order = FALSE,  
           type = "lower",   
           lab = TRUE,
           lab_size = 3,
           p.mat = p_values_adj_ordered,  # Use FDR-corrected p-values
           sig.level = sig_threshold,  # Set significance threshold
           insig = "blank")  # Hide non-significant correlations


#################################
#################################
# Significant factors
#################################
#################################

# Define the variables you want to subset
# variables_old <- c("WBscore", "Age", "Sex", "BMI", "Smoking", "Basophils")
variables <- c("Factor16", "Wellbeing Score", "Age", "BMI", "Smoking", "Eosinophils", "Basophils")

# Subset the correlation matrix
subset_cor_matrix <- as.data.frame(cor_matrix[variables, variables])
subset_p_matrix <- as.data.frame(p_values_adj[variables, variables])

ggcorrplot(subset_cor_matrix, 
           hc.order = FALSE,  
           type = "lower",   
           lab = TRUE,
           lab_size = 3,
           p.mat = subset_p_matrix, # p-values for significance
           sig.level = sig_threshold, # Set significance threshold
           insig = "blank",
           colors = c("darkgreen", "white", "orange"))  # Leave non-significant correlations blank

# save the image #To guess the inches, divide the px shown in the url of the external viewer by 96
#Width 508 == 5.3, height 398 == 4.15
ggsave("CorrelationMatrixSubset.tiff",last_plot(), path = "/data/natalia/Multi-omics/MOFA/Plots", width = 5.3, height = 4.15, dpi=300)

###############
#### Top 20 features combined views
###############

# Choose the factor of interest
factor_index <- 16  # Choose significant factor

# Extract factor loadings for all views
factor_loadings <- get_weights(MOFAobject.trained, factors = factor_index)

# Convert to a data frame for easier manipulation
loadings_df <- do.call(rbind, lapply(names(factor_loadings), function(view) {
  data.frame(
    Feature = rownames(factor_loadings[[view]]),
    Loading = abs(factor_loadings[[view]][, 1]),  # Take absolute value
    Loading_raw = factor_loadings[[view]][, 1],  # Take raw value
    View = view
  )
}))

# Aggregate loadings across views (sum of absolute values for each feature)
aggregated_loadings <- aggregate(Loading ~ Feature, data = loadings_df, sum)

which(grepl("^GENE|^trans", aggregated_loadings[order(-aggregated_loadings$Loading), "Feature"]))[1]
# Select the top 10 features
top_10 <- aggregated_loadings[order(-aggregated_loadings$Loading), ][1:20, ]

# Ensure proper feature ordering (largest to smallest)
top_10$Feature <- factor(top_10$Feature, levels = top_10$Feature)  # Reverse order here

top_10 <- merge(top_10, loadings_df[, c("Feature", "Loading_raw")], by = "Feature", all.x = TRUE)

# Add a new column with the sign of the feature loading
top_10$Sign <- ifelse(top_10$Loading_raw >= 0, "+", "-")
top_10 <- subset(top_10, select = -Loading_raw)

top_10 <- top_10[order(top_10$Loading, decreasing = TRUE), ]

# Preserve order of Feature
top_10$Feature <- factor(top_10$Feature, levels = rev(top_10$Feature))

# Plot
ggplot(top_10, aes(x = Feature, y = Loading)) +
  geom_bar(stat = "identity", width = 0.1, fill = "darkgreen") +  # Thin bars
  geom_text(aes(label = Sign), hjust = -0.2, size = 5) +        # Sign at bar end
  coord_flip() +
  labs(
    title = paste("Top 20 Features for Factor", factor_index),
    x = "Feature",
    y = "Aggregated Feature Importance"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# save the image #To guess the inches, divide the px shown in the url of the external viewer by 96
#Width 854 == 8.9, height 592 == 6.2
ggsave("Top20Features.tiff",last_plot(), path = "/data/natalia/Multi-omics/MOFA/Plots", width = 8.9, height = 6.2, dpi=300)


##################
## GEE for variance explained by Factor 16 in WB
##################
#Just for the bootstrap
WBscore <- scale(combined_data_withPCs$`Wellbeing Score`)
Age <- scale(combined_data_withPCs$Age)
BMI <- scale(combined_data_withPCs$BMI)
Smoking <- as.factor(combined_data_withPCs$Smoking)
Eosinophils <- scale(combined_data_withPCs$Eosinophils)
Basophils <- scale(combined_data_withPCs$Basophils)
D <- scale(combined_data_withPCs$D)
ArrayRownum <- as.factor(combined_data_withPCs$ArrayRownum)
SamplePlate <- as.factor(combined_data_withPCs$SamplePlate)
PC1 <- scale(combined_data_withPCs$PC1)
PC2 <- scale(combined_data_withPCs$PC2)
PC3 <- scale(combined_data_withPCs$PC3)
PC4 <- scale(combined_data_withPCs$PC4)
PC5 <- scale(combined_data_withPCs$PC5)
PC6 <- scale(combined_data_withPCs$PC6)
PC7 <- scale(combined_data_withPCs$PC7)
PC8 <- scale(combined_data_withPCs$PC8)
PC9 <- scale(combined_data_withPCs$PC9)
PC10 <- scale(combined_data_withPCs$PC10)
PC11 <- scale(combined_data_withPCs$PC11)
PC12 <- scale(combined_data_withPCs$PC12)
PC13 <- scale(combined_data_withPCs$PC13)
PC14 <- scale(combined_data_withPCs$PC14)
PC15 <- scale(combined_data_withPCs$PC15)
PC16 <- scale(combined_data_withPCs$PC16)
PC17 <- scale(combined_data_withPCs$PC17)
PC18 <- scale(combined_data_withPCs$PC18)
PC19 <- scale(combined_data_withPCs$PC19)
PC20 <- scale(combined_data_withPCs$PC20)
Factor1 <- scale(combined_data_withPCs$Factor1)
Factor2 <- scale(combined_data_withPCs$Factor2)
Factor3 <- scale(combined_data_withPCs$Factor3)
Factor4 <- scale(combined_data_withPCs$Factor4)
Factor5 <- scale(combined_data_withPCs$Factor5)
Factor6 <- scale(combined_data_withPCs$Factor6)
Factor7 <- scale(combined_data_withPCs$Factor7)
Factor8 <- scale(combined_data_withPCs$Factor8)
Factor9 <- scale(combined_data_withPCs$Factor9)
Factor10 <- scale(combined_data_withPCs$Factor10)
Factor11 <- scale(combined_data_withPCs$Factor11)
Factor12 <- scale(combined_data_withPCs$Factor12)
Factor13 <- scale(combined_data_withPCs$Factor13)
Factor14 <- scale(combined_data_withPCs$Factor14)
Factor15 <- scale(combined_data_withPCs$Factor15)
Factor16 <- scale(combined_data_withPCs$Factor16)
Factor17 <- scale(combined_data_withPCs$Factor17)
Factor18 <- scale(combined_data_withPCs$Factor18)
Factor19 <- scale(combined_data_withPCs$Factor19)
Factor20 <- scale(combined_data_withPCs$Factor20)

#For the rest
combined_data_withPCs$`Wellbeing Score` <- scale(combined_data_withPCs$`Wellbeing Score`)
combined_data_withPCs$Factor16 <- scale(combined_data_withPCs$Factor16)
combined_data_withPCs$Age <- scale(combined_data_withPCs$Age)
# combined_data_withPCs$Sex <- combined_data_withPCs$Sex -1
# combined_data_withPCs$Sex <- as.factor(combined_data_withPCs$Sex)
combined_data_withPCs$BMI <- scale(combined_data_withPCs$BMI)
combined_data_withPCs$Smoking <- as.factor(combined_data_withPCs$Smoking)
combined_data_withPCs$Eosinophils <- scale(combined_data_withPCs$Eosinophils)
combined_data_withPCs$Basophils <- scale(combined_data_withPCs$Basophils)
combined_data_withPCs$D <- scale(combined_data_withPCs$D)
combined_data_withPCs$ArrayRownum <- as.factor(combined_data_withPCs$ArrayRownum)
combined_data_withPCs$SamplePlate <- as.factor(combined_data_withPCs$SamplePlate)
combined_data_withPCs$PC1 <- scale(combined_data_withPCs$PC1)
combined_data_withPCs$PC2 <- scale(combined_data_withPCs$PC2)
combined_data_withPCs$PC3 <- scale(combined_data_withPCs$PC3)
combined_data_withPCs$PC4 <- scale(combined_data_withPCs$PC4)
combined_data_withPCs$PC5 <- scale(combined_data_withPCs$PC5)
combined_data_withPCs$PC6 <- scale(combined_data_withPCs$PC6)
combined_data_withPCs$PC7 <- scale(combined_data_withPCs$PC7)
combined_data_withPCs$PC8 <- scale(combined_data_withPCs$PC8)
combined_data_withPCs$PC9 <- scale(combined_data_withPCs$PC9)
combined_data_withPCs$PC10 <- scale(combined_data_withPCs$PC10)
combined_data_withPCs$PC11 <- scale(combined_data_withPCs$PC11)
combined_data_withPCs$PC12 <- scale(combined_data_withPCs$PC12)
combined_data_withPCs$PC13 <- scale(combined_data_withPCs$PC13)
combined_data_withPCs$PC14 <- scale(combined_data_withPCs$PC14)
combined_data_withPCs$PC15 <- scale(combined_data_withPCs$PC15)
combined_data_withPCs$PC16 <- scale(combined_data_withPCs$PC16)
combined_data_withPCs$PC17 <- scale(combined_data_withPCs$PC17)
combined_data_withPCs$PC18 <- scale(combined_data_withPCs$PC18)
combined_data_withPCs$PC19 <- scale(combined_data_withPCs$PC19)
combined_data_withPCs$PC20 <- scale(combined_data_withPCs$PC20)

# Fit the GEE model 
Factor16corrected <- gee(`Wellbeing Score` ~ Factor16 + Age + BMI + Smoking + Eosinophils + Basophils + D + ArrayRownum + SamplePlate + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, 
                         id = FamilyNumber, 
                         data = combined_data_withPCs, 
                         family = gaussian())

# Extract and calculate p-values in one step, then combine with coefficient table
coef_table <- as.data.frame(summary(Factor16corrected)$coefficients)
coef_table$p_value <- 2 * pnorm(abs(coef_table[, 1] / coef_table[, 2]), lower.tail = FALSE)

# View the combined table
coef_table

betasF16 <- coef_table[2, 1]
(R2_total <- betasF16^2)


# # Fit the GEE model without F14 ####NEED TO ADD PCS
# Factor14corrected_noF <- gee(WBscore ~ Age + Sex + BMI + Smoking + Basophils + ArrayRownum + SamplePlate, 
#                          id = FamilyNumber, 
#                          data = combined_data, 
#                          family = gaussian())
# 
# # Extract and calculate p-values in one step, then combine with coefficient table
# coef_table_noF <- as.data.frame(summary(Factor14corrected_noF)$coefficients)
# coef_table_noF$p_value <- 2 * pnorm(abs(coef_table_noF[, 1] / coef_table_noF[, 2]), lower.tail = FALSE)
# 
# # View the combined table
# coef_table_noF

Fullmodel_r2 <- summary(lm(`Wellbeing Score` ~ Factor16 + Age + BMI + Smoking + Eosinophils + Basophils + D + ArrayRownum + SamplePlate + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = combined_data_withPCs))$r.squared
NoFactor_r2 <- summary(lm(`Wellbeing Score` ~ Age + BMI + Smoking + Eosinophils + Basophils + D + ArrayRownum + SamplePlate + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = combined_data_withPCs))$r.squared

cat("R2_fullmodel:", Fullmodel_r2, 
    "\nR2_covariates:", NoFactor_r2, 
    "\nR2_difference:", Fullmodel_r2 - NoFactor_r2,
    "\nR2_Factor16 (beta2_gee):", R2_total, "\n")


##################
## GEE for variance explained by covariates in Factor 14
##################

# Fit the GEE model 
Factor16_base <- gee(Factor16 ~ Age + BMI + Smoking, 
                         id = FamilyNumber, 
                         data = combined_data_withPCs, 
                         family = gaussian())

# Extract and calculate p-values in one step, then combine with coefficient table
coef_table_base <- as.data.frame(summary(Factor16_base)$coefficients)
coef_table_base$p_value <- 2 * pnorm(abs(coef_table_base[, 1] / coef_table_base[, 2]), lower.tail = FALSE)

# View the combined table
coef_table_base
R2_base <- summary(lm(Factor16 ~ Age + BMI + Smoking, data = combined_data_withPCs))$r.squared


# Fit the GEE model 
Factor16_WBC <- gee(Factor16 ~ Age + BMI + Smoking + Eosinophils + Basophils, 
                             id = FamilyNumber, 
                             data = combined_data_withPCs, 
                             family = gaussian())

# Extract and calculate p-values in one step, then combine with coefficient table
coef_table_WBC <- as.data.frame(summary(Factor16_WBC)$coefficients)
coef_table_WBC$p_value <- 2 * pnorm(abs(coef_table_WBC[, 1] / coef_table_WBC[, 2]), lower.tail = FALSE)

# View the combined table
coef_table_WBC
R2_WBC <- summary(lm(Factor16 ~ Age + BMI + Smoking + Eosinophils + Basophils, data = combined_data_withPCs))$r.squared
  
# Fit the GEE model 
Factor16_EpiBatch <- gee(Factor16 ~ Age + BMI + Smoking + Eosinophils + Basophils + D + ArrayRownum + SamplePlate, 
                    id = FamilyNumber, 
                    data = combined_data_withPCs, 
                    family = gaussian())

# Extract and calculate p-values in one step, then combine with coefficient table
coef_table_EpiBatch <- as.data.frame(summary(Factor16_EpiBatch)$coefficients)
coef_table_EpiBatch$p_value <- 2 * pnorm(abs(coef_table_EpiBatch[, 1] / coef_table_EpiBatch[, 2]), lower.tail = FALSE)

# View the combined table
coef_table_EpiBatch
R2_EpiBatch <- summary(lm(Factor16 ~ Age + BMI + Smoking + Eosinophils + Basophils + D + ArrayRownum + SamplePlate, data = combined_data_withPCs))$r.squared

# Fit the GEE model 
Factor16_PCs <- gee(Factor16 ~ Age + BMI + Smoking + Eosinophils + Basophils + D + ArrayRownum + SamplePlate + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, 
                         id = FamilyNumber, 
                         data = combined_data_withPCs, 
                         family = gaussian())

# Extract and calculate p-values in one step, then combine with coefficient table
coef_table_PCs <- as.data.frame(summary(Factor16_PCs)$coefficients)
coef_table_PCs$p_value <- 2 * pnorm(abs(coef_table_PCs[, 1] / coef_table_PCs[, 2]), lower.tail = FALSE)

# View the combined table
coef_table_PCs
R2_PCs <- summary(lm(Factor16 ~ Age + BMI + Smoking + Eosinophils + Basophils + D + ArrayRownum + SamplePlate + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = combined_data_withPCs))$r.squared

## Difference of R2 per group

cat("R2_base:", R2_base,
    "\nR2_WBC:", R2_WBC,
    "\nR2_EpiBatch:", R2_EpiBatch,
    "\nR2_PCs:", R2_PCs, "\n")

######################
# Variance of WB explained by all factors together
######################

summary(lm(`Wellbeing Score` ~ Factor1 + Factor2 + Factor3 + Factor4 + Factor5 + Factor6 + Factor7 + Factor8 + Factor9 + Factor10 + Factor11 + Factor12 + Factor13 + Factor14 + Factor15 + Factor16 + Factor17 + Factor18 + Factor19 + Factor20, data = combined_data_withPCs))

#Multiple R-squared:  0.02013


#Bootstrapping R2
B <- 5000  # Number of bootstrap samples
rsqWBbyFactors <- numeric(B) # for model 1
rsqWBbyF16_cov<- numeric(B) # for model 2
rsqWBbycov <- numeric(B) # for model 1
rsq_demcov <- numeric(B) # for model 1
rsq_bloodcov <- numeric(B) # for model 1
rsq_batchcov <- numeric(B) # for model 1
rsq_PCcov <- numeric(B) # for model 1
adj.rsqWBbyFactors <- numeric(B) # for model 1
adj.rsqWBbyF16_cov<- numeric(B) # for model 2
adj.rsqWBbycov <- numeric(B) # for model 1
adj.rsq_demcov <- numeric(B) # for model 1
adj.rsq_bloodcov <- numeric(B) # for model 1
adj.rsq_batchcov <- numeric(B) # for model 1
adj.rsq_PCcov <- numeric(B) # for model 1
rsq_diff_boot_F16 <- numeric(B) # for differences in R2s
adj.rsq_diff_boot_F16 <- numeric(B)
rsq_diff_boot_bloodcov <- numeric(B) # for differences in R2s
adj.rsq_diff_boot_bloodcov <- numeric(B)
rsq_diff_boot_batchcov <- numeric(B) # for differences in R2s
adj.rsq_diff_boot_batchcov <- numeric(B)
rsq_diff_boot_PCcov <- numeric(B) # for differences in R2s
adj.rsq_diff_boot_PCcov <- numeric(B)


# with group
group <- combined_data_withPCs$FamilyNumber
# groups should be a vector of family ids (FamilyNumber)
for (i in 1:B) {
  # Resample indices
  n <- length(group)
  #n.orig <- length(unique(group))
  #cl.samp <- split(1 : n, group)
  train_indices <- sample(n, replace=TRUE) # groups
  #train_indices <- unlist(cl.samp[train_indices]) # people
  # Fit models to bootstrap sample
  boot_model <- lm(WBscore[train_indices] ~ Factor1[train_indices] + Factor2[train_indices] + Factor3[train_indices] + Factor4[train_indices] + Factor5[train_indices] + Factor6[train_indices] + Factor7[train_indices] + Factor8[train_indices] + Factor9[train_indices] + Factor10[train_indices] + Factor11[train_indices] + Factor12[train_indices] + Factor13[train_indices] + Factor14[train_indices] + Factor15[train_indices] + Factor16[train_indices] + Factor17[train_indices] + Factor18[train_indices] + Factor19[train_indices] + Factor20[train_indices])
  
  boot_model1 <- lm(WBscore[train_indices] ~ Age[train_indices] + BMI[train_indices] + Smoking[train_indices] + Eosinophils[train_indices] + Basophils[train_indices] + D[train_indices] + ArrayRownum[train_indices] + SamplePlate[train_indices] + PC1[train_indices] + PC2[train_indices] + PC3[train_indices] + PC4[train_indices] + PC5[train_indices] + PC6[train_indices] + PC7[train_indices] + PC8[train_indices] + PC9[train_indices] + PC10[train_indices] + PC11[train_indices] + PC12[train_indices] + PC13[train_indices] + PC14[train_indices] + PC15[train_indices] + PC16[train_indices] + PC17[train_indices] + PC18[train_indices] + PC19[train_indices] + PC20[train_indices])
  boot_model2 <- lm(WBscore[train_indices] ~ Factor16[train_indices] + Age[train_indices] + BMI[train_indices] + Smoking[train_indices] + Eosinophils[train_indices] + Basophils[train_indices] + D[train_indices] + ArrayRownum[train_indices] + SamplePlate[train_indices] + PC1[train_indices] + PC2[train_indices] + PC3[train_indices] + PC4[train_indices] + PC5[train_indices] + PC6[train_indices] + PC7[train_indices] + PC8[train_indices] + PC9[train_indices] + PC10[train_indices] + PC11[train_indices] + PC12[train_indices] + PC13[train_indices] + PC14[train_indices] + PC15[train_indices] + PC16[train_indices] + PC17[train_indices] + PC18[train_indices] + PC19[train_indices] + PC20[train_indices])
  
  boot_model3 <- lm(Factor16[train_indices] ~ Age[train_indices] + BMI[train_indices] + Smoking[train_indices])
  boot_model4 <- lm(Factor16[train_indices] ~ Age[train_indices] + BMI[train_indices] + Smoking[train_indices] + Eosinophils[train_indices] + Basophils[train_indices])
  boot_model5 <- lm(Factor16[train_indices] ~ Age[train_indices] + BMI[train_indices] + Smoking[train_indices] + Eosinophils[train_indices] + Basophils[train_indices] + D[train_indices] + ArrayRownum[train_indices] + SamplePlate[train_indices])
  boot_model6 <- lm(Factor16[train_indices] ~ Age[train_indices] + BMI[train_indices] + Smoking[train_indices] + Eosinophils[train_indices] + Basophils[train_indices] + D[train_indices] + ArrayRownum[train_indices] + SamplePlate[train_indices] + PC1[train_indices] + PC2[train_indices] + PC3[train_indices] + PC4[train_indices] + PC5[train_indices] + PC6[train_indices] + PC7[train_indices] + PC8[train_indices] + PC9[train_indices] + PC10[train_indices] + PC11[train_indices] + PC12[train_indices] + PC13[train_indices] + PC14[train_indices] + PC15[train_indices] + PC16[train_indices] + PC17[train_indices] + PC18[train_indices] + PC19[train_indices] + PC20[train_indices])
  
  # Get R-squared
  rsqWBbyFactors[i] <- summary(boot_model)$r.squared
  
  rsqWBbycov[i] <- summary(boot_model1)$r.squared
  rsqWBbyF16_cov[i] <- summary(boot_model2)$r.squared
  
  rsq_demcov[i] <- summary(boot_model3)$r.squared
  rsq_bloodcov[i] <- summary(boot_model4)$r.squared
  rsq_batchcov[i] <- summary(boot_model5)$r.squared
  rsq_PCcov[i] <- summary(boot_model6)$r.squared
  
  # Get adjusted R-squared
  adj.rsqWBbyFactors[i] <- summary(boot_model)$adj.r.squared
  
  adj.rsqWBbycov[i] <- summary(boot_model1)$adj.r.squared
  adj.rsqWBbyF16_cov[i] <- summary(boot_model2)$adj.r.squared
  
  adj.rsq_demcov[i] <- summary(boot_model3)$adj.r.squared
  adj.rsq_bloodcov[i] <- summary(boot_model4)$adj.r.squared
  adj.rsq_batchcov[i] <- summary(boot_model5)$adj.r.squared
  adj.rsq_PCcov[i] <- summary(boot_model6)$adj.r.squared
  
  # Compute R-squared difference
  rsq_diff_boot_F16[i] <- summary(boot_model2)$r.squared - summary(boot_model1)$r.squared
  
  rsq_diff_boot_bloodcov[i] <- summary(boot_model4)$r.squared - summary(boot_model3)$r.squared
  rsq_diff_boot_batchcov[i] <- summary(boot_model5)$r.squared - summary(boot_model4)$r.squared
  rsq_diff_boot_PCcov[i] <- summary(boot_model6)$r.squared - summary(boot_model5)$r.squared
  
  # Compute adjusted R-squared difference
  adj.rsq_diff_boot_F16[i] <- summary(boot_model2)$adj.r.squared - summary(boot_model1)$adj.r.squared
  
  adj.rsq_diff_boot_bloodcov[i] <- summary(boot_model4)$adj.r.squared - summary(boot_model3)$adj.r.squared
  adj.rsq_diff_boot_batchcov[i] <- summary(boot_model5)$adj.r.squared - summary(boot_model4)$adj.r.squared
  adj.rsq_diff_boot_PCcov[i] <- summary(boot_model6)$adj.r.squared - summary(boot_model5)$adj.r.squared
  
}
# get CI
quantile(rsqWBbyFactors, c(.025,.975))

quantile(rsqWBbycov, c(.025,.975))
quantile(rsqWBbyF16_cov, c(.025,.975))
quantile(rsq_diff_boot_F16, c(.025,.975))

quantile(rsq_demcov, c(.025,.975))
quantile(rsq_bloodcov, c(.025,.975))
quantile(rsq_batchcov, c(.025,.975))
quantile(rsq_PCcov, c(.025,.975))
quantile(rsq_diff_boot_bloodcov, c(.025,.975))
quantile(rsq_diff_boot_batchcov, c(.025,.975))
quantile(rsq_diff_boot_PCcov, c(.025,.975))

#CI for adjusted
quantile(adj.rsqWBbyFactors, c(.025,.975))

quantile(adj.rsqWBbycov, c(.025,.975))
quantile(adj.rsqWBbyF16_cov, c(.025,.975))
quantile(adj.rsq_diff_boot_F16, c(.025,.975))

quantile(adj.rsq_demcov, c(.025,.975))
quantile(adj.rsq_bloodcov, c(.025,.975))
quantile(adj.rsq_batchcov, c(.025,.975))
quantile(adj.rsq_PCcov, c(.025,.975))
quantile(adj.rsq_diff_boot_bloodcov, c(.025,.975))
quantile(adj.rsq_diff_boot_batchcov, c(.025,.975))
quantile(adj.rsq_diff_boot_PCcov, c(.025,.975))

######################
# Visualization of factors
######################
# Load sample metadata
metadata <- samples_metadata(MOFAobject.trained)

# Rename "WBscore" column to "Wellbeing Score"
colnames(metadata)[colnames(metadata) == "WBscore"] <- "Wellbeing Score"

# Set the updated metadata back into the MOFA object
samples_metadata(MOFAobject.trained) <- metadata
#Scatterplot (visualise if the relationship between (top) features and factor values is linear)
#Color them by covariate (sig16: Age, BMI, Smoking, Eosinophils, Basophils)


p1 <- plot_data_scatter(
  MOFAobject.trained,
  view = "Epigenome",
  factor = 16,
  features = 3,
  add_lm = TRUE,
  color_by = "Wellbeing Score"
)

p1 <- p1 + annotate("text", x=9, y=9, label=" ")

p2 <- plot_data_scatter(
  MOFAobject.trained,
  view = "Epigenome",
  factor = 16,
  features = 3,
  add_lm = TRUE,
  color_by = "Age"
)

p2 <- p2 + annotate("text", x=9, y=9, label=" ")

p3 <- plot_data_scatter(
  MOFAobject.trained,
  view = "Epigenome",
  factor = 16,
  features = 3,
  add_lm = TRUE,
  color_by = "BMI"
)

p3 <- p3 + annotate("text", x=9, y=9, label=" ")

combined_plot <- p1 / p2 / p3
combined_plot

# save the image #To guess the inches, divide the px shown in the url of the external viewer by 96
#Width 813 == 8.47, height 965 == 10.05
ggsave("Wellbeing_Age_BMI.tiff",last_plot(), path = "/data/natalia/Multi-omics/MOFA/Plots", width = 8.5, height = 11, dpi=300)


p4 <- plot_data_scatter(
  MOFAobject.trained,
  view = "Epigenome",
  factor = 16,
  features = 3,
  add_lm = TRUE,
  color_by = "Smoking"
)

p4 <- p4 + annotate("text", x=9, y=9, label=" ")

p5 <- plot_data_scatter(
  MOFAobject.trained,
  view = "Epigenome",
  factor = 16,
  features = 3,
  add_lm = TRUE,
  color_by = "Eosinophils"
)

p5 <- p5 + annotate("text", x=9, y=9, label=" ")

p6 <- plot_data_scatter(
  MOFAobject.trained,
  view = "Epigenome",
  factor = 16,
  features = 3,
  add_lm = TRUE,
  color_by = "Basophils"
)

p6 <- p6 + annotate("text", x=9, y=9, label=" ")

combined_plot2 <- p4 / p5 / p6
combined_plot2

# save the image #To guess the inches, divide the px shown in the url of the external viewer by 96
#Width 813 == 8.47, height 965 == 10.05
ggsave("Smoking_Eosinophils_Basophils.tiff",last_plot(), path = "/data/natalia/Multi-omics/MOFA/Plots", width = 8.5, height = 11, dpi=300)

