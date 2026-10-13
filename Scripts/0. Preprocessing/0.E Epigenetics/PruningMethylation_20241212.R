rm(list=ls(all=TRUE)) #Delete every previous object in the environment
setwd("/data/natalia/Multi-omics/MOFA")

# I used the dataset without any missings, because it allows for faster computation of correlations.

load("/data/PUBLIC/Methylation/Blood_450karray_freeze2024/Blood.450k.BetavaluesFunNorm.052024.missingvaluesimputed.RData")

dim(Beta.FunNorm) # 411169   3081

Beta.FunNorm <- t(Beta.FunNorm)  # Transpose, because CpGs should be colums

#######################################################################################################################################################################################

sd <- apply(Beta.FunNorm,2,sd)

summary(sd)

#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.

# 0.003375 0.012764 0.018439 0.023611 0.029966 0.341228

length(which(sd > 0.018439)) # median split  205578

#######################################################################################################################################################################################

#  SPECIFY standard deviation cut-off

orderedsd <- sd[order(sd,decreasing=T)] # order standard deviations from high to low. names(orderedsd) contain the CpGIDs

top50 <- names(orderedsd)[1:205578]   # here, you can also change to 200.000 instead if that’s preferred

head(top50) # CpGids of the variable CpGs

#######################################################################################################################################################################################

# Subsetting the dataset to retain only the selected CpGs
beta2_pruned <- Beta.FunNorm[, top50]

dim(beta2_pruned) # 3081 205578 (ppl x CpGs)


save(beta2_pruned,file="MethylationData_pruned_20241212.RData")
