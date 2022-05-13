###############################
## Celiac standard PRS model ##
###############################

library(bigsnpr)


##################################################
## Define important variables and load in data  ##
##################################################


# Relevant data and variables
base.bigSNP <- snp_attach("merged_train.rds")
G <- base.bigSNP$genotypes # 676,783 SNPs across 22 chromosomes
CHR <- base.bigSNP$map$chromosome #Chr of each SNP (vector)
POS <- base.bigSNP$map$physical.pos #Pos of each SNP (vector)
NCORES <- nb_cores() 


## Load in summary data (already processed in "celiac_pre.R")
summary.data <- read.table("summary_stats.txt", header=TRUE)
## Get column index (G matrix) of SNPs which have a summary statistic
summary.ind <- which(base.bigSNP$map$marker.ID %in% summary.data$ID)


## Load in covariate and phenotype
## PCA covariates already performed in celiac_PRS_model.R using plink


## Load eigenvector covariates and separate to training and testing sets
eigen.all <- read.table("all.eigenvec", header=F)
colnames(eigen.all) <- c("FID", "IID", paste0("PC", 1:6))
## Load age and sex covariates
age.sex <- read.table("age_sex.csv", header=T)


# Extract training data subset
train.covar <- read.table("train_id.txt")
colnames(train.covar) <- c("FID", "IID")
train.covar <- merge(train.covar, eigen.all, by=c("FID", "IID"))
train.covar <- merge(train.covar, age.sex, by=c("FID", "IID"))


# Extract test data subset
test.covar <- read.table("test_id.txt")
colnames(test.covar) <- c("FID", "IID")
test.covar <- merge(test.covar, eigen.all, by=c("FID", "IID"))
test.covar <- merge(test.covar, age.sex, by=c("FID", "IID"))


## Merge with phenotype data 
fam.cname <- c("FID", "IID", "FID.F", "FID.M", "SEX", "PHENO")

# Training set
train.fam <- read.table("cd_train_samples.fam", header=F, col.names=fam.cname)
train.data <- merge(train.covar, train.fam[,c(2,6)], by=c("IID"))
train.data$PHENO <- train.data$PHENO - 1 # Convert to 0 and 1


# Testing set
test.fam <- read.table("cd_test_samples.fam", header=F, col.names=fam.cname)
test.data <- merge(test.covar, test.fam[,c(2,6)], by=c("IID"))
test.data$PHENO <- test.data$PHENO - 1 # Convert to 0 and 1




#####################
## Stacking model ##
####################


## Create BETA and log-pvalue vectors
## Vectors need to match column length of genome FBM. Only fill in summary stat
## SNPs
## summary.ind is the index of summary SNPs in genomic FBM
beta <- rep(NA, ncol(G))
beta[summary.ind] <- summary.data$BETA
lpval <- rep(NA, ncol(G))
lpval[summary.ind] <- -log10(summary.data$P)


## Impute missing data since we can't have missing data
## Imputation argument "random" will impute by sampling
## according to allele frequencies
G2 <- snp_fastImputeSimple(G, method="random", ncores=NCORES)


## LONGEST STEP
## Will take about 1.2 hours at 10400 individuals at almost 600K SNPs
## Perform clumping
# 22 chromosomes x 28 clumping parameter combinations
all_keep <- snp_grid_clumping(G2, CHR, POS, lpS = lpval, 
                              exclude = which(is.na(lpval)), ncores = NCORES)
attr(all_keep, "grid")


## Generate PRS using different combinations of hyperparameters
loc <- paste0(getwd(),"/tmp") # Backing file
multi_PRS <- snp_grid_PRS(G2, all_keep, beta, lpval,
                          backingfile=loc, n_thr_lpS = 50, ncores = NCORES)


## Perform linear regression to learn PRS weights
## Include covariates
## covar.train will include the covariates along with the PRS in the model
final.mod.train <- snp_grid_stacking(multi_PRS,
                                     train.data$PHENO,
                                     ncores=NCORES, K=4,
                                     covar.train=covar_from_df(
                                       train.data[,c(3:10)])
                                     )

# Update with new beta scores from derived weights
train.beta <- final.mod.train$beta.G # PRS weights
train.beta.covar <- final.mod.train$beta.covar # Covariate weightings
train.ind <- which(train.beta != 0) #Skip any beta scores with value of zero


# See performance on training data
train.data.sct <- as_FBM(cbind(G2[,train.ind], train.data[,c(3:10)]))
train.pred.sct <- final.mod.train$intercept + big_prodVec(train.data.sct,
                                                          c(train.beta[train.ind],
                                                            train.beta.covar))
AUCBoot(train.pred.sct, train.data$PHENO)


# See performance on test data
# Load in test subset genome
test.bigSNP <- snp_attach("merged_test.rds")
G.test <- test.bigSNP$genotypes
# Impute
G2.test <- snp_fastImputeSimple(G.test, method="random", ncores=NCORES)

test.data.sct <- as_FBM(cbind(G2.test[,train.ind], test.data[,c(3:10)]))
test.pred.sct <- final.mod.train$intercept + big_prodVec(test.data.sct,
                                                          c(train.beta[train.ind],
                                                            train.beta.covar))
AUCBoot(test.pred.sct, test.data$PHENO)


# Training data performance

# AUCBoot(train.pred.sct, train.data$PHENO)
# Mean       2.5%      97.5%         Sd 
# 0.85247328 0.83382450 0.87033202 0.00936552 


# Test data performance

# AUCBoot(test.pred.sct, test.data$PHENO)
# Mean      2.5%     97.5%        Sd 
# 0.8192766 0.7967727 0.8407858 0.0111474 


# Manhattan plot of effect size or coefficients
# Chromosome 6

## Which parameters produce more variance in scores








