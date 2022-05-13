###############################
## Celiac standard PRS model ##
###############################

setwd("/Users/robryan/Desktop/masters_project/celiac_model")

###########################################
## Account for LD by performing clumping ##
###########################################

## We will only consider SNPs that have a summary statistic --extract argument
# Extract SNPs with summary stats
system("awk '{print $3}' summary_stats.txt > valid.SNP")

# Standard clumping parameters
#
# All p-values
# Reject SNPs in LD (r2) above 0.1
# Clumping window of 250kb

# Perform clumping using plink
system(paste0("./plink \\",
              "--bfile merged_train \\",
              "--clump-p1 1 \\",
              "--clump-r2 0.1 \\",
              "--clump-kb 250 \\",
              "--clump summary_stats.txt \\",
              "--clump-snp-field ID \\",
              "--clump-field P \\",
              "--extract valid.SNP \\", # SNPs with stats
              "--out celiac.train"))
# --clump: 236859 clumps formed from 588601 top variants.
# Results written to celiac.train.clumped .

# Extract SNP ID
system("awk 'NR!=1{print $3}' celiac.train.clumped >  clumped.valid.snp")
# Extract p-value
system("awk '{print $3,$12}' summary_stats.txt > SNP.pvalue")




######################################
## Generate covariate data for both ##
######################################


## Perform pruning before PCA to remove LD correlation
# Sliding "window" = 200 variants
# Eliminate LD above 0.25. 
# Iterate this process by a length of 50 variants per step
system(paste0("./plink \\",
              "--bfile merged_all \\",
              "--indep-pairwise 200 50 0.25 \\",
              "--out all"))
# Pruning complete. 332891 of 676783 variants removed.
# Marker lists written to all.prune.in and all.prune.out .


## PCA
system(paste0("./plink \\",
              "--bfile merged_all \\",
              "--extract all.prune.in \\",
              "--pca 6 \\",
              "--out all"))
# --pca: Results saved to all.eigenval and all.eigenvec


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




###############################
## Merge with phenotype data ##
###############################


fam.cname <- c("FID", "IID", "FID.F", "FID.M", "SEX", "PHENOTYPE")

# Training set
train.fam <- read.table("cd_train_samples.fam", header=F, col.names=fam.cname)
train.data <- merge(train.covar, train.fam[,c(2,6)], by=c("IID"))
train.data$PHENOTYPE <- train.data$PHENOTYPE - 1 # Convert to 0 and 1


# Testing set
test.fam <- read.table("cd_test_samples.fam", header=F, col.names=fam.cname)
test.data <- merge(test.covar, test.fam[,c(2,6)], by=c("IID"))
test.data$PHENOTYPE <- test.data$PHENOTYPE - 1 # Convert to 0 and 1



###################
## Calculate PRS ## 
###################


# Calculate PRS for separate p-value thresholds since optimum
# p-value is unknown
# Create separate files for p-value thresholds
system("echo \"0.001 0 0.001\" > range_list")
system("echo \"0.05 0 0.05\" >> range_list")
system("echo \"0.1 0 0.1\" >> range_list")
system("echo \"0.2 0 0.2\" >> range_list")
system("echo \"0.3 0 0.3\" >> range_list")
system("echo \"0.4 0 0.4\" >> range_list")
system("echo \"0.5 0 0.5\" >> range_list")


## Calculate PRS for each threshold using training dataset
system(paste0("./plink \\",
              "--bfile merged_train \\",
              "--score summary_stats.txt 3 6 13 header \\", # ID, A1, BETA
              "--q-score-range range_list SNP.pvalue \\",
              "--extract clumped.valid.snp \\",
              "--out PRS.train"))


## Get baseline performance of PRS power by first excluding it
# Logistic regression
null.lm <- glm(PHENOTYPE~., family=binomial, data=train.data[,-c(1,2)])
null.pred <- predict.glm(null.lm, newdata=train.data[,3:10], type="response")
# Get mean AUC
null.auc <- as.numeric(AUCBoot(null.pred, train.data$PHENOTYPE)[1])


## Get best p-value threshold
p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5) # P-value range
prs.result <- NULL
for (i in p.threshold){
  score <-  read.table(paste0("PRS.train.",i,".profile"), header=T)[,c(1,2,6)]
  pthres.data <- merge(train.data, score, by=c("FID", "IID"))
  model.p <- glm(PHENOTYPE~., family=binomial, data=pthres.data[,-c(1,2)])
  model.pred <- predict.glm(model.p, newdata=pthres.data[,-c(1,2,11)], 
                            type="response")
  model.auc <- as.numeric(AUCBoot(model.pred, train.data$PHENOTYPE)[1])
  prs.result <- rbind(prs.result, 
                      data.frame("Threshold" = i,
                                 "Total.AUC" = model.auc,
                                 "PRS.AUC" = model.auc - null.auc))
}
prs.result


## Optimal P-value threshold is the lowest one
## This makes intuitive sense since Celiac disease is rare and a subset
## of rare SNPs would be more likely to be causal rather than large numbers
## of common SNPs

## Add optimum p-value PRS to training and testing data
train.data <- merge(train.data, read.table("PRS.train.0.001.profile", 
                                           header=T)[,c(1,2,6)], 
                    by=c("FID", "IID"))




###########################
## Test data performance ##
###########################


# Get PRS for testing data using optimum P-value threshold 
system("echo \"0.001 0 0.001\" >> test_list")
system(paste0("./plink \\",
              "--bfile merged_test \\", # Test .bed file
              "--score summary_stats.txt 3 6 13 header \\", # ID, A1, BETA
              "--q-score-range test_list SNP.pvalue \\", # P: [0, 0.001]
              "--extract clumped.valid.snp \\",
              "--out PRS.test"))
# Results written to PRS.test.0.001.profile.

# Add PRS score to test data
test.data <- merge(test.data, read.table("PRS.test.0.001.profile", 
                                           header=T)[,c(1,2,6)], 
                    by=c("FID", "IID"))


# Train logisitic regression model using training data and optimum P-value PRS
PRS.lm <- glm(PHENOTYPE~., family=binomial, data=train.data[,-c(1,2)])

# Predict on testing data
test.pred <- predict.glm(PRS.lm, 
                         newdata=test.data[,-c(1,2,11)], 
                         type="response")

# Get performance metrics
PRS.metric <- AUCBoot(test.pred, test.data$PHENOTYPE)

# Interestingly, it performs better on test data
# Mean       2.5%      97.5%         Sd 
# 0.84975683 0.82959420 0.86905143 0.01004892 









