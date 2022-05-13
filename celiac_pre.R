###################################
## Celiac dataset pre-processing ##
###################################


##################
## Training set ##
##################

## Generate a merged bed/bim/fam file with the training set
## bed_list.txt is a text file that contains names of all separate
## chromosome .bed file

## Read training fam file
train.fam <- read.table("cd_train_samples.fam")
## Need to read the cd_chr1.fam to get the both the FID and IID
chr1.fam <- read.table("cd_chr1.fam")
## Matching training IID to chr1.fam IID and also extract FID
train.id <- chr1.fam[which(chr1.fam$V2 %in% train.fam$V2),
                     c(1,2)]
## Write to file for --keep argument
write.table(train.id , "train_id.txt", row.names=F, col.names=F, quote=F)

## Merge all the bed files containing only the training samples
system(paste0("./plink \\",
              "--bfile cd_chr1 \\",
              "--merge-list bed_list.txt \\",
              "--keep train_id.txt \\",
              "--make-bed \\",
              "--out merged_train"))


## Create rds file of merged bed file
snp_readBed2("merged_train.bed")



##############
## Test set ##
##############

## Generate a merged bed/bim/fam file with the test set

## Read test fam file
test.fam <- read.table("cd_test_samples.fam")
## Matching test IID to chr1.fam IID and also extract FID
test.id <- chr1.fam[which(chr1.fam$V2 %in% test.fam$V2),
                    c(1,2)]
## Write to file for --keep argument
write.table(test.id , "test_id.txt", row.names=F, col.names=F, quote=F)

## Merge all the bed files containing only the training samples
system(paste0("./plink \\",
              "--bfile cd_chr1 \\",
              "--merge-list bed_list.txt \\",
              "--keep test_id.txt \\",
              "--make-bed \\",
              "--out merged_test"))

## Create rds file of merged_test bed file
snp_readBed2("merged_test.bed")



##################################
## Pre-processing summary files ##
##################################


## Load in the 22 separate summary files and merge
c.names <- c("CHROM",	"POS", "ID", "REF",	"ALT", "A1", "TEST", "OBS_CT", 
             "OR", "LOG_OR_SE", "Z_STAT", "P")
summary.all <- NULL
for (i in 1:22){
  file = paste0("cd_chr",i,"_summary.PHENO2.glm.logistic")
  raw.file = read.table(file, header=FALSE)
  colnames(raw.file) <- c.names
  raw.file = raw.file[which(raw.file$TEST=="ADD"),]
  summary.all <- rbind(summary.all, raw.file)
}

## Add a beta score field (ln(OR))
summary.all$BETA <- log(summary.all$OR)

## Write into a text file
write.table(summary.all, "summary_stats.txt", quote=F, row.names=F)



#######################
## Combined BED file ##
#######################

## Merge all the bed files
system(paste0("./plink \\",
              "--bfile cd_chr1 \\",
              "--merge-list bed_list.txt \\",
              "--make-bed \\",
              "--out merged_all"))
# --make-bed to merged_all.bed + merged_all.bim + merged_all.fam ... done.


## Create rds file of merged_all bed file
snp_readBed2("merged_all.bed")



