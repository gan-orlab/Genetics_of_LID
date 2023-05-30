library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

FILE <- args[1] #your prefix from step1
OUT <- args[2] 
PHENO <- args[3] #A file containing the IID, FID and LID status

prs <- as.data.frame(fread(paste0(FILE,".best"), header = T))
pheno_covar <- as.data.frame(fread(PHENO, header = T))

names(prs)[4] <- "PRS"

#Normalize PRS scores
prs$Zscore <- (prs$PRS - mean(prs$PRS))/sd(prs$PRS)

#If LID status is coded 1/2 instead of 0/1 uncomment this and modify the header names below accorging to your pheno file 
#pheno_covar <- pheno_covar%>%
#mutate(LID = ifelse(LID == 2, 1, 0))

prs_pheno_covar <- inner_join(prs, pheno_covar, by = c("FID", "IID"))

write.csv(prs_pheno_covar, paste0(OUT, "_normalized"), quote = F, row.names = F)
