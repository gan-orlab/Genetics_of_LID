#!/bin/bash

#Run all the steps together 

TARGET=$1
BASE=$2
PHENO_glm=$3
PHENO_lm=$4
PHENO_cox=$5
OUT=$6

bash step1_run.prs_glm.sh $TARGET $BASE ${OUT}_glm
Rscript step2_prs2zscore_glm.R ${OUT}_glm ${OUT}_glm $PHENO_glm --no-save
Rscript step3_quartiles_glm.R ${OUT}_glm --no-save

bash step1_run.prs_cox.sh $TARGET $BASE ${OUT}_cox
Rscript step2_prs2zscore_cox.R ${OUT}_cox ${OUT}_cox $PHENO_cox --no-save
Rscript step3_quartiles_cox.R ${OUT}_cox --no-save
