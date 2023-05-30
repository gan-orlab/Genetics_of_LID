#!/bin/bash

module load r/4.0.0

TARGET=$1
BASE=$2
PREFIX=$3

Rscript PRSice.R \
--prsice PRSice_linux \
--target $TARGET \
--base $BASE \
--beta \
--stat b \
--snp SNP \
--A1 A1 \
--A2 A2 \
--pvalue p \
--pheno covar_glm.txt \
--pheno-col LID \
--cov covar_glm.txt \
--cov-col @PC[1-5],Sex,PD_AAO,LDOPA_cum,HY \
--score avg \
--binary-target T \
--prevalence 0.4 \
--out $PREFIX \
--print-snp \
--perm 10000 \
