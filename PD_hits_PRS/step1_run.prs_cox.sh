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
--pheno covar_cox.txt \
--pheno-col LID \
--cov covar_cox.txt \
--cov-col @PC[1-5],Sex,PD_AAO,LDOPA_total,HY \
--score avg \
--out $PREFIX \
--print-snp \
--fastscore \
--bar-level 1 \
--no-regress 
