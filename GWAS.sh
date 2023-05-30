#!/bin/bash

#LARGE-SCALE GWAS

#############################################Logistic regression#############################################

#Select the samples to include in logisitc regression
plink --bfile Cohort --keep keep_glm.adj.txt --make-bed --out Cohort_glm.adj

plink --bfile Cohort_glm.adj \
--logistic --ci 0.95 --freq \
--hide-covar --covar covar_glm.adj.txt \
--covar-name PD_AAO,Sex,LDOPA_cum,LEDD_cum,DA,HY,Disease_duration,BMI,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--out temp.ADJ

#Order by significance and remove NAs
grep 'BP\|ADD' temp.ADJ.assoc.logistic | sort -gk12 > RESULTS.ADJ.assoc.logistic
sed -i '/NA/d' RESULTS.ADJ.assoc.logistic


# Cox regression is performed with the software SurvialGWAS_SV for internal computational reasons, but it can be performed also with coxph as with the other analyses
# Details on SurvivalGWAS_SV: https://www.liverpool.ac.uk/population-health/research/groups/statistical-genetics/survival-gwas-sv

#######################################Cox regression###########################################

mkdir Cox

#Select just samples with clinical data
cut -f 1,2 covar_cox.adj.txt > keep_cox.adj.txt
plink --bfile Cohort --keep keep_cox.adj.txt --make-bed --out Cohort_cox.adj

#Convert bfile into .gen
plink --bfile Cohort_cox.adj --recode oxford --out Cohort_cox.adj


#Analysis
str1=0 #Start position in genotype file
str=320000 #Number of SNPs/lines in genotype file
no_of_jobs=100 #Number of cores
inc=`expr \( $str - $str1 \) \/ $no_of_jobs` #Increment


for j in {1..100}
do
        nstart=`expr \( $j - 1 \) \* $inc`
        nstop=`expr $nstart + $inc - 1`
        SurvivalGWAS_SV -gf=Cohort_cox.adj.gen -sf=covar_cox.adj.txt -threads=40 -censor=lid -time=LID_onset.yrs -cov=PD_AAO,Sex,LDOPA_cum,LEDD_cum,DA,HY,Disease_duration,BMI,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 -lstart=$nstart -lstop=$nstop -method=cox -p=onlysnp -o=Cox/cox.adj.${j}.output
done


#Concatenate 100 files into one
for j in {1..100}
do
        cat Cox/cox.adj.${j}.output >> Cox/cox.adj.MERGED.output
done

#Remove headers
sed -i '/InputName/d' Cox/cox.adj.MERGED.output

#Remove duplicates and add just header at the beginning of the file
sort -k2 Cox/cox.adj.MERGED.output | uniq > Cox/temp
sort -gk11 Cox/temp > Cox/temp2
mv Cox/temp2 Cox/cox.adj.MERGED.output


head -1 Cox/cox.adj.1.output > Cox/header.txt
cat Cox/header.txt Cox/cox.adj.MERGED.output > Cox/temp
mv Cox/temp Cox/cox.adj.MERGED.output

mkdir Cox/adj
mv Cox/cox.adj.{1..100}.output Cox/adj
