# Genetics_of_LID

## Analyses on GBA1 and LRRK2 risk variants.

The carrier status of the single GBA1 and LRRK2 risk variants were collapsed into single variables, namely the carrier status of GBA1/LRRK2 of any of the included GBA1/LRRK2 risk variants. Such variants include T369M, E326K, N370S and L444P for GBA1; G2019S, M1646T and R1441C/G for LRRK2.

We evaluated the association between GBA1/LRRK2 carrier status and the risk of LID/time to LID

Time to LID (“LID_onset.yrs”) is defined as the time between the start of L-dopa therapy and LID onset or, when LID is not present, last follow-up.

In the adjusted analyses we used the cumulative levodopa dosage/LEDD for logistic regression (LDOPA_total/LEDD_total) and the last time point levodopa dosage/LEDD (LDOPA_cum/LEDD_cum) for Cox regression. 

## PRS using PD hits
To run the PRS analyses (step 1 below) you will need to download the PRSice software:
https://www.prsice.info/ and have PRSice.R and PRSice_linux in the same folder
of your analyses

The PRS includes the PD GWAS hits from Nalls et al.,2019 (PMID: 31701892)

For the logistic and Cox regression, you'll need to run the following
scripts:
- step1_run.prs.sh-->to obtain PRS for each individual
- step2_prs2zscore.R-->to normalize PRS into Zscores
- step3_quartiles.R-->to run the analyses between the individual
Zscores and LID
The all_steps.sh script will run all the steps at once for you, you just need
to specify the arguments.

Input files for the step1
- base file: the PD GWAS summary statistics (taken from Nalls et al.)
- target file: your cohort bfiles recoded by sex and LID status/time to LID
for logistic regression and Cox regression, respectively
Input file for the step2:
- covariate file (used also as the argument of --pheno) with LID status, LID_onset.yrs, sex, PD_AAO...

For target bfiles make sure they have been QC'd, otherwise please follow the 
PRSice instructions on the necessary QC before running the analyses
(https://www.prsice.info/quick_start/#input)

The results from these analyses will be:
- The results table of continuous Zscores
- The results table of quartiles
- The quartiles plot

## Dopamine pathway PRS

To run the PRS analyses (step 1 below) you will need to download the PRSice software:
https://choishingwan.github.io/PRSice/ and have PRSice.R and PRSice_linux in the same folder
of your analyses

The pathways PRS refers to the dopaminergic transmission pathway specified in the "dopamine_pathways_and_LIDgenes.gmt" file, which refers to the gencode.v40lift37.annotation.gtf for annotation (http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh37_mapping/gencode.v40lift37.annotation.gtf.gz).

For each analysis you'll need to run three steps:
- step1-->to obtain PRS for each individual
- step2-->to normalize PRS into Zscores
- step3-->to run the regression between the individual
Zscores and the manifestation of RBD

Steps 2 and 3 are incorporated into a single script.

Input files for the step1
- base file: the meta-analysis summary statistics I provided to the link
- target file: your cohort bfiles recoded by phenotype (LID status) and sex
Input file for the step2:
- covariate file (used also as the argument of --pheno) with LID status, sex, PD_AAO
- gmt and gtf file, explained above

For step1 the only argument is the name of your cohort (which will be used to retrieve multiple files like your target files and also to name your final output, so make sure the name is coherent). 

For step2-3 your first argument is the same as the one of step1, but you will add a suffix depending on your analysis output (specified in step1 if you need a reference). For example, if your cohort is mcgill and your analysis is logistic regression your first argument will be mcgill_glm. Your second argument is your file name containing covariates (Sex, PD_AAO...) and phenotype (i.e., LID).

For target bfiles make sure they have been QC'd, otherwise please follow the 
PRSice instructions on the necessary QC before running the analyses
(https://www.prsice.info/quick_start/#input)

The results from these analyses will be:
- The results csv table of continuous Zscores
- The results csv table of quartiles
