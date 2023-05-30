library(dplyr)
library(tibble)
library(ggplot2)
library(factoextra)
library(ggpubr)
library(MASS)
library(survival)

#########################################Database preparation########################################

#Clinical data
orig_cohort<-read.csv("cohort_whole.csv", fileEncoding = "UTF-8-BOM")

cohort_clin.orig<-orig.cohort%>%
  mutate(Disease_duration = Age_last_follow.up - PD_AAO)%>%
  mutate(LID_onset.yrs = ifelse(LID==0, Age_last_follow.up - LDOPA_start, LID_onset - LDOPA_start))

#Evaluate how many PCs have to be included in the analyses out of top 20 
cohort_PC<-read.csv("plink.eigenvec.csv", fileEncoding="UTF-8-BOM")

#Scree plot to select PCs
eigen<- read.csv("plink.eigenval.csv")
eigen<-slice(eigen, 1:10)
colnames(eigen)="eigenvals"

plot(
    x = seq(nrow(eigen)), y = eigen$eigenvals,
    type = "o",
    xlab = "Principal Component", ylab = "Variance")

#Select relevant columns for adjusted analyses (in this example PC1-8 were selected, based on the scree plot)
PCs_headers<-paste0("PC", 1:8)
COVs_headers<-c("PD_AAO","Sex","LDOPA_cum","LEDD_cum", "BMI","DA", "HY", "Disease_duration")
headers<-c("FID","IID","LID",COVs_headers,PCs_headers)

cohort_clin.adj<-cohort_clin[,headers]

#Genetic data
cohort_snps<-read.csv("cohort_genotypes.csv", fileEncoding = "UTF-8-BOM")

SNPs<-cohort_snps[,-c(1,2)]%>%colnames

#Merge genetic and clinical data
cohort_noadj<-inner_join(cohort_clin,cohort_snps, by = c("IID","FID"))

cohort_adj<-cohort_clin.adj%>%
  na.omit%>%
  inner_join(.,cohort_snps, by = c("FID","IID"))


###################################Logistic regression without adjustments#####################################

glm.noadj<-function(x){
  SNPs.df<-cohort_noadj[,SNPs[x]]
  fit<-glm(LID ~ SNPs.df + PC1 + PC2 + PC3, data = cohort_noadj, family="binomial")
  OR<-exp(cbind("OR" = coef(fit), confint.default(fit, level = 0.95)))
  OR[2,]<-round(as.numeric(OR[2,]),4)
  sfit<-summary(fit)
  sfit$coef[2,]<-round(as.numeric(sfit$coef[2,]),4)
  #N of samples without missing data
  n_samples<-SNPs.df%>%
    na.omit%>%
    length
  #Minor allele frequency
  maf<-sum(SNPs.df, na.rm = T)/(n_samples*2)
  maf<-round(maf,4)
  #Results
  coef<-c(SNPs[x],sfit$coef[2,],OR[2,],n_samples,maf)
}

#Results
res<-data.frame(lapply(1:length(SNPs), glm.noadj))
colnames(res)<-res[1,]
res<-res[-1,]
rownames(res)<-c("Beta","SE","z","p-value", "OR","2.5%","97.5%","N","MAF")
res<-as.data.frame(t(res))

#Order by smallest p-value
res<-res[order(res[,4]),]

write.csv(res, "cohort_log.reg.noadj_RESULTS.csv")


#########################################Logistic regression with adjustments##################################

#Covariates included in the adjusted analyses, when available, are the following:

#    PD age at onset (AAO)
#    PD duration*
#    Sex
#    Principal components
#    BMI
#    Cumulative L-dopa dosage* (if not available, L-dopa dose at last follow-up)
#    Cumulative Levodopa Equivalent Daily Dose* (LEDD) (if not available, LEDD dose at last follow-up)
#    Dopamine agonist use
#    Hoehn & Yahr score at last follow-up

#*PD duration is defined as the time between PD AAO and LID onset or last follow-up if LID is not present.

#*Cumulative L-dopa/LEDD is defined as the L-dopa dosage/LEDD multiplied by the time between the start of L-dopa/PD therapy and LID onset or last follow-up if LID is not present.

#Covariates passed through stepAIC
PCs<-paste0("PC", 4:8, collapse = " + ")
COVs<-paste0("PD_AAO + Sex + LDOPA_cum + LEDD_cum + DA + HY + Disease_duration + BMI + ",PCs, collapse = "") 

glm.adj<-function(x){
  SNPs.df<-cohort_adj[,SNPs[x]]
  base<-glm(LID ~ SNPs.df + PC1 + PC2 + PC3, data = cohort_adj, family = "binomial")
  res.AIC = stepAIC(base, scope =list(lower=formula(base), upper=paste("~. +", COVs), direction="both"), trace = F) #Select the covariates to include in the model with a stepwise regression based on Akaike Information Criterion
  fit<-glm(formula(res.AIC), data = cohort_adj, family = "binomial")
  OR<-exp(cbind("OR" = coef(fit), confint.default(fit, level = 0.95)))
  OR[2,]<-round(as.numeric(OR[2,]),4)
  sfit<-summary(fit)
  sfit$coef[2,]<-round(as.numeric(sfit$coef[2,]),4)
  #N of samples without missing data
  n_samples<-SNPs.df%>%
    na.omit%>%
    length
  #Minor allele frequency
  maf<-sum(SNPs.df, na.rm = T)/(n_samples*2)
  maf<-round(maf,4)
  #Results
  coef<-c(SNPs[x],sfit$coef[2,],OR[2,],n_samples,maf)
}


#Results and final table
res<-data.frame(lapply(1:length(SNPs), glm.adj))
colnames(res)<-res[1,]
res<-res[-1,]
rownames(res)<-c("Beta","SE","z","p-value", "OR","2.5%","97.5%", "N","MAF")

res<-as.data.frame(t(res))

#Order by smallest p-value
res<-res[order(res[,4]),]

write.csv(res, "cohort_log.reg.adj_RESULTS.csv")


###############################################Cox regression without adjustments###############################################
 
COVs<-paste0("PD_AAO + Sex + LDOPA_total + LEDD_total + DA + HY + Disease_duration + BMI + ",PCs, collapse = "") 

cohort_cox.noadj<-cohort_noadj%>%
  filter(complete.cases(LID_onset.yrs))

cox.noadj<-function(x){
  SNPs.df<-cohort_cox.noadj[,SNPs[x]]
  fit<-coxph(Surv(LID_onset.yrs, LID) ~ SNPs.df + PC1 + PC2 + PC3, data = cohort_cox.noadj)
  sfit<-summary(fit)
  coef_CI<-c(sfit$coef[1,], sfit$conf.int[1,])
  coef_CI<-round(as.numeric(coef_CI),4)
  #N of samples without missing data
  n_samples<-SNPs.df%>%
    na.omit%>%
    length
  #Minor allele frequency
  maf<-sum(SNPs.df, na.rm = T)/(n_samples*2)
  maf<-round(maf,4)
  #Results
  coef<-c(SNPs[x],coef_CI,n_samples,maf)
}

#Results and final table
res<-data.frame(lapply(1:length(SNPs), cox.noadj))
colnames(res)<-res[1,]
res<-res[-c(1,7,8),]
rownames(res)<-c("Beta","HR","SE","z", "p-value","2.5%","97.5%","N","MAF")

#Order by smallest p-value
res<-res[order(res[5,])]

res<-as.data.frame(t(res))

write.csv(res, "cohort_cox.noadj_RESULTS.csv")

############################################Cox regression with adjustments##############################################

cohort_cox.adj<-cohort_adj%>%
  filter(LDOPA_start>0)%>%   #Filter out missing LDOPA_start (0)
  #Calculate LID_onset in years from LDOPA start
  #For LID- just use the L.DOPA_duration 
  mutate(LID_onset.yrs = case_when(LID_onset>1 ~ as.integer(LID_onset - LDOPA_start), 
                                   TRUE ~ as.integer(LDOPA_duration) ))%>%
  filter(complete.cases(LID_onset.yrs))

cox.adj<-function(x){
  SNPs.df<-cohort_cox.adj[,SNPs[x]]
  surv.obj<-with(cohort_cox.adj, Surv(LID_onset.yrs, LID))
  base<-coxph(surv.obj ~ SNPs.df + PC1 + PC2 + PC3, data = cohort_cox.adj)
  res.AIC = stepAIC(base, scope = list(lower=formula(base), upper=paste("~. +", COVs), direction="both"), trace = F)
  fit<-coxph(formula(res.AIC), data = cohort_cox.adj)
  sfit<-summary(fit)
  coef_CI<-c(sfit$coef[1,], sfit$conf.int[1,])
  coef_CI<-round(as.numeric(coef_CI),4)
  #N of samples without missing data
  n_samples<-SNPs.df%>%
    na.omit%>%
    length
  #Minor allele frequency
  maf<-sum(SNPs.df, na.rm = T)/(n_samples*2)
  maf<-round(maf,4)
  #Results
  coef<-c(SNPs[x],coef_CI,n_samples,maf)
}


#Results and final table
res<-data.frame(lapply(1:length(SNPs), cox.adj))
colnames(res)<-res[1,]
res<-res[-c(1,7,8),]
rownames(res)<-c("Beta","HR","SE","z", "p-value","2.5%","97.5%", "N","MAF")

#Order by smallest p-value
res<-res[order(as.numeric(res[5,]))]

res<-as.data.frame(t(res))

write.csv(res, "cohort_cox.adj_RESULTS.csv")



