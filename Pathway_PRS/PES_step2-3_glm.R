###################### STEP2 ###########################

library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

FILE <- args[1] #your prefix from step1
PHENO <- args[2] #A file containing the IID, FID, LID status and covariates

prs <- as.data.frame(fread(paste0(FILE,".best"), header = T))
pheno_covar <- as.data.frame(fread(PHENO, header = T))

#If LID status is coded 1/2 instead of 0/1 uncomment this and modify the header names below accorging to your pheno file
#pheno_covar <- pheno_covar%>%
#mutate(LID = ifelse(LID == 2, 1, 0)) %>%
#mutate(Sex = ifelse(Sex == 2, 1, 0))

#Create empty list to store the final results
res_cont <- list()
res_quart <- list()

for (set in 4:ncol(prs)){
	prs$Zscore <- (prs[,set] - mean(prs[,set]))/sd(prs[,set])

	prs_clin <- inner_join(prs, pheno_covar, by = c("FID", "IID"))


####################### STEP3 ###########################

	#Change the covariates depending on the ones available between sex, PD_AAO, LDOPA dosage (cumulative or from single time point), HY
	prs_clin <- prs_clin %>%
	dplyr::filter(complete.cases(Sex))%>%
	dplyr::filter(complete.cases(PD_AAO))%>%
	dplyr::filter(complete.cases(HY))%>%
	dplyr::filter(complete.cases(LDOPA_cum))

#CONTINUOUS
#Aggiungi le covariate a seconda della coorte
	sfit <-  summary(glm(LID ~ Zscore + Sex + PD_AAO + LDOPA_cum + HY +  PC1 + PC2 + PC3 + PC4 + PC5, data = prs_clin, family ="binomial"))

	res_cont[[set-3]] <- as.data.frame(sfit$coef[2,])
	names(res_cont[[set-3]]) <- paste0(FILE, "_", names(prs)[set], "_PRS.score")


#QUARTILES
	prs_clin <- prs_clin[!duplicated(prs_clin$Zscore),] #Remove duplicate Zscores to allow equal breaks

	prs_clin$quartile <- with(prs_clin, cut(Zscore, 
                                breaks=quantile(Zscore, probs=seq(0,1, by=0.25), na.rm=TRUE), 
                                include.lowest=TRUE, labels = 1:4))


	fit_quar <- glm(LID ~ as.factor(quartile) + Sex + PD_AAO + LDOPA_cum + HY + PC1 + PC2 + PC3, data = prs_clin, family = "binomial")
	sfit<-summary(fit_quar)

	lid_quartile_fit <- sfit$coefficients
	if (nrow(lid_quartile_fit) >= 4) {
	  lid_quartile_fit <- lid_quartile_fit[2:4,]
	  } else {
	    next
	  }
	lid_quartile_fit <- as.data.frame(lid_quartile_fit)
	names(lid_quartile_fit) <- c("estimate", "std.error", "statistic", "p.value")

	lid_quartile_fit$beta.UB = lid_quartile_fit$estimate+1.96*lid_quartile_fit$std.error
	lid_quartile_fit$beta.LB = lid_quartile_fit$estimate-1.96*lid_quartile_fit$std.error
	lid_quartile_fit$OR = exp(lid_quartile_fit$estimate)
	lid_quartile_fit$OR.UB = exp(lid_quartile_fit$beta.UB)
	lid_quartile_fit$OR.LB = exp(lid_quartile_fit$beta.LB)
	
	lid_quartile_fit$pheno = "LID risk"
	lid_quartile_fit$quartile <- c(2,3,4)
	lid_quartile_fit$set = names(prs)[set]

	res_quart[[set-3]]<-lid_quartile_fit[,c(12,10,11,7,9,8,2,4)] #Final results should have: set,pheno,quartile,OR,OR.LB,OR.UB,std.error,p.value
}

res_quart_bind<-list()
for (x in 1:length(res_quart)){
  res_quart_bind<-rbind(res_quart_bind,res_quart[[x]])
}

res_cont_bind <- list()
for (x in 1:length(res_cont)){
  res_cont_bind<-rbind(res_cont_bind,t(res_cont[[x]]))
}


write.csv(res_cont_bind, paste0("RESULTS_adj_", FILE, ".csv"), quote = F)
write.csv(res_quart_bind,paste0("RESULTS_adj_QUARTILES_",FILE,".csv"), row.names=F, quote=F)
