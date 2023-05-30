library(data.table)
library(survival)
require(broom)

args <- commandArgs(trailingOnly = TRUE)

FILE <- args[1]

prs_clin <- as.data.frame(fread(paste0(FILE, "_normalized")))

fit <-  coxph(Surv(LID_onset.yrs, LID) ~ Zscore + PD_AAO + Sex + HY + LDOPA_total + PC1 + PC2 + PC3 + PC4 + PC5, data = prs_clin)
sfit <- summary(fit)

res <- as.data.frame(sfit$coef[1,])
names(res) <- paste0(FILE, "_PRS.score")
res<-t(res)

head(res)

write.csv(res, paste0("RESULTS_", FILE), quote = F)


#QUARTILES

prs_clin$quartile <- with(prs_clin, cut(Zscore, 
                                breaks=quantile(Zscore, probs=seq(0,1, by=0.25), na.rm=TRUE), 
                                include.lowest=TRUE, labels = 1:4))

head(prs_clin)

require(ggplot2)

fit_quar <- coxph(Surv(prs_clin$LID_onset.yrs,LID) ~ as.factor(prs_clin$quartile) + PD_AAO + Sex + HY + LDOPA_total + PC1 + PC2 + PC3 + PC4 + PC5, data = prs_clin)
tidy_fitq = tidy(fit_quar, exponentiate = T, conf.int = T)
tidy_fitq$quartile = c(2,3,4)
colnames(tidy_fitq)[2] <- "HR"
tidy_fitq$term = "LID progression"
print("tidy_fitq")
tidy_fitq
lid_quartile_fit = tidy_fitq


print("Results")
res<-lid_quartile_fit[,c(1,8,2,6,7,3,5)]
res
write.csv(res,paste0("RESULTS_QUARTILES_",FILE,".csv"), row.names=F, quote=F)

png(paste0(FILE,"_LID.Progression_PRS_Quartiles.png"), width = 5.5, height = 4, units = "in", res = 300)

x = ggplot(lid_quartile_fit, aes(x=as.factor(quartile), y=HR)) + 
  geom_point(size=3, color="Darkred") + 
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high),position=position_dodge(.2), width = 0.1, color="Darkred")
 
x + labs(title = "LID progression in PD by PD PRS Quartile", x="PD PRS Quartile", y="Hazard ratio & 95% CI") +
  theme_bw() + scale_color_manual(values=c('red4', '#999999', 'darkblue')) + ylim(0,5.5) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

dev.off()
