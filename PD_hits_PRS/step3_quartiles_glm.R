library(data.table)

args <- commandArgs(trailingOnly = TRUE)

FILE <- args[1]

prs_clin <- as.data.frame(fread(paste0(FILE, "_normalized")))

sfit <-  summary(glm(LID ~ Zscore + PD_AAO + Sex + HY + LDOPA_total + PC1 + PC2 + PC3 + PC4 + PC5, data = prs_clin))

res <- as.data.frame(sfit$coef[2,])
names(res) <- paste0(FILE, "_PRS.score")

res <- t(res)

write.csv(res, paste0("RESULTS_", FILE), quote = F)


#QUARTILES

prs_clin$quartile <- with(prs_clin, cut(Zscore, 
                                breaks=quantile(Zscore, probs=seq(0,1, by=0.25), na.rm=TRUE), 
                                include.lowest=TRUE, labels = 1:4))

head(prs_clin)

require(broom)
require(ggplot2)

fit_quar <- glm(prs_clin$LID ~ as.factor(prs_clin$quartile) + PD_AAO + Sex + HY + LDOPA_cum + PC1 + PC2 + PC3 + PC4 + PC5, family = "binomial")
tidy_fitq = tidy(fit_quar)
tidy_fitq$quartile = c(NA,2,3,4)
tidy_fitq = tidy_fitq[2:4,]
tidy_fitq$beta.UB = tidy_fitq$estimate+1.96*tidy_fitq$std.error
tidy_fitq$beta.LB = tidy_fitq$estimate-1.96*tidy_fitq$std.error
tidy_fitq$pheno = "LID risk"


lid_quartile_fit = tidy_fitq
lid_quartile_fit$OR = exp(lid_quartile_fit$estimate)
lid_quartile_fit$OR.UB = exp(lid_quartile_fit$beta.UB)
lid_quartile_fit$OR.LB = exp(lid_quartile_fit$beta.LB)

print("Results")
lid_quartile_fit[,c(9,6,10:12,3,5)]
res<-lid_quartile_fit[,c(9,6,10:12,3,5)]
write.csv(res,paste0("RESULTS_QUARTILES_",FILE,".csv"), row.names=F, quote=F)

png(paste0(FILE,"_LID.risk_PRS_Quartiles.png"), width = 5.5, height = 4, units = "in", res = 300)

x = ggplot(lid_quartile_fit, aes(x=as.factor(quartile), y=exp(estimate))) + 
  geom_point(size=3, color="Darkred") + 
  geom_errorbar(aes(ymin=OR.LB, ymax=OR.UB),position=position_dodge(.2), width = 0.1, color="Darkred")
 
x + labs(title = "Risk for LID in PD by PRS Quartile", x="PRS Quartile", y="Odds ratio & 95% CI") +
  theme_bw() + scale_color_manual(values=c('red4', '#999999', 'darkblue')) + ylim(0,5.5) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

dev.off()
