
library(ggplot2)
library(survival)
library(survminer)
library(svglite)

#Upload csv file with the information regarding survival and decompensation events as survival_decomp, where clinical information is presented as columns and patients as rows. 
#In this case, this data frame contains a column called kleiner_fibrosis_combined that indicates the fibrotic stage of each patient, that will be used to study how the evens develop as the fibrosis increases.
# 'survival_' shows time and status for death outcomes. 'decomp_time_days_compRisk' indicates the time where the event was registered. 'decompStatus' indicates the status for decompesantion events for each patient. 'HCC' shows time and status for HCC events.

##Survival
svglite("./Mortality.svg", width = 9, height = 7)
fit_mortality <- survfit(Surv(survival_time_days, survival_status) ~ kleiner_fibrosis_combined, data = survival_decomp)
ggsurvplot(fit_mortality, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
           cumevents = TRUE, ggtheme = theme_minimal(), 
           legend = "right", pval.coord = c(500, 0.6), 
           risk.table.title = "Number at risk", 
           cumevents.title = "Cumulative number of events", 
           risk.table.y.text = FALSE, cumevents.y.text = FALSE, 
           risk.table.height = 0.18, cumevents.height = 0.18, 
           xlab = "Time to death (days)", 
           ylab = 'Cumulative probability', 
           font.x = c(16), font.y = c(16), font.legend = c(14), 
           fontsize = 3, tables.theme = theme_cleantable(base_family = "Calibri"))
dev.off()

#Performs multivariate Cox regression to ensure there is a significant diffence.
cox_mortality <- coxph(Surv(survival_time_days, survival_status) ~ kleiner_fibrosis_combined, data = survival_decomp)

##Decompensation events

survival_decomp_2 <- survival_decomp[which(survival_decomp$decomp_time_days_compRisk >= 0),] #Remove negative events 
fit_decomp <- survfit(Surv(decomp_time_days_compRisk, decompStatus) ~ kleiner_fibrosis_combined, data = survival_decomp_2)
svglite("./Decomp.svg", width = 9, height = 7)
ggsurvplot(fit_decomp, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
           cumevents = TRUE, ggtheme = theme_minimal(), 
           legend = "right", pval.coord = c(500, 0.6), 
           risk.table.title = "Number at risk", 
           cumevents.title = "Cumulative number of events", 
           risk.table.y.text = FALSE, cumevents.y.text = FALSE, 
           risk.table.height = 0.18, cumevents.height = 0.18, 
           xlab = "Time to decompensation event (days)", 
           ylab = 'Cumulative probability', 
           font.x = c(16), font.y = c(16), font.legend = c(14), 
           fontsize = 3, tables.theme = theme_cleantable(base_family = "Calibri"))
dev.off()
cox_decomp <- coxph(Surv(survival_time_days, survival_status) ~ kleiner_fibrosis_combined, data = survival_decomp_2)

#HCC

survival_decomp_3 <- survival_decomp[which(survival_decomp$HCCtime >= 0), ] #Remove negative events
fit_hcc <- survfit(Surv(HCCtime, HCCstatus) ~ rs58542926, data = survival_decomp_3)
svglite("./HCC.svg", width = 9, height = 7)
ggsurvplot(fit_hcc, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
           cumevents = TRUE, ggtheme = theme_minimal(), 
           legend = "right", pval.coord = c(500, 0.6), 
           risk.table.title = "Number at risk", 
           cumevents.title = "Cumulative number of events", 
           risk.table.y.text = FALSE, cumevents.y.text = FALSE, 
           risk.table.height = 0.18, cumevents.height = 0.18, 
           xlab = "Time to event (days)", 
           ylab = 'Cumulative probability', 
           font.x = c(16), font.y = c(16), font.legend = c(14), 
           fontsize = 3, tables.theme = theme_cleantable(base_family = "Calibri"))
dev.off()
