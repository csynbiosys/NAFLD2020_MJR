##K-M curves for NAFLD genetic variants

##Survival
svglite("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results paper/Genotype/TM6SF2_mortality.svg", width = 9, height = 7)
fit_tm6sf2 <- survfit(Surv(survival_time_days, survival_status) ~ rs58542926, data = survival_decomp)
ggsurvplot(fit_tm6sf2, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
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
svglite("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results paper/Genotype/HSD17B13_mortality.svg", width = 9, height = 7)
fit_hsd17b13_1 <- survfit(Surv(survival_time_days, survival_status) ~ rs62305723, data = survival_decomp)
ggsurvplot(fit_hsd17b13_1, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
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
svglite("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results paper/Genotype/MBOAT7_mortality.svg", width = 9, height = 7)
fit_mboat7 <- survfit(Surv(survival_time_days, survival_status) ~ rs641738, data = survival_decomp)
ggsurvplot(fit_mboat7, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
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
svglite("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results paper/Genotype/HSD17B13_mortality.svg", width = 9, height = 7)
fit_hsd17b13_2 <- survfit(Surv(survival_time_days, survival_status) ~ rs72613567, data = survival_decomp)
ggsurvplot(fit_hsd17b13_2, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
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
cox_hsd17b13 <- coxph(Surv(survival_time_days, survival_status) ~ rs72613567, data = survival_decomp) #Fits a cox model p=0.0126
svglite("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results paper/Genotype/PNPLA3_mortality.svg", width = 9, height = 7)
fit_pnpla3 <- survfit(Surv(survival_time_days, survival_status) ~ rs738409, data = survival_decomp)
ggsurvplot(fit_pnpla3, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
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
cox_pnpla3 <- coxph(Surv(survival_time_days, survival_status) ~ rs738409, data = survival_decomp) #Fits cox model too p=0.004367

##Decompensation events

survival_decomp_2 <- survival_decomp[which(survival_decomp$decomp_time_days_compRisk >= 0),] #Remove negative events 
fit_tm6sf2_d <- survfit(Surv(decomp_time_days_compRisk, decompStatus) ~ rs58542926, data = survival_decomp_2)
svglite("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results paper/Genotype/TM6SF2_decomp.svg", width = 9, height = 7)
ggsurvplot(fit_tm6sf2_d, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
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
svglite("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results paper/Genotype/HSD17B13_decomp.svg", width = 9, height = 7)
fit_hsd17b13_1_d <- survfit(Surv(decomp_time_days_compRisk, decompStatus) ~ rs62305723, data = survival_decomp_2)
ggsurvplot(fit_hsd17b13_1_d, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
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
svglite("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results paper/Genotype/MBOAT7_decomp.svg", width = 9, height = 7)
fit_mboat7_d <- survfit(Surv(decomp_time_days_compRisk, decompStatus) ~ rs641738, data = survival_decomp_2)
ggsurvplot(fit_mboat7_d, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
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
svglite("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results paper/Genotype/HS17B13_2_decomp.svg", width = 9, height = 7)
fit_hsd17b13_2_d <- survfit(Surv(decomp_time_days_compRisk, decompStatus) ~ rs72613567, data = survival_decomp_2)
ggsurvplot(fit_hsd17b13_2_d, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
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
svglite("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results paper/Genotype/PNPLA3_decomp.svg", width = 9, height = 7)
fit_pnpla3_d <- survfit(Surv(decomp_time_days_compRisk, decompStatus) ~ rs738409, data = survival_decomp_2)
ggsurvplot(fit_pnpla3_d, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
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

#HCC

survival_decomp_3 <- survival_decomp[which(survival_decomp$HCCtime >= 0), ] #Remove negative events
fit_tm6sf2_h <- survfit(Surv(HCCtime, HCCstatus) ~ rs58542926, data = survival_decomp_3)
svglite("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results paper/Genotype/TM6SF2_hcc.svg", width = 9, height = 7)
ggsurvplot(fit_tm6sf2_h, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
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
svglite("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results paper/Genotype/HSD17B13_hcc.svg", width = 9, height = 7)
fit_hsd17b13_1_h <- survfit(Surv(HCCtime, HCCstatus) ~ rs62305723, data = survival_decomp_3)
ggsurvplot(fit_hsd17b13_1_h, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
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
svglite("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results paper/Genotype/MBOAT7_hcc.svg", width = 9, height = 7)
fit_mboat7_h <- survfit(Surv(HCCtime, HCCstatus) ~ rs641738, data = survival_decomp_3)
ggsurvplot(fit_mboat7_h, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
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
svglite("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results paper/Genotype/HS17B13_2_hcc.svg", width = 9, height = 7)
fit_hsd17b13_2_h <- survfit(Surv(HCCtime, HCCstatus) ~ rs72613567, data = survival_decomp_3)
ggsurvplot(fit_hsd17b13_2_h, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
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
svglite("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results paper/Genotype/PNPLA3_hcc.svg", width = 9, height = 7)
fit_pnpla3_h <- survfit(Surv(HCCstatus, HCCtime) ~ rs738409, data = survival_decomp_3)
ggsurvplot(fit_pnpla3_h, conf.int = FALSE, pval = TRUE, risk.table = TRUE, fun = "event", 
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