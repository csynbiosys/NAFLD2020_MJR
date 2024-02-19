#Script with the different plots you can obtain

source("./TRS_model.R")
library(ggalluvial)
library(tidyr)

#Plot for cross-validation of lambda
plot(fit, xvar = "lambda")
abline(v = log(lambda.min), col = "red")

#Plot for coefficients
plot(fit, xvar = "lambda")
abline(v = log(lambda.min), col = "red")
abline(v = log(lambda_1se), col = "red")

#Survplot for decompensation events when dividing patients with new score
ggsurvplot(fit2, conf.int = TRUE, pval = TRUE, risk.table = TRUE, fun = "event", 
           cumevents = TRUE, ggtheme = theme_minimal(), 
           legend = "right", pval.coord = c(500, 0.6), 
           risk.table.title = "Number at risk", 
           cumevents.title = "Cumulative number of events", 
           risk.table.y.text = FALSE, cumevents.y.text = FALSE, 
           risk.table.height = 0.18, cumevents.height = 0.18, 
           xlab = "Time to decompensation events (days)", 
           ylab = 'Cumulative probability', 
           font.x = c(16), font.y = c(16), font.legend = c(14), 
           fontsize = 3, tables.theme = theme_cleantable(base_family = "Calibri"))

#Plots for gene expression of each gene across disease stage
ggplot(genes_expr, aes(groups, counts, fill = groups)) + 
  geom_boxplot(alpha = 0.6, outlier.shape = NA) + theme_bw() +
  stat_compare_means(comparisons = comparisons, label = "p.signif") + facet_wrap(gene ~.)

#Alluvial plot 
df <- decomp[, c(20,18, 3)]
names(df) <- c("Risk", "Fibrotic group", "Hepatocyte ballooning")
df <- df %>% mutate(Risk = case_when(Risk == "high"  ~ "High", TRUE ~ "Low"))
df <- df %>% mutate(`Fibrotic group` = case_when(`Fibrotic group` == "F0F1" ~ "F0/F1", TRUE ~ "F3/F4"))
freq <- to_lodes_form(df)
freq$x <- factor(freq$x, levels= c("Fibrotic group", "Risk", "Hepatocyte ballooning"))

ggplot(data = freq,
      aes(x = x, stratum = stratum, alluvium = alluvium)) +
  geom_flow(aes(fill = stratum)) +
  geom_stratum(aes(fill = stratum)) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  xlab("") +
  ylab("Number of patients") + 
  theme_minimal() +
  theme(legend.position = "none")

#Plot for AUC
plot(roc_surv, time = 1825, title = FALSE, col = "red")
plot(roc_surv, time = 1095, add = TRUE, col = "green")
plot(roc_surv, time = 365, add = TRUE, col = "blue")
legend("bottomright",c("AUC in 5 years = 0.83", "AUC in 3 years = 0.81", "AUC in 1 year = 0.86"),
       col=c("red","green", "blue"),lty= 1, cex = 0.65)

#Heatmap for the genes in the TRS highlighting the expression across patients
heatmap_counts <- counts_multi[active.index]
risk_scores <- risk_scores[match(decomp$STUDY_NO, rownames(risk_scores)),]
decomp$risks <- risk_scores$risk 
rownames(decomp) <- decomp$STUDY_NO #It needs sample_id in rownames!
annotations <- decomp[, c(3,5,7,10,20,21)]
colnames(annotations) <- c("NAS ballooning", "NAS steatosis", "NAS inflammation", "NASH-CRN fibrosis", "TRS-defined risk status", "Sex")
annotations <- annotations %>% mutate(`NASH-CRN fibrosis` = case_when(`NASH-CRN fibrosis` == "1a" | `NASH-CRN fibrosis` == "1b" | `NASH-CRN fibrosis` == "1c" ~ "1", TRUE ~ `NASH-CRN fibrosis`))
annotations$`TRS-defined risk status` <- factor(annotations$`TRS-defined risk status`, levels = c("low", "high"))
annotations$Sex <- factor(annotations$Sex, levels = c("0", "1"))
annotations <- annotations[order(annotations$`TRS-defined risk status`),, drop = FALSE]
heatmap_counts <- heatmap_counts[match(rownames(annotations), rownames(heatmap_counts)),]
annotcolfibrosis <- brewer.pal(5, "Greens")
names(annotcolfibrosis) <- c("0", "1", "2", "3", "4")
annotcolsteatosis <- brewer.pal(4, "Oranges")
names(annotcolsteatosis) <- c("0","1","2","3")
annotcolLI <- brewer.pal(4, "Purples")
names(annotcolLI) <- c("0","1","2","3")
annotcolHB <- brewer.pal(3, "Greys")
names(annotcolHB) <- c("0","1","2")
annotcolrisk <- c("#CC0000", "#3399FF")
names(annotcolrisk) <- c("high", "low")
annotcolsex <- c("chartreuse3", "darkorange2")
names(annotcolsex) <- c("0", "1")
ann_colors <- list('TRS-defined risk status' = annotcolrisk, 'NAS ballooning' = annotcolHB, 'NAS steatosis' = annotcolLI, 'NAS inflammation' = annotcolsteatosis, 'NASH-CRN fibrosis' = annotcolfibrosis, 'Sex' = annotcolsex)
heatmap_counts <- t(heatmap_counts)
pheatmap(heatmap_counts, annotation_col = annotations, annotation_colors = ann_colors, cluster_rows = TRUE, scale = "row", show_colnames  = FALSE, cluster_cols = FALSE)

#Waterfall plot
risk_scores <- risk_scores %>% distinct(sums, .keep_all = TRUE) %>% mutate(rank = rank(sums, ties.method = "random")) %>% arrange(rank)
risk_scores$status <- as.character(risk_scores$status)
ggplot(risk_scores, aes(x = rank, y = sums)) + 
  geom_point(aes(color = risk)) + 
  geom_hline(yintercept = median(risk_scores$sums), linetype = "dashed", color = "black") + 
  geom_vline(xintercept = median(risk_scores$rank), linetype = "dashed", color = "black") + 
  xlab("Patients (increasing risk score)") + ylab("Risk scores") + theme_bw()

#Correlation plot
ggplot(risk_scores, aes(x = rank, y = time)) + 
  geom_point(aes(color = status)) + 
  geom_vline(xintercept = median(risk_scores$rank), linetype = "dashed", color = "black") + 
  xlab("Patients (increasing risk score)") + ylab("Survival time (days)") + scale_color_discrete(labels= c("No events", "Events")) + 
  labs(color = "Decompensation") + theme_bw()
