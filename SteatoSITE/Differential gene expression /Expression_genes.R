library(ggplot2)
library(tidyr)
library(ggcorrplot)

#Get pheatmap with the DGE for each stage vs control for the hub genes
mcode_expression <- read.csv("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results_project/Without non-prot coding genes/Hub genes/DGE_hub_genes.csv", row.names = 1)
pheatmap(mcode_expression[, c(1,3,5,7,9)], cluster_cols = FALSE, angle = 45)

#Get the expression of the hub genes
mcode_normalised_expression <- cpm(dgList, log = TRUE , prior.count = 1)
mcode_normalised_expression <- as.data.frame(mcode_normalised_expression[which(rownames(mcode_normalised_expression)%in%mcode_genes$display.name),])
mcode_normalised_expression$genes <- rownames(mcode_normalised_expression)

#Modify the document so we can group each patient according to the disease stage and plot boxplots for log expression
mcode_normalised_longer <- tidyr::pivot_longer(mcode_normalised_expression, 1:ncol(mcode_normalised_expression)-1, names_to = "sample_id", values_to = "counts")
mcode_normalised_longer <- mcode_normalised_longer %>% mutate(category = ifelse(sample_id%in%simple_steatosis, "Simple steatosis", "Control"))
mcode_normalised_longer <- mcode_normalised_longer %>% mutate(category = case_when(sample_id%in%f0f1 == TRUE ~ "F0/F1", TRUE~ category))
mcode_normalised_longer <- mcode_normalised_longer %>% mutate(category = case_when(sample_id%in%f2 == TRUE ~ "F2", TRUE~ category))
mcode_normalised_longer <- mcode_normalised_longer %>% mutate(category = case_when(sample_id%in%f3 == TRUE ~ "F3", TRUE~ category))
mcode_normalised_longer <- mcode_normalised_longer %>% mutate(category = case_when(sample_id%in%f4 == TRUE ~ "F4", TRUE~ category))
mcode_normalised_longer$category <- factor(mcode_normalised_longer$category, levels =c("Control", "Simple steatosis", "F0/F1", "F2", "F3", "F4"))

for (i in unique(mcode_genes$display.name)){
  ggplot(mcode_normalised_longer[which(mcode_normalised_longer$genes == i),], aes(x = category, y = counts)) +
    geom_boxplot() + 
    xlab("") + 
    ylab("Log counts") +
    theme_bw()
  ggsave(paste("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results_project/Without non-prot coding genes/Hub genes/", i, ".pdf", sep = ""), width = 7, height = 5, units = "in")
  ggsave(paste("/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results_project/Without non-prot coding genes/Hub genes/", i, ".svg", sep = ""), width = 7, height = 5, units = "in")
  
}

#Perform correlation of the hub genes
mcode_correlation <- cor(t(mcode_normalised_expression), method = "speaman")
mcode_cor_pmat <- cor_pmat(t(mcode_normalised_expression))
ggcorrplot(mcode_correlation, type = 'lower', lab = TRUE, p.mat = mcode_cor_pmat, insig = 'blank', method = "circle", lab_size = 3, outline.color = "white")
