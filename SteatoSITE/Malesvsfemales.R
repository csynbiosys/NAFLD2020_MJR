##Script to get DGE for each fibrotic stage in males and females
library(clusterProfiler)
library(org.Hs.eg.db)
library(VennDiagram)
library(pheatmap)

#Split dgList into male and female
dgList_female <- dgList_prot[, dgList_prot$samples$sex == 0]
dgList_male <- dgList_prot[, dgList_prot$samples$sex == 1]

#Perform DGE analysis for both males and females and each fib group
f0f1_male <- dea_limma(dgList_male, "NAFLD_Tim", contrast_name = "NAFLD_TimF0F1 - NAFLD_TimControl")
f2_male <- dea_limma(dgList_male, "NAFLD_Tim", contrast_name = "NAFLD_TimF2 - NAFLD_TimControl")
f3_male <- dea_limma(dgList_male, "NAFLD_Tim", contrast_name = "NAFLD_TimF3 - NAFLD_TimControl")
f4_male <- dea_limma(dgList_male, "NAFLD_Tim", contrast_name = "NAFLD_TimF4 - NAFLD_TimControl")

f0f1_female <- dea_limma(dgList_female, "NAFLD_Tim", contrast_name = "NAFLD_TimF0F1 - NAFLD_TimControl")
f2_female <- dea_limma(dgList_female, "NAFLD_Tim", contrast_name = "NAFLD_TimF2 - NAFLD_TimControl")
f3_female <- dea_limma(dgList_female, "NAFLD_Tim", contrast_name = "NAFLD_TimF3 - NAFLD_TimControl")
f4_female <- dea_limma(dgList_female, "NAFLD_Tim", contrast_name = "NAFLD_TimF4 - NAFLD_TimControl")

#Get pathway from each group
f0f1_male_path <- pathway_analysis(f0f1_male, go, reactome, kegg, lv = TRUE, permut = 10000, eps = 0)
f2_male_path <- pathway_analysis(f2_male, go, reactome, kegg, lv = TRUE, permut = 10000, eps = 0)
f3_male_path <- pathway_analysis(f3_male, go, reactome, kegg, lv = TRUE, permut = 10000, eps = 0)
f4_male_path <- pathway_analysis(f4_male, go, reactome, kegg, lv = TRUE, permut = 10000, eps = 0)

f0f1_female_path <- pathway_analysis(f0f1_male, go, reactome, kegg, lv = TRUE, permut = 10000, eps = 0)
f2_female_path <- pathway_analysis(f2_male, go, reactome, kegg, lv = TRUE, permut = 10000, eps = 0)
f3_female_path <- pathway_analysis(f3_male, go, reactome, kegg, lv = TRUE, permut = 10000, eps = 0)
f4_female_path <- pathway_analysis(f4_male, go, reactome, kegg, lv = TRUE, permut = 10000, eps = 0)

#Compare clusters between groups
f0f1_list <- list(male = bitr(f0f1_female$results$genes[f0f1_female$results$adj.P.Val < 0.05 & abs(f0f1_female$results$logFC) >= 1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID, 
                  female =bitr(f0f1_male$results$genes[f0f1_male$results$adj.P.Val < 0.05 & abs(f0f1_male$results$logFC) >= 1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID) 
f0f1_comparison <- compareCluster(geneClusters = f0f1_list, fun = "enrichGO", OrgDb = org.Hs.eg.db)
f0f1_comparison_reactome <- compareCluster(geneClusters = f0f1_list, fun = "enrichPathway")

f2_list <- list(male = bitr(f2_female$results$genes[f0f1_female$results$adj.P.Val < 0.05 & abs(f2_female$results$logFC) >= 1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID, 
                  female =bitr(f2_male$results$genes[f0f1_male$results$adj.P.Val < 0.05 & abs(f2_male$results$logFC) >= 1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID) 
f2_comparison <- compareCluster(geneClusters = f2_list, fun = "enrichGO", OrgDb = org.Hs.eg.db)
f2_comparison_reactome <- compareCluster(geneClusters = f2_list, fun = "enrichPathway")

f3_list <- list(male = bitr(f3_female$results$genes[f3_female$results$adj.P.Val < 0.05 & abs(f3_female$results$logFC) >= 1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID, 
                  female =bitr(f3_male$results$genes[f2_male$results$adj.P.Val < 0.05 & abs(f3_male$results$logFC) >= 1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID) 
f3_comparison <- compareCluster(geneClusters = f3_list, fun = "enrichGO", OrgDb = org.Hs.eg.db)
f3_comparison_reactome <- compareCluster(geneClusters = f3_list, fun = "enrichPathway")

f4_list <- list(male = bitr(f4_female$results$genes[f0f1_female$results$adj.P.Val < 0.05 & abs(f4_female$results$logFC) >= 1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID, 
                  female =bitr(f4_male$results$genes[f0f1_male$results$adj.P.Val < 0.05 & abs(f4_male$results$logFC) >= 1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID) 
f4_comparison <- compareCluster(geneClusters = f4_list, fun = "enrichGO", OrgDb = org.Hs.eg.db)
f4_comparison_reactome <- compareCluster(geneClusters = f4_list, fun = "enrichPathway")

#Common genes between females and males per stage
list_f0f1 <- list(males = f0f1_male$results$genes[1:100], females = f0f1_female$results$genes[1:100])
venn.diagram(list_f0f1, filename = "/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results_project/Without non-prot coding genes/Sex differences/f0f1.png", imagetype = "png")

list_f2 <- list(males = f2_male$results$genes[1:100], females = f2_female$results$genes[1:100])
venn.diagram(list_f2, filename = "/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results_project/Without non-prot coding genes/Sex differences/f2.png", imagetype = "png")

list_f3 <- list(males = f3_male$results$genes[1:100], females = f3_female$results$genes[1:100])
venn.diagram(list_f3, filename = "/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results_project/Without non-prot coding genes/Sex differences/f3.png", imagetype = "png")

list_f4 <- list(males = f4_male$results$genes[1:100], females = f4_female$results$genes[1:100])
venn.diagram(list_f4, filename = "/Users/mariajimenezramos/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Third year/SteatoSITE/Results_project/Without non-prot coding genes/Sex differences/f4.png", imagetype = "png")

all_top_genes <- c(list_f0f1$males, list_f0f1$females, list_f2$males, list_f2$females, list_f3$males, list_f4$males, list_f4$females)
all_top_genes <- unique(all_top_genes)

#Pheatmap with the top genes 
counts_sex_diff <- logcpm_norm_prot[rownames(logcpm_norm_prot)%in%all_top_genes,]
dgList_prot$samples$sex <- factor(dgList_prot$samples$sex, levels = c("0", "1"))
dgList_prot$samples <- dgList_prot$samples[order(dgList_prot$samples$sex),]
counts_sex_diff <- counts_sex_diff[, match(dgList_prot$samples$sample_id, colnames(counts_sex_diff))]

pheatmap(counts_sex_diff, annotation = dgList_prot$samples[,c(6,7,8,10,12,33)], cluster_cols = FALSE,  show_rownames = FALSE, show_colnames = FALSE)

fold_change_malevsfemale <- f0f1_male$results[which(f0f1_male$results$adj.P.Val < 0.05)[1:20], c(1,2,6)]
fold_change_malevsfemale <- merge(fold_change_malevsfemale, f0f1_female$results[which(f0f1_female$results$adj.P.Val < 0.05)[1:20], c(1,2,6)], all = TRUE, by = "genes")
fold_change_malevsfemale <- fold_change_malevsfemale %>%
  left_join(f0f1_male$results[which(f0f1_male$results$adj.P.Val < 0.05), c(1,2,6)], by = "genes") %>% mutate(logFC_male_f0f1 = coalesce(logFC.x, logFC))
fold_change_malevsfemale <- fold_change_malevsfemale %>% mutate(adj.P.Val_male_f0f1 = coalesce(adj.P.Val.x, adj.P.Val))
fold_change_malevsfemale <- fold_change_malevsfemale %>%
  left_join(f0f1_female$results[which(f0f1_female$results$adj.P.Val < 0.05), c(1,2,6)], by = "genes") %>% mutate(logFC_female_f0f1 = coalesce(logFC.y, logFC.y.y))
fold_change_malevsfemale <- fold_change_malevsfemale %>% mutate(adj.P.Val_female_f0f1 = coalesce(adj.P.Val.y, adj.P.Val.y.y))
fold_change_malevsfemale <- fold_change_malevsfemale[, -c(2:7, 10:11)]

fold_change_malevsfemale <- merge(fold_change_malevsfemale, f2_male$results[which(f2_male$results$adj.P.Val < 0.05)[1:20], c(1,2,6)], all = TRUE, by = "genes")
fold_change_malevsfemale <- merge(fold_change_malevsfemale, f2_female$results[which(f2_female$results$adj.P.Val < 0.05)[1:20], c(1,2,6)], all = TRUE, by = "genes")
fold_change_malevsfemale <- fold_change_malevsfemale %>%
  left_join(f2_male$results[which(f2_male$results$adj.P.Val < 0.05), c(1,2,6)], by = "genes") %>% mutate(logFC_male_f2 = coalesce(logFC.x, logFC))
fold_change_malevsfemale <- fold_change_malevsfemale %>% mutate(adj.P.Val_male_f2 = coalesce(adj.P.Val.x, adj.P.Val))
fold_change_malevsfemale <- fold_change_malevsfemale %>%
  left_join(f2_female$results[which(f2_female$results$adj.P.Val < 0.05), c(1,2,6)], by = "genes") %>% mutate(logFC_female_f2 = coalesce(logFC.y, logFC.y.y))
fold_change_malevsfemale <- fold_change_malevsfemale %>% mutate(adj.P.Val_female_f2 = coalesce(adj.P.Val.y, adj.P.Val.y.y))
fold_change_malevsfemale <- fold_change_malevsfemale[, -c(6:11, 14:15)]

fold_change_malevsfemale <- merge(fold_change_malevsfemale, f3_male$results[which(f3_male$results$adj.P.Val < 0.05)[1:20], c(1,2,6)], all = TRUE, by = "genes")
fold_change_malevsfemale <- merge(fold_change_malevsfemale, f3_female$results[which(f3_female$results$adj.P.Val < 0.05)[1:20], c(1,2,6)], all = TRUE, by = "genes")
fold_change_malevsfemale <- fold_change_malevsfemale %>%
  left_join(f3_male$results[which(f3_male$results$adj.P.Val < 0.05), c(1,2,6)], by = "genes") %>% mutate(logFC_male_f3 = coalesce(logFC.x, logFC))
fold_change_malevsfemale <- fold_change_malevsfemale %>% mutate(adj.P.Val_male_f3 = coalesce(adj.P.Val.x, adj.P.Val))
fold_change_malevsfemale <- fold_change_malevsfemale %>%
  left_join(f3_female$results[which(f3_female$results$adj.P.Val < 0.05), c(1,2,6)], by = "genes") %>% mutate(logFC_female_f3 = coalesce(logFC.y, logFC.y.y))
fold_change_malevsfemale <- fold_change_malevsfemale %>% mutate(adj.P.Val_female_f3 = coalesce(adj.P.Val.y, adj.P.Val.y.y))
fold_change_malevsfemale <- fold_change_malevsfemale[, -c(10:15, 18:19)]

fold_change_malevsfemale <- merge(fold_change_malevsfemale, f4_male$results[which(f4_male$results$adj.P.Val < 0.05)[1:20], c(1,2,6)], all = TRUE, by = "genes")
fold_change_malevsfemale <- merge(fold_change_malevsfemale, f4_female$results[which(f4_female$results$adj.P.Val < 0.05)[1:20], c(1,2,6)], all = TRUE, by = "genes")
fold_change_malevsfemale <- fold_change_malevsfemale %>%
  left_join(f4_male$results[which(f4_male$results$adj.P.Val < 0.05), c(1,2,6)], by = "genes") %>% mutate(logFC_male_f4 = coalesce(logFC.x, logFC))
fold_change_malevsfemale <- fold_change_malevsfemale %>% mutate(adj.P.Val_male_f4 = coalesce(adj.P.Val.x, adj.P.Val))
fold_change_malevsfemale <- fold_change_malevsfemale %>%
  left_join(f4_female$results[which(f4_female$results$adj.P.Val < 0.05), c(1,2,6)], by = "genes") %>% mutate(logFC_female_f4 = coalesce(logFC.y, logFC.y.y))
fold_change_malevsfemale <- fold_change_malevsfemale %>% mutate(adj.P.Val_female_f4 = coalesce(adj.P.Val.y, adj.P.Val.y.y))
fold_change_malevsfemale <- fold_change_malevsfemale[, -c(14:19, 22:23)]
rownames(fold_change_malevsfemale) <- fold_change_malevsfemale$genes
fold_change_malevsfemale <- fold_change_malevsfemale[, c(1,2,6,10,14,4,8,12,16,3,5,7,9,11,13, 15,17)]
colnames(fold_change_malevsfemale) <- c("Genes", "Male F0/F1", "Male F2", "Male F3", "Male F4", "Female F0/F1", "Female F2", "Female F3", "Female F4", "pval_male_f0f1", "pval_female_f0f1", "pval_male_f2", "pval_female_f2", "pval_male_f3", "pval_male_f3", "pval_male_f4", "pval_female_f4")
pheatmap(fold_change_malevsfemale[, c(2:9)], na_col = "lightgrey", cluster_rows = FALSE, cluster_cols = FALSE, border_color = "NA", angle_col = 315)
