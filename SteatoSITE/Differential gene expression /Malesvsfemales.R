##Script to get DGE for each fibrotic stage in men and women

library(clusterProfiler)
library(org.Hs.eg.db)
library(VennDiagram)
library(pheatmap)

#DGList prot corresponds to the DGList created with the metadata available

#Split dgList into man and women
dgList_women <- dgList[, dgList$samples$sex == 0]
dgList_man <- dgList[, dgList$samples$sex == 1]

#Perform DGE analysis for both men and women and each fib group
fc = 1
fdr = 0.05
f0f1_man <- dea_limma(dgList_man, fc, fdr, "NAFLD", contrast_name = "NAFLDF0F1 - NAFLDControl")
f2_man <- dea_limma(dgList_man, fc, fdr, "NAFLD", contrast_name = "NAFLDF2 - NAFLDControl")
f3_man <- dea_limma(dgList_man, fc, fdr, "NAFLD", contrast_name = "NAFLDF3 - NAFLDControl")
f4_man <- dea_limma(dgList_man, fc, fdr, "NAFLD", contrast_name = "NAFLDF4 - NAFLDControl")

f0f1_woman <- dea_limma(dgList_woman, fc, fdr, "NAFLD", contrast_name = "NAFLDF0F1 - NAFLDControl")
f2_woman <- dea_limma(dgList_woman, fc, fdr, "NAFLD", contrast_name = "NAFLDF2 - NAFLDControl")
f3_woman <- dea_limma(dgList_woman, fc, fdr, "NAFLD", contrast_name = "NAFLDF3 - NAFLDControl")
f4_woman <- dea_limma(dgList_woman, fc, fdr, "NAFLD", contrast_name = "NAFLDF4 - NAFLDControl")

#Get pathway from each group.
kegg <- read.gmt("../Data/c2.cp.kegg.v7.4.symbols.gmt") #Upload gmt files that contain the different terms and pathways to be used for analysis. This is continuously updated in gsea-msigdb.
reactome <- read.gmt("../Data/c2.cp.reactome.v7.4.symbols.gmt")
go <- read.gmt("../Data/c5.go.v7.4.symbols.gmt")

f0f1_man_path <- pathway_analysis(f0f1_man, go, reactome, kegg, lv = TRUE, permut = 10000, eps = 0)
f2_man_path <- pathway_analysis(f2_man, go, reactome, kegg, lv = TRUE, permut = 10000, eps = 0)
f3_man_path <- pathway_analysis(f3_man, go, reactome, kegg, lv = TRUE, permut = 10000, eps = 0)
f4_man_path <- pathway_analysis(f4_man, go, reactome, kegg, lv = TRUE, permut = 10000, eps = 0)

f0f1_woman_path <- pathway_analysis(f0f1_man, go, reactome, kegg, lv = TRUE, permut = 10000, eps = 0)
f2_woman_path <- pathway_analysis(f2_man, go, reactome, kegg, lv = TRUE, permut = 10000, eps = 0)
f3_woman_path <- pathway_analysis(f3_man, go, reactome, kegg, lv = TRUE, permut = 10000, eps = 0)
f4_woman_path <- pathway_analysis(f4_man, go, reactome, kegg, lv = TRUE, permut = 10000, eps = 0)

#Compare clusters between groups. compareCluster requires entrezid form insteado of symbols, so it needs to be changed too. A list with the genes for both men and women need to be added. GO terms and reactome pathways are retrieved.
f0f1_list <- list(man = bitr(f0f1_woman$results$genes[f0f1_woman$results$adj.P.Val < fdr & abs(f0f1_woman$results$logFC) >= fc], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID, 
                  woman =bitr(f0f1_man$results$genes[f0f1_man$results$adj.P.Val < fdr & abs(f0f1_man$results$logFC) >= fc], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID) 
f0f1_comparison_go <- compareCluster(geneClusters = f0f1_list, fun = "enrichGO", OrgDb = org.Hs.eg.db)
f0f1_comparison_reactome <- compareCluster(geneClusters = f0f1_list, fun = "enrichPathway")

f2_list <- list(man = bitr(f2_woman$results$genes[f0f1_woman$results$adj.P.Val < fdr & abs(f2_woman$results$logFC) >= fc], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID, 
                  woman =bitr(f2_man$results$genes[f0f1_man$results$adj.P.Val < fdr & abs(f2_man$results$logFC) >= fc], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID) 
f2_comparison_go <- compareCluster(geneClusters = f2_list, fun = "enrichGO", OrgDb = org.Hs.eg.db)
f2_comparison_reactome <- compareCluster(geneClusters = f2_list, fun = "enrichPathway")

f3_list <- list(man = bitr(f3_woman$results$genes[f3_woman$results$adj.P.Val < fdr & abs(f3_woman$results$logFC) >= fc], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID, 
                  woman =bitr(f3_man$results$genes[f2_man$results$adj.P.Val < fdr & abs(f3_man$results$logFC) >= fc], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID) 
f3_comparison_go <- compareCluster(geneClusters = f3_list, fun = "enrichGO", OrgDb = org.Hs.eg.db)
f3_comparison_reactome <- compareCluster(geneClusters = f3_list, fun = "enrichPathway")

f4_list <- list(man = bitr(f4_woman$results$genes[f0f1_woman$results$adj.P.Val < fdr & abs(f4_woman$results$logFC) >= fc], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID, 
                  woman =bitr(f4_man$results$genes[f0f1_man$results$adj.P.Val < fdr & abs(f4_man$results$logFC) >= fc], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID) 
f4_comparison_go <- compareCluster(geneClusters = f4_list, fun = "enrichGO", OrgDb = org.Hs.eg.db)
f4_comparison_reactome <- compareCluster(geneClusters = f4_list, fun = "enrichPathway")

#Common genes between women and men per stage in the first 100 genes.
list_f0f1 <- list(men = f0f1_man$results$genes[1:100], women = f0f1_woman$results$genes[1:100])
venn.diagram(list_f0f1, filename = "./f0f1.png", imagetype = "png")

list_f2 <- list(men = f2_man$results$genes[1:100], women = f2_woman$results$genes[1:100])
venn.diagram(list_f2, filename = "./f2.png", imagetype = "png")

list_f3 <- list(men = f3_man$results$genes[1:100], women = f3_woman$results$genes[1:100])
venn.diagram(list_f3, filename = "./f3.png", imagetype = "png")

list_f4 <- list(men = f4_man$results$genes[1:100], women = f4_woman$results$genes[1:100])
venn.diagram(list_f4, filename = "./f4.png", imagetype = "png")

all_top_genes <- c(list_f0f1$men, list_f0f1$women, list_f2$men, list_f2$women, list_f3$men, list_f4$men, list_f4$women)
all_top_genes <- unique(all_top_genes)

#Pheatmap with the top genes. A data frame the top 20 genes (ranked according to |logFC|) that have adj. pval < fdr were selected for each stage for men and women. Gene names, logFC and adj. pval are the chosen columns.

fold_change_manvswoman <- f0f1_man$results[which(f0f1_man$results$adj.P.Val < fdr)[1:20], c(1,2,6)]
fold_change_manvswoman <- merge(fold_change_manvswoman, f0f1_woman$results[which(f0f1_woman$results$adj.P.Val < fdr)[1:20], c(1,2,6)], all = TRUE, by = "genes") #Merge men and women, as well as different stages to the data frame 
fold_change_manvswoman <- fold_change_manvswoman %>%
  left_join(f0f1_man$results[which(f0f1_man$results$adj.P.Val < fdr), c(1,2,6)], by = "genes") %>% mutate(logFC_man_f0f1 = coalesce(logFC.x, logFC))
fold_change_manvswoman <- fold_change_manvswoman %>% mutate(adj.P.Val_man_f0f1 = coalesce(adj.P.Val.x, adj.P.Val))
fold_change_manvswoman <- fold_change_manvswoman %>%
  left_join(f0f1_woman$results[which(f0f1_woman$results$adj.P.Val < fdr), c(1,2,6)], by = "genes") %>% mutate(logFC_woman_f0f1 = coalesce(logFC.y, logFC.y.y))
fold_change_manvswoman <- fold_change_manvswoman %>% mutate(adj.P.Val_woman_f0f1 = coalesce(adj.P.Val.y, adj.P.Val.y.y))
fold_change_manvswoman <- fold_change_manvswoman[, -c(2:7, 10:11)]

fold_change_manvswoman <- merge(fold_change_manvswoman, f2_man$results[which(f2_man$results$adj.P.Val < fdr)[1:20], c(1,2,6)], all = TRUE, by = "genes")
fold_change_manvswoman <- merge(fold_change_manvswoman, f2_woman$results[which(f2_woman$results$adj.P.Val < fdr)[1:20], c(1,2,6)], all = TRUE, by = "genes")
fold_change_manvswoman <- fold_change_manvswoman %>%
  left_join(f2_man$results[which(f2_man$results$adj.P.Val < fdr), c(1,2,6)], by = "genes") %>% mutate(logFC_man_f2 = coalesce(logFC.x, logFC))
fold_change_manvswoman <- fold_change_manvswoman %>% mutate(adj.P.Val_man_f2 = coalesce(adj.P.Val.x, adj.P.Val))
fold_change_manvswoman <- fold_change_manvswoman %>%
  left_join(f2_woman$results[which(f2_woman$results$adj.P.Val < fdr), c(1,2,6)], by = "genes") %>% mutate(logFC_woman_f2 = coalesce(logFC.y, logFC.y.y))
fold_change_manvswoman <- fold_change_manvswoman %>% mutate(adj.P.Val_woman_f2 = coalesce(adj.P.Val.y, adj.P.Val.y.y))
fold_change_manvswoman <- fold_change_manvswoman[, -c(6:11, 14:15)]

fold_change_manvswoman <- merge(fold_change_manvswoman, f3_man$results[which(f3_man$results$adj.P.Val < fdr)[1:20], c(1,2,6)], all = TRUE, by = "genes")
fold_change_manvswoman <- merge(fold_change_manvswoman, f3_woman$results[which(f3_woman$results$adj.P.Val < fdr)[1:20], c(1,2,6)], all = TRUE, by = "genes")
fold_change_manvswoman <- fold_change_manvswoman %>%
  left_join(f3_man$results[which(f3_man$results$adj.P.Val < fdr), c(1,2,6)], by = "genes") %>% mutate(logFC_man_f3 = coalesce(logFC.x, logFC))
fold_change_manvswoman <- fold_change_manvswoman %>% mutate(adj.P.Val_man_f3 = coalesce(adj.P.Val.x, adj.P.Val))
fold_change_manvswoman <- fold_change_manvswoman %>%
  left_join(f3_woman$results[which(f3_woman$results$adj.P.Val < fdr), c(1,2,6)], by = "genes") %>% mutate(logFC_woman_f3 = coalesce(logFC.y, logFC.y.y))
fold_change_manvswoman <- fold_change_manvswoman %>% mutate(adj.P.Val_woman_f3 = coalesce(adj.P.Val.y, adj.P.Val.y.y))
fold_change_manvswoman <- fold_change_manvswoman[, -c(10:15, 18:19)]

fold_change_manvswoman <- merge(fold_change_manvswoman, f4_man$results[which(f4_man$results$adj.P.Val < fdr)[1:20], c(1,2,6)], all = TRUE, by = "genes")
fold_change_manvswoman <- merge(fold_change_manvswoman, f4_woman$results[which(f4_woman$results$adj.P.Val < fdr)[1:20], c(1,2,6)], all = TRUE, by = "genes")
fold_change_manvswoman <- fold_change_manvswoman %>%
  left_join(f4_man$results[which(f4_man$results$adj.P.Val < fdr), c(1,2,6)], by = "genes") %>% mutate(logFC_man_f4 = coalesce(logFC.x, logFC))
fold_change_manvswoman <- fold_change_manvswoman %>% mutate(adj.P.Val_man_f4 = coalesce(adj.P.Val.x, adj.P.Val))
fold_change_manvswoman <- fold_change_manvswoman %>%
  left_join(f4_woman$results[which(f4_woman$results$adj.P.Val < fdr), c(1,2,6)], by = "genes") %>% mutate(logFC_woman_f4 = coalesce(logFC.y, logFC.y.y))
fold_change_manvswoman <- fold_change_manvswoman %>% mutate(adj.P.Val_woman_f4 = coalesce(adj.P.Val.y, adj.P.Val.y.y))
fold_change_manvswoman <- fold_change_manvswoman[, -c(14:19, 22:23)]
rownames(fold_change_manvswoman) <- fold_change_manvswoman$genes
fold_change_manvswoman <- fold_change_manvswoman[, c(1,2,6,10,14,4,8,12,16,3,5,7,9,11,13, 15,17)] #Reorder columns for illustration purposes so men are first and then women. 
colnames(fold_change_manvswoman) <- c("Genes", "Man F0/F1", "Man F2", "Man F3", "Man F4", "Woman F0/F1", "Woman F2", "Woman F3", "Woman F4", "pval_man_f0f1", "pval_woman_f0f1", "pval_man_f2", "pval_woman_f2", "pval_man_f3", "pval_man_f3", "pval_man_f4", "pval_woman_f4")
pheatmap(fold_change_manvswoman[, c(2:9)], na_col = "lightgrey", cluster_rows = FALSE, cluster_cols = FALSE, border_color = "NA", angle_col = 315)
