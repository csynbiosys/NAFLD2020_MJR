library(UpSetR)

load("./DGE_and_GSEA.R")

#UpSet plot for the stages vs control
##Match genes so the order is the same in all the stages of the disease.
F3[[1]]$results <- F3[[1]]$results[match(F4[[1]]$results$genes, F3[[1]]$results$genes),]
F2[[1]]$results <- F2[[1]]$results[match(F4[[1]]$results$genes, F2[[1]]$results$genes),]
F0F1[[1]]$results <- F0F1[[1]]$results[match(F4[[1]]$results$genes, F0F1[[1]]$results$genes),]
Simple_steatosis[[1]]$results <- Simple_steatosis[[1]]$results[match(F4[[1]]$results$genes, Simple_steatosis[[1]]$results$genes),]

##Create a data frame that contains the genes as rows and the stages as columns.
de <- as.data.frame(F4[[1]]$results$logFC)
colnames(de) <- "F4"
de$F3 <- F3[[1]]$results$logFC
de$F2 <- F2[[1]]$results$logFC
de$`F0/F1` <- F0F1[[1]]$results$logFC
de$`Simple steatosis` <- Simple_steatosis[[1]]$results$logFC
de$F4_fdr <- F4[[1]]$results$adj.P.Val
de$F3_fdr <- F3[[1]]$results$adj.P.Val
de$F2_fdr <- F2[[1]]$results$adj.P.Val
de$`F0/F1_fdr` <- F0F1[[1]]$results$adj.P.Val
de$`Simple steatosis_fdr` <- Simple_steatosis[[1]]$results$adj.P.Val
rownames(de) <- F4$genes

##Make a list to use as input in the UpSet plot.
list_genes <- lapply(c(1:5), function(i){
  rownames(de)[(abs(de[[i]]) >= 1 & de[[i+5]] < 0.05)]
})
names(list_genes) <- c("F4", "F3", "F2", "F0_F1", "Simple_steatosis")
list_genes <- rev(list_genes)

##Create plot
upset(fromList(list_genes), intersections = list(list("Simple_steatosis", "F0_F1", "F2", "F3", "F4"), list("Simple_steatosis"), list("F0_F1"), list("F2"), list("F3"), list("F4")), keep.order = T)

#UpSet plot for the fibrotic stages vs simple steatosis
F3vsss <- F3vsss[match(F4vsss$genes, F3vsss$genes),]
F2vsss <- F2vsss[match(F4vsss$genes, F2vsss$genes),]
F0F1vsss <- F0F1vsss[match(F4vsss$genes, F0F1vsss$genes),]

##Create a data frame that contains the genes as rows and the stages as columns.
de <- as.data.frame(F4vsss$logFC)
colnames(de) <- "F4"
de$F3 <- F3vsss$logFC
de$F2 <- F2vsss$logFC
de$`F0/F1` <- F0F1vsss$logFC
de$F4_fdr <- F4vsss$FDR
de$F3_fdr <- F3vsss$FDR
de$F2_fdr <- F2vsss$FDR
de$`F0/F1_fdr` <- F0F1vsss$FDR
rownames(de) <- F4$genes

##Make a list to use as input in the UpSet plot.
list_genes <- lapply(c(1:4), function(i){
  rownames(de)[(abs(de[[i]]) >= 1 & de[[i+4]] < 0.05)]
})
names(list_genes) <- c("F4", "F3", "F2", "F0_F1")
list_genes <- rev(list_genes)

##Create plot
upset(fromList(list_genes), intersections = list(list("F0_F1", "F2", "F3", "F4"), list("F0_F1"), list("F2"), list("F3"), list("F4")), keep.order = T)