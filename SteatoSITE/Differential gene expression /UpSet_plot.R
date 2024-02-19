library(UpSetR)

load("./DGE_and_GSEA.R")

#UpSet plot for the stages vs control
##Match genes so the order is the same in all the stages of the disease.
F3[[1]]$results <- F3[[1]]$results[match(F4[[1]]$results$genes, F3[[1]]$results$genes),]
F2[[1]]$results <- F2[[1]]$results[match(F4[[1]]$results$genes, F2[[1]]$results$genes),]
F0F1[[1]]$results <- F0F1[[1]]$results[match(F4[[1]]$results$genes, F0F1[[1]]$results$genes),]
Simple_steatosis[[1]]$results <- Simple_steatosis[[1]]$results[match(F4[[1]]$results$genes, Simple_steatosis[[1]]$results$genes),]

##Create a data frame that contains the genes as rows and the stages as columns. logFC and FDR for each stage are added in this data frame.
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

##Create plot. In this case, I am only interested in the in the intersection between all the stages together, and each stage separately. This can be modified to user's preferences.
upset(fromList(list_genes), intersections = list(list("Simple_steatosis", "F0_F1", "F2", "F3", "F4"), list("Simple_steatosis"), list("F0_F1"), list("F2"), list("F3"), list("F4")), keep.order = T)

