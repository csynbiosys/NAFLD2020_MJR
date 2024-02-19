library(clusterProfiler)

#This function performs GSEA analysis to obtain GO terms, reactome and KEGG pathways.
#Results correspond to object obtained with eirht DGE_edgeR or DGE_limma_voom
#go, reactome and kegg corresponds to the gmt objects processed with read.gmt
#eps sets boundaries for calculating p-values
#permut is the number of permutations the GSEA function will perform for each analysis.
#lv can be either TRUE or FALSE. If TRUE, results provided were obtained with DGE_limma_voom function. If not, it was obtained with DGE_edgeR.
pathway_analysis <- function(results, go, reactome, kegg, eps = 1e-10, permut = 1000, lv = TRUE){
  if (lv == FALSE){
    genes <- results$results[which(results$results$FDR < 0.05),]
  }else{
    genes <- results$results[which(results$results$adj.P.Val < 0.05),]
  }
  genes <- genes[which(!duplicated(genes$genes)),] #Ensure duplicated or genes without a symbols name are eliminated
  genes <- genes[which(!genes$genes == "<NA>"),] 
  geneList <- genes$logFC
  names(geneList) <- genes$genes
  geneList <- sort(geneList, decreasing = TRUE) #Sort list according to logFC
  gsea_go <- GSEA(geneList, TERM2GENE = go, eps = eps, nPermSimple = permut)
  go_results <- gsea_go@result
  gsea_kegg <- GSEA(geneList, TERM2GENE = kegg, eps = eps, nPermSimple = permut)
  kegg_results <- gsea_kegg@result
  gsea_reactome <- GSEA(geneList, TERM2GENE = reactome, eps = eps, nPermSimple = permut)
  reactome_results <- gsea_reactome@result
  list = list(GO = go_results, KEGG = kegg_results, Reactome = reactome_results, GSEA_GO = gsea_go, GSEA_KEGG = gsea_kegg, GSEA_reactome = gsea_reactome)
  return(list)
}
