library(clusterProfiler)

pathway_analysis <- function(results, go, reactome, kegg, eps = 1e-10, permut = 1000, lv = TRUE){
  if (lv == FALSE){
    genes <- results$results[which(results$results$FDR < 0.05),]
  }else{
    genes <- results$results[which(results$results$adj.P.Val < 0.05),]
  }
  genes <- genes[which(!duplicated(genes$genes)),]
  genes <- genes[which(!genes$genes == "<NA>"),]
  geneList <- genes$logFC
  names(geneList) <- genes$genes
  geneList <- sort(geneList, decreasing = TRUE)
  gsea_go <- GSEA(geneList, TERM2GENE = go, eps = eps, nPermSimple = permut)
  go_results <- gsea_go@result
  gsea_kegg <- GSEA(geneList, TERM2GENE = kegg, eps = eps, nPermSimple = permut)
  kegg_results <- gsea_kegg@result
  gsea_reactome <- GSEA(geneList, TERM2GENE = reactome, eps = eps, nPermSimple = permut)
  reactome_results <- gsea_reactome@result
  list = list(GO = go_results, KEGG = kegg_results, Reactome = reactome_results, GSEA_GO = gsea_go, GSEA_KEGG = gsea_kegg, GSEA_reactome = gsea_reactome)
  return(list)
}
