library(clusterProfiler)
library(enrichplot)

Pathway_analysis <- function(dir, type, fc_threshold){
  file <- read.csv(dir, row.names = 1)
  file <- subset(file, (file$FDR < 0.05))
  file <- subset(file, (file$logFC > fc_threshold|file$logFC < fc_threshold))
  file <- subset(file, file$entrez!= "<NA>")
  file <- file[which(duplicated(file$entrez) == FALSE), ]
  geneList <- file$logFC
  names(geneList) <- as.character(file$entrez)
  geneList <- sort(geneList, decreasing =TRUE)
  
  reactome <- read.gmt("./c2.cp.reactome.v7.2.entrez.gmt")  #gmt file that contains the reactome sets (downloaded from Braod Institute)
  reactome_enrich <- enricher(names(geneList), TERM2GENE= reactome)
  reactome_enrich_results <- reactome_enrich@result
  if (reactome_enrich@result$p.adjust[1] >= 0.05){
    print("No enrich reactome pathway with p.adjust value < 0.05!")
  } else {
    write.csv(reactome_enrich_results, paste("./Pathways/Reactome", type, "enrich.csv", sep = "_"))
    pdf(paste("./Pathways/Barplot_reactome", type, "enrich.pdf", sep = "_"), width = 15, height = 12)
    p <- barplot(reactome_enrich)
    print(p)
    dev.off()
    ##The following part cannot be done for now as org.Hs.eg.db has some issues and cannot be used with R at the moment.
    #pdf(paste("./Pathways/Network_reactome", type, "enrich.pdf"), width = 30, height = 20)
    #network <- setReadable(reactome_enrich, org.Hs.eg.db, "ENTREZID")  #This function transforms the entrez id into symbols for the networks. 
    #p <- cnetplot(network, foldChange = geneList, showCategory = 4)  #It creates a network with the top 4 pathways obtained and the genes that have been found from that pathway in the DEG list. Those genes will have different colour depending on the fold-change. 
    #print(p)
    #dev.off()
  }
  
  go <- read.gmt("./c5.go.v7.3.entrez.gmt")  #gmt file that contains GO terms (CC, BP and MF).
  go_enrich <- enricher(names(geneList), TERM2GENE = go, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  go_enrich_results <- go_enrich@result
  if (go_enrich@result$p.adjust[1] >= 0.05){
    print("No enrich GO term with p.adjust value < 0.05!")
  } else {
    write.csv(go_enrich_results, paste("./Pathways/GO", type, "enrich.csv", sep = "_"))
    pdf(paste("./Pathways/Dotplot", type, "go_enrich.pdf", sep = "_"), width = 15, height = 14)
    p <- dotplot(go_enrich, font.size = 8)  #It creates a dotplot with the top 10 GO terms linked with the DEG list. Each dot will have a different colour, depending on their p.adjust value, and a different size, depending on the number of counts on each GO term.
    print(p)
    dev.off()
  }
  
  kegg <- read.gmt("./c2.cp.kegg.v7.2.entrez.gmt")  #gmt file that contains the KEGG pathways (downloaded from Broad Institute)
  kegg_enrich <- enricher(names(geneList), TERM2GENE = kegg, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  kegg_enrich_results <- kegg_enrich@result
  if (kegg_enrich@result$p.adjust[1] >= 0.05){
    print("No enrich KEGG pathway with p.adjust value < 0.05!")
  } else {
    write.csv(kegg_enrich_results, paste("./Pathways/Kegg", type, "enrich.csv", sep = "_"))
    pdf(paste("./Pathways/Dotplot", type, "kegg_enrich.pdf", sep = "_"), width = 15, height = 12)
    p <- dotplot(kegg_enrich)
    print(p)
    dev.off()
  }
}
