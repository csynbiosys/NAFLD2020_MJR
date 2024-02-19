library(ggplot2)

load("./DGE_and_GSEA.R")

#Create a data frame that contains the different pathways for each stage of the disease. The pathway objects correspond to those obtained in DGE_and_GSEA.R with pathway_analysis funtion. At the end the data frame will contain all the information retrieved from performing GSEA and and column called "type" that will indicate the stage of the disease.
pathways <- F4_path_lv$KEGG
pathways$type <- "F4"
F0F1_path_lv$KEGG$type <- "F0/F1"
F2_path_lv$KEGG$type <- "F2"
F3_path_lv$KEGG$type <- "F3"
simple_steatosis_path_lv$KEGG$type <- "Isolated steatosis"
pathways <- rbind(pathways, F0F1_path_lv$KEGG)
pathways <- rbind(pathways, F2_path_lv$KEGG)
pathways <- rbind(pathways, F3_path_lv$KEGG)
pathways <- rbind(pathways, simple_steatosis_path_lv$KEGG)

#Clean the KEGG terms a bit so it looks nicer in the plot.
pathways$ID <- gsub("_", " ", pathways$ID)
pathways$ID <- gsub("KEGG", "", pathways$ID)

#Factor the different stages of the disease according to the order you want to see them in the plot. Also, an additional column is created to differentiate upregulated and downregulated pathways.
pathways$type <- factor(pathways$type, levels = c("Isolated steatosis", "F0/F1", "F2", "F3", "F4"))
pathways <- pathways %>% mutate(Dysregulation = case_when(enrichmentScore > 0 ~ "Upregulated", TRUE ~ "Downregulated"))

#Create a final data frame that contains the first X -in this case 6- KEGG terms from each stage, and check if they are present in each stage. Then eliminate duplicates pathways.
n_cases = 6
plot_pathways <- data.frame()
for (i in pathways$type){
  plot_pathways <- rbind(plot_pathways, pathways[which(pathways$ID%in%subset(pathways[pathways$type == i & pathways$Dysregulation == "Downregulated",], ID%in%ID[1:n_cases])$ID),])
  plot_pathways <- rbind(plot_pathways, pathways[which(pathways$ID%in%subset(pathways[pathways$type == i & pathways$Dysregulation == "Upregulated",], ID%in%ID[1:n_cases])$ID),])
}

plot_pathways <- unique(plot_pathways) #Remove repeated pathways

#Generate the final plot for either downregulated or upregulated pathways. The x axis is going to indicate the disease stage, the y axis the different KEGG terms, the color of the dots the adjusted p-value and their size the number of genes present in the pathway.
ggplot(plot_pathways[which(plot_pathways$Dysregulation == "Upregulated"),])+
  geom_point(mapping = aes(x= reorder(ID,NES), y=type, color= p.adjust, size = setSize))+
  scale_color_gradient(low = "blue", high = "red") +
  coord_fixed(ratio = 0.5)+
  labs(y = "", x = "", fill = "p.adjust") +
  theme(axis.text=element_text(size=8)) +
  ggtitle('') +
  theme(plot.title = element_text(2, face = "bold", hjust = 1), legend.key.size = unit(2, "line")) +
  coord_flip() + 
  theme_bw()

#Barplots for GO terms. This can be perform with each object obtained in DGE_and_GSEA.R with pathway_analysis function.
simple_steatosis_path_lv$GO$qscore <- - log(simple_steatosis_path_lv$GO$p.adjust, base = 10)
ggplot(F4_path_lv$GO[1:10,]) +
  geom_bar(aes(y = ID, x= qscore, fill = p.adjust), stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(y = "", x = "qscore", fill = "p.adjust") +
  theme_bw()

#Cnetplot for different pathways. In this case, F3 has been used.
fdr = 0.05
genes <- F3[[1]]$results[which(F3[[1]]$results$adj.P.Val < fdr),]
genes <- genes[which(!genes$genes == "<NA>"),] #Remove those genes that have no symbol name.
geneList <- genes$logFC
names(geneList) <- genes$genes
geneList <- sort(geneList, decreasing = TRUE) #Rank genes according to logFC.
cnetplot(F3_path_lv$GSEA_KEGG, categorySize = "pvalue", foldChange =geneList, circular = TRUE, colorEdge = TRUE, showCategory = 3) #Show only first 3 categories.


