source('../SteatoSITE/Pheatmap_function.R')
source('../SteatoSITE/MA_volcano_plot.R')
source('../SteatoSITE/de_analysis.R')
source('../SteatoSITE/Pathway_analysis.R')

#Information about the pathological of each patient.
NAFLD_scores <- read.csv("./NAFLD_scores_pilot.csv", row.names = 1, stringsAsFactors = FALSE)
#We combine the fibrosis stages 1a, 1b and 1c into 1, and transformed the categories Hepatocyte.Ballooning, Steatosis and Lobular.inflammation into characters.
NAFLD_scores <- NAFLD_scores %>% mutate(fibrosis_combined = case_when(Fibrosis == '1a' ~ '1', Fibrosis == '1b' ~ '1', Fibrosis == '1c' ~ '1', TRUE ~ Fibrosis))
NAFLD_scores$Hepatocyte.ballooning <- as.character(NAFLD_scores$Hepatocyte.ballooning)
NAFLD_scores$Lobular.inflammation <- as.character(NAFLD_scores$Lobular.inflammation)
NAFLD_scores$Steatosis <- as.character(NAFLD_scores$Steatosis)
#Upload of the counts matrix. #We eliminate the beginning of the patients ID number so there are only 4 numbers (that matches with the pathological metadata uploaded before).
#DGEList creates DGE list that contains the counts for each sample. The output is a counts matrix with the patients as rownames and the counts of each sample on each column.
feature_counts <- read.csv("./featurecounts+coordsorted+gene_pilot.csv", row.names = 1, stringsAsFactors = FALSE)
names(feature_counts)<- gsub(pattern = "SMS.IC.014.", replacement = "", names(feature_counts))
names(feature_counts) <- gsub(pattern = "SMSIC.014.", replacement = "",names(feature_counts))
counts_matrix <- DGEList(counts = feature_counts, genes = rownames(feature_counts), samples = NAFLD_scores, group = NAFLD_scores$fibrosis_combined)

#A data frame with different colours for each group from the NAFLD score (the pathological data of each patient) is created.
annotcolNAS<- brewer.pal(9, "Paired")
names(annotcolNAS) <- c("0","1","2","3","4","5","6","7","8")
annotcolfibrosis <- brewer.pal(7, "Set2")
names(annotcolfibrosis) <- c("0", "1a", "1b", "1c", "2", "3", "4")
annotcolsteatosis <- brewer.pal(4, "Oranges")
names(annotcolsteatosis) <- c("0","1","2","3")
annotcolLI <- brewer.pal(4, "Purples")
names(annotcolLI) <- c("0","1","2","3")
annotcolHB <- brewer.pal(3, "Greys")
names(annotcolHB) <- c("0","1","2")
annotCol <- list("NAS.score" = annotcolNAS, "Fibrosis" = annotcolfibrosis, "Steatosis" = annotcolsteatosis, "Lobular.inflammation" = annotcolLI, "Hepatocyte.ballooning" = annotcolHB)

#######################
# Pre-processing data #
#######################

#The log2_boxplot function is used to create boxplots with the log2-transformed raw counts.
logcounts <- log2(counts_matrix$counts + 1)
boxplot(logcounts, main = "", xlab = "", ylab = "Raw read counts per gene (log2)", axes = FALSE)
axis(2)
axis(1,  at=c(1:length(rownames(counts_matrix$samples))),labels=colnames(counts_matrix$counts),las=2,cex.axis=0.6)

#We filter the counts to keep genes that have count-per-million (CPM) above k in at least n samples.
#k is the minimum count and n is determined by the design matrix/group (the smallest  number of samples in a group).
keep <- rowSums(cpm(counts_matrix)> 0.1) >= 10
summary(keep) #Check which counts are kept 
counts_filtered <- counts_matrix[keep, , keep.lib.sizes = FALSE]

#Another graph was obtained for the filtered counts.
logcounts_filtered <- log2(counts_filtered$counts + 1)
boxplot(logcounts_filtered, main = "", xlab = "", ylab = "Filtered read counts per gene (log2)", axes = FALSE)
axis(2)
axis(1,  at=c(1:length(rownames(counts_filtered$samples))),labels=colnames(counts_mfiltered$counts),las=2,cex.axis=0.6)

#We perform TMM normalisation to eliminate composition biases between the libraries. It assumes that the majority of the genes are not DE. 
#It estimates scale factors between samples that can be incorporated into statistical analysis for DE.
counts_normalised <- calcNormFactors(counts_filtered, method = "TMM")
counts_normalised$samples  #Check the norm factors have been modified (it does not have a value of 1 anymore)
log_cpm <- cpm(counts_normalised$counts, log = TRUE, prior.count = 1)

#Boxplots with log2-transformed filtered and normalised counts are obtained.
boxplot(log_cpm, main = "", xlab = "", ylab = "Normalised cpm read counts per gene (log2)", axes = FALSE)
axis(2)
axis(1,  at=c(1:length(rownames(counts_normalised$samples))),labels=colnames(counts_normalised$counts),las=2,cex.axis=0.6)

#The top 50 most variant genes are obtained. This is calculated by performing the log2 transformation to the filtered and normalised counts.
rownames(log_cpm) <- counts_normalised$genes[,1]
select <- order(rowMeans(log_cpm), decreasing =TRUE)[1:50]
highexprgenes_counts <- log_cpm[select,] #Create a data frame with the scores of the pathological data.
gene_symbols <- read.csv('../SteatoSITE/Gene_symbols_ensembl.csv', row.names = 1) #Upload the information regarding the gene symbols for each ensembl ID.
symbols <- gene_symbols[,2]
names(symbols) <- gene_symbols[,1]
symbols2 <- symbols[rownames(highexprgenes_counts)]
rownames(highexprgenes_counts) <- paste(rownames(highexprgenes_counts), symbols2, sep = '-')
rownames(highexprgenes_counts)[rownames(highexprgenes_counts) == "ENSG00000281454-NA"] <- "ENSG00000281454-AC099314.1"
pdf(file = "./Pre_processed/Pheatmap_highexpressed.pdf", width = 13, height = 14)
pheatmap(highexprgenes_counts, legend = FALSE, annotation_col = counts_normalised$samples[4:8], scale = "row", border_color = NA, show_colnames = FALSE, annotation_colors = annotCol)
dev.off()

#Clustering of the samples with the filtered and normalised values for all the samples. The Euclidean distance is performed.
sampleDists <- dist(t(log_cpm))
sampleDistMatrix <- as.matrix(sampleDists)
dataframe <- data.frame(row.names = c(colnames(log_cpm)))
rownames(sampleDistMatrix) <- rownames(dataframe)
colnames(sampleDistMatrix) <- rownames(dataframe)
pdf(file = "./Pre_processed/HClust_sampledist.pdf")
plot(hclust(sampleDists))
dev.off()
pdf(file = "./Pre_processed/Pheatmap_sampledist.pdf", width = 13, height = 14)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, annotation_col = counts_normalised$samples[, c(4,6,8,10)], legend = FALSE, annotation_colors = annotCol, border_color = NA)
dev.off()

#PCA analysis: this is going to extract the fundamental structure of the data. This will provide information about 
#the proportion of explained variance by calculating Eigenvalues. Any outliers or possible batch effects can be explained by this PCA plot.
PCA <- pca(log_cpm, metadata = counts_normalised$samples)
##The SCREE plot will explain the number of components that can explain the variability and how many components we should be looking at.
pdf(file = "./Pre_processed/Screeplot.pdf")
screeplot(PCA)
dev.off()

##The bi-plot will show two components against each other. This will indicate how the data is distributed.
pdf(file = "./Pre_processed/PC1-PC2_fibrosis-steatosis.pdf", width = 15, height = 14)
biplot(pcaobj = PCA, x = "PC1", y ="PC2", colby = "Fibrosis", shape = "Steatosis", legendPosition = "bottom", title= "Biplot of the logCPM", colLegendTitle = "Fibrosis", shapeLegendTitle = "Steatosis")
dev.off()
pdf(file = "./Pre_processed/PC1-PC2_LI-HB.pdf", width = 15, height = 14)
biplot(pcaobj = PCA, x = "PC1", y = "PC2", colby = "Lobular.inflammation",shape = "Hepatocyte.ballooning", legendPosition = "bottom", title= "Biplot of the logCPM", colLegendTitle ="Lobular inflammation", shapeLegendTitle = "Hepatocyte ballooning")
dev.off()

pdf(file = "./Pre_processed/PC2-PC3_fibrosis-steatosis.pdf", width = 15, height = 14)
biplot(pcaobj = PCA, x = "PC2", y = "PC3", colby = "Fibrosis", shape = "Steatosis", legendPosition = "bottom", title= "Biplot of the logCPM", colLegendTitle = "Fibrosis", shapeLegendTitle = "Steatosis")
dev.off()
pdf(file = "./Pre_processed/PC2-PC3_LI-HB.pdf", width = 15, height = 14)
biplot(pcaobj = PCA, x = "PC2", y = "PC3", colby = "Lobular.inflammation", shape = "Hepatocyte.ballooning", legendPosition = "bottom", title= "Biplot of the logCPM", colLegendTitle = "Lobular inflammation", shapeLegendTitle = "Hepatocyte ballooning")
dev.off()

##############################
#   Differential expression  #
##############################

##Fibrosis high - fibrosis low
Cirrhosis_results <- de_analysis(counts_normalised, group_var = 'Fibrosis_level', blocking_vars = 'Lobular.inflammation', contrast_name = 'Fibrosis_levelFibrosis_high-Fibrosis_levelFibrosis_low')
symbols2 <- symbols[rownames(Cirrhosis_results$results)] #Include the gene symbols in the results table.
Cirrhosis_results$results$symbol <- symbols2
entrez_id <- read.csv('./Entrez_id.txt')
entrez <- entrez_id[,1] #Include the entrez id in the results table.
names(entrez) <- entrez_id[,2]
entrez2 <- entrez[Cirrhosis_results$results$genes]
Cirrhosis_results$results$entrez <- entrez
write.csv(Cirrhosis_results$results, './New_results/DE/Cirrhosis.csv')
top_genes <- Cirrhosis_results$results[1:50,]
pdf('./New_results/DE/Cirrhosis_maplot.pdf', width = 15, height = 10)
ma_plot(Cirrhosis_results$results, 2, 0.05)
dev.off()
pdf('./New_results/DE/Cirrhosis_pheatmap.pdf', width = 15, height = 10)
plot_pheatmap(rownames(top_genes), 'Top 50 expressed genes in Cirrhosis-control', Cirrhosis_results, Cirrhosis_results$dgList[, which(Cirrhosis_results$dgList$samples$Fibrosis == "4"|Cirrhosis_results$dgList$samples$Fibrosis == "0")], c('Fibrosis', 'Hepatocyte.ballooning', 'Lobular.inflammation', 'Steatosis'))
dev.off()
pdf('./New_results/DE/Cirrhosis_volcanoplot.pdf', width = 15, height = 10)
volcano_plot(Cirrhosis_results$results)
dev.off()

##Fibrosis med- fibrosis low
Fibrosis_med_low <- de_analysis(counts_normalised, group_var = 'Fibrosis_level', blocking_vars = 'Lobular.inflammation', contrast_name = 'Fibrosis_levelFibrosis_med-Fibrosis_levelFibrosis_low')
symbols2 <- symbols[rownames(Fibrosis_med_low$results)]
Fibrosis_med_low$results$symbol <- symbols2
entrez2 <- entrez[Fibrosis_med_low$results$genes]
Fibrosis_med_low$results$entrez <- entrez2
write.csv(Fibrosis_med_low$results, './New_results/DE/Fibrosis_med-low.csv')
top_genes <- Fibrosis_med_low$results[1:50,]
pdf('./New_results/DE/Fibrosis_med-low_maplot.pdf', width = 15, height = 10)
ma_plot(Fibrosis_med_low$results, 2, 0.05)
dev.off()
pdf('./New_results/DE/Fibrosis_med-low_volcanoplot.pdf', width = 15, height = 10)
volcano_plot(Fibrosis_med_low$results)
dev.off()
pdf('./New_results/DE/Fibrosis_med-low_pheatmap.pdf', width = 15, height = 10)
plot_pheatmap(rownames(top_genes), 'Top 50 most expressed genes in Fibrosis med-low', Fibrosis_med_low, Fibrosis_med_low$dgList[, which(Fibrosis_med_low$dgList$samples$Fibrosis == '2'| Fibrosis_med_low$dgList$samples$Fibrosis == '0'| Fibrosis_med_low$dgList$samples$Fibrosis == '1a'|Fibrosis_med_low$dgList$samples$Fibrosis == '1b'|Fibrosis_med_low$dgList$samples$Fibrosis == '1c'| Fibrosis_med_low$dgList$samples$Fibrosis== '3')], c('Fibrosis', 'Hepatocyte.ballooning', 'Steatosis', 'Lobular.inflammation'))
dev.off()

##Fibrosis high - med
Fibrosis_high_med <- de_analysis(counts_normalised, 'Fibrosis_level', 'Lobular.inflammation', contrast_name = 'Fibrosis_levelFibrosis_high-Fibrosis_levelFibrosis_med')
symbols2 <- symbols[rownames(Fibrosis_high_med$results)]
Fibrosis_high_med$results$symbol <- symbols2
entrez2 <- entrez[Fibrosis_high_med$results$genes]
Fibrosis_high_med$results$entrez <- entrez2
write.csv(Fibrosis_high_med$results, './New_results/DE/High-med.csv')
top_genes <- Fibrosis_high_med$results[1:50,]
pdf('./New_results/DE/Fibrosis_high-med_maplot.pdf', width = 15, height = 10)
ma_plot(Fibrosis_high_med$results, 2, 0.05)
dev.off()
pdf('./New_results/DE/Fibrosis_high-med_volcanoplot.pdf', width = 15, height = 10)
volcano_plot(Fibrosis_high_med$results)
dev.off()
pdf('./New_results/DE/Fibrosis_high-med_pheatmap.pdf', width = 15, height = 10)
plot_pheatmap(rownames(top_genes), 'Top 50 most expressed genes in Fibrosis mhigh - med', Fibrosis_high_med, Fibrosis_high_med$dgList[, which(Fibrosis_high_med$dgList$samples$Fibrosis == '2'| Fibrosis_high_med$dgList$samples$Fibrosis == '4'| Fibrosis_high_med$dgList$samples$Fibrosis == '1a'|Fibrosis_high_med$dgList$samples$Fibrosis == '1b'|Fibrosis_high_med$dgList$samples$Fibrosis == '1c'| Fibrosis_high_med$dgList$samples$Fibrosis== '3')], c('Fibrosis', 'Hepatocyte.ballooning', 'Steatosis', 'Lobular.inflammation'))
dev.off()

##Steatosis high-low
variant <- c("CG","CG","C","C", "CG", "C", "C", "C", "GG", "C", "C", "C", "CG", "C", "CG", "GT", "C" ,"CG", "CG", "C", "CG", "C",  "GG", "GG", "CG" ,"C", "GG", "CG", "GG", "C", "C", "GG", "C", "C", "CG", "C", "C", "C", "C", "C", "C", "unknown", "C", "C", "CG", "CG", "CG", "C", "C", "CG", "CG", "C", "CG", "C")
Steatosis <- de_analysis(counts_normalised, group_var = 'Steatosis_level', blocking_vars = 'variant', 'Steatosis_levelSteatosis_high-Steatosis_levelSteatosis_low')
symbols2 <- symbols[rownames(Steatosis$results)]
Steatosis$results$symbol <- symbols2
entrez2 <- entrez[Steatosis$results$genes]
Steatosis$results$entrez <- entrez2
write.csv(Steatosis$results, './New_results/DE/Steatosis.csv')
top_genes <- Steatosis$results[1:50,]
pdf('./New_results/DE/Steatosis_maplot.pdf', width = 15, height = 10)
ma_plot(Steatosis$results, 2, 0.05)
dev.off()
pdf('./New_results/DE/Steatosis_volcanoplot.pdf', width = 15, height = 10)
volcano_plot(Steatosis$results)
dev.off()
pdf('./New_results/DE/Steatosis_pheatmap.pdf', width = 15, height = 10)
plot_pheatmap(rownames(top_genes), 'Top 50 most expressed genes in Steatosis high -low', Steatosis, Steatosis$dgList, c('Steatosis', 'Hepatocyte.ballooning', 'Fibrosis', 'Lobular.inflammation'))
dev.off()

#Hepatocyte ballooning high -low
Hepatocyte_ballooning <- de_analysis(counts_normalised, 'HB', contrast_name = 'HBHB_high-HBHB_low')
symbols2 <- symbols[rownames(Hepatocyte_ballooning$results)]
Hepatocyte_ballooning$results$symbol <- symbols2
entrez2 <- entrez[Hepatocyte_ballooning$results$genes]
Hepatocyte_ballooning$results$entrez <- entrez2
write.csv(Hepatocyte_ballooning$results, './New_results/DE/Hepatocyte_ballooning.csv')
top_genes <- Hepatocyte_ballooning$results[1:50,]
pdf('./New_results/DE/Hepatocyte_ballooning_maplot.pdf', width = 15, height = 10)
ma_plot(Hepatocyte_ballooning$results, 2, 0.05)
dev.off()
pdf('./New_results/DE/Hepatocyte_ballooning_volcanoplot.pdf', width = 15, height = 10)
volcano_plot(Hepatocyte_ballooning$results)
dev.off()
pdf('./New_results/DE/Hepatocyte_ballooning_pheatmap.pdf', width = 15, height = 10)
plot_pheatmap(rownames(top_genes), 'Top 50 most expressed genes in hepatocyte ballooning high-low', Hepatocyte_ballooning, Hepatocyte_ballooning$dgList, c('Hepatocyte.ballooning', 'Steatosis', 'Fibrosis', 'Lobular.inflammation'))
dev.off()

#Loublar inflammation high-low
Lobular_inflammation <- de_analysis(counts_normalised, 'LI', c('Fibrosis', 'variant'), 'LILI_high-LILI_low')
symbols2 <- symbols[rownames(Lobular_inflammation$results)]
Lobular_inflammation$results$symbol <- symbols2
entrez2 <- entrez[Lobular_inflammation$results$genes]
Lobular_inflammation$results$entrez <- entrez2
write.csv(Lobular_inflammation$results, './New_results/DE/Lobular_inflammation.csv')
top_genes <- Lobular_inflammation$results[1:50,]
pdf('./New_results/DE/Lobular_inflammation_maplot.pdf', width = 15, height = 10)
ma_plot(Lobular_inflammation$results, 2, 0.05)
dev.off()
pdf('./New_results/DE/Lobular_inflammation_pheatmap.pdf', width = 15, height = 10)
volcano_plot(Lobular_inflammation$results)
dev.off()
pdf('./New_results/DE/Lobular_inflammation_pheatmap.pdf', width = 15, height = 10)
plot_pheatmap(rownames(top_genes), 'Top 50 most expressed genes in lobular inflammation high-low', Lobular_inflammation, Lobular_inflammation$dgList, c('Lobular.inflammation', 'Hepatocyte.ballooning', 'Steatosis', 'Fibrosis'))
dev.off()


#######################
# Functional analysis #
#######################

#Cirrhosis
Pathway_analysis('./DE/Cirrhosis.csv', 'cirrhosis', 1)
Pathway_analysis('./DE/Fibrosis_med-low.csv', 'fibrosis_med-low', 1)
Pathway_analysis('./DE/Fibrosis_high-med.csv', 'fibrosis_high_med', 1)
Pathway_analysis('./DE/Steatosis.csv', 'steatosis',1)
Pathway_analysis('./DE/Lobular_inflammation.csv', 'lobular_inflammation', 1)
Pathway_analysis('./DE/Hepatocyte_ballooning.csv', 'hepatocyte_ballooning', 1)
