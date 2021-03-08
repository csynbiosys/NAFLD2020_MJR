library("edgeR")
library("pheatmap")
library("RColorBrewer")
library("PCAtools")
library("matrixStats")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("genefilter")
library("ggplot2")
library("pathview")
library("clusterProfiler")
library("enrichplot")

#Information about the pathological of each patient that would be used in the model matrix for DE analysis.
HB_score <- factor(c("score_0", "score_1", "score_1", "score_1", "score_1", "score_0", "score_2", "score_2", "score_2", "score_1", "score_0", "score_0", "score_1", "score_1", "score_2", "score_2", "score_2", "score_1", "score_2", "score_0", "score_1", "score_0", "score_1", "score_2", "score_0", "score_2", "score_2","score_2", "score_2", "score_0", "score_1", "score_2", "score_0", "score_0", "score_2", "score_2", "score_2", "score_1", "score_1", "score_1", "score_0", "score_2", "score_0", "score_1", "score_0", "score_2", "score_0", "score_2", "score_0", "score_1", "score_2", "score_0", "score_2", "score_1"))
HB <- factor(c("HB_low", "HB_high", "HB_high", "HB_high", "HB_high", "HB_low", "HB_high", "HB_high", "HB_high", "HB_high", "HB_low", "HB_low", "HB_high", "HB_high", "HB_high", "HB_high", "HB_high", "HB_high", "HB_high", "HB_low", "HB_high", "HB_low","HB_high", "HB_high", "HB_low", "HB_high", "HB_high", "HB_high", "HB_high", "HB_low", "HB_high", "HB_high", "HB_low", "HB_low", "HB_high", "HB_high", "HB_high", "HB_high", "HB_high", "HB_high", "HB_low", "HB_high", "HB_low", "HB_high", "HB_low", "HB_high", "HB_low", "HB_high", "HB_low", "HB_high", "HB_high", "HB_low", "HB_high", "HB_high"))
LI <- factor(c("LI_low", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_low", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_low", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_low", "LI_low", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_low", "LI_high", "LI_low", "LI_low", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high", "LI_high"))
LI_score <- factor(c("score_0", "score_1", "score_1", "score_1", "score_2", "score_1", "score_1", "score_2", "score_1", "score_1", "score_0", "score_1", "score_1", "score_1", "score_2", "score_2", "score_1", "score_2", "score_1", "score_1", "score_2", "score_1", "score_1", "score_3", "score_0", "score_2", "score_3", "score_3", "score_1", "score_1", "score_2", "score_3", "score_0", "score_0", "score_1", "score_1", "score_2", "score_3", "score_3", "score_2", "score_1", "score_3", "score_0", "score_2", "score_0", "score_0", "score_1", "score_3", "score_1", "score_3", "score_2", "score_2", "score_2", "score_2"))
steatosis_score <- factor(c("score_0", "score_2", "score_1", "score_1", "score_2", "score_1", "score_1", "score_1", "score_2", "score_2", "score_0", "score_1", "score_1", "score_1", "score_1", "score_1", "score_2", "score_2", "score_3", "score_1", "score_2", "score_2", "score_1", "score_2", "score_0", "score_3", "score_3", "score_3", "score_3", "score_2", "score_2", "score_2", "score_1", "score_1", "score_2", "score_2", "score_1", "score_1", "score_3", "score_1", "score_2", "score_3", "score_0", "score_1", "score_2", "score_0", "score_1", "score_2", "score_0", "score_2", "score_1", "score_3", "score_1", "score_1"))
steatosis <- factor(c("steatosis_low", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_low", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_low", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_low", "steatosis_high", "steatosis_high", "steatosis_low", "steatosis_high", "steatosis_high", "steatosis_low", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high", "steatosis_high"))
fibrosis_score <- c("score_0", "score_1c", "score_0", "score_4", "score_1c", "score_1a", "score_1c", "score_3", "score_3", "score_1a", "score_4", "score_1a", "score_3", "score_2", "score_1b", "score_2", "score_4", "score_2", "score_4", "score_0", "score_3", "score_1a", "score_4", "score_4", "score_4", "score_3", "score_3", "score_3", "score_0", "score_0", "score_3", "score_3", "score_0", "score_3", "score_2", "score_0", "score_1b", "score_0", "score_4", "score_0", "score_0", "score_2", "score_1c", "score_3", "score_3", "score_4", "score_1a", "score_1a", "score_4", "score_3", "score_3", "score_2", "score_4", "score_1c")
fibrosis <- factor(c("fibrosis_low", "fibrosis_med", "fibrosis_low", "fibrosis_high", "fibrosis_med", "fibrosis_med", "fibrosis_med", "fibrosis_med", "fibrosis_med", "fibrosis_med", "fibrosis_high", "fibrosis_med", "fibrosis_med", "fibrosis_med", "fibrosis_med", "fibrosis_med", "fibrosis_high", "fibrosis_med", "fibrosis_high", "fibrosis_low", "fibrosis_med", "fibrosis_med", "fibrosis_high", "fibrosis_high", "fibrosis_high", "fibrosis_med", "fibrosis_med", "fibrosis_med", "fibrosis_low", "fibrosis_low", "fibrosis_med", "fibrosis_med", "fibrosis_low", "fibrosis_med", "fibrosis_med", "fibrosis_low", "fibrosis_med", "fibrosis_low", "fibrosis_high", "fibrosis_low", "fibrosis_low", "fibrosis_med", "fibrosis_med", "fibrosis_med", "fibrosis_med", "fibrosis_high", "fibrosis_med", "fibrosis_med", "fibrosis_high", "fibrosis_med", "fibrosis_med", "fibrosis_med", "fibrosis_high", "fibrosis_med"))
samples <- data.frame(HB, HB_score, LI,LI_score, steatosis, steatosis_score, fibrosis, fibrosis_score)
samples

#Upload of the counts matrix. We create a DGE list that contains the counts for each sample. 
#The output is a counts matrix with the patients as rownames and the counts of each sample in each column.
#The names are going to be shortened with the gsub fuction for simplicity.
feature_counts <- read.table("../featurecounts+coordsorted+gene_pilot.csv", header =TRUE, sep = ",")
counts_matrix <- DGEList(feature_counts[,2:55], genes = feature_counts[,1], group = factor(fibrosis))
colnames(counts_matrix) <- gsub(pattern = "SMS.IC.014.", replacement = "", colnames(counts_matrix))
colnames(counts_matrix) <- gsub(pattern = "SMSIC.014.", replacement = "", colnames(counts_matrix))
counts_matrix$samples
head(counts_matrix$counts)
rownames(samples) <- colnames(counts_matrix)

#Density plots: this can give an idea of log-intensity distribution of each library. It is expected to be similar, but not identical. 
pdf(file = "./Pre_processed/Density_plot.pdf")
logcounts <- log(counts_matrix$counts[,1],10)
d <- density(logcounts)
plot(d, xlim = c(1,9), main = "", ylim= c(0, 0.5), xlab = "Raw read counts per gene (log10)", ylab = "Density")
for (s in 2:length(rownames(counts_matrix$samples))){
  logcounts <- log(counts_matrix$counts[,s], 10)
  d <- density(logcounts)
  lines(d)
}
dev.off()

#Boxplot for each sample with the raw counts. For visualisation the log2 transformation was calculated.
logcounts_raw <- log2(counts_matrix$counts + 1)
pdf(file = "./Pre_processed/Boxplot_notnormalised.pdf")
boxplot(logcounts_raw, main = "", xlab = "", ylab = "Raw read counts per gene (log2)", axes = FALSE)
axis(2)
axis(1,  at=c(1:length(rownames(counts_matrix$samples))),labels=colnames(counts_matrix$counts),las = 2,cex.axis = 0.6)
dev.off()

#We filter the counts to keep genes that have count-per-million (CPM) above k in at least n samples.
#k is the minimum count and n is determined by the design matrix/group (the smallest  number of samples in a group).
keep <- rowSums(cpm(counts_matrix) > 0.1) >= 10
summary(keep) #Check which counts are kept -> 37315
counts_filtered <- counts_matrix[keep, , keep.lib.sizes = FALSE]
counts_filtered$samples #The library sizes are smaller now, as it contains less reads.
dim(counts_filtered) #Final amount of counts kept and number of samples.

#Boxplot of each sample for the filtered read counts after a log2 transformation. The number of outliers should have decreased. 
logcounts_filtered <- log2(counts_filtered$counts +1)
pdf(file = "./Pre_processed/Boxplot_filtered.pdf")
boxplot(logcounts_filtered, main = "", xlab = "", ylab = "Filtered read counts per gene (log2)", axes = FALSE)
axis(2)
axis(1,  at=c(1:length(rownames(counts_filtered$samples))),labels=colnames(counts_filtered$counts),las=2,cex.axis=0.6)
dev.off()

#We perform TMM normalisation to eliminate composition biases between the libraries. It assumes that the majority of the genes are not DE. 
#It estimates scale factors between samples that can be incorporated into statistical analysis for DE.
counts_normalised <- calcNormFactors(counts_filtered, method = "TMM")
counts_normalised$samples  #Check the norm factors have been modified (it does not have a value of 1 anymore)
head(counts_normalised$counts) 

#Boxplot after normalisation of the samples. For visualisation purposes the normalised counts are represented in cpm with a log2 transformation.
#The reason why cpm is performed now is because the counts are now normalised, the log2 transformation have to be scaled to the library size and 
#the normalised library sizes are taken into account for the cpm values.
logCPMs <- cpm(counts_normalised, log = TRUE)
pdf(file = "./Pre_processed/Boxplot_normalised.pdf")
boxplot(logCPMs, main = "", xlab = "", ylab = "Normalised read counts per gene (log2)", axes = FALSE)
axis(2)
axis(1,  at=c(1:length(rownames(counts_normalised$samples))),labels=colnames(counts_normalised$counts),las=2,cex.axis=0.6)
dev.off()

#The top 50 most variant genes are obtained. This is calculated by performing the log2 transformation to the filtered and normalised counts.
#The euclidean distances are obtained to cluster the samples, and 1- Pearson correlation to obtain the cluster of the genes.
rownames(logCPMs) <- counts_normalised$genes[,1]
select <- order(rowMeans(logCPMs), decreasing =TRUE)[1:50]
highexprgenes_counts <- logCPMs[select,]
annot <- as.data.frame(samples[,c(2,4,6,8)])  #Create a data frame with the scores of the pathological data.
annotCol <- list("fibrosis_score" = c("score_0" = "#8DD3C7", "score_1a" = "#FFFFB3", "score_1b" = "#BEBADA", "score_1c" = "#FB8072", "score_2" = "#FDB462", "score_3" = "#80B1D3", "score_4" = "#B3DE69"), "steatosis_score" = c("score_0" = "#E41A1C", "score_1" = "#377EB8", "score_2" = "#4DAF4A", "score_3" = "#984EA3"), "LI_score" = c("score_0" = "#1B9E77", "score_1" = "#D95F02", "score_2" = "#7570B3", "score_3" = "#E7298A"))
rownames(annot) <- colnames(highexprgenes_counts)
rownames(highexprgenes_counts) <- paste(rownames(highexprgenes_counts), (mapIds(org.Hs.eg.db, keys = rownames(highexprgenes_counts), keytype = "ENSEMBL", column = "SYMBOL")), sep = "-") #Each ensembl ID would be link to its Symbol annotation.
pdf(file = "./Pre_processed/Pheatmap_highexpressed.pdf", width = 12, height = 14)
pheatmap(highexprgenes_counts, legend = FALSE, annotation_col = annot, annotation_colors = annotCol, scale = "row", border_color = NA, clustering_distance_rows = "correlation")
dev.off()

#Clustering of the samples with the filtered and normalised values for all the samples. The Euclidean distance is performed.
sampleDists <- dist(t(logCPMs))
sampleDistMatrix <- as.matrix(sampleDists)
dataframe <- data.frame(row.names = c(colnames(logCPMs)))
rownames(sampleDistMatrix) <- rownames(dataframe)
colnames(sampleDistMatrix) <- rownames(dataframe)
rownames(annot) <- colnames(sampleDistMatrix)  #The rownames of the annot data frame are modified to match the current data.
pdf(file = "./Pre_processed/HClust_sampledist.pdf")
plot(hclust(sampleDists))
dev.off()
pdf(file = "./Pre_processed/Pheatmap_sampledist.pdf", width = 12, height = 14)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, annotation_col = annot, legend = FALSE, annotation_colors = annotCol, border_color = NA)
dev.off()

#PCA analysis: this is going to extract the fundamental structure of the data. This will provide information about 
#the proportion of explained variance by calculating Eigenvalues. Any outliers or possible batch effects can be explained by this PCA plot.
select2 <- order(rowVars(logCPMs), decreasing =TRUE)
PCA <- pca(logCPMs[select2,])
PCA$metadata <- data.frame(samples$fibrosis_score, samples$steatosis_score, samples$LI_score, samples$HB_score)
##The SCREE plot will explain the number of components that can explain the variability and how many components we should be looking at.
pdf(file = "./Pre_processed/Screeplot.pdf")
screeplot(PCA)
dev.off()
##The bi-plot will show the two components with the biggest variance against each other. This will indicate how the data is distributed.
pdf(file = "./edgeR/Pre_processed/Biplot_fibrosis-steatosis.pdf")
biplot(pcaobj = PCA, x = "PC1", y = "PC2", colby = "samples.fibrosis_score", shape = "samples.steatosis_score", legendPosition = "bottom", title= "Biplot of the logCPM")
dev.off()
pdf(file = "./edgeR/Pre_processed/Biplot_LI-HB.pdf")
biplot(pcaobj = PCA, x = "PC1", y = "PC2", colby = "samples.LI_score", shape = "samples.HB_score", legendPosition = "bottom", title= "Biplot of the logCPM")
dev.off()

pdf(file = "./edgeR/Pre_processed/Pairsplot_fibrosis-steatosis.pdf")
biplot(pcaobj = PCA, x = "PC2", y = "PC3", colby = "samples.fibrosis_score", shape = "samples.steatosis_score", legendPosition = "bottom", title= "Biplot of the logCPM")
dev.off()
pdf(file = "./edgeR/Pre_processed/Pairsplot_fibrosis-steatosis.pdf")
biplot(pcaobj = PCA, x = "PC2", y = "PC3", colby = "samples.LI_score", shape = "samples.HB_score", legendPosition = "bottom", title= "Biplot of the logCPM")
dev.off()

#Mean-difference plot of all log fold change against average count size of the first sample to verify the TMM normalisation. They have to be around the expression
#log ratio 0.
png(file = "./Pre_processed/MDplot_normalisedcounts.png")
plotMD(logCPMs)
abline(h = 0, col="red", lty=2, lwd=2)
dev.off()

##############################
#   Differential expression  #
##############################

##HB DEG

#We create the design matrix that is going to contain all the groups and distribution of the samples within the groups. 
#First we need to relevel the HB group, so the comparison is high-low. The columns of the design matrix are renamed for simplicity.
HB <- relevel(HB, ref = "HB_low")
design <- model.matrix(~HB, data = counts_normalised$samples)
colnames(design) <- c("HB_low", "HB_high")
design

#We estimate the dispersion of the normalised counts.This function calculates a matrix of likelihoods for each tag at a set of dispersion
#grid points, and then applies weighted likelihood empirical Bayes method to obtain dispersion estimates. 
estimated_dispersion <- estimateDisp(counts_normalised, design, robust = TRUE)
estimated_dispersion$common.dispersion
sqrt(estimated_dispersion$common.dispersion)

#To visualise the estimated dispersion, we generate a plot for the Biological Coefficient of Variation (BCV). 
#This is the square root of the negative binomial dispersion. It is going to represent the common, tagwise and trended BCV estimates.
#When the gene abundance is small, the BCV is bigger, and viceversa.
pdf(file = "./DE/BCV_plot_HB.pdf")
plotBCV(estimated_dispersion)
dev.off()

#Once dispersion estimates are obtained, the glmQLFit functions fits a quasi-likelihood negative binomial generalised log-linear 
#model (GLM) to count data. Then, the following step is to determine the differential expression through the quasi-likelihood (QL) 
#F-test or likelihood ratio test. Depending of the number of replicates and samples available, it is better to use one test or 
#another. In here, the QLFTest is performed. 
fit <- glmQLFit(estimated_dispersion, design, robust = TRUE)
head(fit$coefficients)

#The next step  is to conduct the likelihood ratio test for one or more coefficients in the linear model.
#The null hypothesis is that all the coefficients are equal to zero. This will give back a table with the logFC, logCPM, LR, p-values and FDR.
qlf_hb <- glmQLFTest(fit)
de <- decideTests(qlf_hb) #53 down
pdf(file = "./DE/SmearPlot_HB.pdf")
plotSmear(qlf_hb, de.tags = rownames(qlf_hb)[as.logical(de)])
dev.off()

#The most differential expressed genes can be extracted with topTags. In this case, those genes with a p-value < 0.05 and FDR < 0.05 -order by logFC- have been arranged in a table.
#An histogram to observe the frequency of genes with p-value < 0.05 is recommended as a diagnostic plot. 
#For visual aid, gene symbols and entrez ID have been added to the final table.
results <- topTags(qlf_hb, n= Inf)$table
pdf(file ="./DE/Hist_pvalue_HB.pdf")
hist(results$PValue, breaks = 0:30/30,
     col = "grey50", border = "white", main = "p-value histogram", xlab = "p-values")
dev.off()

sum(results$FDR < 0.05, na.rm = TRUE)
results_FDR <- subset(results, FDR < 0.05)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_annot <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", mart = ensembl, values = results_FDR$genes, uniqueRows = FALSE)
results_FDR <- results_FDR[order(results_FDR$genes), ] 
results_FDR$entrez.id <- ensembl_annot[,1]
results_FDR$symbol <- ensembl_annot[,3]

results_FDR <- results_FDR[order(results_FDR$logFC), ]
resOrderedDF <- as.data.frame(results_FDR)
write.csv(resOrderedDF, file = "./DE/QLF_HB.csv")

#Volcano plot illustrating differential expression between groups HB_low and HB_high as defined by HB.
results$threshold <- as.factor(results$FDR < 0.05)
pdf(file = "./DE/Volcanoplot_HB.pdf")
ggplot(data = results, aes(x = logFC, y= -log10(FDR), colour = threshold)) + 
  geom_point(size = 0.8) + 
  scale_color_manual(values = c("Grey", "Red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = -log2(2), linetype = "dashed", size = 0.3) +
  xlim(c(-6, 6)) +
  xlab("log2 fold change") + ylab("-log10 FDR") +
  theme_bw() +
  theme(legend.position= "none")
dev.off()

#Pheatmap clustering of the top 50 most DEG between HB_low and HB_high as defined in HB. This is ordered by FDR.
rownames(logCPMs) <- results$genes
colnames(logCPMs) <- paste(colnames(logCPMs), HB, sep = "-")
o <- order(results_FDR$FDR)
deg_HB <- results_FDR[o[1:50],]
logCPMs_HB <- logCPMs[deg_HB$genes, ]
rownames(logCPMs_HB) <- paste(rownames(logCPMs_HB), deg_HB$symbol, sep = "-")
rownames(annot) <- colnames(deg_HB)
pdf("./DE/Pheatmap_HB.pdf", width = 12, height = 14)
pheatmap(logCPMs_HB, legend = FALSE, scale = "row", border_color = NA, clustering_distance_cols = "euclidean", clustering_distance_rows = "correlation", annotation_col = annot, annotation_colors = annotCol, clustering_method = "complete")
dev.off()

#GSEA analysis
#First we eliminate those genes that have no entrez ID or are duplicated. Then we rank the list by logFC. This ranked list is 
#introduced in the GSEA function. The gene sets that are going to be used as input have been obtained from the Broad Institute.
results_FDR <- subset(results_FDR, results_FDR$entrez!="NA")
results_FDR <- results_FDR[which(duplicated(results_FDR$entrez) == F), ]

#Another way to rank the genes is with the following lines, there is no preference for ranking.
##trial$sign <- sign(trial$logFC)
##trial$logP <- -log10(trial$PValue)
##trial$metric <- trial$sign * trial$logP
##trial[,c("metric")]

geneList <- results_FDR$logFC
names(geneList) <- as.character(results_FDR$entrez.id)
geneList <- sort(geneList, decreasing = TRUE)

#Reactome analysis. This will indicate those pathways that have been linked to some genes that appear in the geneList.
#Due to the small number of genes in the geneList, no term enriched under the specific pvalueCutoff has been obtained.
reactome <- read.gmt("./c2.cp.reactome.v7.2.entrez.gmt")
reactome_hb <- GSEA(geneList = geneList, TERM2GENE = reactome, nPermSimple = 1000000, eps = 0)
reactome_hb_results <- reactome_hb@result
write.csv(reactome_hb_results, "./DE/Reactome_hb.csv")
pdf("./DE/Ridgeplot_reactome_hb.pdf", width = 15, height = 20)
ridgeplot(reactome_hb, label_format = 20, showCategory = 20)
dev.off()
pdf("./DE/GSEAplot_reactome_hb_1.pdf")
gseaplot2(reactome_hb, geneSetID = 1)
dev.off()

#GO analysis. This will indicate which GO terms have been linked with genes present in the geneList.
#Due to the small number of genes in the geneList, no term enriched under the specific pvalueCutoff has been obtained.
go_bp <- read.gmt("./c5.go.bp.v7.2.entrez.gmt")
go_cc <- read.gmt("./c5.go.cc.v7.2.entrez.gmt")
go_mf <- read.gmt("./c5.go.mf.v7.2.entrez.gmt")
gobp_hb <- GSEA(geneList = geneList, TERM2GENE = go_bp, nPermSimple = 1000000, eps =0)
gobp_hb_results <- gobp_fib_high_low@result
write.csv(gobp_hb_results, "./DE/GObp_hb.csv")
gocc_hb <- GSEA(geneList = geneList, TERM2GENE = go_cc, nPermSimple = 1000000, eps = 0)
gocc_hb_results <- gocc_hb@result
write.csv(gocc_fib_high_low_results, "./DE/GOcc_hb.csv") 
gomf_hb <- GSEA(geneList = geneList, TERM2GENE = go_mf, nPermSimple = 1000000, eps = 0)
gomf_hb_results <- gomf_hb@result
write.csv(gomf_hb_results, "./DE/GOmf_hb.csv")
pdf("./DE/DotplotGObp_hb.pdf", width = 12, height = 14)
dotplot(gobp_hb)
dev.off()
pdf("./DE/DotplotGOcc_hb.pdf", width = 12, height = 14)
dotplot(gocc_hb)
dev.off()
pdf("./DE/DotplotGOmf_hb.pdf", width = 12, height = 14)
dotplot(gomf_hb)
dev.off()

#KEGG pathways. This will give any KEGG pathway that has been linked with specific genes that are present in the geneList. 
#The main difference with the reactome pathways is that the terms linked to the KEGG are less broad, hence there are not as many pathways as in the reactome list.
#Due to the small number of genes in the geneList, no term enriched under the specific pvalueCutoff has been obtained.
kegg <- read.gmt("./c2.cp.kegg.v7.2.entrez.gmt")
kegg_hb <- GSEA(geneList = geneList, TERM2GENE = kegg, nPermSimple = 1000000, eps = 0)
kegg_hb_results <- kegg_hb@result
write.csv(kegg_hb_results, "./DE/KEGG_fib_high_med.csv")
pathview(gene.data  = geneList, pathway.id = "", species = "hsa", limit = list(gene=max(abs(geneList)), cpd=1))

## FIBROSIS HIGH-LOW

#First, the patients classified with high level of fibrosis is taken out, so only fibrosis high and fibrosis low are compared.
#Since LI_score will be used as a blocking effect in the model matrix, the LI_score of these specific patients are selected too.
fib_high <- factor(subset(fibrosis, fibrosis!= "fibrosis_med"))
fib_high <- relevel(fib_high, ref = "fibrosis_low")
select_fib <- subset(counts_normalised$samples, counts_normalised$samples$group!= "fibrosis_med")
counts_normalised_fib_high <- counts_normalised[, rownames(select_fib)]
samples_fib_high <- subset(samples,rownames(samples) %in% colnames(counts_normalised_fib_high))
LI_score <- factor(c(samples_fib_high$LI_score))

#We create the design matrix that is going to contain all the groups and distribution of the samples within the groups.
design <- model.matrix(~LI_score + fib_high, data = counts_normalised_fib_high$samples)
colnames(design) <- gsub(pattern = "LI_score", replacement = "", colnames(design))
colnames(design) <- gsub(pattern = "fib_high", replacement = "", colnames(design))
design

#We estimate the dispersion of the normalised counts.
estimated_dispersion <- estimateDisp(counts_normalised_fib_high, design, robust = TRUE)
estimated_dispersion$common.dispersion
sqrt(estimated_dispersion$common.dispersion)

#To visualise the estimated dispersion, we generate a plot for the Biological Coefficient of Variation (BCV). 
pdf(file = "./DE/BCV_plot_Fib_high-Fib_low.pdf")
plotBCV(estimated_dispersion)
dev.off()

#Next we fit a quasi-likelihood negative binomial generalized log-linear model to count data and conduct genewise statistical 
#tests for a given coefficient or contrast.
fit <- glmQLFit(estimated_dispersion, design, robust =TRUE)
head(fit$coefficients)
plotQLDisp(fit)
qlf_fib_high <- glmQLFTest(fit)
de <- decideTests(qlf_fib_high) #3480 down, 9971 up
pdf(file = "./DE/SmearPlot_Fibrosis_high-Fibrosis_low.pdf")
plotSmear(qlf_fib_high, de.tags = rownames(qlf_fib_high)[as.logical(de)])
dev.off()

#We perform the diagnostic plot and extract the DEG.
results <- topTags(qlf_fib_high, n= Inf)$table
pdf(file ="./DE/Hist_pvalueFibrosis_high-Fibrosis_low.pdf")
hist(results$PValue, breaks = 0:30/30,
     col = "grey50", border = "white", main = "P-value histogram", xlab = "p-values")
dev.off()

sum(results$FDR < 0.05, na.rm = TRUE) #13451
results$symbol <- mapIds(org.Hs.eg.db, keys=results$genes, column="SYMBOL", keytype="ENSEMBL")
results$entrez <- mapIds(org.Hs.eg.db, keys= results$genes,column="ENTREZID", keytype="ENSEMBL")

results_FDR <- subset(results, FDR < 0.05)
results_FDR <- results_FDR[order(results_FDR$logFC), ]
resOrderedDF <- as.data.frame(results_FDR)
write.csv(resOrderedDF, file = "./DE/QLF_Fibrosis_high-Fibrosis_low.csv")

#Volcano plot illustrating differential expression between groups fibrosis_low and fibrosis_high as defined by fibrosis.
results$threshold <- as.factor(results$FDR < 0.05)
pdf("./DE/Volcano_plot_Fibrosis-high_Fibrosis-low.pdf")
ggplot(data = results, aes(x =logFC, y= -log10(FDR), colour = threshold)) + 
  geom_point(size = 0.8) + 
  scale_color_manual(values = c("Grey", "Red")) +
  geom_hline(yintercept = 1.30103, linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = c(log2(2), -log2(2)), linetype = "dashed", size = 0.3) +
  xlim(c(-10, 10)) +
  xlab("log2 fold change") + ylab("-log10 FDR") +
  theme_bw() +
  theme(legend.position= "none")
dev.off()

#Pheatmap clustering of the top 50 most DEG between fibrosis_low and fibrosis_high as defined by fibrosis. This is ordered by FDR.
logCPMs_fib_high <- cpm(counts_normalised_fib_high, log = TRUE)
rownames(logCPMs_fib_high) <- counts_normalised_fib_high$genes[,1]
colnames(logCPMs_fib_high) <- paste(colnames(logCPMs_fib_high), fib_high, sep = "-")
o <- order(results_FDR$FDR)
results_fib_high <- results_FDR[o[1:50],]
deg_fib_high <- logCPMs_fib_high[results_fib_high$genes,]
annot_fib_high <- as.data.frame(samples_fib_high[, c(2,4,6,8)])
rownames(deg_fib_high) <- paste(rownames(deg_fib_high), mapIds(org.Hs.eg.db, keys = rownames(deg_fib_high), keytype = "ENSEMBL", column = "SYMBOL"), sep = "-")
rownames(annot_fib_high) <- colnames(deg_fib_high)
pdf("./DE/Pheatmap_Fibrosis-high_Fibrosis-low.pdf", width = 12, height = 14)
pheatmap(deg_fib_high, scale = "row", legend = FALSE, border_color = NA, clustering_distance_rows = "correlation", clustering_distance_cols = "euclidean", annotation_col = annot_fib_high, clustering_method = "complete")
dev.off()

#GSEA analysis
#First we eliminate those genes that have no entrez ID or are duplicated. This ranked list is introduced in the GSEA function. 
#The gene sets that are going to be used as input have been obtained from the Broad Institute.
results_FDR <- subset(results_FDR, results_FDR$entrez!="<NA>")
results_FDR <- results_FDR[which(duplicated(results_FDR$entrez) == F), ]

geneList <- results_FDR$logFC
names(geneList) <- as.character(results_FDR$entrez)
geneList <- sort(geneList, decreasing = TRUE)

#Reactome analysis
reactome <- read.gmt("./c2.cp.reactome.v7.2.entrez.gmt")
reactome_fib_high_low <- GSEA(geneList = geneList, TERM2GENE = reactome, nPermSimple = 1000000, eps = 0)
reactome_fib_high_low_results <- reactome_fib_high_low@result
write.csv(reactome_fib_high_low_results, "./DE/Reactome_fib_high_low.csv")
pdf("./DE/Ridgeplot_reactome_fib_high_low.pdf", width = 15, height = 20)
ridgeplot(reactome_fib_high_low, label_format = 20, showCategory = 20)
dev.off()
pdf("./DE/GSEAplot_reactome_fib_high_low_1.pdf")
gseaplot2(reactome_fib_high_low, geneSetID = 1)
dev.off()

#GO analysis
go_bp <- read.gmt("./c5.go.bp.v7.2.entrez.gmt")
go_cc <- read.gmt("./c5.go.cc.v7.2.entrez.gmt")
go_mf <- read.gmt("./c5.go.mf.v7.2.entrez.gmt")
gobp_fib_high_low <- GSEA(geneList = geneList, TERM2GENE = go_bp, nPermSimple = 1000000, eps =0)
gobp_fib_high_low_results <- gobp_fib_high_low@result
write.csv(gobp_fib_high_low_results, "./DE/GObp_fib_high_low.csv")
gocc_fib_high_low <- GSEA(geneList = geneList, TERM2GENE = go_cc, nPermSimple = 1000000, eps = 0)
gocc_fib_high_low_results <- gocc_fib_high_low@result
write.csv(gocc_fib_high_low_results, "./DE/GOcc_fib_high_low.csv") 
gomf_fib_high_low <- GSEA(geneList = geneList, TERM2GENE = go_mf, nPermSimple = 1000000, eps = 0)
gomf_fib_high_low_results <- gomf_fib_high_low@result
write.csv(gomf_fib_high_low_results, "./DE/GOmf_fib_high_low.csv")
pdf("./DE/DotplotGObp_fib_high_low.pdf", width = 12, height = 14)
dotplot(gobp_fib_high_low)
dev.off()
pdf("./DE/DotplotGOcc_fib_high_low.pdf", width = 12, height = 14)
dotplot(gocc_fib_high_low)
dev.off()
pdf("./DE/DotplotGOmf_fib_high_low.pdf", width = 12, height = 14)
dotplot(gomf_fib_high_low)
dev.off()

#KEGG pathways
kegg <- read.gmt("./c2.cp.kegg.v7.2.entrez.gmt")
kegg_fib_high_low <- GSEA(geneList = geneList, TERM2GENE = kegg, nPermSimple = 1000000, eps = 0)
kegg_fib_high_low_results <- kegg_fib_high_low@result
write.csv(kegg_fib_high_low_results, "./DE/KEGG_fib_high_med.csv")
pathview(gene.data  = geneList, pathway.id = "", species = "hsa", limit = list(gene=max(abs(geneList)), cpd=1))

## FIBROSIS MED - LOW

#First, the patients classified with high level of fibrosis is taken out, so only fibrosis med and fibrosis low are compared.
#Since LI_score will be used as a blocking effect in the model matrix, the LI_score of these specific patients are selected too.
fib_med <- factor(subset(fibrosis, fibrosis!= "fibrosis_high"))
fib_med <- relevel(fib_med, ref = "fibrosis_low")
select_fib <- subset(counts_normalised$samples, counts_normalised$samples$group!= "fibrosis_high")
counts_normalised_fib_med <- counts_normalised[, rownames(select_fib)]
samples_fib_med <- subset(samples, rownames(samples) %in% colnames(counts_normalised_fib_med))
LI_score_fib_med <- samples_fib_med$LI_score

#We create the design matrix that is going to contain all the groups and distribution of the samples within the groups.
design <- model.matrix(~LI_score_fib_med + fib_med, data = counts_normalised_fib_med$samples)
colnames(design) <- gsub(pattern = "LI_score_fib_med", replacement = "", colnames(design))
colnames(design) <- gsub(pattern = "fib_med", replacement = "", colnames(design))
design

#We estimate the dispersion of the normalised counts.
estimated_dispersion <- estimateDisp(counts_normalised_fib_med, design, robust = TRUE)
estimated_dispersion$common.dispersion
sqrt(estimated_dispersion$common.dispersion)

#To visualise the estimated dispersion, we generate a plot for the Biological Coefficient of Variation (BCV). 
pdf(file = "./DE/BCV_plot_fib_med.pdf")
plotBCV(estimated_dispersion)
dev.off()

#Next we fit a quasi-likelihood negative binomial generalized log-linear model to count data and conduct genewise statistical 
#tests for a given coefficient or contrast.
fit <- glmQLFit(estimated_dispersion, design, robust =TRUE)
head(fit$coefficients)
plotQLDisp(fit)
qlf_fib_med <- glmQLFTest(fit)
de <- decideTests(qlf_fib_med) #497 down 5186 up
pdf(file = "./DE/MDPlot_Fibrosis_med-Fibrosis_low.pdf")
plotSmear(qlf_fib_med, de.tags = rownames(qlf_fib_med)[as.logical(de)])
dev.off()

#We perform a diagnostic plot and extract the DEG.
results <- topTags(qlf_fib_med, n= Inf)$table
pdf(file ="./DE/Hist_pvalueFibrosis_med-Fibrosis-low.pdf")
hist(results$PValue, breaks = 0:30/30, col = "grey50", border = "white", main = "P-value histogram", xlab = "p-values")
dev.off()

sum(results$FDR < 0.05, na.rm = TRUE) #5683
results$symbol <- mapIds(org.Hs.eg.db, keys=results$genes, column="SYMBOL", keytype="ENSEMBL")
results$entrez <- mapIds(org.Hs.eg.db, keys= results$genes,column="ENTREZID", keytype="ENSEMBL")
results_FDR <- subset(results, FDR < 0.05)
results_FDR <- results_FDR[order(results_FDR$logFC), ]
resOrderedDF <- as.data.frame(results_FDR)
write.csv(resOrderedDF, file = "./DE/QLF_Fibrosis_med-Fibrosis_low.csv")

#Volcano plot illustrating differential expression between groups fibrosis_low and fibrosis_med as defined by fibrosis.
results$threshold <- as.factor(results$FDR < 0.05)
pdf("./DE/Volcano_plot_Fibrosis_med-Fibrosis_low.pdf")
ggplot(data = results, aes(x =logFC, y= -log10(FDR), colour = threshold)) + 
  geom_point(size = 0.8) + 
  scale_color_manual(values = c("Grey", "Red")) +
  geom_hline(yintercept = 1.30103, linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = c(log2(2), -log2(2)), linetype = "dashed", size = 0.3) +
  xlim(c(-10, 10)) +
  xlab("log2 fold change") + ylab("-log10 FDR") +
  theme_bw() +
  theme(legend.position= "none")
dev.off()

#Pheatmap clustering of the top 50 most DEG between fibrosis_low and fibrosis_med as defined by fibrosis. This is ordered by FDR.
logCPMs_fib_med <- cpm(counts_normalised_fib_med, log =TRUE)
rownames(logCPMs_fib_med) <- counts_normalised_fib_med$genes[,1]
colnames(logCPMs_fib_med) <- paste(colnames(logCPMs_fib_med), fib_med, sep = "-")
o <- order(results_FDR$FDR)
results_fib_med <- results_FDR[o[1:50],]
deg_fib_med <- logCPMs_fib_med[results_fib_med$genes,]
annot_fib_med <- as.data.frame(samples_fib_med[, c(2,4,6,8)])
rownames(deg_fib_med) <- paste(rownames(deg_fib_med), mapIds(org.Hs.eg.db, keys = rownames(deg_fib_med), keytype = "ENSEMBL", column = "SYMBOL"), sep = "-")
rownames(annot_fib_med) <- colnames(deg_fib_med)
pdf("./DE/Pheatmap_Fibrosis-med_Fibrosis-low.pdf", width = 12, height = 14)
pheatmap(deg_fib_med, scale = "row", legend = FALSE, border_color = NA, annotation_color = annotCol, clustering_distance_rows = "correlation", clustering_distance_cols = "euclidean", annotation_col = annot_fib_med, clustering_method = "complete")
dev.off()

#GSEA analysis
#First we eliminate those genes that have no entrez ID or are duplicated. This ranked list is introduced in the GSEA function. 
#The gene sets that are going to be used as input have been obtained from the Broad Institute.
results_FDR <- subset(results_FDR, results_FDR$entrez!="<NA>")
results_FDR <- results_FDR[which(duplicated(results_FDR$entrez) == F), ]

geneList <- results_FDR$logFC
names(geneList) <- as.character(results_FDR$entrez)
geneList <- sort(geneList, decreasing = TRUE)

#Reactome analysis
reactome <- read.gmt("./c2.cp.reactome.v7.2.entrez.gmt")
reactome_fib_med_low <- GSEA(geneList = geneList, TERM2GENE = reactome, nPermSimple = 1000000, eps = 0)
reactome_fib_med_low_results <- reactome_fib_med_low@result
write.csv(reactome_fib_med_low_results, "./DE/Reactome_fib_med_low.csv")
pdf("./DE/Ridgeplot_reactome_fib_med_low.pdf", width = 15, height = 20)
ridgeplot(reactome_fib_med_low, label_format = 20, showCategory = 20)
dev.off()
pdf("./DE/GSEAplot_reactome_fib_med_low_1.pdf")
gseaplot2(reactome_fib_med_low, geneSetID = 1)
dev.off()

#GO analysis
go_bp <- read.gmt("./c5.go.bp.v7.2.entrez.gmt")
go_cc <- read.gmt("./c5.go.cc.v7.2.entrez.gmt")
go_mf <- read.gmt("./c5.go.mf.v7.2.entrez.gmt")
gobp_fib_med_low <- GSEA(geneList = geneList, TERM2GENE = go_bp, nPermSimple = 1000000, eps =0)
gobp_fib_med_low_results <- gobp_fib_med_low@result
write.csv(gobp_fib_med_low_results, "./DE/GObp_fib_med_low.csv")
gocc_fib_med_low <- GSEA(geneList = geneList, TERM2GENE = go_cc, nPermSimple = 1000000, eps = 0)
gocc_fib_med_low_results <- gocc_fib_med_low@result
write.csv(gocc_fib_med_low_results, "./DE/GOcc_fib_med_low.csv") 
gomf_fib_med_low <- GSEA(geneList = geneList, TERM2GENE = go_mf, nPermSimple = 1000000, eps = 0)
gomf_fib_med_low_results <- gomf_fib_med_low@result
write.csv(gomf_fib_med_low_results, "./DE/GOmf_fib_med_low.csv")
pdf("./DE/DotplotGObp_fib_med_low.pdf", width = 12, height = 14)
dotplot(gobp_fib_med_low)
dev.off()
pdf("./DE/DotplotGOcc_fib_med_low.pdf", width = 12, height = 14)
dotplot(gocc_fib_med_low)
dev.off()
pdf("./DE/DotplotGOmf_fib_med_low.pdf", width = 12, height = 14)
dotplot(gomf_fib_med_low)
dev.off()

#KEGG pathways
kegg <- read.gmt("./c2.cp.kegg.v7.2.entrez.gmt")
kegg_fib_med_low <- GSEA(geneList = geneList, TERM2GENE = kegg, nPermSimple = 1000000, eps = 0)
kegg_fib_med_low_results <- kegg_fib_med_low@result
write.csv(kegg_fib_med_low_results, "./DE/KEGG_fib_high_med.csv")
pathview(gene.data  = geneList, pathway.id = "hsa04910", species = "hsa", limit = list(gene=max(abs(geneList)), cpd=1))

## FIBROSIS MED - HIGH

#First, the patients classified with high level of fibrosis is taken out, so only fibrosis med and fibrosis high are compared.
#Since LI_score will be used as a blocking effect in the model matrix, the LI_score of these specific patients are selected too.
fib_high_med <- factor(subset(fibrosis, fibrosis!= "fibrosis_low"))
fib_high_med <- relevel(fib_high_med, ref = "fibrosis_med")
select_fib <- subset(counts_normalised$samples, counts_normalised$samples$group!= "fibrosis_low")
counts_normalised_fib_high_med <- counts_normalised[, rownames(select_fib)]
samples_fib_high_med <- subset(samples, rownames(samples) %in% colnames(counts_normalised_fib_high_med))
LI_score_fib_high_med <- factor(samples_fib_high_med$LI_score)

#We create the design matrix that is going to contain all the groups and distribution of the samples within the groups.
design <- model.matrix(~LI_score_fib_high_med + fib_high_med, data = counts_normalised_fib_high_med$samples)
colnames(design) <- gsub(pattern = "LI_score_fib_high_med", replacement = "", colnames(design))
colnames(design) <- gsub(pattern = "fib_high_med", replacement = "", colnames(design))
design

#We estimate the dispersion of the normalised counts. 
estimated_dispersion <- estimateDisp(counts_normalised_fib_high_med, design, robust = TRUE)
estimated_dispersion$common.dispersion
sqrt(estimated_dispersion$common.dispersion)

#To visualise the estimated dispersion, we generate a plot for the Biological Coefficient of Variation (BCV). 
pdf(file = "./DE/BCV_plot_fib_high_med.pdf")
plotBCV(estimated_dispersion)
dev.off()

#Next we fit a quasi-likelihood negative binomial generalized log-linear model to count data and conduct genewise statistical 
#tests for a given coefficient or contrast.
fit <- glmQLFit(estimated_dispersion, design, robust =TRUE)
head(fit$coefficients)
qlf_fib_high_med <- glmQLFTest(fit)
de <- decideTests(qlf_fib_high_med) #1186 down 1308 up
pdf(file = "./DE/MDPlot_Fibrosis_high-Fibrosis_med.pdf")
plotSmear(qlf_fib_high_med, de.tags = rownames(qlf_fib_high_med)[as.logical(de)])
dev.off()

#We perform a diagnostic plot and extract the DEG.
results <- topTags(qlf_fib_high_med, n= Inf)$table
pdf(file ="./DE/Hist_pvalueFibrosis_high-Fibrosis-med.pdf")
hist(results$PValue, breaks = 0:30/30, col = "grey50", border = "white", main = "P-value histogram", xlab = "p-values")
dev.off()

sum(results$FDR < 0.05, na.rm = TRUE) #2494
results$symbol <- mapIds(org.Hs.eg.db, keys=results$genes, column="SYMBOL", keytype="ENSEMBL")
results$entrez <- mapIds(org.Hs.eg.db, keys= results$genes,column="ENTREZID", keytype="ENSEMBL")

results_FDR <- subset(results, FDR < 0.05)
results_FDR <- results_FDR[order(results_FDR$logFC), ]
resOrderedDF <- as.data.frame(results_FDR)
write.csv(resOrderedDF, file = "./DE/QLF_Fibrosis_high-Fibrosis_med.csv")

#Volcano plot illustrating differential expression between groups fibrosis_high and fibrosis_med as defined by fibrosis.
results$threshold <- as.factor(results$FDR < 0.05)
pdf("./DE/Volcano_plot_Fibrosis_high-Fibrosis_med.pdf")
ggplot(data = results, aes(x =logFC, y= -log10(FDR), colour = threshold)) + 
  geom_point(size = 0.8) + 
  scale_color_manual(values = c("Grey", "Red")) +
  geom_hline(yintercept = 1.30103, linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = c(log2(2), -log2(2)), linetype = "dashed", size = 0.3) +
  xlim(c(-5, 5)) +
  xlab("log2 fold change") + ylab("-log10 FDR") +
  theme_bw() +
  theme(legend.position= "none")
dev.off()

#Pheatmap clustering of the top 50 most DEG between fibrosis_high and fibrosis_med as defined by fibrosis. This is ordered by FDR.
logCPMs_fib_high_med <- cpm(counts_normalised_fib_high_med, log =TRUE)
rownames(logCPMs_fib_high_med) <- counts_normalised_fib_high_med$genes[,1]
colnames(logCPMs_fib_high_med) <- paste(colnames(logCPMs_fib_high_med), fib_high_med, sep = "-")
o <- order(results_FDR$FDR)
results_fib_high_med <- results_FDR[o[1:50],]
deg_fib_high_med <- logCPMs_fib_high_med[results_fib_high_med$genes,]
annot_fib_high_med <- as.data.frame(samples_fib_high_med[, c(2,4,6,8)])
rownames(deg_fib_high_med) <- paste(rownames(deg_fib_high_med), mapIds(org.Hs.eg.db, keys = rownames(deg_fib_high_med), keytype = "ENSEMBL", column = "SYMBOL"), sep = "-")
rownames(annot_fib_high_med) <- colnames(deg_fib_high_med)
pdf("./DE/Pheatmap_Fibrosis-high_Fibrosis-med.pdf", width = 12, height = 14)
pheatmap(deg_fib_high_med, scale = "row", legend = FALSE, border_color = NA, clustering_distance_rows = "correlation", clustering_distance_cols = "euclidean", annotation_col = annot_fib_high_med, annotation_color = annotCol, clustering_method = "complete")
dev.off()

#GSEA analysis
#First we eliminate those genes that have no entrez ID or are duplicated. This ranked list is introduced in the GSEA function. 
#The gene sets that are going to be used as input have been obtained from the Broad Institute.
results_FDR <- subset(results_FDR, results_FDR$entrez!="<NA>")
results_FDR <- results_FDR[which(duplicated(results_FDR$entrez) == F), ]

geneList <- results_FDR$logFC
names(geneList) <- as.character(results_FDR$entrez)
geneList <- sort(geneList, decreasing = TRUE)

#Reactome analysis
reactome <- read.gmt("./c2.cp.reactome.v7.2.entrez.gmt")
reactome_fib_high_med <- GSEA(geneList = geneList, TERM2GENE = reactome, nPermSimple = 1000000, eps = 0)
reactome_fib_high_med_results <- reactome_fib_high_med@result
write.csv(reactome_fib_high_med_results, "./DE/Reactome_fib_high_med.csv")
pdf("./DE/Ridgeplot_reactome_fib_high_med.pdf", width = 15, height = 20)
ridgeplot(reactome_fib_high_med, label_format = 20, showCategory = 20)
dev.off()
pdf("./DE/GSEAplot_reactome_fib_high_med_1.pdf")
gseaplot2(reactome_fib_high_med, geneSetID = 1)
dev.off()

#GO analysis
go_bp <- read.gmt("./c5.go.bp.v7.2.entrez.gmt")
go_cc <- read.gmt("./c5.go.cc.v7.2.entrez.gmt")
go_mf <- read.gmt("./c5.go.mf.v7.2.entrez.gmt")
gobp_fib_high_med <- GSEA(geneList = geneList, TERM2GENE = go_bp, nPermSimple = 1000000, eps =0)
gobp_fib_high_med_results <- gobp_fib_high_med@result
write.csv(gobp_fib_high_med_results, "./DE/GObp_fib_high_med.csv")
gocc_fib_high_med <- GSEA(geneList = geneList, TERM2GENE = go_cc, nPermSimple = 1000000, eps = 0)
gocc_fib_high_med_results <- gocc_fib_high_med@result
write.csv(gocc_fib_high_med_results, "./DE/GOcc_fib_high_med.csv") 
gomf_fib_high_med <- GSEA(geneList = geneList, TERM2GENE = go_mf, nPermSimple = 1000000, eps = 0)
gomf_fib_high_med_results <- gomf_fib_high_med@result
write.csv(gomf_fib_high_med_results, "./DE/GOmf_fib_high_med.csv")
pdf("./DE/DotplotGObp_fib_high_med.pdf", width = 12, height = 14)
dotplot(gobp_fib_high_med)
dev.off()
pdf("./DE/DotplotGOcc_fib_high_med.pdf", width = 12, height = 14)
dotplot(gocc_fib_high_med)
dev.off()
pdf("./DE/DotplotGOmf_fib_high_med.pdf", width = 12, height = 14)
dotplot(gomf_fib_high_med)
dev.off()

#KEGG pathways
kegg <- read.gmt("./c2.cp.kegg.v7.2.entrez.gmt")
kegg_fib_high_med <- GSEA(geneList = geneList, TERM2GENE = kegg, nPermSimple = 1000000, eps = 0)
kegg_fib_high_med_results <- kegg_fib_high_med@result
write.csv(kegg_fib_high_med_results, "./DE/KEGG_fib_high_med.csv")
pathview(gene.data  = geneList, pathway.id = "hsa00260", species = "hsa", limit = list(gene=max(abs(geneList)), cpd=1))

## STEATOSIS

#Steatosis levels (high and low) are compared.
#Since rs738409 variant will be used as a blocking effect in the model matrix, this is selected  for the patients too.
rs738409 <- factor(c("M", "M", "N", "N", "M", "N", "N", "N", "M", "N", "N", "N", "M", "N", "M", "M", "N", "M", "M", "N", "M", "N", "M", "M", "M", "N", "M", "M", "M", "N", "N", "M", "N", "N", "M", "N", "N", "N", "N", "N", "N", "N", "N", "N", "M", "M", "M", "N", "N", "M", "M", "N", "N", "N"))
steatosis_deg <- relevel(steatosis, ref = "steatosis_low")
design <- model.matrix(~rs738409 + steatosis_deg, data = counts_normalised$samples)
colnames(design) <- gsub(pattern = "N", replacement = "", colnames(design))
colnames(design) <- gsub(pattern = "steatosis_deg", replacement = "", colnames(design))
design

#We estimate the dispersion of the normalised counts.
estimated_dispersion <- estimateDisp(counts_normalised, design, robust = TRUE)
estimated_dispersion$common.dispersion
sqrt(estimated_dispersion$common.dispersion)

#To visualise the estimated dispersion, we generate a plot for the Biological Coefficient of Variation (BCV). 
pdf(file = "./DE/BCV_plot_fib_steatosis.pdf")
plotBCV(estimated_dispersion)
dev.off()

#Next we fit a quasi-likelihood negative binomial generalized log-linear model to count data and conduct genewise statistical 
#tests for a given coefficient or contrast.
fit <- glmQLFit(estimated_dispersion, design, robust =TRUE)
head(fit$coefficients)
qlf_fib_steatosis <- glmQLFTest(fit)
de <- decideTests(qlf_fib_steatosis) #1932 down 687 up
pdf(file = "./DE/MDPlot_Fibrosis_steatosis.pdf")
plotSmear(qlf_fib_steatosis, de.tags = rownames(qlf_fib_high_med)[as.logical(de)])
dev.off()

#We perform a diagnostic plot and extract the DEG.
results <- topTags(qlf_fib_steatosis, n= Inf)$table
pdf(file ="./DE/Hist_pvalueSteatosis.pdf")
hist(results$PValue, breaks = 0:30/30, col = "grey50", border = "white", main = "P-value histogram", xlab = "p-values")
dev.off()

sum(results$FDR < 0.05, na.rm = TRUE) #2619
results$symbol <- mapIds(org.Hs.eg.db, keys=results$genes, column="SYMBOL", keytype="ENSEMBL")
results$entrez <- mapIds(org.Hs.eg.db, keys= results$genes,column="ENTREZID", keytype="ENSEMBL")

results_FDR <- subset(results, FDR < 0.05)
results_FDR <- results_FDR[order(results_FDR$logFC), ]
resOrderedDF <- as.data.frame(results_FDR)
write.csv(resOrderedDF, file = "./DE/QLF_Steatosis.csv")

#Volcano plot illustrating differential expression between groups steatosis_low and steatosis_high as defined by steatosis.
results$threshold <- as.factor(results$FDR < 0.05)
pdf("./DE/Volcano_plot_Steatosis.pdf")
ggplot(data = results, aes(x =logFC, y= -log10(FDR), colour = threshold)) + 
  geom_point(size = 0.8) + 
  scale_color_manual(values = c("Grey", "Red")) +
  geom_hline(yintercept = 1.30103, linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = c(log2(2), -log2(2)), linetype = "dashed", size = 0.3) +
  xlim(c(-9, 9)) +
  xlab("log2 fold change") + ylab("-log10 FDR") +
  theme_bw() +
  theme(legend.position= "none")
dev.off()

#Pheatmap clustering of the top 50 most DEG between steatosis_low and steatosis_high as defined by steatosis. This is ordered by FDR.
logCPMs_steatosis <- logCPMs
colnames(logCPMs_steatosis) <- paste(colnames(logCPMs), steatosis, sep = "-")
o <- order(results_FDR$FDR)
results_steatosis <- results_FDR[o[1:50],]
deg_steatosis <- logCPMs_steatosis[results_steatosis$genes,]
annot_steatosis <- as.data.frame(samples[, c(2,4,6,8)])
rownames(deg_steatosis) <- paste(rownames(deg_steatosis), mapIds(org.Hs.eg.db, keys = rownames(deg_steatosis), keytype = "ENSEMBL", column = "SYMBOL"), sep = "-")
rownames(annot_steatosis) <- colnames(deg_steatosis)
pdf("./DE/Pheatmap_Steatosis.pdf", width = 12, height = 14)
pheatmap(deg_steatosis, scale = "row", legend = FALSE, border_color = NA, clustering_distance_rows = "correlation", clustering_distance_cols = "euclidean", annotation_col = annot_steatosis, annotation_colors = annotCol, clustering_method = "complete")
dev.off()

#GSEA analysis
#First we eliminate those genes that have no entrez ID or are duplicated. This ranked list is introduced in the GSEA function. 
#The gene sets that are going to be used as input have been obtained from the Broad Institute.
results_FDR <- subset(results_FDR, results$entrez!="<NA>")
results_FDR <- results_FDR[which(duplicated(results_FDR$entrez) == F), ]

geneList <- results_FDR$logFC
names(geneList) <- as.character(results_FDR$entrez)
geneList <- sort(geneList, decreasing = TRUE)

#Reactome analysis
reactome <- read.gmt("./c2.cp.reactome.v7.2.entrez.gmt")
reactome_steatosis <- GSEA(geneList = geneList, TERM2GENE = reactome, nPermSimple = 1000000, eps = 0)
reactome_steatosis_results <- reactome_steatosis@result
write.csv(reactome_steatosis_results, "./DE/Reactome_steatosis.csv")
pdf("./DE/Ridgeplot_reactome_steatosis.pdf", width = 12, height = 16)
ridgeplot(reactome_steatosis, label_format = 20, showCategory = 20)
dev.off()
pdf("./DE/GSEAplot_reactome_steatosis_1.pdf")
gseaplot2(reactome_steatosis, geneSetID = 1)
dev.off()

#GO analysis
go_bp <- read.gmt("./c5.go.bp.v7.2.entrez.gmt")
go_cc <- read.gmt("./c5.go.cc.v7.2.entrez.gmt")
go_mf <- read.gmt("./c5.go.mf.v7.2.entrez.gmt")
gobp_steatosis <- GSEA(geneList = geneList, TERM2GENE = go_bp, nPermSimple = 1000000, eps =0)
gobp_steatosis_results <- gobp_steatosis@result
write.csv(gobp_steatosis_results, "./DE/GObp_steatosis.csv")
gocc_steatosis <- GSEA(geneList = geneList, TERM2GENE = go_cc, nPermSimple = 1000000, eps = 0)
gocc_steatosis_results <- gocc_steatosis@result
write.csv(gocc_steatosis_results, "./DE(GOcc_steatosis.csv") 
gomf_steatosis <- GSEA(geneList = geneList, TERM2GENE = go_mf, nPermSimple = 1000000, eps = 0)
gomf_steatosis_results <- gomf_steatosis@result
write.csv(gomf_steatosis_results, "./DE/GOmf_steatosis.csv")
pdf("./DE/DotplotGObp_steatosis.pdf", width = 12, height = 14)
dotplot(gobp_steatosis)
dev.off()
pdf("./DE/DotplotGOcc_steatosis.pdf", width = 12, height = 14)
dotplot(gocc_steatosis)
dev.off()
pdf("./DE/DotplotGOmf_steatosis.pdf", width = 12, height = 14)
dotplot(gomf_steatosis)
dev.off()

#KEGG pathways
kegg <- read.gmt("./c2.cp.kegg.v7.2.entrez.gmt")
kegg_steatosis <- GSEA(geneList = geneList, TERM2GENE = kegg, nPermSimple = 1000000, eps = 0)
kegg_steatosis_results <- kegg_steatosis@result
write.csv(kegg_steatosis_results, "./DE/KEGG_steatosis.csv")
pathview(gene.data  = geneList, pathway.id = "hsa00071", species = "hsa", limit = list(gene=max(abs(geneList)), cpd=1))

## LOBULAR INFLAMMATION

#Lobular inflammation levels (high and low) are compared.
#rs738409 variant and fibrosis score will be used as a blocking effect in the model matrix.
LI_deg <- relevel(LI, ref = "LI_low")
design <- model.matrix(~ fibrosis_score + rs738409 + LI_deg, data = counts_normalised$samples)
colnames(design) <- gsub(pattern = "fibrosis_score", replacement = "", colnames(design))
colnames(design) <- gsub(pattern = "N", replacement = "", colnames(design))
colnames(design) <- gsub(pattern = "LI_deg", replacement = "", colnames(design))
design

#We estimate the dispersion of the normalised counts.
estimated_dispersion <- estimateDisp(counts_normalised, design, robust = TRUE)
estimated_dispersion$common.dispersion
sqrt(estimated_dispersion$common.dispersion)

#To visualise the estimated dispersion, we generate a plot for the Biological Coefficient of Variation (BCV). 
pdf(file = "./DE/BCV_plot_fib_LI.pdf")
plotBCV(estimated_dispersion)
dev.off()

#Next we fit a quasi-likelihood negative binomial generalized log-linear model to count data and conduct genewise statistical 
#tests for a given coefficient or contrast.
fit <- glmQLFit(estimated_dispersion, design, robust =TRUE)
head(fit$coefficients)
qlf_LI <- glmQLFTest(fit)
de <- decideTests(qlf_LI) #205 down 4 up
pdf(file = "./DE/SmearPlot_LI.pdf")
plotSmear(qlf_LI, de.tags = rownames(qlf_LI)[as.logical(de)])
dev.off()

#We perform a diagnostic plot and extract the DEG.
results <- topTags(qlf_LI, n= Inf)$table
pdf(file ="./DE/Hist_pvalueLI.pdf")
hist(results$PValue, breaks = 0:30/30, col = "grey50", border = "white", main = "P-value histogram", xlab = "p-values")
dev.off()

sum(results$FDR < 0.05, na.rm = TRUE) #209
results$symbol <- mapIds(org.Hs.eg.db, keys=results$genes, column="SYMBOL", keytype="ENSEMBL")
results$entrez <- mapIds(org.Hs.eg.db, keys= results$genes,column="ENTREZID", keytype="ENSEMBL")

results_FDR <- subset(results, FDR < 0.05)
results_FDR <- results_FDR[order(results_FDR$logFC), ]
resOrderedDF <- as.data.frame(results_FDR)
write.csv(resOrderedDF, file = "./DE/QLF_LI.csv")

#Volcano plot illustrating differential expression between groups.
results$threshold <- as.factor(results$FDR < 0.05)
pdf("./DE/Volcano_plot_LI.pdf")
ggplot(data = results, aes(x =logFC, y= -log10(FDR), colour = threshold)) + 
  geom_point(size = 0.8) + 
  scale_color_manual(values = c("Grey", "Red")) +
  geom_hline(yintercept = 1.30103, linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = c(log2(2), -log2(2)), linetype = "dashed", size = 0.3) +
  xlim(c(-6, 6)) +
  xlab("log2 fold change") + ylab("-log10 FDR") +
  theme_bw() +
  theme(legend.position= "none")
dev.off()

#Pheatmap clustering of the top 50 most DEG between LI_low and LI_high as defined by LI. This is ordered by FDR.
logCPMs_LI <- logCPMs
colnames(logCPMs_LI) <- paste(colnames(logCPMs), LI, sep = "-")
o <- order(results_FDR$FDR)
results_LI <- results_FDR[o[1:50],]
deg_LI <- logCPMs_LI[results_LI$genes,]
annot_LI <- as.data.frame(samples[, c(2,4,6,8)])
rownames(deg_LI) <- paste(rownames(deg_LI), mapIds(org.Hs.eg.db, keys = rownames(deg_LI), keytype = "ENSEMBL", column = "SYMBOL"), sep = "-")
rownames(annot_LI) <- colnames(deg_LI)
pdf("./DE/Pheatmap_LI.pdf", width = 12, height = 14)
pheatmap(deg_LI, scale = "row", legend = FALSE, border_color = NA, clustering_distance_rows = "correlation", clustering_distance_cols = "euclidean", annotation_col = annot_LI, annotation_colors = annotCol, clustering_method = "complete")
dev.off()

#GSEA analysis
#First we eliminate those genes that have no entrez ID or are duplicated. This ranked list is introduced in the GSEA function. 
#The gene sets that are going to be used as input have been obtained from the Broad Institute.
results_FDR <- subset(results_FDR, results$entrez!="<NA>")
results_FDR <- results_FDR[which(duplicated(results_FDR$entrez) == F), ]

geneList <- results_FDR$logFC
names(geneList) <- as.character(results_FDR$entrez)
geneList <- sort(geneList, decreasing = TRUE)

#Reactome analysis
reactome <- read.gmt("./c2.cp.reactome.v7.2.entrez.gmt")
reactome_LI <- GSEA(geneList = geneList, TERM2GENE = reactome, nPermSimple = 1000000, eps = 0)
reactome_LI_results <- reactome_LI@result
write.csv(reactome_LI_results, "./DE/Reactome_LI.csv")
pdf("./DE/Ridgeplot_reactome_LI.pdf", width = 12, height = 16)
ridgeplot(reactome_LI, label_format = 20, showCategory = 20)
dev.off()
pdf("./DE/GSEAplot_reactome_LIs_1.pdf")
gseaplot2(reactome_LI, geneSetID = 1)
dev.off()

#GO analysis
go_bp <- read.gmt("./c5.go.bp.v7.2.entrez.gmt")
go_cc <- read.gmt("./c5.go.cc.v7.2.entrez.gmt")
go_mf <- read.gmt("./c5.go.mf.v7.2.entrez.gmt")
gobp_LI <- GSEA(geneList = geneList, TERM2GENE = go_bp, nPermSimple = 1000000, eps =0)
gobp_LI_results <- gobp_LI@result
write.csv(gobp_LI_results, "./DE/GObp_LI.csv")
gocc_LI <- GSEA(geneList = geneList, TERM2GENE = go_cc, nPermSimple = 1000000, eps = 0)
gocc_LI_results <- gocc_LI@result
write.csv(gocc_LI_results, "./DE(GOcc_LI.csv") 
gomf_LI <- GSEA(geneList = geneList, TERM2GENE = go_mf, nPermSimple = 1000000, eps = 0)
gomf_LI_results <- gomf_LI@result
write.csv(gomf_LI_results, "./DE/GOmf_LI.csv")
pdf("./DE/DotplotGObp_LI.pdf", width = 12, height = 14)
dotplot(gobp_LI)
dev.off()
pdf("./DE/DotplotGOcc_LI.pdf", width = 12, height = 14)
dotplot(gocc_LI)
dev.off()
pdf("./DE/DotplotGOmf_LI.pdf", width = 12, height = 14)
dotplot(gomf_LI)
dev.off()

#KEGG pathways
kegg <- read.gmt("./c2.cp.kegg.v7.2.entrez.gmt")
kegg_LI <- GSEA(geneList = geneList, TERM2GENE = kegg, nPermSimple = 1000000, eps = 0)

