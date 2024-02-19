library(PCAtools)
library(ggcorrplot)
library(rgl)
library(dplyr)
library(RColorBrewer)

load("./SteatoSITE/DGE_and_GSEA.R")
#Upload final csv file with counts with genes as rows and patient ID as columns called feature_counts
#Upload final csv file with pathological data with patient ID as rows and clinical information as columns called path_metadata

#CORRELATIONS BETWEEN VARIABLES
##Create copy of the path metadata to transform some of the variables into numeric to perform correlations later.
path_copy <- path_metadata
path_copy <- as.data.frame(lapply(path_copy, as.factor))
path_copy <- as.data.frame(lapply(path_copy, as.numeric))

##Change colnames to a presentable format
colnames(path_copy) <- c("sample_id", "nash_diagnosis", "Hepatocyte ballooning", "Steatosis", "Lobular inflammation", "Modified Ishak score", "Kleiner fibrosis sep", "Age", "Sex", "Status", "Diabetes","BMI", "Kleiner fibrosis", "NAFLD", "NAS", "NAS scores", "Mapped reads", "Dup PCR rate per gene", "% mapped reads", "rs62305723",  "rs738409", "rRNA", "Specimen type", "Specimen date", "Health board", "Sections cut", "Curls cut", "Batch no.")

##Obtain correlation matrix with pvalues for the different variables.
correlation_matrix <- cor(path_copy[, c(3:5, 8:9, 11:13, 16:ncol(path_copy))], use = "complete.obs") #NA values are not taken into account when performing correlation.
cor_pmat <- cor_pmat(path_copy[, c(3:5, 8:9, 11:13, 16:ncol(path_copy))])
#Plot correlations with p-values.
ggcorrplot(correlation_matrix, type = 'lower', lab = TRUE, p.mat = cor_pmat, insig = 'blank', method = "circle", lab_size = 3, outline.color = "white")

#PCA. Normalised log2 CPM counts are used for this. 
PCA <- pca(logcpm_norm, metadata = dgList$samples)

screeplot <- screeplot(PCA) #Plot screeplot to observe the variances against the number of PCs.

##Plot biplots for the first three PCs. It saves the plots in PDF version
for(i in names(PCA$metadata[6:ncol(PCA$metadata)])){
  pdf(paste("./Biplot_PC1_PC2_", i, ".pdf", sep = ""), width = 10, height = 10)
  print(biplot(pcaobj = PCA, x = "PC1", y = "PC2", colby = i, title= "", legendPosition = "bottom", lab = NULL))
  dev.off()
}

for(i in names(PCA$metadata[6:ncol(PCA$metadata)])){
  pdf(paste("./Biplot_PC2_PC3_", i, ".pdf", sep = ""), width = 10, height = 10)
  print(biplot(pcaobj = PCA, x = "PC2", y = "PC3", colby = i, title= "", legendPosition = "bottom", lab = NULL))
  dev.off()
}

##Plot component loadings to observe those variables that drives the variation in the chosen PCs.
pca_loadings <- plotloadings(PCA, components = getComponents(PCA, c(2)),
                             rangeRetain = 0.10, absolute = TRUE,
                             col = c('black', 'pink', 'red4'), labSize = 4, drawConnectors = TRUE) + coord_flip()

##Create correlation plot 
PCA$metadata <- as.data.frame(lapply(PCA$metadata, as.factor))
PCA$metadata <- as.data.frame(lapply(PCA$metadata, as.numeric))
names(PCA$metadata) <- c("group", "lib.sizes", "norm.factors", "sample_id", "NASH diagnosis", "Hepatocyte ballooning", "Steatosis", "Lobular inflammation", "Modified Ishak Score", "Kleiner fibrosis sep", "Age", "Sex", "Status", "Diabetes", "BMI", "Kleiner fibrosis", "NAFLD", "NAS", "NAS_scores", "Mapped reads", "Dup PCR rate per gene", "% mapped reads", "rs62305723", "rs738409", "rRNA", "Specimen type", "Specimen date", "Health board", "Sections cut", "Curls cut", "Year",  "Batch no.")
PCA$metadata <- PCA$metadata %>% dyplr::select(c(1:10, 16:19, 11:15, 20:32)) #Re-order names of metadata
#Correlate PCs to variables
eigencorplot(PCA, components = getComponents(PCA, 1:15), metavars = c("NASH diagnosis", "Hepatocyte ballooning", "Steatosis", "Lobular inflammation", "Modified Ishak Score", "Kleiner fibrosis", "Age", "Sex", "Status", "Diabetes", "BMI", "Health board", "Specimen type", "Sections cut", "Curls cut", "Batch no.", "rs738409", "rs62305723", "Mapped reads", "% mapped reads", "Dup PCR rate per gene", "rRNA", "Year") , col= c("#2171b5", "#6baed6", "#bdd7e7", "#eff3ff", "white", "#fee5d9", "#fcae91", "#fb6a4a", "#cb181d"), corFUN = 'spearman',
             corUSE = 'pairwise.complete.obs',
             corMultipleTestCorrection = 'BH',
             signifSymbols = c('****', '***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
             rotLabX = 45)
##R-squared values
r_squared <- eigencorplot(PCA, components = getComponents(PCA, 1:15), metavars = c("Hepatocyte ballooning", "Steatosis", "Lobular inflammation", "Kleiner fibrosis", "Mapped reads", "Dup PCR rate per gene", "% mapped reads", "rs62305723", "rs738409", "rRNA", "Specimen type", "Health board", "Sections cut", "Curls cut", "Batch no.", "Sex", "Year") , col= c("white", "#bdd7e7", "#6baed6", "#2171b5"), corFUN = 'spearman',
                          corUSE = 'pairwise.complete.obs',
                          corMultipleTestCorrection = 'BH',
                          signifSymbols = c('****', '***', '**', '*', ''),
                          signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                          rotLabX = 45, plotRsquared = TRUE)

##Plot first 3 PCs in 3D
get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}

brewer.pal(n = 5, name = "Paired") 
cols <- get_colors(as.numeric(PCA$metadata$`Kleiner fibrosis`), c("#A6CEE3", "#B2DF8A", "#E31A1C", "#FF7F00", "#CAB2D6"))
cols2 <- get_colors(as.numeric(as.factor(dgList$samples$specimen_type)), c("#A6CEE3", "#B2DF8A", "#E31A1C","#FF7F00","#CAB2D6"))
plot3d(PCA$rotated, col = cols, size = 6, type = "p")
legend3d("topright", legend = c("Biopsy", "Explant", "Resection"), col = c("#A6CEE3", "#B2DF8A", "#E31A1C", "#FF7F00", "#CAB2D6"), pch = 10)
grid3d(c("x", "y", "z"))

#Calculate distance centroids from different groups in PCA. In this case, explants between F4 and rest of samples were studied.
pca_rotated <- PCA$rotated
pca_rotated$groups <- dgList$samples$NAFLD
pca_centroids <- aggregate(pca_rotated[,1:(ncol(pca_rotated)-1)], list(Type = pca_rotated$groups), mean)
dist(rbind(pca_centroids[pca_centroids$Type == "F4_explant", 2:3], pca_centroids[pca_centroids$Type == "F4", 2:3]), method = "euclidean")
dist(rbind(pca_centroids[pca_centroids$Type == "F4_explant", 2:3], pca_centroids[pca_centroids$Type == "Control", 2:3]), method = "euclidean")
dist(rbind(pca_centroids[pca_centroids$Type == "F4", 2:3], pca_centroids[pca_centroids$Type == "Control", 2:3]), method = "euclidean")
