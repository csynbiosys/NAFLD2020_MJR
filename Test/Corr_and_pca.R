library(PCAtools)
library(ggcorrplot)
library(rgl)

load("./SteatoSITE/DGE_and_GSEA.R")
feature_counts <- read.csv("./Data/Clean_counts.csv", stringsAsFactors = FALSE, row.names = 1)
path_metadata <- read.csv("./Data/Path_metadata.csv", stringsAsFactors = FALSE, row.names =1, header = TRUE)


#CORRELATIONS BETWEEN VARIABLES
##Create copy of the path metadata to transform some of the variables into numeric to perform later the correlation.
path_copy <- path_metadata
path_copy <- lapply(path_copy, as.factor)
path_copy <- lapply(path_copy, as.numeric)

##Change colnames to a presentable format
colnames(path_copy) <- c("sample_id", "nash_diagnosis", "Hepatocyte ballooning", "Steatosis", "Lobular inflammation", "Modified Ishak score", "Kleiner fibrosis sep", "Kleiner fibrosis", "Mapped reads", "Dup PCR rate per gene", "% mapped reads", "rs62305723", "rs738409", "rRNA", "Specimen type", "Specimen date", "Health board", "Sections cut", "Curls cut", "Batch no", "Sex", "NAS score", "Steatosis level", "Inflammation level", "Ballooning level", "Fibrosis level", "nafld", "nash_fibrosis", "fib", "nas_groups", "Year", "month" )

##Obtain correlation matrix with pvaluesfor the different variables.
correlation_matrix <- cor(path_copy[, c("Kleiner fibrosis", 'Lobular inflammation', 'Steatosis', 'Hepatocyte ballooning', 'rRNA', 'rs738409', 'rs62305723', 'Mapped reads', 'Dup PCR rate per gene', '% mapped reads','Specimen type', 'Health board', 'Sections cut','Curls cut', 'Batch no', 'Sex', "Year")], method = 'spearman')
cor_pmat <- cor_pmat(path_copy[, c("Kleiner fibrosis", 'Lobular inflammation', 'Steatosis', 'Hepatocyte ballooning', 'rRNA', 'rs738409', 'rs62305723', 'Mapped reads', 'Dup PCR rate per gene', '% mapped reads','Specimen type', 'Health board', 'Sections cut','Curls cut', 'Batch no', 'Sex', "Year")], method = 'spearman')
#Plot correlations with p-values.
ggcorrplot(correlation_matrix, type = 'lower', lab = TRUE, p.mat = cor_pmat, insig = 'blank', method = "circle", lab_size = 3, outline.color = "white")

#PCA 
PCA <- pca(logcpm_norm, metadata = dgList$samples)

screeplot <- screeplot(PCA) #Plot screeplot to observe the variances against the number of PCs.

##Plot biplots for the first three PCs
for(i in names(PCA$metadata[6:ncol(PCA$metadata)])){
  pdf(paste("./Results/Biplot_PC1_PC2_", i, ".pdf", sep = ""), width = 10, height = 10)
  print(biplot(pcaobj = PCA, x = "PC1", y = "PC2", colby = i, title= "", legendPosition = "bottom", lab = NULL))
  dev.off()
}
for(i in names(PCA$metadata[6:ncol(PCA$metadata)])){
  pdf(paste("./Results/Biplot_PC2_PC3_", i, ".pdf", sep = ""), width = 10, height = 10)
  print(biplot(pcaobj = PCA, x = "PC2", y = "PC3", colby = i, title= "", legendPosition = "bottom", lab = NULL))
  dev.off()
}

##Plot component loadings to observe those variables that drives the variation in the chosen PCs.
pca_loadings <- plotloadings(PCA, components = getComponents(PCA, c(2)),
                             rangeRetain = 0.10, absolute = TRUE,
                             col = c('black', 'pink', 'red4'), labSize = 4, drawConnectors = TRUE) + coord_flip()

##Create correlation plot 
PCA$metadata <- lapply(PCA$metadata, as.numeric)
names(PCA$metadata) <- c("group", "lib.size", "norm.factors", "sample_id", "NASH diagnosis", "Hepatocyte ballooning", "Steatosis", "Lobular inflammation", "Modified Ishak score", "Kleiner fibrosis sep", "Kleiner fibrosis", "Mapped reads", "Dup PCR rate per gene", "% mapped reads", "rs62305723", "rs738409", "rRNA", "Specimen type", "Specimen date", "Health board", "Sections cut", "Curls cut", "Batch no", "Sex", "NAS score", "Steatosis level", "Inflammation level", "Ballooning level", "Fibrosis level", "NAFLD", "NASH fibrosis", "Fib", "NAS groups", "Year", "Month") 
#Correlate PCs to variables
eigencorplot(PCA, components = getComponents(PCA, 1:15), metavars = c("Hepatocyte ballooning", "Steatosis", "Lobular inflammation", "Kleiner fibrosis", "Mapped reads", "Dup PCR rate per gene", "% mapped reads", "rs62305723", "rs738409", "rRNA", "Specimen type", "Health board", "Sections cut", "Curls cut", "Batch no", "Sex", "Year") , col= c("#2171b5", "#6baed6", "#bdd7e7", "#eff3ff", "white", "#fee5d9", "#fcae91", "#fb6a4a", "#cb181d"), corFUN = 'spearman',
             corUSE = 'pairwise.complete.obs',
             corMultipleTestCorrection = 'BH',
             signifSymbols = c('****', '***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
             rotLabX = 45)
##R-squared values
r_squared <- eigencorplot(PCA, components = getComponents(PCA, 1:15), metavars = c("Hepatocyte ballooning", "Steatosis", "Lobular inflammation", "Kleiner fibrosis", "Mapped reads", "Dup PCR rate per gene", "% mapped reads", "rs62305723", "rs738409", "rRNA", "Specimen type", "Health board", "Sections cut", "Curls cut", "Batch no", "Sex", "Year") , col= c("white", "#bdd7e7", "#6baed6", "#2171b5"), corFUN = 'spearman',
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
cols <- get_colors(as.numeric(fib), c("#A6CEE3", "#B2DF8A", "#E31A1C", "#FF7F00", "#CAB2D6"))
cols2 <- get_colors(as.numeric(as.factor(dgList$samples$specimen_type)), c("#A6CEE3", "#B2DF8A", "#E31A1C","#FF7F00","#CAB2D6"))
plot3d(PCA$rotated, col = cols2, size = 6, type = "p")
legend3d("topright", legend = c("Biopsy", "Explant", "Resection"), col = c("#A6CEE3", "#B2DF8A", "#E31A1C", "#FF7F00", "#CAB2D6"), pch = 10)
grid3d(c("x", "y", "z"))

#ANOVA between PCs and experiment
pc_experiment <- function(pca, pcameta, last_pc){

  # Create an empty matrix to store ANOVA p values associating components with experimental factors

pvals <- matrix(
  data = NA, 
  ncol = ncol(pcameta), 
  nrow = last_pc, 
  dimnames = list(as.character(1:last_pc), colnames(pcameta))
)

  # Fill the matrix with anova p values

for (i in 1:ncol(pcameta)) {
  for (j in 1:last_pc) {
    fit <- aov(pca$rotated[, j] ~ factor(t(pcameta[, i])))
    if ("Pr(>F)" %in% names(summary(fit)[[1]])) {
      pvals[j, i] <- summary(fit)[[1]][["Pr(>F)"]][[1]]
    }
  }
}

  # Calculate the percent variance explained by each component

fraction_explained <- round((pca$sdev)^2/sum(pca$sdev^2), 3) * 100
names(fraction_explained) <- colnames(pca$rotated)
rownames(pvals) <- paste(paste("PC", 1:last_pc, sep = ""), " (", fraction_explained[1:last_pc],  "%)", sep = "" )

pvals
}

  # Factor plot to illustrate association of components with categorical variables

factor_plot <- function(pca, pvals){
  
  pvals <- pvals[rev(rownames(pvals)), , drop = FALSE]
  pvals <- -log10(pvals)
  # Make a heatmap of the p values
  
  ggplot(data = reshape2::melt(pvals), aes(Var2, Var1, fill = value)) +
    ylab("Principal components") +
    xlab("Experimental factor") +
    geom_tile(color = "white") +
    scale_fill_gradientn(
      colors = c('white', 'Salmon', 'FireBrick'),
      values = c(0,0.1,1),
      space = "Lab",
      name = "ANOVA\n -log10(p value)"
    ) +
    geom_text(aes(label = signif(value, 3)), color = "white", size = 4) +
    theme_minimal() +
    theme(axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      size = 12,
      hjust = 1
    )) +
    coord_equal(ratio=.3) +
    theme(
      axis.text.x=element_text(size=14),
      axis.text.y=element_text(size=14),
      axis.title.x=element_text(size=16),
      axis.title.y=element_text(size=16),
      legend.text = element_text(size=12),
      legend.title = element_text(size=14)
    )
}

pca_anova <- pc_experiment(PCA, PCA$metadata[, c("Hepatocyte ballooning", "Steatosis", "Lobular inflammation", "Kleiner fibrosis", "Mapped reads", "Dup PCR rate per gene", "% mapped reads", "rs62305723", "rs738409", "rRNA", "Specimen type", "Health board", "Sections cut", "Curls cut", "Batch no", "Sex", "Year") , col= c("#2171b5", "#6baed6", "#bdd7e7", "#eff3ff", "white", "#fee5d9", "#fcae91", "#fb6a4a", "#cb181d")])
fact_plot <- factor_plot(PCA, pca_anova)
