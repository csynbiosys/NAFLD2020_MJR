library(ggplot2)

geneinfo_matrix <- data.table::fread("../Data/geneinfo_beta.txt") #Uplaod gene info matrix

#Determine number of false positives when creating random gene signatures.
#Set the number of signatures and genes
num_signatures <- 100
num_genes <- 800

# Load a sample gene list (you can replace this with your own gene list)
library(biomaRt)

# Use the Ensembl database to retrieve gene names for humans
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(attributes = c("hgnc_symbol"), mart = ensembl)

# Extract the gene names from the gene list
gene_names <- genes$hgnc_symbol

# Create a list to store the signatures
signatures <- list()

# Generate the signatures
set.seed(123)
for (i in 1:num_signatures) {
  # Sample 800 genes from the gene list without replacement
  signature_genes <- sample(gene_names, num_genes, replace=FALSE)
  # Add the signature genes to the signatures list
  signatures[[i]] <- signature_genes
}

drugs <- list()

# Match genes from lists created to those available in the clue.io database and retrieve the gene names stored. Create data frames with the compounds retrieved for each random gene list.
for (i in 1:length(signatures)){
  genes_list <- geneinfo_matrix[match(signatures[[i]], geneinfo_matrix$gene_symbol),]
  genes_list$logFc <-  nrow(genes_list):1
  genes_list <- genes_list[which(!is.na(genes_list$gene_id)),]
  signature <- as.matrix(genes_list$logFc)
  rownames(signature) <- genes_list$gene_id
  drug <- classifyprofile(signature,PRLs = PRLs, type = "fixed", uplength = 100, downlength = 100, signif.fdr = 0.05, nperm = 1000, avgstat = "mean", no.signif = 30, stat = "KS", drugClusters = ap_communities, cluster = "single")
  if (is.null(drug)) {
    drugs[[i]] <- NA
  } else {
    drugs[[i]] <- as.data.frame(drug[[1]])
  }
}

# Extract drugs of each dataframe and store as a list
col_list <- lapply(drugs, function(df) df[[1]][1])

# Create a named list of unique values
unique_vals <- list()
for (i in 1:length(col_list)) {
  unique_vals[[i]] <- unique(col_list[[i]])
  names(unique_vals)[i] <- paste0("df", i)
}

# Create a frequency table of unique values
freq_table <- as.data.frame(table(unlist(unique_vals)))

# plot the frequency table
barplot(freq_table, main = "Frequency of unique values in first column of dataframes")
ggplot(freq_table, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity") + 
  xlab("Drugs") + 
  ylab("Frequency") + 
  theme(axis.text.x = element_blank())