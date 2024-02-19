library(edgeR)
library(dplyr)
source("./SteatoSITE/DGE_edgeR.R")
source("./SteatoSITE/DGE_limma_voom.R")
source("./SteatoSITE/Pathway_analysis.R")

#Upload final csv file with counts with genes as rows and patient ID as columns called feature_counts
#Upload final csv file with pathological data with patient ID as rows and clinical information as columns called path_metadata

feature_counts <- feature_counts[, path_metadata$sample_id] #Match order of the samples

#Create DGElist object.
dgList_raw <- DGEList(feature_counts, 
                      genes = rownames(feature_counts), 
                      samples = path_metadata, 
                      group = path_metadata$NAFLD)
head(dgList_raw$samples)

dgList_unfiltered <- dgList_raw
outliers <- c("S1071", "S1211", "S1215", "S1216", "S1220", "S1221", "S1223", "S1225", "S1257", "S1281", "S1407", "S1410", "S1412", "S1414", "S1631") #Remove outliers observed in PCA plot
outliers_2 <- path_metadata$sample_id[which(path_metadata$percent_dup_gene >= 80)]
dgList_unfiltered <- dgList_unfiltered[, which(!rownames(dgList_unfiltered$samples) %in% outliers)]
dgList_unfiltered <- dgList_unfiltered[, which(!rownames(dgList_unfiltered$samples) %in% outliers_2)]

#Filter data to eliminate the minimum number of genes according to the smaller group of the dataset.
keep <- rowSums(cpm(dgList_unfiltered) > 1) >= 11 #The smallest group (1b) has 11 samples
summary(keep) 
dgList <- dgList_unfiltered[keep, , keep.lib.sizes = FALSE]

#Normalise data with the trimmed mean of M-values method and obtain the logCPM normalised counts for future uses.
dgList <- calcNormFactors(dgList)
logcpm_norm <- cpm(dgList, log = TRUE, prior.count = 1) 

#Obtain the DEG for the different stages of the disease. Both limma+voom and edgeR can be used.
dea_functions <- function(dgList, fc, fdr, group_var = "group", blocking_vars = c(), contrast_name, print =TRUE){
  dge_lv <- dea_limma(dgList, fc, fdr, group_var, blocking_vars, contrast_name, print =TRUE)
  dge_edgeR <- dea_edgeR(dgList, fc, fdr, group_var, blocking_vars, contrast_name, print =TRUE)
  return(list(dge_lv, dge_edgeR))
}

####Obtain comparisons against controls. In this case, the category used for DGE is NAFLD, and an additional covariate, sex, has been added.
fc = 1 #This can be modified according to user's preferences
fdr = 0.05 #This can be modified according to user's preferences
F4 <- dea_functions(dgList, fc, fdr, group_var = 'NAFLD', blocking_vars = 'sex', contrast_name = 'NAFLDF4-NAFLDControl')
F3 <- dea_functions(dgList, fc, fdr, 'NAFLD', 'sex', 'NAFLDF3-NAFLDControl')
F2 <- dea_functions(dgList, fc, fdr, 'NAFLD', 'sex', 'NAFLDF2-NAFLDControl')
F0F1 <- dea_functions(dgList, fc, fdr, 'NAFLD', "sex", 'NAFLDF0F1-NAFLDControl')
simple_steatosis <- dea_functions(dgList, 'NAFLD', 'sex', 'NAFLDSimple_steatosis-NAFLDControl')

#Perform GSEA to obtain enriched GO terms, reactome and KEGG pathways
kegg <- read.gmt("../Data/c2.cp.kegg.v7.4.symbols.gmt") #Upload gmt files that contain the different terms and pathways to be used for analysis. This is continuously updated in gsea-msigdb.
reactome <- read.gmt("../Data/c2.cp.reactome.v7.4.symbols.gmt")
go <- read.gmt("../Data/c5.go.v7.4.symbols.gmt")
###NAFLD stages vs controls
F4path_lv <- pathway_analysis(F4[[1]], go, reactome, kegg, eps = 0, permut = 100000, lv = TRUE)
F4path_edgeR <- pathway_analysis(F4[[2]], go, reactome, kegg, eps = 0, permut = 1000, lv = FALSE)
F3path_lv <- pathway_analysis(F3[1], go, reactome, kegg, eps = 0, permut = 1000, lv = TRUE)
F3path_edgeR <- pathway_analysis(F3[[2]], go, reactome, kegg, eps = 0, permut = 1000, lv = FALSE)
F2path_lv <- pathway_analysis(F2[[1]], go, reactome, kegg, eps = 0, permut = 1000, lv = TRUE)
F2path_edgeR <- pathway_analysis(F2[[2]], go, reactome, kegg, eps = 0, permut = 1000, lv = FALSE)
F0F1path_lv <- pathway_analysis(F0F1[[1]], go, reactome, kegg, eps = 0, permut = 1000, lv = TRUE)
F0F1path_edgeR <- pathway_analysis(F0F1[[2]], go, reactome, kegg, eps = 0, permut = 1000, lv = FALSE)
simple_steatosispath_lv <- pathway_analysis(simple_steatosis[[1]], go, reactome, kegg, eps = 0, permut = 1000, lv = TRUE)
simple_steatosispath_edgeR <- pathway_analysis(simple_steatosis[[2]], go, reactome, kegg, eps = 0, permut = 1000, lv = FALSE)

