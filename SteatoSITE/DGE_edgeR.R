library(edgeR)

#Create function to perform differential gene expression where the chosen differentially expressed genes have a FDR < 0.05 and |fold-change| >=2. The output is a list that contains the main objects of the analysis, as well as the genes that meet the criteria.
dea_edgeR <- function(dgList, group_var = "group", blocking_vars = c(), contrast_name, print =TRUE){
  fc_threshold <- 2
  fdr_threshold <- 0.05
  formula <- paste(c('~0', c(group_var, blocking_vars)), collapse = ' + ' ) # create the formula for he experimental design
  design <- model.matrix(as.formula(formula), data = dgList$samples)
  dgGlm <- estimateDisp(dgList[,rownames(design)], design, robust = TRUE) # estimate the dispersion
  fit <- glmQLFit(dgGlm, design, robust = TRUE) # fit the data to the design
  contrast_matrix <- makeContrasts(contrasts =contrast_name, levels=design) #create the contrast
  de <- glmQLFTest(fit, contrast=contrast_matrix) # test for differentially expressed genes
  results <- topTags(de, n=nrow(dgList) )$table # extract the results table
  de_genes <- rownames(results)[abs(results$logFC) >= log2(fc_threshold) & results$FDR <= fdr_threshold ] # extract genes passing FDR and fold change thresholds.
  if (print){ 
    print(paste(length(de_genes), 'genes are differentially expressed for contrast', contrast_name, 'at a fold change of at least', fc_threshold, 'and a maximum FDR of', fdr_threshold, '.', length(results$logFC[results$logFC >= log2(fc_threshold) & results$FDR <= fdr_threshold]), 'genes are upregulated and', length(results$logFC[results$logFC <= log2(fc_threshold)*-1 & results$FDR <= fdr_threshold]), 'are downregulated.'))
  }
  list(  #Create a list with several information it might be needed for downstream analyses.
    dgGlm = dgGlm,
    dgList = dgList,
    fit =fit,
    results = results,
    diffexp_genes = de_genes, 
    contrast_matrix = contrast_matrix, 
    design = design,
    de = de_genes
  ) 
}