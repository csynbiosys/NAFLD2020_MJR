library(limma)
library(edgeR)

#Differential gene expression with limma+voom. It chooses those genes according to chosen values for FDR < fc and |FC| >= fc. The output returns a list of objects obtained during the analysis, including those genes that meet the criteria.

dea_limma <- function(dgList, fc, fdr, group_var = "group", blocking_vars = c(), contrast_name, print =TRUE){
  fc_threshold <- fc
  fdr_threshold <- fdr
  formula <- paste(c('~0', c(group_var, blocking_vars)), collapse = ' + ' ) # create the formula for he experimental design
  design <- model.matrix(as.formula(formula), data = dgList$samples)
  v <- voom(dgList, design) #Transform count data to logCPM, estimate the mean-variance relationship and use this to compute appropriate observation-level weights. 
  fit <- lmFit(v, design) #Fit linear model for each gene.
  contrast_matrix <- makeContrasts(contrasts =contrast_name, levels=design) #create the contrast
  fit_tmp <- contrasts.fit(fit, contrast_matrix)
  fit_tmp <- eBayes(fit_tmp, robust = TRUE) #compute moderated t-statistics, moderated F-statistic, and log-odds of de by empirical Bayes moderation of the standard errors towards a global value.
  results <- topTable(fit_tmp, sort.by = "logFC", number = "Inf")# extract the results table
  de_genes <- rownames(results)[abs(results$logFC) >= log2(fc_threshold) & results$adj.P.Val <= fdr_threshold] # extract genes passing FDR and fold change thresholds.
  if (print){ 
    print(paste(length(de_genes), 'genes are differentially expressed for contrast', contrast_name, 'at a fold change of at least', fc_threshold, 'and a maximum FDR of', fdr_threshold, '.', length(results$logFC[results$logFC >= log2(fc_threshold) & results$adj.P.Val <= fdr_threshold]), 'genes are upregulated and', length(results$logFC[results$logFC <= log2(fc_threshold)*-1 & results$adj.P.Val <= fdr_threshold]), 'are downregulated.'))
  }
  list(  #Create a list with several information it might be needed for downstream analyses.
    v = v,
    dgList = dgList,
    fit =fit,
    fit_tmp = fit_tmp,
    results = results,
    diffexp_genes = de_genes, 
    contrast_matrix = contrast_matrix, 
    design = design,
    de = de_genes
  ) 
}
