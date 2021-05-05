library(ggplot2)

ma_plot<- function(results_table, fc_threshold, fdr_threshold){
  results_table$significant <- 'no'
  results_table$significant[ results_table$FDR <= fdr_threshold ] <- 'yes'
  ggplot(results_table, aes(logCPM, logFC, color=significant )) + 
    geom_point(size = 0.8) + # This alters the transparency of the points 
    scale_colour_manual(name = 'significant', # provides the label
                        values = setNames(c('red','grey'),c('yes', 'no'))) +
    geom_hline(yintercept=log2(fc_threshold), linetype= "dashed") +
    geom_hline(yintercept=-1*log2(fc_threshold), linetype= "dashed") + 
    theme_bw()
}

volcano_plot<- function(results_table, fc_threshold=2, fdr_threshold=0.05,label_threshold=9){
  genes<-results_table
  genes$Significant<-ifelse((genes$FDR < fdr_threshold & (genes$logFC>=log2(fc_threshold)|genes$logFC<=log2(fc_threshold)*-1)), 'Yes', 'No') #adds a column in the `genes` table to show whether each gene is significantly expressed
  ggplot(data=genes, aes(logFC, -log10(FDR))) + # the values for the x and y axis
    geom_point(aes(col=Significant),size=0.8)+# shows how the points will look
    geom_text_repel(aes(label=ifelse(abs(logFC)>=label_threshold,as.character(row.names(genes)),''))) + # labels all genes with a log fold change greater than the specified threshold
    scale_colour_manual(name = 'Significant', values =
                          setNames(c('red','grey'),c('Yes', 'No'))) +# shows how to decide the colour
    geom_vline(xintercept=c(log2(fc_threshold),log2(fc_threshold)*-1), linetype = "dashed") + # adds the vertical indicating our chosen fold change threshold
    theme_bw()
}