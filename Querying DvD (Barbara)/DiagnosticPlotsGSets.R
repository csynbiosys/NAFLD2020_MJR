#libraries needed
library(affy) 
library(simpleaffy)
library(affyPLM)
library(GEOquery)

#This function generates the diagnostic plots for the datasets used.
generateDiagPlots<-function(accession){
  dir.accession=paste("../GEO data/",accession,sep="")
  #### Loading the Affy Batch object ####
  #the text file is created in generateProfiletoDVD available at ProcessingGeneSetsFeedintoDvD.R
  cat("Loading the Affy Batch object")
  celfiles = read.affy(covdesc="phenoData.txt", path=paste(dir.accession,"untareddata",sep="/"))
  ph = celfiles@phenoData
  
  cat("\nCreating diagnostic plots")
  #### Looking at the raw data and identify severe problems with diagnostic plots####
  cat("\nPlotting the boxplot with log expression before normalization")
  png(filename = paste("../Plots",accession,"boxplot_beforenorm.png",sep="/"),width=800,height=484)
  boxplot(celfiles,which="pm",names=rep("",length(ph@data$Names)),main=paste("log expression (before normalization) ",accession,sep=""))
  dev.off()  #shuts down the specified (by default the current) device.
  # Note: boxplot(affybatch) creates a boxplots of log base 2 intensities by default!So it's directly comparable with the boxplot of normalized data
  cat("\nComputing quality control")
  celfiles.qc <- fitPLM(celfiles)  #fitplm: This function converts an AffyBatch into an PLMset by fitting a specified robust linear model to the probe level data.
  cat("\nPlotting RLE")
  png(filename = paste("../Plots",accession,"rle_beforenorm.png",sep="/"),width=800,height=484)
  RLE(celfiles.qc, main=paste("RLE ",accession,sep=""),names=rep("",length(ph@data$Names)))  #Run Length Encoding: compute the lengths and values of runs of equal values in a vector -- or the reverse operation
  dev.off()  
  cat("\nPlotting NUSE")
  png(filename = paste("../Plots",accession,"nuse_beforenorm.png",sep="/"),width=800,height=484)
  NUSE(celfiles.qc, main=paste("NUSE ",accession,sep=""),names=rep("",length(ph@data$Names)))  #Produce boxplot of Normalized Unscaled Standard Errors (NUSE) for the set of arrays. 
  dev.off()
  
  cat("\nComputing quality control measures")
  celfiles.qcmet = qc(celfiles)  #returns quoted array of character items
  scfact=sfs(celfiles.qcmet)  #Sequential Forward Selection: Applies the Sequential Forward Selection algorithm for Feature Selection.
  if ((max(scfact)-min(scfact))>3){print("Scale factors are not within 3 fold of each other")}
  perc=percent.present(celfiles.qcmet)
  p=which(perc<20)
  if(length(p)>0){print(paste("Sample", names(perc)[p],"has a very low percent present calls",sep=" "))}
  degrat=ratios(celfiles.qcmet)
  d=which(degrat>3)
  if(length(d)>0){print(paste(rownames(degrat)[d]," shows more RNA degradation than desired",sep=""))}
  
  ####Normalization with rma - boxplot####
  cat("\nNormalization with rma")
  celfiles.rma = rma(celfiles)  #This function converts an AffyBatch object into an ExpressionSet object using the robust multi-array average (RMA) expression measure.
  cat("\nPlotting the boxplot with log expression after rma normalization\n")
  png(filename = paste("../Plots",accession,"boxplot_afternormRMA.png",sep="/"),width=800,height=484)
  boxplot(celfiles.rma,names=rep("",length(ph@data$Names)),main=paste("log expression (after rma normalization) - ",accession,sep=""))
  dev.off()
  
  #### Normalization available in GEO -boxplot ####
  geoobject = getGEO(filename=paste("../GEO R objects/",accession,"/",accession,"_series_matrix.txt.gz",sep=""), GSEMatrix =TRUE)
  cat("\nPlotting the boxplot with log expression after their normalization\n")
  png(filename = paste("../Plots",accession,"boxplot_afternormTheirChoice.png",sep="/"),width=800,height=484)
  boxplot(celfiles.rma,names=rep("",length(geoobject$data_row_count)),main=paste("log expression (after their normalization) - ",accession,sep=""))
  dev.off()
}

accessions=c("GSE6764","GSE28619","GSE49541","GSE14323")
for (accession in accessions){
generateDiagPlots(accession)
}

#In the last accession there is an attaching package hgu133a2cdf, but it doesn't match the hgu133plus2cdf. In every accession after plotting the boxplot with log expression after rma normalisation, there are 62 parsing failures. All the plots are identical to the ones in Barbara's dissertation.
