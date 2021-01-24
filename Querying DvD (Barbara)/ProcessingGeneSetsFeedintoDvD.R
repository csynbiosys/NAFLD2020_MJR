#Querying DvD with CLDs datasets from GEO
#Processing GEO series (GSE) so that they can be fed into DvD function classifyprofiles

#libraries needed
library(affy) 
library(simpleaffy)
library(affyPLM)
library(DrugVsDisease)
library(GEOquery)
library(hgu133plus2hsentrezg.db)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pd.hg.u133.plus.2")

library(statmod)
library("illuminaRatv1.db")

#using drugvsdisease data
data(annotationlist)
data(genelist)

#load scripts with useful functions
source("InternalFunctionsGenerateProfiles.R")

#GEO accessions of the datasets used in the analysis
accessions=c("GSE6764","GSE63726","GSE49541","GSE14323","GSE63067","GSE28619")
#these are the regex expressions needed to dived the samples into groups for the differential gene expression analysis
reg_ex_gr=list("_.+","\\_[1-9]","^Stage: (.*?) \\(.*$",",.*","[\\s+|,|-].*$","Liver sample from|hepatitis group|group|\\s|\\(.*\\)")
names(reg_ex_gr)=accessions
#attempt to include patients as random effect - datasets that I end up using don't have several samples from the same patient so it's not needed
reg_ex_repl=list(NULL,NULL,NULL,NULL,NULL,NULL,NULL)
names(reg_ex_repl)=accessions

#This function is the equivalent of DvD's function generateprofiles. I couldn't use their function because I was using GEO series datasets (GSE) and not GDS. DvD doesn't process GSE datsets
#Input: accession: GEO accession identifier
#       reg_ex_gr: vector with the regex expression to retrieve the classification of the samples according to the group they belong to (e.g. disease vs healthy). Each element corresponds to the regex for a GSE dataset
#       reg_ex_repl:vector with the regex expression to retrieve the classification of the samples according to the patient they belong to. Each element corresponds to the regex for a GSE. NULL indicates that patient is not consider a random effect for the linear model for differential expression
#       selgroups: vector with the two groups used for the comparison of gene expression profiles
#       contrast : string with the cotnrast to make "disease-control" where disease should be replaced with the group label for disease samples and control with the group label for healthy samples
generateProfiletoDVD=function(accession,reg_ex_gr,reg_ex_repl,selgroups,contrast){
#### Loading the respective geo expression set object ####
cat(paste("\nLoading the geo expression set object for ",accession," \n",sep=""))
geoobject = getGEO(filename=paste("../GEO R objects/",accession,"/",accession,"_series_matrix.txt.gz",sep=""), GSEMatrix =TRUE)
samplegroups=sample_groupNames(accession,geoobject,reg_ex_gr[[accession]],reg_ex_repl[[accession]])
sample_names=samplegroups$samples
groups=samplegroups$groups
if(accession=="GSE6764"){groups=tolower(groups)}
replicates=samplegroups$replicates
platform = annotation(geoobject)

#### Organizing the information to be used to create the Affy Batch object (R object that described the gene set)####
if(platform=="GPL570"|platform=="GPL14877"|platform=="GPL571"){
cat("\nOrganizing the information to be used to create the Affy Batch object \n")
createphenoDataFile(accession,geoobject,reg_ex_gr[[accession]])}

#### Updating groups in geoobject ####
phdata=data.frame(row.names = colnames(exprs(geoobject)),Groups=sample_names)
phdata=as(phdata,"AnnotatedDataFrame")
geoobject=ExpressionSet(assayData = exprs(geoobject),phenoData = phdata)

####Clustering####
cat("\nClustering samples \n")
clusteringSamples(accession,geoobject,sample_names)

#### Log2 transform, if needed (reused from GEO2R) ####
cat("\nLog2 transform, if needed \n")
geoobject=logtransformIfNeeded(geoobject)


####Mapping probes to genes, using Biomart####
cat("\nMapping probes to genes, using Biomart \n")
if(platform =="GPL570"){
  platform = "hgu133plus2"
  annotation=retrieveAnnotation(geoobject,platform)
}else if (platform == "GPL571"){
  platform = "hgu133a2"
  annotation=retrieveAnnotation(geoobject,platform)
}else if (platform == "GPL14877"){
  #this platform has its own annotation, cannot be retrieved from BioMart
  x <- hgu133plus2hsentrezgSYMBOL
  mapped_probes <- mappedkeys(x)
  xx <- as.list(x[mapped_probes])
  annotation=matrix(c(as.character(names(xx)),as.character(unname(xx))),ncol=2)
  probes=rownames(exprs(geoobject))
  annotation=annotation[which(annotation[,1]%in%probes),]
  
}else if (platform=="GPL6101"){
  #The samples of this platform come from rattus norvegicus, so human orthologous genes have to be found first
  x <- illuminaRatv1SYMBOL
  mapped_probes <- mappedkeys(x)
  xx <- as.list(x[mapped_probes])
  annotationrat=matrix(c(as.character(names(xx)),as.character(unname(xx))),ncol=2)
  probes=rownames(exprs(geoobject))
  annotationrat=annotationrat[which(annotationrat[,1]%in%probes),]
  rgenes=annotationrat[,2]
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  rattus = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  
  genes_RandH= getLDS(attributes = c("rgd_symbol"), filters = "rgd_symbol", values =  annotationrat[,2] , mart = rattus, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  hgene.symbol=genes_RandH$HGNC.symbol[match(rgenes,genes_RandH$RGD.symbol)]
  annotation=matrix(c(as.character(annotationrat[,1]),hgene.symbol),ncol=2)
  annotation[,2][is.na(annotation[,2])]=""
  }

    
#Removing probes that map to multiple genes and when there isn't a 1-1 mapping (several probes to the same gene), the average of the expression values of those probes is computed
annotationProcess=combineEnsembl(exprs(geoobject),annotation,genelist,type="average")
datagenes=annotationProcess[[1]]
statsAnnot=annotationProcess[[2]]

####Ranking the genes using limma package####
cat("\nRanking the genes using limma package \n")
results=rankingGenes(geoobject,groups,datagenes,selgroups, contrast,replicates)
png(paste("../Plots/",accession,"/pvalues",contrast,".png",sep=""))
hist(results$P.Value)
dev.off()

####Tidying up the results so that they can be fed into DvD####
cat("\nTidying up the results so that they can be fed into DvD\n")
results=results[order(results$logFC,decreasing=TRUE),]
RankedGenes=list(ranklist=as.matrix(results$logFC),pvalues=as.matrix(results$adj.P.Val))
attr(RankedGenes$ranklist,"dimnames")=list(rownames(results),contrast)
attr(RankedGenes$pvalues,"dimnames")=list(rownames(results),contrast)
cat("\nDone")
return(RankedGenes)}

#This function uses DvD to find drug repurposing candidates based on the disease gene expression profiles. 
#A table with the drugs that match the gene expression profile is produced and saved in DvDResults folder
findDrugRepurposingCandidates=function(accession,reg_ex_gr,reg_ex_repl,selgroups,contrast,fdr=0.05){
#generate disease gene expression profile
RankedGenes=generateProfiletoDVD(accession,reg_ex_gr,reg_ex_repl,selgroups,contrast)
#Matching disease gene expression profile to CMap drugs expression profile
cat("\nRetrieving drug repurposing candidates using DvD\n")
classify=tryCatch(classifyprofile(RankedGenes$ranklist,adj="qvalue",signif.fdr = fdr,case="disease"),error=function(e){classify=classifyprofile(RankedGenes$ranklist,adj="qvalue",signif.fdr = fdr,case="disease")
return(classify)})
#Saving the results
filename=paste("../DvDResults/",accession,"/",contrast,".csv",sep="")
cat(paste("\nSaving the results in",filename))
write.csv(classify,filename,row.names = FALSE)
return(classify)
}


#### Finding Drug Repurposing candidates for the set of GEO datasets ####

classify6764_cirvsc=findDrugRepurposingCandidates("GSE6764",reg_ex_gr,reg_ex_repl,c("c","ci"),"ci-c",fdr=0.05)
classify49541=findDrugRepurposingCandidates("GSE49541",reg_ex_gr,reg_ex_repl,c("mild","advanced"),"advanced-mild",fdr=0.05)
classify14323=findDrugRepurposingCandidates("GSE14323",reg_ex_gr,reg_ex_repl,c("Normal","cirrhosis"),"cirrhosis-Normal",fdr=0.05)
classify28619=findDrugRepurposingCandidates("GSE28619",reg_ex_gr,reg_ex_repl,c("Control","Alcoholic"),"Alcoholic-Control",fdr=0.05)

#Datasets that were removed
#classify63726=tryCatch(findDrugRepurposingCandidates("GSE63726",reg_ex_gr,reg_ex_repl,c("DEN_HC","PBS_HC"),"DEN_HC-PBS_HC",fdr=0.05),error=function(e){classify63726=findDrugRepurposingCandidates("GSE63726",reg_ex_gr,reg_ex_repl,c("DEN_HC","PBS_HC"),"DEN_HC-PBS_HC",fdr=0.05)
#return(classify63726)}) # little overlap due to probe conversions
#classify63067=findDrugRepurposingCandidates("GSE63067",reg_ex_gr,reg_ex_repl,c("Healthy","Non"),"Non-Healthy",fdr=0.05) # no signficant results
#classify33650=findDrugRepurposingCandidates("GSE33650",reg_ex_gr,reg_ex_repl,c("LowHP","HighHP"),"HighHP-LowHP",fdr=0.2) # problems with clustering. few differentially expressed genes


##Test - check that I could replicate DvD example accession=gse17906, reg_ex_gr = ".*_" replicates=null selgroups=c("CaP","NP"), contrast="CaP-NP"
## Top negative match result was estradiol, the same drug they indicated!

