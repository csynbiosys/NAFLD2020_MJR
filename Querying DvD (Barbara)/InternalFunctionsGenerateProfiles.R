#### Internal functions to process gene sets to feed into DvD ####

#### Organizing the information to be used to create the Affy Batch object (R object that describes the gene set)####
#Retrieve human readable sample names from GEO object - this is almost done manually, because a regex has to be found for each dataset
sample_groupNames=function(accession,geoobject,reg_ex_gr,reg_ex_repl=NULL){
  d=pData(geoobject)  #Combine data, type, comments, and metadata information to create a new pdata object, or check such an object for consistency
  if (accession=="GSE49541")
  {sample_names=d$characteristics_ch1
    sep="\\1"
  }else{sample_names=d$title
  sep=""}
  groups=gsub(reg_ex_gr,sep,sample_names)
  if(is.null(reg_ex_repl)){
    replicates=NULL
  }else{
  replicates=gsub(reg_ex_repl,"",sample_names)}
  return(list(samples=sample_names,groups=groups,replicates=replicates))}

#Create the samples/groups information in a file
createphenoDataFile=function(accession,geoobject,reg_ex){
  #cel file names
  dir.accession=paste("../GEO data/",accession,sep="")
  celspath = paste(dir.accession,"/","untareddata/",sep="")
  cels = list.files(celspath, pattern = ".CEL$")
  #sample and group names
  samplegroups=sample_groupNames(accession,geoobject,reg_ex)
  sample_names=samplegroups$samples
  sample_names=gsub(" ","_",sample_names)
  groups=samplegroups$groups
  #create the file and save it
  expInfo=data.frame(FileNames=cels,Names=sample_names,Groups=groups)
  write.table(expInfo,file = paste(dir.accession,"untareddata/phenoData.txt",sep="/"),sep="\t",quote = FALSE,row.names = FALSE)
  cat(paste("Finished writing phenoData.txt for",accession))
}

####Clustering####
#To check that samples are clustered in a logical way
clusteringSamples=function(accession,geoobject,sample_names){
  eset = exprs(geoobject)
  LogC=notlogarithmized(eset)
  if(!LogC){eset=2^eset}
  cat("Computing Pearson distances between samples")
  distance = as.dist(1-cor(eset))  #cor: correlation
  cat("\nGenerating the clusters")
  clusters = hclust(distance)
  cat("\nProducing the clustering plot\n")
  png(filename = paste("../Plots",accession,"clusteringHC_complete_samplefile.png",sep="/"),width=800,height=484)
  plot(clusters,main=paste("Hierarchical clustering of samples - ",accession,sep=""))
  dev.off()
  
  colnames(eset)=as.character(sample_names)
  cat("Computing Pearson distances between samples")
  distance = as.dist(1-cor(eset))
  cat("\nGenerating the clusters")
  clusters = hclust(distance)
  cat("\nProducing the clustering plot\n")
  png(filename = paste("../Plots",accession,"clusteringHC_complete.png",sep="/"),width=800,height=484)
  plot(clusters,main=paste("Hierarchical clustering of samples - ",accession,sep=""))
  dev.off()}

#### Log transform if needed####
#limma package deals with logarithmized data
logtransformIfNeeded=function(geoobject){
ex <- exprs(geoobject)
LogC=notlogarithmized(ex)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(geoobject) <- log2(ex) 
cat("The data was log transformed to be suited for limma package")}
return(geoobject)}

notlogarithmized=function(ex){
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  return(LogC)
}

#### Retrieve and refine the annotation - Averaging the expression values of probe ids that correspond to the same gene#####
retrieveAnnotation=function(geoobject,platform){
  probenames = rownames(exprs(geoobject))
  annref=annotationlist[which(annotationlist[,1]==platform),2]
  
  #Accessing biomart - matching probes to gene names
  cat("\nAccessing Biomart...\n")
  human_mart<-useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="useast.ensembl.org") #Possible ensembl index (update): http://uswest.ensembl.org/index.html  http://useast.ensembl.org/index.html http://asia.ensembl.org/index.html
  annotation_ensembl<-getBM(attributes=c(annref,"hgnc_symbol"),filters=annref,values=probenames,mart=human_mart)
  annotation=cbind(annotation_ensembl[[1]],annotation_ensembl[[2]])
  return(annotation)}

#reused internal functions from DvD with some alterations
combineEnsembl <-
  function(data,annotation,listgenes,type="average"){
    stats=list()
    #remove probes which map to more than one gene...
    countprobe<-table(annotation[,1])
    stats$initNrProbes=nrow(data)
    print(paste("Initial number of probes: ",nrow(data),sep=""))
    selp<-names(which(countprobe==1))
    stats$percTotProbesRemoved=round((nrow(data)-length(selp))*100/nrow(data),2)
    print(paste(round((nrow(data)-length(selp))*100/nrow(data),2),"percentage of total probes were removed because they mapped to more than one gene",sep=" "))
    data<-as.matrix(data[selp,])
    annotation<-annotation[which(annotation[,1]%in%selp),]
    #find out how many probes map to a given gene
    countgene<-table(annotation[,2])
    #find all those which are 1-1 probe-gene mappings
    sel<-names(which(countgene==1))
    print(paste(length(sel),"probes match to exactly one gene",sep=" "))
    probes<-annotation[which(annotation[,2]%in%sel),1]
    temp<-as.matrix(data[probes,])
    rownames(temp)<-annotation[which(annotation[,2]%in%sel),2]
    #find those which are not a 1-1 mapping
    combinegenes<-names(which(countgene!=1))
    if(length(combinegenes)>0){
    print(paste("The expression values of", length(combinegenes)-1,"genes are the result of averaging more than 1 probe",sep=" ")) #-1 because "" is not a gene
    subdata<-as.matrix(data[-(which(rownames(data)%in%probes)),]) #subset of the matrix for the probe cases where multiple probes match to a gene
    subannot<-annotation[which(annotation[,1]%in%rownames(subdata)),] # subset of the annotation matrix for these cases
    #sort data according to the genes they map to (so that probes that map to the same gene are in adjacent rows in annotation and data matrices)
    sorder<-subannot[order(subannot[,2]),]
    sdorder<-subdata[sorder[,1],]
    sdorder<-as.matrix(sdorder)
    tabann<-table(sorder[,2])
    #combine results where there is not a 1-1 mapping
    if(type=="average"){
      other<-combined(tabann,sdorder)}
    rownames(other)<-combinegenes
    d<-rbind(temp,other)}
    else {
      print(paste("All probes match exactly one gene",sep=" "))
      d=temp}
    
    l1<-length(which(rownames(d)%in%listgenes))
    if(l1!=length(listgenes)){cat('Note: Ensembl genes do not match listgenes in reference data. Consider uploading pre-processed lists to classifyprofiles',fill=TRUE)}
    print(paste("Original number of genes",length(unique(setdiff(rownames(d),""))), "of which",length(sel)*100/(length(unique(setdiff(rownames(d),"")))),"percent were mapped by a unique probe"))
    stats$annotatedGenes=length(unique(rownames(d)))-1
    stats$percNonPromiscProbesOneGene=length(sel)*100/(length(unique(rownames(d)))-1)
    stats$percNonPromiscProbesMultipleGenes= (length(setdiff(combinegenes,"")))*100/(length(unique(setdiff(rownames(d),""))))
    
    print(paste(length(which(!setdiff(rownames(d),"")%in%listgenes)),"genes were left out because they are not in the reference data",sep=" "))
    d=d[which(rownames(d)%in%listgenes),]
    print(paste("Number of genes common to the data and the reference",length(unique(rownames(d)))))
    stats$intersectionWithRef=length(unique(rownames(d)))
    return(list(d,stats))
  }
#This internal function combines the values of the ranks of probes belonging to the same gene, by averaging them
combined <-
  function(tabann,sdorder){
    
    refs<-cumsum(tabann)  #Returns a vector whose elements are the cumulative sums of the elements of the argument.
    other2<-matrix(nrow=(length(tabann)-1),ncol=ncol(sdorder))
    if(ncol(sdorder)==1){
      other<-mean(sdorder[1:refs[1],])
      for(i in 2:length(tabann)){
        other2[i-1,]<-mean(sdorder[(refs[i-1]+1):refs[i],]) #fixed bug in DvD
      }
    }else{
      other<-colMeans(sdorder[1:refs[1],])
      for(i in 2:length(tabann)){
        other2[i-1,]<-colMeans(sdorder[(refs[i-1]+1):refs[i],]) #fixed bug in DvD
      }
    }
    colnames(other2)<-colnames(sdorder)
    
    rbind(other,other2) # other corresponds to the "" entries. These end up not being considered because "" is not part of listgeness (list of genes in CMap drugs expression profiles)
  }


#This function performs the differential gene expression analysis with limma package
rankingGenes=function(geoobject,groups,datagenes,selgroups,contrast,replicates=NULL){
  #Eliminate samples that aren't relevant for the contrast
  sel=which(groups%in%selgroups)
  geoobject=geoobject[,sel]
  samples=as.factor(groups[sel])
  geoobject$Groups=samples

  #confirm that there are at least two replicates for each group
  noreplicates=which(table(groups[sel])==1)
  if(length(noreplicates)>0){print(paste(names(table(groups[sel]))[noreplicates],"has no replicates"))}
  
  #define a matrix with the sample labels
  design = model.matrix(~0+Groups,geoobject)
  colnames(design) = as.character(levels(samples))
  
  #fit a linear model
  subset_genes=datagenes[,sel]
  if(!is.null(replicates)){
    replicates=as.factor(replicates[sel])
    dc=duplicateCorrelation(subset_genes,design,block=replicates)
    cat("\nFitting a multi-level model, with replicates as random effect")
    fit=lmFit(subset_genes,design,block=replicates,correlation=dc$consensus.correlation)
    
  }else{
    cat("\nFitting a linear model, with no random effects")
  fit = lmFit(subset_genes, design)}
  
  #define the comparisons to make between groups
  cmd <- paste("contrast.matrix <- makeContrasts(DisvsCont=",contrast, ", levels =design)", sep ='"')
  
  eval(parse(text = cmd))
  #produce the results
  fits = contrasts.fit(fit, contrast.matrix)
  res = eBayes(fits)
  
  results = topTable(res, coef=1, number=Inf)
  return(results)}
