####Downloading the gene sets and tidying up####
#This takes more than 30 min, it depends on your internet connection. If the internet connection for some reason breaks, the script fails too
accessions=c("GSE6764","GSE14323","GSE28619","GSE49541")
library(GEOquery)

#retrieve the raw data
dir.create(path="../GEO data")
dir.create(path="../GEO R objects")
dir.create(path="../Plots")
dir.create(path="../DvDResults")
for (accession in accessions){
  dir.accession=paste("../GEO data/",accession,sep="")
  if(!dir.exists(dir.accession)){dir.create(path=dir.accession)}
  if(!dir.exists(paste("../Plots",accession,sep="/"))){
    dir.create(path=paste("../Plots",accession,sep="/"))}
  #download and untar raw files
  geofiles=getGEOSuppFiles(accession,baseDir = "../GEO data/")
  untar(paste(dir.accession,"/",accession,"_RAW.tar",sep=""), exdir=paste(dir.accession,"untareddata",sep="/")) #untar and put everything in a file
  
  #separate the cel files and unzip them
  celspath = paste(dir.accession,"/","untareddata/",sep="")
  cels = list.files(celspath, pattern = ".CEL.gz") # only the cel files are relevant  #list.files: These functions produce a character vector of the names of files or directories in the named directory.
  sapply(paste(celspath, cels, sep="/"), gunzip)
  cels = list.files(celspath, pattern = ".CEL$") # cel files names
  if (accession=="GSE14323"){
    removecels=cels[c(1:9)]
    for(c in removecels){
    file.remove(paste(celspath,c,sep=""))
  }
  }
}
#download the normalized data
for (accession in accessions){
  dir.accession=paste("../GEO R objects/",accession,sep="")
  if(!dir.exists(dir.accession)){dir.create(path=dir.accession)}
  getGEO(GEO=accession, GSEMatrix =TRUE,destdir=dir.accession)
}

#create the folders to store the results
for (accession in accessions){
  dir.accession=paste("../DvDResults/",accession,sep="")
  if(!dir.exists(dir.accession)){dir.create(path=dir.accession)}
}

file.remove("../GEO R objects/GSE14323/GSE14323-GPL96_series_matrix.txt.gz")
file.rename("../GEO R objects/GSE14323/GSE14323-GPL571_series_matrix.txt.gz","../GEO R objects/GSE14323/GSE14323_series_matrix.txt.gz")
