#Reproducing the clusters obtained used in DvD - Ran on EDDIE output: DvD_distances/DistanceMatrix100avg.txt
## Step 2: Computing the distances, using the PRLs built previously

#args <- commandArgs(TRUE)
#file_prls=args[1]
#DRUG_PRLs=read.table(file_prls,header=TRUE,check.names = FALSE)

packages <- c("parallel", "fastmatch")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
library(parallel)

#' retrieveSignature returns the top 100 and bottom 100 genes. These are the signature of each drug
#'
#' @param PRLa The PRL (prototype ranekd list) of a given drug
#' @return signature - a list where Top contains top 100 genes in the PRL and Bottom the bottom 100 ones.
retrieveSignature=function(PRLa){
  top100=PRLa[1:100]
  len=length(PRLa)
  bottom100=PRLa[(len-100):len]
  signature=list("Top"=top100,"Bottom"=bottom100)
  return(signature)
}

library(fastmatch)

#' enrichmentScore measures how much a set of genes (belonging to signature of drug A) is at the top of the PRL of drug B .
#' computed as described in Subramanian 2005 paper
#'
#' @param p - vector with the genes belonging to the signature of drug A
#' @param B - PRL of drug B (vector with genes ordered by their rankings)
#' @return ES - enrichement score
enrichmentScore=function(p,B){
  Nh=length(p)
  N=length(B)
  Nr=length(intersect(p,B))
  sums=vector(mode="numeric",length=N)
  sum=0
  ES=0
  matches=B%fin%p
  for (i in 1:N){
    if(matches[i]){
      sum=sum+1
    }
    sums[i]=sum
  }
  Phit=sums/Nr
  nmiss=seq(1, N, by=1)-sums
  Pmiss=1/(N-Nh)*nmiss
  dif=Phit-Pmiss
  ES=dif[which.max(abs(dif))]
  return(ES)
  
}

#' totalEnrichmentScore computes the inverse total enrichment score of the drug signature of drug A regarding PRL of drug B
#'
#' @param enUp - result of enrichmentScore(p,B) where p is the top signature of A and B the PRL of drug B
#' @param enDown - result of enrichmentScore(p,B) where p is the bottom signature of A and B the PRL of drug B
#' @return TES - inverse total enrichement score
totalEnrichmentScore=function(enUp,enDown){
  TES=1-(enUp-enDown)/2
  return(TES)
}

#' maxESDistance computes the maximum enrichment score distance between two drugs
#' This is the one Iorio uses for the computation of distances
#' 
#' @param TESa - result of totalEnrichmentScore of signature of drug A regarding PRL of drug B
#' @param TESb - result of totalEnrichmentScore of signature of drug B regarding PRL of drug A
#' @return D - maximum enrichment score distance
maxESDistance=function(TESa,TESb){
  D=min(c(TESa,TESb))
  return(D)
}

#' avgESDistance computes the average enrichment score distance between two drugs
#' 
#' @param TESa - result of totalEnrichmentScore of signature of drug A regarding PRL of drug B
#' @param TESb - result of totalEnrichmentScore of signature of drug B regarding PRL of drug A
#' @return D - average enrichment score distance
avgESDistance=function(TESa,TESb){
  D=(TESa+TESb)/2
  return(D)
}


#' computeDist computes the distance between two drugs A and B
#' 
#' @param PRLa - PRL of drug A
#' @param PRLb - PRL of drug B
#' @return D - distance between drugs A and B
computeDist=function(PRLa,PRLb){
  signaturesA=retrieveSignature(PRLa)
  signaturesB=retrieveSignature(PRLb)
  ESupA=enrichmentScore(signaturesA[[1]],PRLb)
  ESdownA=enrichmentScore(signaturesA[[2]],PRLb)
  ESupB=enrichmentScore(signaturesB[[1]],PRLa)
  ESdownB=enrichmentScore(signaturesB[[2]],PRLa)
  
  TESa=totalEnrichmentScore(ESupA,ESdownA)
  TESb=totalEnrichmentScore(ESupB,ESdownB)
  dist=avgESDistance(TESa,TESb)
  
  return(dist)
}

#' buildDistMatrix computes the distance matrix corresponding to the drugs whose PRLs are in DRUG_PRLs
#' 
#' @param DRUG_PRLs - dataframe where each column represents a drug, each row a position (rank) and each entry is the gene in that rank
#' @return distances - matrix with the distances between the drugs
buildDistMatrix=function(DRUG_PRLs){
  cl=makeCluster(8,type = "FORK")
  clusterEvalQ(cl, library(fastmatch))
  clusterExport(cl=cl,varlist=c("retrieveSignature", "enrichmentScore", "totalEnrichmentScore", "maxESDistance","computeDist"))
  
  combin=combn(ncol(DRUG_PRLs),2)
  distances= parApply(cl,combin,2,function(x){round(computeDist(DRUG_PRLs[,x[1]],DRUG_PRLs[,x[2]]),5)})
  n=length(DRUG_PRLs)
  distances_mat=matrix(ncol=n,nrow=n)
  distances_mat[,]=0
  for( i in 1:length(combin[1,])){
    a=combin[1,i] 
    b=combin[2,i]
    distances_mat[a,b]=distances_mat[b,a]=distances[i]
  }
  
  rownames(distances_mat)=colnames(DRUG_PRLs)
  colnames(distances_mat)=colnames(DRUG_PRLs)
  stopCluster(cl)
  return(distances_mat)
}


DRUG_PRLs=read.csv("DrugPRLsgenes.csv",check.names = FALSE)  #missing file, not able to find it and Barbara doesn't have access to it.
distMatrix=buildDistMatrix(DRUG_PRLs)
write.table(distMatrix,"/exports/eddie/scratch/s1669835/DistanceMatrix100avg_plus1.txt",sep="\t")
