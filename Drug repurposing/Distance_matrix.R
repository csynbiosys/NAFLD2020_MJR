library(parallel)
library(fastmatch)
library(apcluster)
source("./Prototype_ranked_lists.R")

#This script will get the distance matrix for the PRLs generated.

#Returns first and last 250 genes of the gene signature. This number was determined by F.Iorio and colleagues.
retrieve_signature <- function(PRLa){
  top250 <- PRLa[1:250] #head(PRLa, 250)
  len <- length(PRLa)
  bottom250 <- PRLa[(len-249):len] #tail(PRLa, 250)
  signature <- list("Top"= top250,"Bottom"= bottom250)
  return(signature)
}

#Enrichment Score (ES) is based on the Kolmogorov-Smirnov statistics, which can determine how much a set of genes is at the top of the list. It ranges from -1 to 1, the closer to 1 the closer the genes are to the top of the list; the closer to -1, the closer the genes are to the bottom of the list.
#The function will take to inputs:`input_sig`, a vector with the gene signature that is used as input and the PRL that is compared with the gene signature, with genes ordered by their rankings.
ES <- function(input_sig,PRL){
  NH <- length(input_sig)
  N <- length(PRL)
  NR <- length(intersect(as.character(input_sig), as.character(PRL)))
  hits <- PRL%fin%input_sig     
  hitCases <- cumsum(hits)
  missCases <- cumsum(1-hits)
  Phit <- hitCases/NR
  Pmiss <- missCases/(N-NH)
  dif <- Phit-Pmiss
  ES <- dif[which.max(abs(dif))]
  return(ES)
}

#Then, Inverse Total Enrichment Score (TES) is calculated and determines how many genes from the bottom are placed at the bottom and viceversa.
TES <- function(ES_top, ES_bottom){
  TES <- 1 - (ES_top - ES_bottom)/2
  return(TES)
}

#The next function computes the Maximum Enrichment Score distance between two drugs.
maxESDist <- function(TESa, TESb){
  D <- min(c(TESa, TESb))
  return(D)
}

#Final function that will compute the final distance between the signatures.
compute_dist <- function(PRLa,PRLb){
  sigA <- retrieve_signature(PRLa)
  sigB <- retrieve_signature(PRLb)
  EStopA <- ES(sigA[[1]],PRLb)
  ESbottomA <- ES(sigA[[2]],PRLb)
  EStopB <- ES(sigB[[1]],PRLa)
  ESbottomB <- ES(sigB[[2]],PRLa)
  
  TESa <- TES(EStopA,ESbottomA)
  TESb <- TES(EStopB,ESbottomB)
  dist <- maxESDist(TESa,TESb)
  
  return(dist)
}

dist_matrix <- function(PRLs){
  ncores <- detectCores() - 1
  cl<- makeCluster(ncores, type = "FORK")
  clusterEvalQ(cl, library(fastmatch))
  clusterExport(cl=cl,varlist=c("retrieve_signature", "ES", "TES", "maxESDist","compute_dist"))
  
  combin <- combn(ncol(PRLs),2)
  distances <- parApply(cl,combin,2,function(x){
    round(compute_dist(PRLs[,x[1]],PRLs[,x[2]]),5)
  })
  n <- length(PRLs)
  distances_mat <- matrix(ncol=n,nrow=n)
  distances_mat[,] <- 0
  for (i in 1:length(combin[1,])){
    a <- combin[1,i] 
    b <- combin[2,i]
    distances_mat[a,b] = distances_mat[b,a] = distances[i]
  }
  
  rownames(distances_mat) <- colnames(PRLs)
  colnames(distances_mat) <- colnames(PRLs)
  stopCluster(cl)
  return(distances_mat)
}

PRLS_withgenes <- apply(PRLs,2,function(x){sortedcol=sort(x) 
return(names(sortedcol))}) #Returns matrix where the genes are ranked for each drug.

distmat <- dist_matrix(PRLS_withgenes)
write.table(distmat, '../Data/Distance_matrix.txt', sep = '\t')

#Generate ap cluster for the sensitivity analysis and create data frame that contains the name of the drug as one column and the cluster they belong to as another.
ap_result <- apcluster(1-as.matrix(distmat), details = TRUE, seed = 10) #In this case a similarity matrix is needed. Since distmat is measures the distance between the drugs (the farther away, the most different), to obtain the similarity matrix we need to do 1-distmat.
ap_communities <- function(ap_result, distmat){
  AP_communities <- data.frame(Cluster = 1:length(colnames(distmat)), Drugs = colnames(distmat))
  for (i in 1:length(ap_result@clusters)){
    AP_communities$Cluster[which(AP_communities$Drugs%in%names(ap_res@clusters[[i]]))] <- i
  }
  #vect <- as.vector(AP_communities$Cluster)
  #names(vect) <- AP_communities$Drugs
  #stats <- cluster.stats(as.dist(1-distmat), vect) #Returns statistics if necesssary
  return(AP_communities)
}

ap_comm <- ap_communities(distmat)

write.csv(ap_comm, "./Data/AP_clusters.txt", sep = "\t")