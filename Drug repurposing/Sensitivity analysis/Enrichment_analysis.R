library(apcluster)
library(cluster)
library(e1071)
library(parallel)
library(tidyr)

load("./DT_MoA.R")
source("./Internal_functions.R")

ap_vector <- durg_communities$Cluster
names(ap_vector) <- drug_communities$Drug

#AP clustering
obtain_AP <- function(preferences){
  ap_result <- apcluster(1-distmat,details=TRUE,seed=10,p=preferences) 
  Drugs <- colnames(distmat)
  drug_communities <- data.frame("Cluster"=numeric(length(Drugs)), "Drug"=Drugs)
  for (i in 1:length(ap_result)){
    drug_communities$Cluster[which(drug_communities$Drug%in%names(ap_result[[i]]))] <- i
  }
  return(drug_communities)
}

enriched_AP <- function(p,signific,level, soft = FALSE, distmat){
  ap_table <- obtain_AP(p)
  return(enrichedclusters(ap_table,signific,level,soft, distmat))
}

#HC clustering
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")  

agglo_coef <- function(x) {
  agnes(as.dist(distmat), method = x)$ac
}

ac_comparison <- sapply(m, agglo_coef) #Comparison to check which linkage method is the best

obtain_HC <- function(i, type){
  hc_clust <- hclust(as.dist(distmat),method= type) #since it's needed a dissimilarity matrix we just need to use the as.dist function with the distance matrix created before.
  hc_res <- cutree(hc_clust, i)  #With cutree we can cut the dendogram into several groups according to the number of clusters.
  HC_df <- data.frame(Cluster= integer(nrow(distmat)), Drug=character(nrow(distmat)))
  HC_df$Drug <- names(hc_res)
  HC_df$Cluster <- hc_res
  return(HC_df)
}

enriched_HC  <- function(i, type, signific,level, soft = FALSE, distmat){
  hc_table <- obtainCommTableHC(i, type)
  return(enrichedclusters(hc_table,signific, level, soft, distmat))
}

#PAM clustering
obtain_PAM <- function(i){
  pam_res <- pam(as.dist(distmat),i,diss=TRUE)  #i: n of clusters
  pam_df <- data.frame(Cluster= integer(nrow(distmat)), Drug= character(nrow(distmat)))
  pam_df$Drug <- names(pam_res$clustering)
  pam_df$Cluster <- pam_res$clustering
  return(pam_df)
}

enrichedclustersPAM <- function(i,signific,level, soft = FALSE, distmat){
  pam_table <- obtainCommTablePAM(i)
  return(enrichedclusters(pam_table,signific,level, soft, distmat))
}

#ClusterOne clustering

#First we need to perform the clustering with java on the terminal. The code used is the following:
# for dens in (seq 0.0125 0.01 0.7425); do java -jar cluster-one.jar #inputfile.txt -s 1 -d $dens > #outputfile.txt; done
#The distances have to be a similarity matrix with values between 0 and 1, so they were transformed into s = 1/(1+d)

drug_distances <- as.data.frame(distmat)
drug_distances[, 1] <- rownames(distmat)
rownames(drug_distances) <- drug_distances[, 1]
drug_distances <- pivot_longer(drug_distances, cols = 2:(ncol(drug_distances)), names_to = 'id2', values_to = 'weights')
drug_distances <- drug_distances[which(drug_distances$weights != '0'), ]
drug_distances <- drug_distances[drug_distances$weights < quantile(drug_distances$weights, probs = 0.05), ]
drug_distances$weights <- 1/(1+drug_distances$weights)
write.table(drug_distances, '../Data/ClOne/Drug_dist_CLone.txt', row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

obtain_ClONE <- function(dens){
  results.list <- read.table(paste("../Results/ClusterOne/Results_CL_",dens,".txt",sep=""), sep = ",", row.names = 1, header = TRUE)
  return(results.list)
}

enriched_ClONE <- function(dens,signific,level,soft,distmat){     #Computes densities of observations in parameterized MVN mixtures.
  cl_table <- obtainCommTableClONE(dens)
  return(enrichedclusters(cl_table,signific,level,soft, distmat))
}

#Fuzzy c-means
obtainCommTableFCM <- function(i){
  fcm_res <- cmeans(x = 1-distmat, centers = i) #i: n of clusters
  fcm_df <- data.frame(Cluster= integer(nrow(distmat)), Drug= character(nrow(distmat)))
  fcm_df$Drug <- colnames(distmat)
  fcm_df$Cluster <- fcm_res$cluster
  return(list(fcm_df, fcm_res$membership))
}

enrichedclustersFCM <- function(i,signific,level, soft, distmat){
  fcm <- obtainCommTableFCM(i)
  fcm_table <- fcm[[1]]
  fcm_stats <- fcm[[2]]
  return(enrichedclusters(fcm_table,signific,level,soft, fcm_stats, distmat))
}

#We will compare these clusters with random ones
enrichedRandomAP <- function(p,signific,level){
  ap_table <- obtainCommTableAP(p)
  res_random <- randomParallel(ap_table,0.05,4)
  return(res_random)
}

enrichedRandomClONE <- function(d,signific,level){
  clone_table <- obtainCommTableClONE(d)
  res_random <- randomParallel(clone_table,signific,level, soft = soft, distmat = distmat)
  return(res_random)
}

enrichedRandomHC <- function(i,signific,level,soft, distmat){
  hc_tab <- obtainCommTableHC(i)
  res_random <- randomParallel(hc_tab,signific,level, soft = soft, distmat = distmat)
  return(res_random)
}

enrichedRandomPAM <- function(i,signific,level,soft, distmat){
  pam_table <- obtainCommTablePAM(i)
  res_random <- randomParallel(pam_table,signific,level, soft = soft, distmat = distmat)
  return(res_random)
}

