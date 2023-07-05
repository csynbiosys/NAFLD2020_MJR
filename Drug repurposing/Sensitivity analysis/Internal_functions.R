#Internal functions for other scripts

library(parallel)
library(mclust)
library(cluster)

#Load files that we need for these functions
merged_ATC <- read.csv("../Data/List_ACTcodes.csv", stringsAsFactors = FALSE, header = TRUE)
Drug_Targets <- read.csv("../Data/Drug_Targets.csv", stringsAsFactors = FALSE, header = TRUE)
distmat <- read.table("../Data/Distance_matrix.txt", sep = "\t")
drug_communities <- read.table("../Data/AP_clusters.txt", sep = "\t")

#The Fisher test is going to assess the overrepresentation of the ATC code in a specific cluster.
fisher_test_community <- function(atc, atc_community, data, attribute){
  atc_occur_community <- length(atc_community$Drugs[which(atc_community[,attribute]==atc)])
  atc_occur_outside <- length(data$Drugs[which(data[,attribute]==atc)])-atc_occur_community
  other_atc_community <- length(unique(atc_community$Drugs))-atc_occur_community #only drugs with known ATC are counted
  other_atc_outside <- length(unique(data$Drugs))-atc_occur_community- atc_occur_outside - other_atc_community
  cmatrix <- matrix(c(atc_occur_community, other_atc_community, atc_occur_outside,other_atc_outside), ncol=2)
  test=fisher.test(cmatrix, alternative = "greater")
  return(test$p.value)
}

#Retrieves info for Fisher test.
performEnrichTestinCommunity <- function(i, drug_communities, data, attribute){
  #retrieve the drugs in the community
  drugs_community <- drug_communities$Drug[which(drug_communities$Cluster==i)] #1
  
  #retrieve their ATC codes
  atc_community  <- data[data$Drugs%in%drugs_community,]
  
  #enrichment test is only done for atc codes that correspond to at least two drugs
  frequency_codes <- table(atc_community[,attribute])
  freq_codes <- names(frequency_codes)[which(frequency_codes>1)]
  if(length(freq_codes)!=0){
    res <- sapply(freq_codes,function(atc){fisher_test_community(atc, atc_community, data, attribute)})
  }else{
    res <- NA 
    names(res) <- NA
  }
  return(res)
}

#Performs enrichment analysis in a specific cluster.
enrichmentClustering <- function(drug_communities, data, attribute){
  communities <- unique(drug_communities$Cluster)
  list <- lapply(communities,function(x){performEnrichTestinCommunity(x, drug_communities,data,attribute)})
  if (is.list(list)){
    ResultsEnrichment <- data.frame(Cluster=rep(communities,sapply(list, FUN=length)),attribute=as.vector(names(unlist(list))),Pvalues=as.numeric(unlist(list)))
  }else{
    ResultsEnrichment <- data.frame(Cluster=rep(communities,sapply(list, FUN=length)),attribute=as.vector(names(list)),Pvalues=as.numeric(list))
  }
  colnames(ResultsEnrichment)=c("Cluster",attribute,"Pvalues")

  #Correction for multiple testing
  ResultsEnrichment$Pvalues=p.adjust(ResultsEnrichment$Pvalues, method = "BH")
  return(ResultsEnrichment)
}

#Builds data frame that indicates whether every pair of drugs share ATC code or are in the same cluster
buildDrugPairsTable <- function(distmat, merged_ATC, Drug_Targets){
  DrugPairs <- as.matrix(distmat) #retrieve all drug pairs
  DrugPairs[upper.tri(DrugPairs)] <- NA
  diag(DrugPairs) <- NA
  DrugPairs <- melt(DrugPairs, na.rm = TRUE)
  DrugPairs$value <- NULL
  colnames(DrugPairs) <- c("Drug_1","Drug_2")
  
  #Identify which drugs have the same ATC
  cl <- makeCluster(detectCores()-1)
  clusterExport(cl,list=ls("merged_ATC", "DrugPairs"),envir=environment())
  DrugPairs$SameATC <- sapply(c(1:nrow(DrugPairs)),function(x){
    match <- intersect(merged_ATC$Codes[which(merged_ATC$Drugs==DrugPairs$Drug_1[x])], merged_ATC$Codes[which(merged_ATC$Drugs==DrugPairs$Drug_2[x])])
    if (length(match)==0){
      return(FALSE)} 
    else{return(TRUE)}})
  stopCluster(cl)
  
  DrugsATC <- unique(merged_ATC$Drugs)
  DrugPairs$SameATC[!DrugPairs$Drug_1%in%DrugsATC] <- "Not known"
  DrugPairs$SameATC[!DrugPairs$Drug_2%in%DrugsATC] <- "Not known"
  
  #Identify which drugs have the same MoA
  cl <- makeCluster(detectCores()-1)
  clusterExport(cl, list=ls("Drug_Targets", "DrugPairs"),envir=environment())
  DrugPairs$SameDrugTarget <- sapply(c(1:nrow(DrugPairs)),function(x){
    match <- intersect(Drug_Targets$MoA[which(as.character(Drug_Targets$Drugs)==DrugPairs$Drug_1[x])],Drug_Targets$MoA[which(as.character(Drug_Targets$Drugs)==DrugPairs$Drug_2[x])])
    if (length(match)==0){
      return(FALSE)} 
    else{return(TRUE)}})
  
  DrugsDT <- unique(as.character(Drug_Targets$Drugs))
  DrugPairs$SameDrugTarget[!DrugPairs$Drug_1%in%DrugsDT] <- "Not known"
  DrugPairs$SameDrugTarget[!DrugPairs$Drug_2%in%DrugsDT] <- "Not known"
  return(DrugPairs)
}

DrugPairs2 <- buildDrugPairsTable(as.matrix(distmat), merged_ATC, DrugTargets = Drug_Targets)
write.table(DrugPairs, "../Data/DrugPairs.txt", sep = "\t")

addCommunitiestoDrugPairs <- function(DrugPairs,drug_communities){
  DrugPairs$community1 <- drug_communities$Cluster[match(DrugPairs$Drug_1,drug_communities$Drug)]
  DrugPairs$community2 <- drug_communities$Cluster[match(DrugPairs$Drug_2,drug_communities$Drug)]
  DrugPairs$SameCommunity <- FALSE
  DrugPairs$SameCommunity[which(DrugPairs$community1==DrugPairs$community2)] <- TRUE
  return(DrugPairs)
}

runFisherTestDrugPairs <- function(DrugPairs,attribute){
  sameCommunity_sameAttribute <- length(which(DrugPairs$SameCommunity==TRUE & DrugPairs[,attribute]==TRUE))
  sameCommunity_difAttribute <- length(which(DrugPairs$SameCommunity==TRUE & DrugPairs[,attribute]==FALSE))
  difCommunity_sameAttribute <- length(which(DrugPairs$SameCommunity==FALSE & DrugPairs[,attribute]==TRUE))
  difCommunity_difAttribute <- length(which(DrugPairs$SameCommunity==FALSE & DrugPairs[,attribute]==FALSE))
  
  cmatrix <- matrix(c(sameCommunity_sameAttribute,sameCommunity_difAttribute,difCommunity_sameAttribute,difCommunity_difAttribute),ncol=2)
  test <- fisher.test(cmatrix,alternative="greater")
  return(test$p.value)
}

#With all the functions designed above, we can finally create one final function that can perform the enrichment analysis of a clustering according to their ATC codes and MoA. 
#This function will give as output a data frame with the number of clusters, the drugs included, the enriched clusters, the percentage of enriched clusters out of the total number of clusters, Fisher's exact test p-value, the Adjusted Rand Index compared to reference communities and the overlap with the communities found.

enrichedclusters <- function(drug_communities,signific,level, soft, cluster_stats, distmat){
  merged_ATC$Codes <- substring(merged_ATC$Codes,1,level)
  merged_ATC <- unique(merged_ATC)
  
  #Performing the enrichment of the communities with AP clustering method -the one chosen by Iorio- to then compared with the alternative clustering results 
  enriched_ATC_AP <- enrichmentClustering(AP_communities,merged_ATC,"Codes")
  enriched_ATC_AP <- as.character(unique(enriched_ATC_AP$Codes[which(enriched_ATC_AP$Pvalues<0.05)]))
  enriched_DT_AP <- enrichmentClustering(AP_communities,Drug_Targets,"MoA")
  enriched_DT_AP <- as.character(unique(enriched_DT_AP$MoA[which(enriched_DT_AP$Pvalues<0.05)]))
  
  #Enrichment analysis of the clustering results
  enrichment <- enrichmentClustering(drug_communities,merged_ATC,"Codes")
  enriched_clusters_ATC <- length(unique(enrichment$Cluster[which(enrichment$Pvalues < signific)]))
  enriched_ATC <- as.character(unique(enrichment$Codes[which(enrichment$Pvalues<0.05)]))
  intersect_AP <- intersect(enriched_ATC_AP,enriched_ATC)
  percent_overlap_ATC <- length(intersect_AP)/length(enriched_ATC_AP)
  
  enrichment_DT <- enrichmentClustering(drug_communities,Drug_Targets,"MoA")
  enriched_clusters_DT <- length(unique(enrichment_DT$Cluster[which(enrichment_DT$Pvalues < signific)]))
  
  enriched_DT <- as.character(unique(enrichment_DT$MoA[which(enrichment_DT$Pvalues<0.05)]))
  intersect_DT_AP <- intersect(enriched_DT_AP,enriched_DT)
  percent_overlap_DT <- length(intersect_DT_AP)/length(enriched_DT_AP)
  
  #Compute other info about the clustering
  RemovedDrugs <- nrow(distmat)-length(unique(drug_communities$Drug))
  IncludedDrugs <- nrow(distmat)-RemovedDrugs-length(which(table(drug_communities$Cluster)==1))
  totalClusters <- length(which(table(drug_communities$Cluster)>1))
  percentage_ATC <- round(enriched_clusters_ATC*100/totalClusters, 2)
  percentage_DT <- round(enriched_clusters_DT*100/totalClusters,2)
  
  #1 Fisher Test
  DrugPairs_k <- addCommunitiestoDrugPairs(DrugPairs,drug_communities)
  pvalue_ATC <- runFisherTestDrugPairs(DrugPairs_k,"SameATC")
  pvalue_DT <- runFisherTestDrugPairs(DrugPairs_k,"SameDrugTarget")
  vect_res <- as.numeric(drug_communities$Cluster)
  names(vect_res) <- drug_communities$Drug
  
  #Adjusted Rand Index
  if (soft == TRUE){
    part_entropy <- PE(cluster_stats)
    part_coeff <- PC(cluster_stats)
    modpart_coeff <- MPC(cluster_stats)
    vect_res <- vect_res[match(names(ap_vector),names(vect_res))]
    adj <- adjustedRandIndex(ap_vector,vect_res)
    res <- data.frame(K=length(unique(drug_communities$Cluster)),"IncDrugs"= IncludedDrugs, "EnrichedClusters"=enriched_clusters_ATC, "Percentage"=percentage_ATC, "EnrichedClustersDT"=enriched_clusters_DT, "PercentageDT"=percentage_DT, "Pval_ATC"=pvalue_ATC,"Pval_DT"=pvalue_DT,"ARIndex"=adj, "Fuzzy Silhouette width" = sil_index, 'Partition Entropy' = part_coeff, 'Partition Coefficient' = part_coeff, 'Modified Partition Coefficient' = modpart_coeff, "Overlap_ATC" = percent_overlap_ATC, "Overlap_DT" = percent_overlap_DT)
  }else{
    vect_res <- vect_res[match(names(ap_vector),names(vect_res))]
    stats <- cluster.stats(as.dist(distmat), vect_res, ap_vector)
    res <- data.frame("K"=length(unique(drug_communities$Cluster)),"IncDrugs"= IncludedDrugs, "EnrichedClusters"=enriched_clusters_ATC, "Percentage"=percentage_ATC, "EnrichedClustersDT"=enriched_clusters_DT, "PercentageDT"=percentage_DT, "Pval_ATC"=pvalue_ATC,"Pval_DT"=pvalue_DT,"ARIndex"=stats$corrected.rand,"Dunn" = stats$dunn, "Silhouette width" = stats$avg.silwidth, "Variation of information" = stats$vi, "Overlap_ATC"= percent_overlap_ATC, "Overlap_DT"= percent_overlap_DT)
  }
  return(res)
}

#Functions for random enrichment analyses
enrichmentRandom <- function(x,tab,signific,level, soft, distmat){
  communities_random <- tab
  communities_random$Cluster <- sample(communities_random$Cluster)
  return(enrichedclusters(communities_random,signific,level, soft = soft, distmat = distmat))
}

randomParallel <- function(tab,signific,level, soft, distmat){
  ncores <- detectCores() - 1
  cl <- makeCluster(ncores, type = "FORK")
  clusterExport(cl,varlist=ls(),envir=environment())
  clusterEvalQ(cl, library(cluster))
  clusterEvalQ(cl, library(mclust))
  res_random <- parSapply(cl,c(1:100),enrichmentRandom,tab=tab,signific=signific,level=level, soft = soft, distmat = distmat)
  stopCluster(cl)
  res <- data.frame(K=quantile(as.numeric(res_random["K",]),0.50),EnrichedClusters=quantile(as.numeric(res_random["EnrichedClusters",]),0.50),Percentage=quantile(as.numeric(res_random["Percentage",]),0.50),EnrichedClustersDT=quantile(as.numeric(res_random["EnrichedClustersDT",]),0.50),PercentageDT=quantile(as.numeric(res_random["PercentageDT",]),0.50),Pval_ATC=quantile(as.numeric(res_random["Pval_ATC",]),0.50),Pval_DT=quantile(as.numeric(res_random["Pval_DT",]),0.50))
  return(res)
}