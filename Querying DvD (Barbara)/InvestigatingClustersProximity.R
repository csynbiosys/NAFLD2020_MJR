#### Investigate the proximity of the clusters these drugs yielded by DvD belong to ####
setwd("../CODE_AND_OBJECTS/QueryingDvD/CODE/")
library(DrugVsDisease)

#DvD results - drugs with significant matches to the profile
classify6764_cirvsc=read.csv("../DvDResults/GSE6764/ci-c.csv")
classify14323=read.csv("../DvDResults/GSE14323/cirrhosis-Normal.csv")
classify49541=read.csv("../DvDResults/GSE49541/advanced-mild.csv")
classify28619=read.csv("../DvDResults/GSE28619/Alcoholic-Control.csv")

classify6764_cirvsc=classify6764_cirvsc[c(2,3,7,8,9),]
classify14323=classify14323[c(1,2,6,7,9),]

#Retrieving the distances between drugs computed almost like DvD
drugdistances100avg=read.table("../DvD_distances/DistanceMatrix100avg.txt",check.names=FALSE)
#Retrieve the exemplars of each cluster
exemplars=read.csv("../DvD_distances/exemplars.csv")
exemplars=as.character(exemplars$x)
exemplars_drugdist=drugdistances100avg[exemplars,exemplars]

listclassify=list(classify14323,classify6764_cirvsc,classify49541,classify14323)

#Function to compare the distributions of distances between clusters of drug repurposed drugs and distances of clusters of random drugs
#Input: listclassify: a list with the outputs of DvD classifyprofile function
#	drugdistances: a dataframe with the distances between drugs
#	exemplars_drugdist: the subset of this dataframe with the distances between exemplars
#	exemplars: vector with the names of the exemplars
#	clusters: drug clusters id
# 	drug.clust : NULL if using DvD clusters or the drug communties data frame
compareDistancesHypot2=function(listclassify,drugdistances,exemplars_drugdist,exemplars,clusters,drug.clust){
  #Obtain the average distance between the exemplars of the clusters of the drugs yielded by DvD
  distances=unlist(sapply(listclassify,drug_res_distHypot2,drugdistances,exemp=TRUE,exemplars=exemplars,drug.clust=drug.clust))
  drugs=sapply(listclassify,nrDrugs)
  
  #Testing the difference between the distances between the exemplar drugs that represent the drugs yielded by DvD and the total exemplar drugs, using a bootstrap approach
  
  nulldist=as.vector(sapply(c(1:100000),distExempl,data=exemplars_drugdist,groups=drugs,clusters=clusters,exemplars=exemplars))

  #pvalue for one tailed test
  test=ks.test(jitter(distances),ecdf(nulldist),alternative="greater") # note that the null hypothesis refers to the cumulative distributions so smaller distances corresponds to greater cumulative distribution
  
  #plotting the distances between exemplar drugs and the exemplar drugs that represent the drugs yielded by DvD
  png("../../../DissertationPlots/ClusterProximity.png",width=483,height=412)
  plot(density(nulldist),xlab="Distances between clusters",main="",lwd=3,ylim=c(0,3),font.lab=2,cex.lab=1.2)
  lines(density(distances),lwd=3,col="darkblue")
  text(x=0.4,y=1.5,paste("pval = ",round(test$p.value,3)))
  dev.off()
   return(test$p.value)
}
distExempl=function(i,data,groups,clusters,exemplars){
  alldist=sapply(groups,function(gr){
    indices=sample(c(1:1309),gr,replace = FALSE)
    drug.clust=clusters[indices]
    exemp=as.character(exemplars[drug.clust])
    d=data[exemp,exemp]
    d=d[upper.tri(d)]
    return(d)})
  return(unlist(alldist))
}

#function that retrieves the number of drugs that were yielded by DvD and present negative correlation
nrDrugs=function(classify){
  classify=as.data.frame(classify)
  return(length(classify[,3][which(classify[,4]==-1)]))}

#function that retrieves the distances between every pair of (exemplar) drugs yielded by DvD (that present negative correlation)
#exempl=TRUE results in distances between every pair of the exemplars of the clusters the drugs yielded by DvD belong to
drug_res_distHypot2=function(classify,drugdistances,exempl=TRUE,exemplars,drug.clust=NULL){
  classify=as.data.frame(classify)
  if(exempl){
    if(is.null(drug.clust)){
      drug.clust=classify[,3][which(classify[,4]==-1)] #clusters of the drugs yielded by DvD
    }
    else{
      drugs=as.character(classify[,1][which(classify[,4]==-1)])
      drug.clust=drug.clust$cID[which(drug.clust$DRUGS%in%drugs)]
      
    }
    exemp=as.character(exemplars[drug.clust])}
  else{exemp=as.character(classify[,1][which(classify[,4]==-1)]) }
  subset=drugdistances[exemp,exemp]
  dist_exemp=subset[upper.tri(subset)] #distances between  every pair of (exemplar) drugs yielded by DvD
  return(dist_exemp)
}

#Fisher tests for communities 53 and 98

mat53=matrix(c(2,12,22,1273),ncol=2)
mat98=matrix(c(2,12,68,1227),ncol=2)
fisher.test(mat53)#significant <0.05
fisher.test(mat98)#not significant pval>0.1

#Compare the distances
compareDistancesHypot2(listclassify,drugdistances100avg,exemplars_drugdist,exemplars,drugClusters$Cluster,NULL)
#0.1003974   ##I know obtain 0.101213.
