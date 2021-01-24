#Script to convert the matrix available in cmap2data package (after matching probes to genes) to a matrix similar to DRUG_PRLs
#that can be fed into my script for computing drug distances
library(cMap2data)
data("drugRL")
drugRL=drugRL

PRLS_withgenes=apply(drugRL,2,function(x){sortedcol=sort(x) 
return(names(sortedcol))})

#Trying to replicate the 103 clusters
library(mclust)

drugdistances100avg=read.table("../DVD_Distances/DistanceMatrix100avg.txt",check.names=FALSE)
apres=apcluster(1-as.matrix(drugdistances100avg),p=0.0081,details=TRUE)
Drugs=as.character(colnames(drugdistances100avg))
my_DRUG_COMMUNITIES=data.frame("cID"=numeric(length(Drugs)),"DRUGS"=Drugs)
for (i in 1:length(apres)){
     my_DRUG_COMMUNITIES$cID[which(my_DRUG_COMMUNITIES$DRUGS%in%names(apres[[i]]))]=i
   }

#This function generates a vector which shows the number of drugs that are different between "equivalent" communities, with equivalent meaning that they are the communities which the exemplars I've generated belong to
equalCommunities=function(apres){
  exemplars = names(apres@exemplars)
  nr=1
  eq=vector()
  for (i in exemplars){
    NR_COM=drugClusters$Cluster[which(drugClusters$Drug==i)] # cID of the examplar in DvD communities
    DRUG_COM=as.character(drugClusters$Drug[which(drugClusters$Cluster==NR_COM)]) # extract the drugs in that community
    DRUG_AP=names(apres[[nr]]) # drugs in my equivalent community (obtained with my code)
    if(!length(union(setdiff(DRUG_COM,DRUG_AP),setdiff(DRUG_AP,DRUG_COM)))==0){
      eq[nr]=length(union(setdiff(DRUG_COM,DRUG_AP),setdiff(DRUG_AP,DRUG_COM))) # nr of drugs that are different between the two "equivalent" communities
    }else{eq[nr]=0}
    nr=nr+1
  }
return(eq)
}
dif=equalCommunities(apres)
barplot(dif,main="Number of drugs that are different between \n DvD and my communities that \nmatch to the exemplars I produced")

adjustedRandIndex(my_DRUG_COMMUNITIES$cID,drugClusters$Cluster) #yields 0.9886966  #this time is 0.9885966
#Note : length(which(dif!=0)) = 16, so 16 communities are different
exemplars=names(apres@exemplars) # these are the exemplars I've generated

minavg=sapply(c(1:103),function(i){
  drugs=as.character(drugClusters$Drug[which(drugClusters$Cluster==i)])
  subset=drugdistances100avg[drugs,drugs]
  subset[subset==0]=NA
  avg=colMeans(subset,na.rm = TRUE)
  return(names(which.min(avg)))})

setdiff(minavg,exemplars) # only one exemplar is different (testosterone). This drug belongs to a community whose drugs do not match to any of the exemplars I produced.


exemplars=minavg #These were set as the exemplars

write.csv(exemplars,"exemplars.csv",row.names = FALSE)

