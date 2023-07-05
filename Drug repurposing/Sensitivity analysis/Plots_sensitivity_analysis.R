library(ggplot2)
library(reshape)
library(pheatmap)
library(cluster)

merged_ATC <- read.csv("./Data/ProcessedATCcodes.csv",stringsAsFactors = FALSE,header=TRUE)
Drug_Targets <- read.csv("./Data/ProcessedDTargets.csv",header=TRUE,stringsAsFactors = FALSE)

colours <- c("cornflowerblue","forestgreen","gray65","red3","gold", "darkorchid1", "burlywood4")
names(colours) <- levels(as.factor(c("AP", "HC_average", "PAM","HC_wardD2", "HC_complete", "FCM", "ClOne")))
colScale <- scale_colour_manual(values = colours,  labels = c("AP", "ClOne", "FCM", "HC average", "HC complete", "HC ward.D2", "PAM"))
colfillScale <- scale_fill_manual(values = colours, labels = c("AP", "ClOne", "FCM", "HC average", "HC complete", "HC ward.D2", "PAM"))

#Function to rearrange the data obtained for easier manipulation for visualisation
tidyclusters <- function(filename,algorithm){
  table <- read.csv(filename,header=TRUE)
  names <- as.character(table$X)
  table <- t(table[,-1])
  colnames(table) <- names
  table <- as.data.frame(table)
  table$algorithm <- algorithm
  return(table)
}

HC_wardD2 <- tidyclusters('/Users/mariajimenezramos/Desktop/Google Drive/Documents/PhD/Second year/Drug network/Cluster analyses/Results_HC_wardD2.csv', 'HC_wardD2')
HC_average <- tidyclusters('/Users/mariajimenezramos/Desktop/Google Drive/Documents/PhD/Second year/Drug network/Cluster analyses/Results_HCaverage.csv', 'HC_average')
PAM <- tidyclusters('/Users/mariajimenezramos/Desktop/Google Drive/Documents/PhD/Second year/Drug network/Cluster analyses/Results_PAM.csv', 'PAM')
FCM <- tidyclusters("/Volumes/GoogleDrive/My Drive/Documents/PhD/Second year/Drug network/Cluster analyses/Results_FCM.csv", "FCM")
ClOne <- tidyclusters("/Volumes/GoogleDrive/My Drive/Documents/PhD/Second year/Drug network/Cluster analyses/Results_CLOne.csv", "ClOne")
HC_complete <- tidyclusters("/Volumes/GoogleDrive/My Drive/Documents/PhD/Second year/Drug network/Cluster analyses/Results_HC_complete.csv", "HC_complete")
AP <- tidyclusters("/Volumes/GoogleDrive/My Drive/Documents/PhD/Second year/Drug network/Cluster analyses/Results_AP.csv", "AP")

algorithm_nf <- rbind(HC_wardD2, HC_average, PAM, HC_complete, AP)
comparison <- data.frame()
list_algorithm <- list(HC_average, PAM, HC_wardD2, HC_complete, AP)
for (list in list_algorithm){
  subset <- list[, c('K', 'EnrichedClusters', 'EnrichedClustersDT', 'algorithm')]
  comparison <- rbind(comparison, subset)
}

k_cluster <- intersect(HC_average$K, PAM$K)
k_cluster <- intersect(HC_wardD2$K, k_cluster)
k_cluster <- intersect(HC_complete$K, k_cluster)
#k_cluster <- intersect(FCM$K, k_cluster) #Not selected due to poor performance
k_cluster <- intersect(PAM$K, k_cluster)
#k_cluster <- intersect(ClOne$K, k_cluster) #Not selected due to poor performance
k_cluster <- intersect(AP$K, k_cluster)

median_enrichment <- function(k){
  enriched_ATC <- median(c(comparison$EnrichedClusters[which(comparison$K==k)], 1))
  enriched_DT <- median(c(comparison$EnrichedClustersDT[which(comparison$K==k)], 1))
  return(list(ATC=enriched_ATC, DT=enriched_DT))
}

values <- as.data.frame(t(as.data.frame(sapply(k_cluster, median_enrichment))))
medians <- data.frame(K = k_cluster, 'Third level ATC' = as.numeric(values$ATC), 'Mode of Action'= as.numeric(values$DT))
subset <- reshape::melt(medians, id = 'K')

ggplot(data=subset, aes(x=K, y=value, colour=variable)) +
  geom_line(size=1.08)+
  labs(x="Number of Clusters",y = "Number of Enriched Communities")+
  xlim(0,400)+
  labs(colour = "")+
  theme_bw()+
  scale_color_manual(values = c("cornflowerblue", "gold"), labels = c("ATC Code (3rd level)", "Mechanism of Action"))

#Function for plotting the ATC clusters 
plotclusters <- function(alg_comparison, comparison, ylab, filename){
  percentage_ATC <- alg_comparison[,c("K",comparison,"algorithm")]
  percentage_ATC <- melt(percentage_ATC,id=c("algorithm","K"))
  ggplot(data=percentage_ATC, aes(x=K, y=value, colour=algorithm)) + 
    geom_line(size=1.08)+ 
    #geom_vline(xintercept=##best number of communities,linetype="dashed")+
    labs(x="Number of Clusters", y = ylab, colour = '')+
    theme_bw()+
    xlim(0,400)+
    colScale
}

plotclusters(algorithm_nf, 'ARIndex', 'Adjusted Rand Index', 'ARIndexNonfiltered')

obtainEnrichedSet <- function(table,attribute,reference_table,level=4,signific=0.05){
  if (attribute == "Codes"){
    reference_table$Codes <- substring(reference_table$Codes,1,level)
    reference_table <- unique(reference_table)
  }
  results <- enrichmentClustering(table,reference_table,attribute)
  enriched <- results[which(results$Pvalues<signific),][,attribute]
  enriched <- as.character(unique(enriched))
  return(enriched)
}

heatmap_function <- function(enriched,enrPAM,enrHC_complete, enrHC_average, enrHC_wardD2, enrClOne, enrAP, Ref,level,filtered, filename){
  df <- data.frame(matrix(0, nrow = 6, ncol = length(enriched)))
  colnames(df) <- enriched
  rownames(df) <- c("PAM","HC_average", "HC_complete", "HC_wardD2", "ClOne", "AP")
  algorithm <- c("PAM","HC_average", "HC_complete", "HC_wardD2", "ClOne", "AP")
  enr <- list(PAM=enrPAM,HC_complete=enrHC_complete, HC_wardD2 = enrHC_wardD2, HC_average = enrHC_average, ClOne=enrClOne, AP= enrAP)
  for (i in c(1:6)){
    df[which(rownames(df)==algorithm[i]),enr[[i]]]=1
  }
  png(paste("/Users/mariajimenezramos/Desktop/Google Drive/Documents/PhD/Second year/Drug network/Cluster analyses/",filename,".png",sep=""),width = 1000, height = 500)
  if (Ref=="ATC"){
    pheatmap(as.matrix(df),color=c("white","grey37"), border_color = "grey", legend = FALSE, show_colnames = FALSE,fontsize_row=16)
    dev.off()
  }else if (Ref=="MoA"){
    pheatmap(as.matrix(df),color=c("white","grey37"), border_color = "grey", legend = FALSE, show_colnames = FALSE,fontsize_row=16)
    dev.off()
  }
}

heatmaps_algorithm <- function(PAM,HC_complete, HC_average, HC_wardD2 ,ClOne, AP, Ref,level,filtered, filename){
  if (Ref=="ATC"){
    enrPAM <- obtainEnrichedSet(PAM,"Codes",merged_ATC,level)
    enrHC_complete <- obtainEnrichedSet(HC_complete,"Codes",merged_ATC,level)
    enrHC_average <- obtainEnrichedSet(HC_average,"Codes",merged_ATC,level)
    enrHC_wardD2 <- obtainEnrichedSet(HC_wardD2,"Codes",merged_ATC,level)
    enrClOne <- obtainEnrichedSet(ClOne,"Codes",merged_ATC,level)
    enrAP <- obtainEnrichedSet(AP,"Codes",merged_ATC,level)
    ATC <- Reduce(union, list(enrHC_complete, enrPAM,enrHC_average, enrHC_wardD2, enrClOne, enrAP))
    heatmap_function(ATC,enrHC_complete, enrPAM, enrHC_average, enrHC_wardD2, enrClOne, enrAP, Ref,level,filtered,filename)
  } else if (Ref=="MoA"){
    enrPAM <- obtainEnrichedSet(PAM,"MoA",Drug_Targets,level)
    enrHC_complete <- obtainEnrichedSet(HC_complete,"MoA",Drug_Targets,level)
    enrHC_average <- obtainEnrichedSet(HC_average, "MoA",Drug_Targets,level)
    enrHC_wardD2 <- obtainEnrichedSet(HC_wardD2,"MoA",Drug_Targets,level)
    enrClOne <- obtainEnrichedSet(ClOne, "MoA",Drug_Targets,level)
    enrAP <- obtainEnrichedSet(AP,"MoA",Drug_Targets,level)
    MoA <- Reduce(union, list(enrHC_complete, enrPAM,enrHC_average, enrHC_wardD2, enrClOne, enrAP))
    heatmap_function(MoA,enrHC_complete, enrPAM,enrHC_average, enrHC_wardD2, enrClOne, enrAP, Ref,level,filtered,filename)
  }
  return(TRUE)
}

pam_alg <- obtainCommTablePAM(260)
hc_wardD2_alg<- obtainCommTableHC(260, "ward.D2")
hc_complete_alg <- obtainCommTableHC(260, "complete")
hc_average_alg <- obtainCommTableHC(260, "average")
clone_alg <- obtainCommTableClONE(0.3425)
#fcm_alg <- obtainCommTableFCM()
ap_alg <- obtainCommTableAP(0.0644)

heatmaps_algorithm(pam_alg,hc_complete_alg, hc_wardD2_alg, hc_average_alg, clone_alg, ap_alg, "ATC",4,FALSE, "Enriched_ATC")
heatmaps_algorithm(pam_alg,hc_complete_alg, hc_wardD2_alg, hc_average_alg, clone_alg, ap_alg, "MoA",4,FALSE, "Enriched_MoA")

distribution_table <- function(communities,algorithm){
  breaks <- c(2,5,10,15,20,25,50,100)
  categories <- c("2-4","5-9","10-14","15-19","20-24","25-49","50-99","+100")
  df <- data.frame(algorithm=rep(algorithm,8),Category=categories,Frequency=numeric(8))
  for (i in c(1:(length(breaks)-1))){
    df$Frequency[i] <- length(which(table(communities$Cluster)>=breaks[i] & table(communities$Cluster)<breaks[i+1]))
  }
  df$Frequency[i+1] <- length(which(table(communities$Cluster)>=breaks[i+1]))
  return(df)
}

distHC_complete <- distribution_table(hc_complete_alg,"HC_complete")
distHC_average <- distribution_table(hc_average_alg,"HC_average")
distHC_wardD2 <- distribution_table(hc_wardD2_alg,"HC_wardD2")
distPAM <- distribution_table(pam_alg,"PAM")
distClOne <- distribution_table(clone_alg,"ClOne")
distAP <- distribution_table(ap_alg, "AP")
distributions <- rbind(distHC_complete,distHC_average, distHC_wardD2, distPAM,distClOne, distAP)

categories=c("2-4","5-9","10-14","15-19","20-24","25-49","50-99","+100")

ggplot(distributions, aes(x=Category,y=Frequency, fill=algorithm)) +
  geom_bar(stat="identity",position=position_dodge())+
  colfillScale+ 
  scale_x_discrete(limits = categories)+
  labs(x="Number of Drugs",y="Frequency of clusters",fill="")+theme_bw()

compareEnrichedClustersRandom <- function(algorithm,comparison,ylab){
  if(comparison =="ATC"){
    subset <- algorithm[,c("K","EnrichedClusters","algorithm")]}
  else if (comparison == "DT"){
    subset<- algorithm[,c("K","EnrichedClustersDT","algorithm")] 
  }
  subset <- melt(subset,id=c("algorithm","K"))
  p <- ggplot(data=subset, aes(x=K, y=value,colour=algorithm))+
    geom_line(size=1.3,aes(color=algorithm))+
    labs(x="Number of Clusters",y = ylab)+
    labs(colour = "")+
    xlim(0,400)+
    ylim(0,130)+
    colScale+ 
    geom_vline(xintercept=260,linetype="dashed") + theme_bw()
  return(p)
}

algorithm_nf <- rbind(HC_wardD2, HC_average, PAM, HC_complete, ClOne, AP, FCM)
compareEnrichedClustersRandom(algorithm_nf,"ATC","Enriched Clusters (ATC -3rd level)") 
compareEnrichedClustersRandom(algorithm_nf,"DT","Enriched Clusters (MoA)")
