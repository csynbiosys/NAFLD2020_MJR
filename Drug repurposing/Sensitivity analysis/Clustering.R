#Enrichemnt_analysis
load("./Enrichment_analysis.R")

##############
# ClusterOne #
##############

dens <-  seq(0.0125,0.7425,by=0.01)

#Random
results_clone_random <- sapply(dens,enrichedRandomClONE,signific=0.05,level=4, soft = FALSE, distmat = distmat) 
results_clone_random <- as.matrix(results_clone_random)
write.csv(results_clone_random,"../Results/Results_CLOne_random.csv",sep="")

#Not random
results_clone <- sapply(dens,enrichedclustersClONE,signific=0.05,level=4, soft = FALSE, distmat= distmat) 
results_clone <- as.matrix(results_clone)
write.csv(results_clone, "../Results/Results_CLOne.csv")

###########################
# Hierarchical clustering #
###########################

ks=seq(2,400,by=2)

#Random
results_HC_random <- sapply(ks,enrichedRandomHC,signific=0.05,level=4, soft = FALSE, distmat = distmat)
results_HC_random <- as.matrix(results_HC_random)
write.csv(results_HC_random, "../Results/Results_HC_random.csv")

#Not random
results_HC <- sapply(c(2,10),enrichedclustersHC,signific=0.05,level=4, soft = FALSE, distmat = distmat)
results_HC <- as.matrix(results_HC)
write.csv(results_HC, "./Results/Results_HC.csv")

#################
# AP clustering #
#################

ps=c(-30,-10,-5,-2.2,-2,-1.5,-1,-0.75,-0.5,seq(-0.45,0.25,by=0.005))

#Random
results_AP_random <- sapply(ps, enrichedRandomAP, signific= 0.05, level = 4, soft =FALSE, dismtat= distmat)
results_AP_random <<- as.matrix(results_AP_random)
write.csv(results_AP_random, "../Results/Results_AP_random.csv")

#Not random
results_AP <- sapply(ps, enrichedclustersAP, signific = 0.05, level = 4,  soft = FALSE, distmat)
results_AP <- as.matrix(results_AP)
write.csv(results_AP, "../Results/Results/Results_AP.csv")

##################
# PAM clustering #
##################

#Random clustering
ks=seq(2,400,by=2)

results_PAM_random <- sapply(ks,enrichedRandomPAM,signific=0.05,level=l, soft = FALSE, distmat = distmat) 
results_PAM_random <- as.matrix(results_PAM_random)
write.csv(results_PAM_random, "../Results/Results_PAM_random.csv")

#Not random
results_PAM <- sapply(ks,enrichedclustersPAM,signific=0.05,level=4, soft = FALSE, distmat) 
results_PAM <- as.matrix(results_PAM)
write.csv(results_PAM, "../Results/Results_PAM",l,".csv")

##################
# FCM clustering #
##################

#Random
results_FCM_random <- sapply(ks, enrichedRandomFCM, signific= 0.05, level = 4, soft = TRUE, distmat =distmat)
results_FCM_random <<- as.matrix(results_FCM_random)
write.csv(results_FCM_random, "../Results/Results_FCM_random.csv")

#Not random
results_FCM <- sapply(ks, enrichedclustersFCM, signific = 0.05, level = 4, soft =TRUE, distmat= distmat)
results_FCM <- as.matrix(results_FCM)
write.csv(results_FCM, "../Results/Results/Results_FCM.csv")
