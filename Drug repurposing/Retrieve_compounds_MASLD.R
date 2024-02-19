

#Upload code from DvD to make certain modifications to make the code more functional and solve issues
geneinfo_matrix <- data.table::fread("../Data/geneinfo_beta.txt") #Upload genetic information 
ap_communities <- read.csv("../Data/AP_clusters.csv")
PRLs <- read.csv("../Data/DrugPRLs.csv")

#Function to rank a matrix. It takes as input the matrix to be ranked, the dimension where the ranking is going to be applied to in the matrix and boolean value -TRUE or FALSE- indicating whether the matrix should be ranked in decreasing or ascending order.
rankmat <- function(data, ref, decreasing){
  if(is.matrix(data)){
    rankmatrix <- apply(data, ref, function(x){names(x)[sort.list(as.integer(x), decreasing = decreasing, method = "radix")]})
  }else{
    rankedmat <- names(sort(data,decreasing = decreasing))
  }
  return(rankmatrix) 
}

#Function to perform the Kolmogorov-Smirnov statistics that is going to take as input a ranked matrix obtained from random profiles for permutation tests, the reference ranked matrix that is going to be compared to and the length of the profile.
KS <- function(refmatrix, inputmatrix, lengthtest = 100){
  refmatrix <- as.matrix(refmatrix)
  inputmatrix <- as.matrix(inputmatrix)
  N <- lengthtest
  Nh <- nrow(refmatrix)
  nototest <- ncol(inputmatrix)
  outmat <- c()
  for (i in 1:nototest){
    toptail <- inputmatrix[1:lengthtest, i]
    toptoref <- refmatrix%in%toptail
    In <- matrix(0, nrow = nrow(refmatrix), ncol = ncol(refmatrix))
    Out <- matrix(0, nrow = nrow(refmatrix), ncol = ncol(refmatrix))
    In[toptoref] <- 1/N
    Out[!toptoref] <- 1/(Nh-N)
    if(ncol(refmatrix)==1){
      In <- as.matrix(In, ncol = 1)
      Out <- as.matrix(Out, ncol = 1)
    }
    In <- apply(In, 2, cumsum)
    Out <- apply(Out, 2, cumsum)
    dif <- In-Out
    refvalues <- apply(dif, 2, function(x){x[which.max(abs(x))]})
    outmat <- rbind(outmat, refvalues)
  }
  rownames(outmat) <- colnames(inputmatrix)
  colnames(outmat) <- colnames(refmatrix)
  return(outmat)
}

#This funtion is going to perform a weight signed statistics. The input is the same as the function before.
WSS <- function(refmatrix, inputmatrix, uplength = 100, downlength = 100){
  refmatrix <- as.matrix(refmatrix)
  inputmatrix <- as.matrix(inputmatrix)
  N <- nrow(refmatrix)
  m <- uplength+downlength
  inputranksignmat <- apply(inputmatrix,2, function(x){rank(abs(x))})
  refranksignmat <- apply(refmatrix, 2, function(x){rank(abs(x))})
  name_input <- rankmat(inputranksignmat,2, TRUE)
  geneset <- as.matrix(name_input[1:m,])
  stat <- matrix(nrow = ncol(inputmatrix), ncol = ncol(refmatrix))
  for (i in 1:ncol(inputmatrix)){
    inputrank <- rank(abs(inputmatrix[geneset[,i],i]))
    inputvals <- inputrank*sign(inputmatrix[geneset[,i],i])
    refvals <- refranksignmat[geneset[,i],]*sign(refmatrix[geneset[,i],])
    statvec <- refvals*inputvals
    stat[i,] <- apply(statvec, 2, sum)
  }
  max_score <- sum((N-1:m+1)*(m-1:m+1))
  stat <- stat/max_score
  return(stat*-1)
}

#This function will provide the ES scores.
getES_scores <- function(diseasematup, diseasematdown, uplength = 100, downlength = 100, PRLs, stat = c("KS", "WSS")){
  if(stat=="KS"){
    PRLs <- abs(PRLs)
    ranknamePRLs <- rankmat(PRLs, 2, TRUE)
    ranknamePRLs_inverse <- rankmat(PRLs, 2, FALSE)
    up_profile <- KS(refmatrix = ranknamePRLs, inputmatrix = diseasematup, lengthtest = uplength)
    up_ref <- KS(refmatrix = diseasematup, inputmatrix= ranknamePRLs, lengthtest= uplength)
    up_ref <- t(up_ref)
    down_profile <- KS(refmatrix = ranknamePRLs, inputmatrix = diseasematdown, lengthtest = downlength)
    down_ref <- KS(refmatrix = diseasematup, inputmatrix = ranknamePRLs_inverse, lengthtest = downlength)
    down_ref <- t(down_ref)
    
    score1 <- (up_profile - down_profile)/2
    score2 <- (up_ref - down_ref)/2
    check1 <- sign(up_profile) == sign(down_profile)
    check2 <- sign(up_ref) == sign(down_ref)
    score1[check1] <- 0
    score2[check2] <- 0
    scores <- score1+score2
    all <- as.matrix(scores/2)
  }else{
    all <- WSS(refmatrix = PRLs, inputmatrix = diseasematup, uplength = uplength, downlength = downlength)
    all <- as.matrix(all)
  }
  if (is.matrix(diseasematup)){
    rownames(all) <- colnames(diseasematup)
  }
  colnames(all) <- colnames(PRLs)
  return(all)
}

#This function identifies the significant compounds according to either KS or WSS function.
findSignifCompounds <- function(scores, rankdata = NULL, rankdata_inverse = NULL, uplength = 100, downlength = 100, stat, nperm = 100, signif.fdr){
  rankdata <- as.matrix(rankdata)
  size <- nrow(rankdata)
  smatrix <- sapply(1:nperm, function(x){sample.int(size)})
  rownames(smatrix) <- rankdata[,1]
  samplemat <- rankmat(smatrix, 2, TRUE)
  samplemat_inverse <- rankmat(smatrix, 2, FALSE)
  if (stat == "KS"){
    permup <- KS(refmatrix = samplemat, inputmatrix = rankdata, lengthtest = uplength)
    permupref <- KS(refmatrix = rankdata, inputmatrix = samplemat, lengthtest = uplength)
    permupref <- t(permupref)
    permdown <- KS(refmatrix = samplemat, inputmatrix = rankdata_inverse, lengthtest = downlength)
    permdownref <- KS(refmatrix = rankdata, inputmatrix = samplemat_inverse, lengthtest = downlength)
    permdownref <- t(permdownref)
    score1 <- (permup-permdown)/2
    score2 <- (permupref-permdownref)/2
    check1 <- sign(permup) == sign(permdown)
    check2 <- sign(permupref) == sign(permdownref)
    score1[check1] <- 0
    score2[check2] <- 0
    Pscores <- score1+score2
    Pscores <- Pscores/2
  }else{
    samplesign <- sapply(1:nperm, function(x){sample(c(-1,1), size = size, replace =TRUE)})
    samplemat <- smatrix*samplesign
    rownames(samplemat) <- rownames(rankdata)
    Pscores <- WSS(refmatrix = samplemat, inputmatrix = rankdata, uplength = uplength, downlength = downlength)
  }
  empcdf <- ecdf(abs(Pscores))
  emppvals <- vector(length=length(scores))
  emppvals <- 1-empcdf(abs(scores))
  pvaladj <- p.adjust(emppvals, "BH")
  sel <- which(pvaladj < signif.fdr) #SI CAMBIO ESTO PUEDO CONSEGUIR TODOS LOS COMPUESTOS INDEPENDIENTEMENTE DEL PVAL?????
  signifcompounds <- scores[sel]
  names(signifcompounds) <- names(scores)[sel]
  return(signifcompounds)
}

#This function gets permutation based empirical p-values for enrichment scores.	
significantES <- function(ESscores, rankdata = NULL, rankdata_inverse = NULL, uplength = 250, downlength = 250, stat, nperm, signif.fdr){
  if (is.list(ESscores)){
    results <- list(length = length(ESscores))
    for (i in 1:length(ESscores)){
      scorematrix <- ESscores[[i]]
      if(is.matrix(scorematrix)){
        signifcompounds <- list(length = nrow(scorematrix))
        for (i in 1:nrow(scorematrix)){
          signifcompounds[[i]] <- findSignifCompounds(scorematrix[i,], rankdata = rankdata[,i], rankdata_inverse = rankdata_inverse[,i], uplength = uplength, downlength = downlength, stat = stat, nperm = nperm, signif.fdr = signif.fdr)
        }
        names(signifcompounds) <- rownames(scorematrix)
      }else{
        signifcompounds <- findSignifCompounds(scorematrix, rankdata = rankdata, rankdata_inverse = rankdata_inverse, uplength = uplength, downlength =downlength, stat = stat, nperm = nperm, signif.fdr = signif.fdr)
      }
      results[[i]] <- signifcompounds
    }
    names(results) <- names(ESscores)
  }else{
    scorematrix <- ESscores
    if (is.matrix(scorematrix)){
      signifcompounds <- list(length = nrow(scorematrix))
      for (i in 1:nrow(scorematrix)){
        signifcompounds[[i]] <- findSignifCompounds(scorematrix[i,], rankdata = rankdata[,i], rankdata_inverse = rankdata_inverse[,i], uplength = uplength, downlength = downlength, nperm = nperm, stat = stat, signif.fdr = signif.fdr)
      }
      names(signifcompounds) <- rownames(scorematrix)
    }else{
      signifcompounds <- list(length = 1)
      signifcompounds[[1]] <- findSignifCompounds(scorematrix, rankdata = rankdata, rankdata_inverse = rankdata_inverse, uplength = uplength, downlength = downlength, stat= stat, signif.fdr = signif.fdr)
      names(signifcompounds) <- names(ESscores)
    }
    results <- signifcompounds
  }
  return(results)
}


#We create another function to calculate the ES score from an input profile. The input is the following:
  
#  *`data`: matrix gene expression profiles, where the rows are the genes and columns are different profiles.
#  *`lengthtest`: integer giving the number of genes that will be used to generate the set size to use for look for drug profiles.
#  *`signif.fdr`: false discovery rate that will be used to determine the significance of enrichment scores.
#  *`nperm`: integer for the number of permutation profiles.
#  *`PRLs`: the PRLs lists that contain the drug profiles  -not ranked, not distance.
#  *`stat`: either Kolmogorov-Smirnov (KS) or Weigthed Signed Ranked Score (WSS). The first one equally weights all elements of the gene set, while the second one uses both the sign and the position in the ranked list to calculate the enrichment scores.

# The output is a list with the ES scores and their significance.

calculateES <- function(data, PRLs, stat=c("WSS", "KS"), uplength = uplength, downlength =downlength, signif.fdr, nperm, pvalues = NULL, type= c("fixed", "dynamic"), dynamic.fdr = 0.05){
  if(stat == "WSS"){
    ngenes <- nrow(PRLs)
    splitn <- round(ngenes/2)+1
    refDB <- splitn - PRLs
  }else{
    refDB <- PRLs
  }
  refDB <- as.matrix(refDB)
  data <- as.matrix(data)
  lt<-sum(rownames(data)%in%rownames(refDB))
  if (lt!=nrow(data)){
    stop('Rownames of input data and reference data do not match')}
  refDB <- refDB[rownames(data),]
  if(stat=="KS"){
    datarl <- apply(data, 2, function(x){as.integer(rank(x))})
    datarl <- as.matrix(datarl)
    rownames(datarl) <- rownames(data)
    colnames(datarl) <- colnames(data)
    rankdata <- rankmat(datarl, 2, TRUE)
    rankdata_inverse <- rankmat(datarl, 2, FALSE)
  }else{
    rankdata <- data
    rankdata_inverse <- data
  }
  if (type == "fixed"){
    outk <- getES_scores(rankdata, rankdata_inverse, PRLs = refDB, uplength = uplength, downlength = downlength, stat = stat)
    signif <- significantES(ESscores = outk, rankdata = rankdata, rankdata_inverse = rankdata_inverse, signif.fdr = signif.fdr, uplength = uplength, downlength = downlength, stat = stat, nperm= nperm)
    names(signif) <- colnames(data)
  }else{
    pvals <- as.matrix(pvalues)
    pvals <- apply(pvals, 2, p.adjust)
    signifgenes <- apply(pvals, 2, function(x){x<dynamic.fdr})
    rankdrugPRLs <- rankmat(refDB, 2, TRUE)
    rankdrugPRLs_inverse <- rankmat(refDB, 2, FALSE)
    if (ncol(pvals)>1){
      signif <- list(length = ncol(pvals))
      outk <- list(length = ncol(pvals))
      postmat <- apply(data,2,function(x){x>0})
      upsignif <- postmat*signifgenes
      for (i in 1:ncol(pvals)){
        if (sum(signifgenes[,i])>0){
          uplength <- sum(upsignif[,i])
          downlength <- sum(signifgenes[,i])-uplength
          if (uplength <50){
            warning("Number of significantly up-regulated genes is less than 50, top 50 used instead")
            uplength<-50
          }
          if (downlength <50){
            warning("Number of significantly down-regulated genes is less than 50, top 50 used instead")
            downlength<-50
          }
          outk[[i]] <- getES_scores(rankdata[,i], rankdata_inverse[,i], uplength = uplength, downlength = downlength, PRLs= refDB, stat = stat)
          signif[[i]] <- significantES(ESscores = outk[[i]], uplength = uplength, downlength= downlength, rankdata = rankdata, rankdata_inverse = rankdata_inverse, nperm =nperm, signif.fdr = signif.fdr, stat = stat)
        }else{
          warning(paste("No significantly differentially expressed genes for experiment",colnames(pvalues)[i]))
          outk[[i]]<-NULL
          signif[[i]]<-NULL
        }
      }
      names(outk) <- colnames(data)
      names(signif) <- colnames(data)
    }else{
      if (sum(signifgenes)>0){
        postmat <- data >0
        upsignif <- postmat*signifgenes
        uplength <- sum(upsignif)
        downlegnth <- sum(signifgenes) - uplength
        if (uplength <50){
          warning("Number of significantly up-regulated genes is less than 50, top 50 used instead")
          uplength<-50
        }
        if (downlength < 50){
          warning("Number of significantly down-regulated genes is less than 50, top 50 used instead")
          downlength<-50
        }
        outk <- getESscores(rankdata, rankdata_inverse, PRLs = refDB, uplength = uplength, downlength = downlength, stat = stat)
        signif <- signficantES(ESscores = ouk, rankdata= rankdata, rankdata_inverse = rankdata_inverse, nperm = nperm, signif.fdr = signif.fdr, uplength = uplength, downlength =downlength, stat = stat)
      }else{
        stop(paste("No significantly differentially expressed genes, consider changing FDR threshold"))
      }
      names(signif) <- colnames(data)
      names(outk) <- colnames(data)
    }
  }
  return(list(scores = outk, significance = signif))
}

#The final drugs will be clustered together by a single linkage cluster following DrugVsDisease package. Two functions have been used to to achieve this. The input is the following:
  
#  *`ESscores`: enrichment scores obtaind by comparing the input profile to the drug network.
#  *`statistics`: whether the median or the mean is used to perform the average cluster. 
#  *`drugClusters`: data frame that contains the different drugs and the cluster where it belongs.
#  *`stat`: whether WSS or KS is used to get the enrichment scores.

#  `classifysingleprofile` output is a data frame that will contained the drugs with the ES scores, the cluster and the RPS.

findCluster <-
  function(SignifScores, drugClusters=NULL, no.signif=30){
    #check that we have some significant results otherwise return na
    if (length(SignifScores)==0){return(NA)}
    if (is.character(SignifScores)){return(NA)}
    #check if there are more significant results than have been asked to assign to a cluster
    if(length(SignifScores) > no.signif){cat(paste('Number of Significant results greater than', no.signif, 'Using top',no.signif, 'hits - consider using average linkage instead'),fill=TRUE)}
    #order by most significant results (i.e. the highest absolute Enrichment scores) and take the top no.signif to classify
    absscores <- abs(SignifScores)
    SortScores <- sort(absscores,TRUE)
    SScores <- SignifScores[names(SortScores)]
    SignifScores <- SScores[1:no.signif]
    signifnames <-names(SignifScores)
    #find the cluster for the significant matches, and output to data frame along with the distance between profiles (1-ES), and name of the profile.
    sel <- sapply(signifnames,function(x){which(drugClusters[,"Drug"]==x)})
    sel <- unlist(sel)
    scores <- SignifScores[names(sel)]
    distscores<-1-abs(scores)
    cors<-ifelse(scores>0,-1,1)
    results.clust<- data.frame(names(sel), distscores, drugClusters[sel,"Cluster"], cors)
    colnames(results.clust)<-c("Drug","ES Distance","Cluster","RPS")
    return(results.clust)
  }

classifysinglelinkage <-
  function(significance,drugClusters=NULL,no.signif=30){
    #significance output from significantES
    #if length(significance)>1 then we have dynamic set sizes for the KS scor
    Clusterlist<-list(length=length(significance))
    for(i in 1:length(significance)){		
      SignifScores <- significance[[i]]
      #see if SignifScores is a list - i.e. more than one disease profile
      if(is.list(SignifScores)){
        Clust<-lapply(SignifScores, function(x){findCluster(x, drugClusters = drugClusters, no.signif = 10)})
        names(Clust) <- names(SignifScores)
      }else{
        #just have one profile:
        Clust<- findCluster(SignifScores,drugClusters= drugClusters,no.signif=no.signif)
      }
      Clusterlist[[i]] <- Clust		
    }
    names(Clusterlist) <- names(significance)
    return(Clusterlist)
  }


# We can finally now create the final function that will give us the drug profiles that could be used as treatment for our ranked gene profile. The input of the function is:
  
#  *`data`: matrix gene expression profiles, where the rows are the genes and columns are different profiles.
#  *`lengthtest`: integer giving the number of genes that will be used to generate the set size to use for look for drug profiles.
#  *`signif.fdr`: false discovery rate that will be used to determine the significance of enrichment scores.
#  *`nperm`: integer for the number of permutation profiles.
#  *`avgstat`: character string that indicates wheter the median or the mean is used to cluster the drugs in the average linked cluster.
#  *`no.signif`: integer giving the number of significant enrichment scores to return. Default is 10.
#  *`stat`: either Kolmogorov-Smirnov (KS) or Weigthed Signed Ranked Score (WSS). The first one equally weights all elements of the gene set, while the second one uses both the sign and the position in the ranked list to calculate the enrichment scores.

# The output is a data frame that will contain the names of the profiles with significant scores, the distance between the input profile and the reference profile (1- ES), the cluster number of the node.

classifyprofile <- function(data, pvalues = NULL, PRLs, type = c("fixed", "dynamic"), uplength = 100, downlength = 100, signif.fdr=0.05, dynamic.fdr = 0.05, nperm=1000, avgstat=c("mean","median"), no.signif=30, stat=c("KS","WSS"), drugClusters){
  if(!is.matrix(data)){
    data<-read.table(data,header=TRUE,row.names=1)
    data<-as.matrix(data)
  }
  if(!is.null(pvalues)&&!is.matrix(pvalues)){
    pvalues<-read.table(pvalues,header=TRUE,row.names=1)
    pvalues<-as.matrix(pvalues)
  }
  if(!is.matrix(PRLs)){PRLs<-read.table(PRLs,header=TRUE,row.names=1)}
  if(!is.data.frame(drugClusters)){drugClusters<-read.table(drugClusters,header=TRUE)}
  refnames<-colnames(PRLs)
  refnames<-make.names(refnames)
  clustnames<-drugClusters[,which(colnames(drugClusters)%in%c("Drug","Disease"))]
  clustnames<-make.names(clustnames)
  intersection<-which(refnames%in%clustnames)
  if(length(refnames)!=length(intersection)){stop('Profiles in custom rank profiles and clustom clusters do not match')}
  #calculate enrichment scores
  ESvals<- calculateES(data,PRLs,type = type, dynamic.fdr = dynamic.fdr, signif.fdr=signif.fdr,nperm=nperm,uplength=uplength, downlength = downlength,stat=stat)
  #assign significant enrichment scores to clusters	
  clusterassignments<- classifysinglelinkage(ESvals$significance, drugClusters= drugClusters,no.signif)
  #check have at least one significant results
  nacheck<-unlist(clusterassignments)
  if(is.list(nacheck)){
    nacheck2<-unlist(nacheck)
    nas<-which(is.na(nacheck2))
    if(length(nas)==length(nacheck2)){warning('No significant matches found. Consider changing the FDR threshold')}
  }else{
    nas<-which(is.na(nacheck))
    if(length(nas)==length(nacheck)){warning('No significant matches found. Consider changing the FDR threshold')}
  }
  return(clusterassignments)
}

# Next the transcriptomic profiles are uploaded as steatosis, f0f1, f2, f3, f4. These are csv files of the results obtained when performing differential gene expression analysis. Those genes according to |logFC| = fc and adj. p-val < fdr in each data frame are selected.
fc = 1
fdr = 0.05

signatures <- list(steatosis, f0f1, f2, f3, f4)

for (i in signatures){
  i <- i[which(abs(i$logFC) >= fc & i$adj.P.Val < fdr),]
  i <- i[order(i$logFC, decreasing = TRUE),]
}

#Select those genes that are present in the gene_info matrix and retrieve the gene name provided in clue.io
test <- list()
for (i in 1:5){
  genes_list <- geneinfo_matrix[match(signatures[[i]]$genes, geneinfo_matrix$gene_symbol),]
  genes_list$logFc <- signatures[[i]]$logFC 
  genes_list <- genes_list[which(!is.na(genes_list$gene_id)),]
  test[[i]] <- as.matrix(genes_list$logFc)
  rownames(test[[i]]) <- genes_list$gene_id
}

#Obtain drugs with the KS signature
steatosis_drugs_ks <- classifyprofile(test[[1]], PRLs = PRLs, type = "fixed", uplength = 44, downlength = 100, signif.fdr = 0.05, nperm = 1000, avgstat = "mean", no.signif = 30, stat = "KS", drugClusters = ap_communities,)
f0f1_drugs_ks <- classifyprofile(test[[2]], PRLs = PRLs, type = "fixed", uplength = 56, downlength = 80, signif.fdr = 0.05, nperm = 1000, avgstat = "mean", no.signif = 30, stat = "KS", drugClusters = ap_communities)
f2_drugs_ks <- classifyprofile(test[[3]], PRLs = PRLs, type = "fixed", uplength = 100, downlength = 91, signif.fdr = 0.05, nperm = 1000, avgstat = "mean", no.signif = 30, stat = "KS", drugClusters = ap_communities)
f3_drugs_ks <- classifyprofile(test[[4]], PRLs = PRLs, type = "fixed", uplength = 100, downlength = 100, signif.fdr = 0.05, nperm = 1000, avgstat = "mean", no.signif = 30, stat = "KS", drugClusters = ap_communities)
f4_drugs_ks <- classifyprofile(test[[5]], PRLs = PRLs, type = "fixed", uplength = 100, downlength = 100, signif.fdr = 0.05, nperm = 1000, avgstat = "mean", no.signif = 30, stat = "KS", drugClusters = ap_communities)

#Obtain drugs with the WSR signature
steatosis_drugs_wsr <- classifyprofile(test[[1]], PRLs = PRLs, type = "fixed", uplength = 44, downlength = 100, signif.fdr = 0.05, nperm = 1000, avgstat = "mean", no.signif = 30, stat = "WSR", drugClusters = ap_communities)
f0f1_drugs_wsr <- classifyprofile(test[[2]], PRLs = PRLs, type = "fixed", uplength = 56, downlength = 80, signif.fdr = 0.05, nperm = 1000, avgstat = "mean", no.signif = 30, stat = "WSR", drugClusters = ap_communities)
f2_drugs_wsr <- classifyprofile(test[[3]], PRLs = PRLs, type = "fixed", uplength = 100, downlength = 91, signif.fdr = 0.05, nperm = 1000, avgstat = "mean", no.signif = 30, stat = "WSR", drugClusters = ap_communities)
f3_drugs_wsr <- classifyprofile(test[[4]], PRLs = PRLs, type = "fixed", uplength = 100, downlength = 100, signif.fdr = 0.05, nperm = 1000, avgstat = "mean", no.signif = 30, stat = "WSR", drugClusters = ap_communities)
f4_drugs_wsr <- classifyprofile(test[[5]], PRLs = PRLs, type = "fixed", uplength = 100, downlength = 100, signif.fdr = 0.05, nperm = 1000, avgstat = "mean", no.signif = 30, stat = "WSR", drugClusters = ap_communities)

#Intersect between the drug signatures
steatosis_drugs_finals <- intersect(steatosis_drugs_ks[[1]]$Drug, steatosis_drugs_wsr[[1]]$Drug)
f0f1_drugs_finals <- intersect(f0f1_drugs_ks[[1]]$Drug, f0f1_drugs_wsr[[1]]$Drug)
f2_drugs_finals <- intersect(f2_drugs_ks[[1]]$Drug, f2_drugs_wsr[[1]]$Drug)
f3_drugs_finals <- intersect(f3_drugs_ks[[1]]$Drug, f3_drugs_wsr[[1]]$Drug)
f4_drugs_finals <- intersect(f4_drugs_ks[[1]]$Drug, f4_drugs_wsr[[1]]$Drug)

