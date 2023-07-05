library(cmapR)

#Upload info from CMap2020
ds_path <- "../Data/level5_beta_trt_cp_n720216x12328.gctx"
siginfo_matrix <- data.table::fread("../Data/siginfo_beta.txt")
geneinfo_matrix <- data.table::fread("../Data/geneinfo_beta.txt")
compoundinfo_matrix <- data.table::fread("../Data/compoundinfo_beta.txt")

#Obtain the gene expression profile (GEP) for each drug observation that has at least three replicates.
siginfo_matrix <- siginfo_matrix[which(siginfo_matrix$pert_type == 'trt_cp')] #We select only those perturbations that are compounds.
##drug_gctx <- parse_gctx(ds_path, cid = test$sig_id, rid = drug_matrix <- mat(drug_gctx)
##rank_matrix <- apply(-drug_matrix, 2, rank, ties.method = 'first')
drugs_nsamples <- siginfo_matrix[siginfo_matrix$nsample >= 3] #We choose only those compounds that have at least 3 replicates
#freq <- table(drugs_nsamples$cmap_name)
#drugs_nsamples <- drugs_nsamples[which(drugs_nsamples$cmap_name %in% names(freq[freq >= 3]))]
drugs_nsamples$cmap_name[drugs_nsamples$cmap_name == "cdk/crk-inhibitor"] <- "cdk_crk-inhibitor" #Change name so it doesn't affect the creation of the files in the loop
#drugs_nsamples <- drugs_nsamples[!grep("^BRD", drugs_nsamples$cmap_name)] #Remove those samples that are BRD compounds
for (drug in unique(siginfo_matrix$cmap_name)){  #33609 drugs in 45.5 GB
  drug_instances <- drugs_nsamples$sig_id[which(drugs_nsamples$cmap_name==drug)] # These are the signature id that are present in the final drug matrix with all the information per drug.
  drug_gctx <- parse_gctx(ds_path, cid = drug_instances) #Retrieve gene signatures for a specific drug
  rank_matrix <- rank_gct(drug_gctx, dim = "col", decreasing = TRUE) #Convert gct object into ranked matrix in decreasing order
  rank_matrix <- mat(rank_matrix)
  rownames(rank_matrix) <- drug_gctx@rid
  colnames(rank_matrix) <- drug_gctx@cid
  if(length(drug_instances)>1){
    drug_GEPs <- rank_matrix[,drug_instances]
  }else {
    name <- drugs_nsamples$sig_id[drugs_nsamples$sig_id%in%drug_instances]
    drug_GEPs <- data.frame(rank_matrix[,drug_instances],row.names=rownames(rank_matrix))
    colnames(drug_GEPs) <- name
  }
  write.table(drug_GEPs, file=paste('../Data/GEPs for each drug', drug, '.out' ,sep=''), row.names=TRUE,col.names=TRUE)
}

#######################
# Spearman's Footrule #
#######################

#This function will take the ranked matrix obtained from the Clue.io database and select to columns, i and j, which would correspond to the GEP of two different drugs. The absolute distance of these two drugs will be calculated and the final result would correspond to its sum.  

spearmansfootrule_ij<-function(i,j,rank_matrix){
  dist <- sum(abs(rank_matrix[,i]- rank_matrix[,j]))
  return(dist)
}
spearmans_footrule <- Vectorize(spearmansfootrule_ij,vectorize.args=list("i","j")) #The function is vectorized by `Vectorize` for further analyses.

##########################
# Borda Merging Function #
##########################

#The two closer ranked lists will be summed and its result is ranked. 

borda_merging <- function(GEPa,GEPb){
  ps <- GEPa+GEPb
  mergedGEP <- rank(ps,ties.method="first") #In case of ties, the gene that is first (following an alphanumeric order)
  return(mergedGEP)
}

#####################
# Kruskal Algorthim #
#####################

#Kruskal algorithm will merge all the observations for each drug until one single list is obtained. The function first create an array that contains the distance of the every pair of GEPs according to the Spearman's Footrule. 
#Then the minimum distance in the array is obtained and the drugs that belong to that distance (which would be column i and j in the ranked matrix) will be merged according to the Borda Merging function. 
#The list of GEPs is modified so the two lists merged are substituted by the merged one. 

kruskal_algorithm <- function(drug_info){
  n <- length(drug_info)
  while (n>1){
    dist <- outer(1:n,1:n, FUN = spearmans_footrule, rank_matrix = drug_info)
    min_drug <- min(dist[upper.tri(dist)]) # the diagonal is always zero and distance is symmetrical
    i <- which(dist == min_drug, arr.ind = TRUE)[[1,1]]   #arr.ind is to return x if it is array
    j <- which(dist == min_drug, arr.ind = TRUE)[[1,2]]
    y <-borda_merging(drug_info[,i],drug_info[,j])
    drug_info[,i] <- NULL
    drug_info[,j] <- NULL
    drug_info <- cbind(drug_info,y)
    n <- length(drug_info)
  }
  #The final PLR list is rearranged.
  drug_info$genes <- rownames(drug_info)
  drug_info <- drug_info[order(drug_info[,1]),]
  drug_info[,1] <- NULL
  rownames(drug_info) <- NULL
  return(drug_info)
}

##################
# Final PRL list #
##################

#We need to reduce the number of observations per drugs (Spearman's Footrule takes up to 5 min for each drug when there is 100 observations). 
#For this reason, we select up to 100 random observations in those drugs that have more than that.
reduce_n_obs <- function(drug_info, drug){
  n_seed <- match(drug, unique(drugs_nsamples$cmap_name))
  set.seed(n_seed)
  drug_info2 <- drug_info[,sample(ncol(drug_info), size = 100)]
  return(drug_info2)
}

generate_PRL<- function(drug){
  drug_info <- read.table(paste('GEPs for each drug/', drug, '.out', sep = ''), header = TRUE)
  if (length(drug_info) > 100){
    drug_info <- reduce_n_obs(drug_info, drug)
  }
  PRL <- kruskal_algorithm(drug_info)
  drug <- gsub(" ", "_", drug)
  colnames(PRL) <- drug
  return(PRL)
}

PRL_matrix <- function(drugs){
  PRLs <- data.frame(matrix(NA, nrow = 12328, ncol = 0)) #nrow corresponds to the total number of genes available in the matrix
  for (drug in drugs){
    PRL<- generate_PRL(drug)
    PRLs<- cbind(PRLs,PRL)
  }
  PRLs <- PRLs[ , order(colnames(PRLs))]
  return(PRLs)
}

PRLs <- PRL_matrix(unique(drugs_nsamples$cmap_name))
write.table(PRLs,"../Data/DrugPRLs.txt")
