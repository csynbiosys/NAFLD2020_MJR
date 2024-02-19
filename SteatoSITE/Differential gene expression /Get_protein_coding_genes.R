##Choose only protein-coding genes and check distributions between both datasets.

library(ensembldb)

#Upload final csv file with counts with genes as rows and patient ID as columns called feature_counts

edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
prot <- genes(edb, filter =GeneBiotypeFilter("protein_coding"), return.type = "data.frame")
prot_counts <- feature_counts[which(rownames(feature_counts)%in%prot$gene_name),]






