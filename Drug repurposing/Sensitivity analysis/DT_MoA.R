library(XML)
load("./Distance_matrix.R")


#We check the number of enriched clusters according to ATC and MoA
#Retrieve ATC from Bioportal and DrugBank
bioportal <- read.csv('../Data/ATC.csv')
bioportal$Class.ID <- gsub('http://purl.bioontology.org/ontology/UATC/', '', bioportal$Class.ID)
bioportal <- bioportal[nchar(bioportal$Class.ID) == 7, ] #We only keep those ATC codes that have 7 characters, which belong to level 5
bioportal_atc <- bioportal[bioportal$Preferred.Label%in%drugs_nsamples$cmap_name, ] 

drugbank <- xmlParse("../Data/full database.xml")
root <- xmlRoot(drugbank)
df <- data.frame()
for (i in 1:xmlSize(root)){
  a <- xmlToList(root[[i]][["atc-codes"]])
  names(a) <- NULL
  b <- xmlValue(root[[i]][["name"]])
  if(!is.null(a)){
    for (i in 1:length(a)){
      row <- c(a[[i]]$.attrs, b)
      df <- rbind(df, row)}
  }
}

names(df) <- c("Class.ID", "Preferred.Label")
matches <- as.vector(sapply(drugs_nsamples$cmap_name, function(x){matches <- grep(x, df$Preferred.Label, ignore.case = TRUE)}))
drugbank_atc <- melt(matches)
drugbank_atc <- drugbank_atc[which(!duplicated(drugbank_atc)),] #Remove duplicated rows
colnames(drugbank_atc) <- c("row_original", "Drugs")
drugbank_atc <- cbind(drugbank_atc, df[drugbank_atc$row_original,])
atcfinal_drugbank <- data.frame(Drugs = drugbank_atc$Preferred.Label, Codes = drugbank_atc$Class.ID)

##Check which drugs are in common with bioportal
Drugs_atc <- atcfinal_drugbank$Drugs[match(atcfinal_drugbank$Codes, bioportal_atc$Class.ID)]
Drugs_atc <- Drugs_atc[!is.na(Drugs_atc)]

##Add additional atc codes to the list of CMap Drugs.
dbDrugs <- atcfinal_drugbank$Drugs[!atcfinal_drugbank$Codes%in%bioportal_atc$Class.ID]
atcfinal_drugbank[atcfinal_drugbank$Drugs%in%as.character(dbDrugs),]
colnames(bioportal_atc) <- c('Codes', 'Drugs')
bioportal_atc <- bioportal_atc[c(2, 1)]
merged_ATC<- rbind(bioportal_atc, atcfinal_drugbank[atcfinal_drugbank$Drugs%in%as.character(dbDrugs),])
merged_ATC <- merged_ATC[which(!duplicated(merged_ATC)), ] #Remove duplicated rows
write.csv(merged_ATC,"../Data/List_ATCcodes.csv",row.names = FALSE)

#We retrieve the drug targets from CMap2020
Drug_Targets <- data.frame('Drugs' = compoundinfo_matrix$cmap_name, 'MoA' = compoundinfo_matrix$moa, 'Drug_Target' = compoundinfo_matrix$target)
Drug_Targets <- Drug_Targets[which(Drug_Targets$MoA != ''),]
Drug_Targets <- Drug_Targets[which(Drug_Targets$Drug_Target != ''),]
Drug_Targets <- unique(Drug_Targets)
write.csv(Drug_Targets, '../Data/DrugTargets_MoA.csv', row.names = FALSE)