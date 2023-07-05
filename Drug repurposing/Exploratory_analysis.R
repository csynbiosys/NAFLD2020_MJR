#Exploratory analyses of the CMap2020 database.

#Upload files from CMap2020
ds_path <- "./Data/level5_beta_trt_cp_n720216x12328.gctx" #Level 5 contains the signatures of each observation 
siginfo_matrix <- data.table::fread("./Data/siginfo_beta.txt") #Matrix with gene signatures
geneinfo_matrix <- data.table::fread("./Data/geneinfo_beta.txt") #Matrix with the genetic information 
compoundinfo_matrix <- data.table::fread("./Data/compoundinfo_beta.txt") #Matrix with the compounds information

#Plot different types of perturbations
pert_type <- ggplot(siginfo_matrix, aes(pert_type))+
  geom_bar(fill = "cornflowerblue") + 
  labs(x = "Pertubation types", y = "Counts")+
  theme_bw()

siginfo_matrix_comp <- siginfo_matrix[siginfo_matrix$pert_type == "trt_cp"] #We are only interested in those pertubations that are compounds, so we will use that for now on.

#Plot number of replicates per compound
cp_nsamples <- ggplot(siginfo_matrix_comp, aes(siginfo_matrix_comp$nsample)) +
  geom_bar(fill = "cornflowerblue") + 
  labs(x = "No of replicates", y  = "No of compounds", title = "No of replicates per compound") + 
  theme_bw() ##36 is the biggest one

#Plot number of cell lines used for each compound
cp_countcells_info <- data.frame(siginfo_matrix_comp$cmap_name, siginfo_matrix_comp$cell_iname, row.names = NULL) #The values we want to summarise is character type, we need to isolate the information to manipulate it later
cp_countcells_info <- lapply(unique(cp_countcells_info), table) #This gives two lists: one summarising compounds per cell line and another one summarising cell lines per compound. We are interested in the latest.
df_countcells <- as.data.frame(cp_countcells_info$siginfo_matrix.siginfo_matrix.pert_type.....trt_cp...cmap_name) 
cp_countcells <- ggplot(df_countcells, aes(Freq)) +
  geom_bar(fill = "cornflowerblue") +
  labs(x = "No of cell lines", y = "No of compounds", title = "No of cell lines per compound")+
  theme_bw()

#Plot different concentrations investigated for each compound
cp_dose_info <- data.frame(siginfo_matrix_comp$cmap_name, siginfo_matrix_comp$pert_idose, row.names = NULL) #There are 43 different concentrations
cp_dose_info <- lapply(unique(cp_dose_info), table)
df_dose <- as.data.frame(cp_dose_info$siginfo_matrix.siginfo_matrix.pert_type.....trt_cp...cmap_name)
cp_dose <- ggplot(df_dose, aes(Freq)) +
  geom_bar(fill = "cornflowerblue") +
  labs(x = "No of concentration", y = "No of compounds", title = "No of concentrations per compound") +
  theme_bw()

#Plot different time doses use for each compound
cp_length_info <- data.frame(siginfo_matrix_comp$cmap_name, siginfo_matrix_comp$pert_itime, row.names = NULL)  #Durations of treatments = 24 h 3 h 48 h 6 h
cp_length_info <- lapply(unique(cp_length_info), table)
df_length <- as.data.frame(cp_length_info$siginfo_matrix.siginfo_matrix.pert_type.....trt_cp...cmap_name)
cp_length <- ggplot(df_length, aes(Freq)) +
  geom_bar(fill = "cornflowerblue") +
  labs(x= "Treatment length", y = "No of compounds", title = "Treatment length per compound") + 
  theme_bw()

#Calculate total number of MoA present in the database for further analyses.
moa_info <- data.frame(compoundinfo_matrix$cmap_name, compoundinfo_matrix$moa, row.names = NULL)
moa_info <- lapply(unique(moa_info),table)
df_moa <- as.data.frame(moa_info$compoundinfo_matrix.moa)
nrow(moa_info$compoundinfo_matrix.moa) #There are 658 moa in the matrix