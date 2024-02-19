library(survival)
library(survminer)
library(glmnet)
library(dplyr)
library(timeROC)
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(svglite)
library(clusterProfiler)

#dgList obtained in ../Differential gene expression/DGE_and_GSEA.R needs to be uploaded as dgList.
#Upload final csv file with pathological data with patient ID as rows and clinical information as columns called path_metadata.
#Upload final csv file with outcomes as decomp.

#For this section, only biopsies are selected.
dgList_biopsy <- dgList[, which(dgList$samples$specimen_type == "Biopsy")]

#Comparison between advanced and early fibrosis is performed.
F3F4vsF0F1 <- dea_limma(dgList_biopsy, group_var = "NAFLD", blocking_vars = "sex", contrast_name = "NAFLDF3F4 - NAFLDF0F1")

#Clean and arrange file.
decomp <- decomp[match(dgList_biopsy$samples$sample_id, decomp$STUDY_NO), ] #Match order between dgList and decomp file.
decomp <- decomp[which(!is.na(decomp$decompStatus)),] #Remove those samples that have not been censored for decompensation events.
decomp <- decomp[which(decomp$nash.crn_kleiner_fibrosis_stage == "3" | decomp$nash.crn_kleiner_fibrosis_stage == "4" | decomp$nash.crn_kleiner_fibrosis_stage == "0" | decomp$nash.crn_kleiner_fibrosis_stage == "1a" | decomp$nash.crn_kleiner_fibrosis_stage== "1b" | decomp$nash.crn_kleiner_fibrosis_stage == "1c"),] #Select those samples that present either stage F3/F4 or F0/F1.
decomp <- decomp[which(decomp$decomp_time_days_compRisk >= 0), ] #Remove those events that took place before the biopsy was taken.
decomp <- decomp[which(decomp$decomp_time_days_compRisk <= 5478.75), ] #Remove those events that took place after 20 years
rownames(decomp) <- decomp$STUDY_NO
decomp <- decomp %>% mutate(Fib = case_when(nash.crn_kleiner_fibrosis_stage == "3" | nash.crn_kleiner_fibrosis_stage == "4" ~ "F3F4", TRUE ~ "F0F1")) #Change variable names for easier reading

#Create matrix with log CPM counts and an additional column of time and status
counts_uni <- logcpm_norm[,which(colnames(logcpm_norm)%in%decomp$STUDY_NO)]
counts_uni <- as.data.frame(t(counts_uni)) #Data frame needs to be transposed so genes are as columns and patients as rows.
counts_uni <- counts_uni[match(decomp$STUDY_NO, rownames(counts_uni)),]
counts_uni <- counts_uni[, which(colnames(counts_uni)%in%F3F4vsF0F1$diffexp_genes)]
counts_uni$time <- decomp$decomp_time_days_compRisk + 1 #Add an extra unit for each day so there is no issues with time point = 0.
counts_uni$status <- decomp$decompStatus
colnames(counts_uni) <- gsub(".", "_", colnames(counts_uni), fixed = TRUE)
colnames(counts_uni) <- gsub("-", "_", colnames(counts_uni), fixed = TRUE)

#Run univariate cox regression model as initial step
univ_formulas <- sapply(colnames(counts_uni), function(x){as.formula(paste("Surv(time, status)~", x))})
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = counts_uni)})

univ_results <- lapply(univ_models, function(x){
  x <- summary(x)
  p_vals <- signif(x$sctest[3], digits = 5)
  logrank_test <- signif(x$sctest[1], digits = 5)
  beta <- signif(x$coef[1], digits = 5)
  HR <- signif(x$coef[2], digits = 5)
  HR_confint_lower <- signif(x$conf.int[, "lower .95"], 5)
  HR_confint_upper <- signif(x$conf.int[, "upper .95"], 5)
  HR <- paste(HR, " (", HR_confint_lower, "-", HR_confint_upper, ") ")
  res <- c(beta, HR, logrank_test, p_vals)
  names(res) <- c("beta", "HR (95% CI for HR", "wald test", "p-value")
  return(res)
})

res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
res$`p-value` <- as.numeric(res$`p-value`)
res_multi <- res[which(res$`p-value` < 0.01),] #Choose those genes with pval < 0.01

#Run regularised cox regression model with 10 runs of 10-fold CV
counts_multi <- counts_uni[, which(colnames(counts_uni)%in%rownames(res_multi))]

nreps = 10
nfold = 10

cv_iterations <- lapply(1:nreps, function(i){
  cv.fit <- cv.glmnet(data.matrix(counts_multi[,1:(ncol(counts_multi)-2)]), Surv(counts_multi$time, counts_multi$status), family = "cox", nfolds = nfold)
  data.frame(MSE = cv.fit$cvm, lambda = cv.fit$lambda, se = cv.fit$cvsd)
})

cv_iterations <- do.call(rbind, cv_iterations)
summarised_cv_iter <- cv_iterations %>% group_by(lambda) %>% summarise(meanMSE = mean(MSE), meanse = mean(se)) %>% arrange(desc(lambda))
idx <- which.min(summarised_cv_iter$meanMSE)
lambda.min <- summarised_cv_iter$lambda[idx] #Lambda.min is saved as an object
index_1se <- with(summarised_cv_iter,which(meanMSE < meanMSE[idx]+meanse[idx])[1])
lambda_1se <- summarised_cv_iter$lambda[index_1se] #Lambda +/- 1 SE has been saved too. This value could be used instead of lambda.min

#I have used lambda min for the analyses
fit <- glmnet(data.matrix(counts_multi[,1:(ncol(counts_multi)-2)]), Surv(counts_multi$time, counts_multi$status), family = "cox", maxit = 100000)
Coefficients <- coef(fit, s = lambda.min)
active.index <- which(Coefficients != 0)
active.coefficients <- Coefficients[active.index]
genes_lambdamin <- names(counts_multi)[active.index] 

risk_scores <- counts_multi[active.index]
for (i in 1:ncol(risk_scores)){
  risk_scores[i] <- risk_scores[i] * active.coefficients[i]
}

#Create new risk score and divide patients into high and low-risk score
risk_scores$sums <- rowSums(risk_scores)
risk_scores <- risk_scores %>% mutate(risk = case_when(sums >= median(sums) ~ "high", TRUE ~ "low")) #Median is -0.7360888
risk_scores$status <- counts_multi$status
risk_scores$time <- counts_multi$time
risk_scores$Fib <- decomp$Fib

fit2 <- survfit(Surv(time, status) ~ risk, data = risk_scores)

#Create file to obtain gene expression of each gene of the TRS according to the fibrotic stage
genes_expr <- data.frame()
for (i in 1:length(active.index)){
  m2 <- list(counts = logcpm_norm[which(rownames(logcpm_norm) == colnames(counts_multi)[active.index[i]]),], groups = as.factor(dgList_biopsy$samples$NAFLD[which(dgList_biopsy$samples$sample_id%in%decomp$STUDY_NO)]))
  m2 <- as_tibble(m2)
  m2 <- m2[which(m2$groups == "F4" | m2$groups == "F0F1" | m2$groups == "F3" | m2$groups == "F2"),]
  m2$gene <- colnames(counts_multi)[active.index[i]]
  genes_expr <- rbind(genes_expr, m2)
}

comparisons <- list(c("F0F1", "F2"), c("F2", "F3"), c("F3", "F4"), c("F0F1", "F3"), c("F2", "F4"), c("F0F1", "F4"))

#Time-dependent AUC for 1, 3 and 5-year from biopsy
roc_surv <- timeROC(T = risk_scores$time, delta = risk_scores$status, marker = risk_scores$sums, cause =1, weighting = "marginal", times = c(365, 1095, 1825))





