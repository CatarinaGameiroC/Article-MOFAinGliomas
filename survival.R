source("setup.R")

## Load data 
mrna <- as.matrix(read.csv("DataSets/processed_assays/assay_rna_coding.csv", row.names = 1))
dna <- as.matrix(read.csv("DataSets/processed_assays/assay_methylation.csv", row.names = 1))

classification <- as.matrix(read.csv("DataSets/SIMPLIFIED_CLASSIFICATION_TCGA_2016_2021.csv", row.names = 1))
survival.complete <- as.matrix(read.csv("DataSets/processed_assays/survival_complete.csv", row.names = 1))


classification <- classification[,"classification.2021"]
classification <- ifelse(classification == "glioblastoma", "GBM", 
                         ifelse(classification  == "astrocytoma", "ASTRO",
                                ifelse(classification  == "oligodendroglioma", "OLIGO",
                                       NA)))


mrna <- mrna[intersect(rownames(mrna), intersect(rownames(survival.complete), names(classification))), ] # 652
dna <- dna[intersect(rownames(dna), intersect(rownames(survival.complete), names(classification))), ] # 624

survival.complete <- survival.complete[rownames(mrna), ]
classification <- classification[rownames(mrna)]
survival.complete <- as.data.frame(survival.complete)

survival.complete_dna <- survival.complete[rownames(dna), ]
classification_dna <- classification[rownames(dna)]
survival.complete_dna <- as.data.frame(survival.complete_dna)


## Functions
surv_group_mrna <- function(vec_genes, omic_dataset, group, survival){
  res <- data.frame(Gene = character(), P_Value = numeric(), stringsAsFactors = FALSE)
  samples <- names(classification)[which(classification %in% group)]
  data_up <- as.data.frame(survival[samples, ])
  for (gene in vec_genes){

    median <- median(omic_dataset[samples, gene], na.rm = TRUE)
    data_up$gene_group <- ifelse(omic_dataset[samples, gene] >= median, "High", "Low")
    
    surv_obj <- Surv(time = data_up[,"Time"], event = data_up[ ,"Status"])
    km_fit <- survfit(surv_obj ~ gene_group, data = data_up)
 
    log_rank_test <- survdiff(surv_obj ~ gene_group, data = data_up)
    p_value <- log_rank_test$pvalue
    
    res <- rbind(res, data.frame(Gene = gene, P_Value = p_value))
    #p <- ggsurvplot(km_fit, data = data_up,   pval = TRUE, isk.table = FALSE)
    #print(p)
  }
  
  rownames(res) <- res[,1]
  res <- res[, -1, drop=F]
  return(res)
}

surv_group_methy <- function(vec_genes, omic_dataset, group, survival){
  res <- data.frame(Gene = character(), P_Value = numeric(), stringsAsFactors = FALSE)
  samples <- names(classification_dna)[which(classification_dna %in% group)]
  data_up <- as.data.frame(survival[samples, ])
  for (gene in vec_genes){
    
    median <- median(omic_dataset[samples, gene], na.rm = TRUE)
    data_up$gene_group <- ifelse(omic_dataset[samples, gene] >= median, "High", "Low")
    
    surv_obj <- Surv(time = data_up[,"Time"], event = data_up[ ,"Status"])
    km_fit <- survfit(surv_obj ~ gene_group, data = data_up)
    
    log_rank_test <- survdiff(surv_obj ~ gene_group, data = data_up)
    p_value <- log_rank_test$pvalue
    
    res <- rbind(res, data.frame(Gene = gene, P_Value = p_value))
    #p <- ggsurvplot(km_fit, data = data_up,   pval = TRUE, isk.table = FALSE)
    #print(p)
  }
  rownames(res) <- res[,1]
  res <- res[, -1, drop=F]
  return(res)
}
#ASTRO   GBM OLIGO 
# 248   207   165 

gbm_mrna <- surv_group_mrna(c(mrna_sel_f1, mrna_sel_f2, mrna_sel_f3), mrna, "GBM", survival.complete)
lgg_mrna <- surv_group_mrna(c(mrna_sel_f1, mrna_sel_f2, mrna_sel_f3), mrna, c("ASTRO", "OLIGO"), survival.complete)

                       
gbm_methy <- surv_group_methy(c(methy_sel_f1, methy_sel_f3), dna, "GBM", survival.complete_dna)
lgg_methy <- surv_group_methy(c(methy_sel_f1, methy_sel_f3), dna, c("ASTRO", "OLIGO"),survival.complete_dna)
                                              
                       
                       

## SAVE
write.csv(gbm_mrna, "survival-results/gbm_mrna.csv", row.names = TRUE)
write.csv(lgg_mrna, "survival-results/lgg_mrna.csv", row.names = TRUE)

write.csv(gbm_methy, "survival-results/gbm_methy.csv", row.names = TRUE)
write.csv(lgg_methy, "survival-results/lgg_methy.csv", row.names = TRUE)
