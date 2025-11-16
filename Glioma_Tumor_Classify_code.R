#####Load libraries#####
options(future.globals.maxSize=6000*1024^2)
library(gplots)
library(biomaRt)
library(tximport)
library(edgeR)
library(tidyr)
library(dplyr)
library(limma)
library(ggplot2)
library(GSVA)
library(reshape2)
library(vcd)
library(fastDummies)
library(BiocParallel)
library(parallel)
library(stringr)
library(cowplot)
library(enrichR)
library(survival)
library(survminer)
library(nestedcv)
library(pROC)
library(ggmosaic)
library(ggvenn)
library(DRIMSeq)
library(immunedeconv)
library(openxlsx)
library(carSurv)
library(limma)
library(caret)
library(Seurat)
library(AUCell)
library(scCustomize)
library(ROCR)
library(variancePartition)
library(umap)
library(STRINGdb)
library(pheatmap)
library(car)
library(enrichR)
library(ggupset)
library(matrixStats)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(gplots)
library(ggnewscale)
library(forcats)
library(ggsankey)
library(ComplexHeatmap)
setwd("C:/Users/HanLab_1/Documents/Mike_Work/GBM_projects/")

############################################################
#####Prepare covariate data (Bulk Glioma)#####
############################################################
factors_gbmlgg_tcga=read.csv("TCGA_GBMLGG-Clinical-Zakharova_IJMS_2022.csv")
factors_gbmlgg_tcga=factors_gbmlgg_tcga[factors_gbmlgg_tcga$TCGA_Diagnosis!="Normal",]
rownames(factors_gbmlgg_tcga)=factors_gbmlgg_tcga[,1]
write.csv(factors_gbmlgg_tcga[c(1,3:4,8,10:12,23:24,15:16)],"Pub_Classify/TCGA_Glioma_Clinical_Info.csv",row.names = F)
factors_gbmlgg_tcga[,c(11,15:22)]=sapply(factors_gbmlgg_tcga[,c(11,15:22)],as.numeric)

factors_gbmlgg_cgga=read.csv("CGGA_GBMLGG-Clinical.csv")
rownames(factors_gbmlgg_cgga)=factors_gbmlgg_cgga[,1]
write.csv(factors_gbmlgg_cgga[c(1:2,4,7:9,12:13,10:11)],"Pub_Classify/CGGA_Glioma_Clinical_Info.csv",row.names = F)
factors_gbmlgg_cgga[,c(8,10:11)]=sapply(factors_gbmlgg_cgga[,c(8,10:11)],as.numeric)

factors_gbmlgg_glass=read.csv("GLASS_GBMLGG-Clinical-cBioPortal.csv")
rownames(factors_gbmlgg_glass)=factors_gbmlgg_glass[,1]
write.csv(factors_gbmlgg_glass[c(1,3,5,8:10,13:14,11:12)],"Pub_Classify/GLASS_Glioma_Clinical_Info.csv",row.names = F)
factors_gbmlgg_glass[,c(9,11:12)]=sapply(factors_gbmlgg_glass[,c(9,11:12)],as.numeric)

gene_labels=read.csv("GENCODEv47_Ensemblv113-Gene_Annotations.csv")
tx2gene=readRDS("tx2gene.rds")
tx2gene$GENEID=gsub("\\..*","",tx2gene$GENEID)
tx2gene$TXNAME=gsub("\\..*","",tx2gene$TXNAME)

############################################################
#####Prepare gene count data (Bulk Glioma)#####
############################################################
data_transcript_gbmlgg_tcga=readRDS("TCGA_GBMLGG-Transcript-Counts-Salmon-Gencode_GRCh38_v47.rds")
rownames(data_transcript_gbmlgg_tcga$abundance)=gsub("\\..*","",rownames(data_transcript_gbmlgg_tcga$abundance))
rownames(data_transcript_gbmlgg_tcga$counts)=gsub("\\..*","",rownames(data_transcript_gbmlgg_tcga$counts))
rownames(data_transcript_gbmlgg_tcga$length)=gsub("\\..*","",rownames(data_transcript_gbmlgg_tcga$length))
data_transcript_gbmlgg_tcga$abundance=data_transcript_gbmlgg_tcga$abundance[,colnames(data_transcript_gbmlgg_tcga$abundance)%in%rownames(factors_gbmlgg_tcga)]
data_transcript_gbmlgg_tcga$counts=data_transcript_gbmlgg_tcga$counts[,colnames(data_transcript_gbmlgg_tcga$counts)%in%rownames(factors_gbmlgg_tcga)]
data_transcript_gbmlgg_tcga$length=data_transcript_gbmlgg_tcga$length[,colnames(data_transcript_gbmlgg_tcga$length)%in%rownames(factors_gbmlgg_tcga)]

data_gene_gbmlgg_tcga=summarizeToGene(data_transcript_gbmlgg_tcga,tx2gene)
cts_gene_gbmlgg_tcga=makeCountsFromAbundance(data_gene_gbmlgg_tcga$counts,data_gene_gbmlgg_tcga$abundance,data_gene_gbmlgg_tcga$length,"scaledTPM")
cts_gene_gbmlgg_tcga=cts_gene_gbmlgg_tcga[rownames(cts_gene_gbmlgg_tcga) %in% gene_labels[!grepl("^RPL|^RPS|^MRPL|^MRPS|N/A",gene_labels$HGNC_Symbol) & grepl("protein|^IG_",gene_labels$Gene_Biotype) & !grepl("pseudo|RNA|vault|artifact",gene_labels$Gene_Biotype) & !grepl("MT",gene_labels$Chromosome),]$Ensembl_Gene_ID,colnames(cts_gene_gbmlgg_tcga) %in% rownames(factors_gbmlgg_tcga)]
keep <- filterByExpr(cts_gene_gbmlgg_tcga,min.prop=0.1,min.count=15)
cts_gene_gbmlgg_tcga=cts_gene_gbmlgg_tcga[keep,]

cts_gene_gbmlgg_cgga=readRDS("CGGA-STAR-Counts-Gencode_GRCh38_v47.rds")
cts_gene_gbmlgg_cgga=cts_gene_gbmlgg_cgga[rownames(cts_gene_gbmlgg_cgga) %in% gene_labels[!grepl("^RPL|^RPS|^MRPL|^MRPS|N/A",gene_labels$HGNC_Symbol) & grepl("protein|^IG_",gene_labels$Gene_Biotype) & !grepl("pseudo|RNA|vault|artifact",gene_labels$Gene_Biotype) & !grepl("MT",gene_labels$Chromosome),]$Ensembl_Gene_ID,colnames(cts_gene_gbmlgg_cgga) %in% rownames(factors_gbmlgg_cgga)]
keep <- filterByExpr(cts_gene_gbmlgg_cgga,min.prop=0.1,min.count=15)
cts_gene_gbmlgg_cgga=cts_gene_gbmlgg_cgga[keep,]

data_transcript_gbmlgg_glass=readRDS("GLASS_GBMLGG-Transcript-Counts-Salmon.rds")
colnames(data_transcript_gbmlgg_glass$abundance)=gsub("\\.","-",colnames(data_transcript_gbmlgg_glass$abundance))
colnames(data_transcript_gbmlgg_glass$counts)=gsub("\\.","-",colnames(data_transcript_gbmlgg_glass$counts))
colnames(data_transcript_gbmlgg_glass$length)=gsub("\\.","-",colnames(data_transcript_gbmlgg_glass$length))
data_transcript_gbmlgg_glass$abundance=data_transcript_gbmlgg_glass$abundance[,colnames(data_transcript_gbmlgg_glass$abundance)%in%rownames(factors_gbmlgg_glass)]
data_transcript_gbmlgg_glass$counts=data_transcript_gbmlgg_glass$counts[,colnames(data_transcript_gbmlgg_glass$counts)%in%rownames(factors_gbmlgg_glass)]
data_transcript_gbmlgg_glass$length=data_transcript_gbmlgg_glass$length[,colnames(data_transcript_gbmlgg_glass$length)%in%rownames(factors_gbmlgg_glass)]

data_gene_gbmlgg_glass=summarizeToGene(data_transcript_gbmlgg_glass,tx2gene)
cts_gene_gbmlgg_glass=makeCountsFromAbundance(data_gene_gbmlgg_glass$counts,data_gene_gbmlgg_glass$abundance,data_gene_gbmlgg_glass$length,"scaledTPM")
cts_gene_gbmlgg_glass=cts_gene_gbmlgg_glass[,!duplicated(colnames(cts_gene_gbmlgg_glass))]
cts_gene_gbmlgg_glass=cts_gene_gbmlgg_glass[rownames(cts_gene_gbmlgg_glass) %in% gene_labels[!grepl("^RPL|^RPS|^MRPL|^MRPS|N/A",gene_labels$HGNC_Symbol) & grepl("protein|^IG_",gene_labels$Gene_Biotype) & !grepl("pseudo|RNA|vault|artifact",gene_labels$Gene_Biotype) & !grepl("MT",gene_labels$Chromosome),]$Ensembl_Gene_ID,colnames(cts_gene_gbmlgg_glass) %in% rownames(factors_gbmlgg_glass)]
keep <- filterByExpr(cts_gene_gbmlgg_glass,min.prop=0.1,min.count=15)
cts_gene_gbmlgg_glass=cts_gene_gbmlgg_glass[keep,]

intersecting_genes=intersect(rownames(cts_gene_gbmlgg_tcga),intersect(rownames(cts_gene_gbmlgg_cgga),rownames(cts_gene_gbmlgg_glass)))

cts_gene_gbmlgg_tcga=cts_gene_gbmlgg_tcga[rownames(cts_gene_gbmlgg_tcga) %in% intersecting_genes,]
cts_gene_gbmlgg_cgga=cts_gene_gbmlgg_cgga[rownames(cts_gene_gbmlgg_cgga) %in% intersecting_genes,]
cts_gene_gbmlgg_glass=cts_gene_gbmlgg_glass[rownames(cts_gene_gbmlgg_glass) %in% intersecting_genes,]

dge_gene_gbmlgg_tcga <- DGEList(cts_gene_gbmlgg_tcga)
dge_gene_gbmlgg_tcga <- calcNormFactors(dge_gene_gbmlgg_tcga)
design <- ~0+Read_Length
v_gene_gbmlgg_tcga <- voomWithDreamWeights(dge_gene_gbmlgg_tcga, design, factors_gbmlgg_tcga, BPPARAM = SnowParam(workers = 10,exportglobals = F),plot = F)
cve_gene_gbmlgg_tcga <- removeBatchEffect(v_gene_gbmlgg_tcga$E,batch=factors_gbmlgg_tcga$Read_Length)

dge_gene_gbmlgg_cgga <- DGEList(cts_gene_gbmlgg_cgga)
dge_gene_gbmlgg_cgga <- calcNormFactors(dge_gene_gbmlgg_cgga)
design <- ~0+CGGA_Dataset_Detailed
v_gene_gbmlgg_cgga <- voomWithDreamWeights(dge_gene_gbmlgg_cgga, design, factors_gbmlgg_cgga, BPPARAM = SnowParam(workers = 10,exportglobals = F),plot = F)
cve_gene_gbmlgg_cgga <- removeBatchEffect(v_gene_gbmlgg_cgga$E,batch=factors_gbmlgg_cgga$CGGA_Dataset_Detailed)

dge_gene_gbmlgg_glass <- DGEList(cts_gene_gbmlgg_glass)
dge_gene_gbmlgg_glass <- calcNormFactors(dge_gene_gbmlgg_glass)
design <- ~0+Tissue_Source_Site_Detailed
v_gene_gbmlgg_glass <- voomWithDreamWeights(dge_gene_gbmlgg_glass, design, factors_gbmlgg_glass, BPPARAM = SnowParam(workers = 10,exportglobals = F),plot = F)
cve_gene_gbmlgg_glass <- removeBatchEffect(v_gene_gbmlgg_glass$E,batch=factors_gbmlgg_glass$Tissue_Source_Site_Detailed)

##################################################################
#####Gene Signature Classification (Bulk Glioma)#####
##################################################################
limma_filter <- function(x, y, nfilter) {
  library(limma)
  library(dplyr)
  y <- as.factor(y)
  design <- model.matrix(~ 0 + y)
  fit <- limma::lmFit(t(x), design)
  contrasts=as.vector(sapply(levels(y), function(cls) { paste0("y",cls, " - (", paste0("y",setdiff(levels(y), cls), collapse = " + "), ")/", length(setdiff(levels(y), cls))) }))
  contrast_matrix =limma::makeContrasts(contrasts=contrasts,levels=design)
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)
  all_top_genes <- unique(unlist(
    lapply(1:ncol(contrast_matrix), function(i) {
        top <- head(limma::topTable(fit2, coef = i, number = ncol(x), sort.by = "p", adjust.method="BH")  %>% .[order(.$adj.P.Val),],nfilter)
      rownames(top)
    })
  ))
  matched_idx <- sort(match(all_top_genes, colnames(x)))
  matched_idx
}

x=as.matrix(t(cve_gene_gbmlgg_tcga[,colnames(cve_gene_gbmlgg_tcga) %in% rownames(factors_gbmlgg_tcga[!grepl("N/A",factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis),])]))
y=factors_gbmlgg_tcga[!grepl("N/A",factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis),]$WHO_CNS5_2021_Diagnosis
cvfit_glm_lasso=nestcv.glmnet(y = y, x = x, family = "multinomial", cv.cores = 10, alphaSet = 1, n_outer_folds = 10, n_inner_folds = 10, finalCV = T, filterFUN = limma_filter, filter_options = list(nfilter=500), verbose = T)
cvfit_glm_ridge=nestcv.glmnet(y = y, x = x, family = "multinomial", cv.cores = 10, alphaSet = 0, n_outer_folds = 10, n_inner_folds = 10, finalCV = T, filterFUN = limma_filter, filter_options = list(nfilter=500), verbose = T)
cvfit_glm_elasticnet=nestcv.glmnet(y = y, x = x, family = "multinomial", cv.cores = 10, alphaSet = 0.25, n_outer_folds = 10, n_inner_folds = 10, finalCV = T, filterFUN = limma_filter, filter_options = list(nfilter=500), verbose = T)
cvfit_rf=nestcv.train(y = y, x = x, method = "rf",n_outer_folds = 10,n_inner_folds = 10, filterFUN = limma_filter, filter_options = list(nfilter=500),finalCV = T,tuneGrid=expand.grid(mtry=sqrt(ncol(x))),cv.cores = 10,savePredictions = "final")
cvfit_xgb=nestcv.train(y = y, x = x, method = "xgbDART",n_outer_folds = 10,n_inner_folds = 10, filterFUN = limma_filter, filter_options = list(nfilter=500),finalCV = T,tuneGrid=expand.grid(nrounds = 10, max_depth = 2, eta = 0.30, rate_drop = 0.10, skip_drop = 0.10, colsample_bytree = 0.90, min_child_weight = 2, subsample = 0.75, gamma = 0.10),cv.cores = 10,savePredictions = "final")


saveRDS(cvfit_glm_lasso,"Pub_Classify/Glioma_GLM_Lasso_Classifier.rds")
saveRDS(cvfit_glm_ridge,"Pub_Classify/Glioma_GLM_Ridge_Classifier.rds")
saveRDS(cvfit_glm_elasticnet,"Pub_Classify/Glioma_GLM_ElasticNet_Classifier.rds")
saveRDS(cvfit_rf,"Pub_Classify/Glioma_RandomForest_Model.rds")
saveRDS(cvfit_xgb,"Pub_Classify/Glioma_XGBoost_Model.rds")


factors_gbmlgg_tcga$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis=predict(object = cvfit_glm_lasso,newdata = as.matrix(t(cve_gene_gbmlgg_tcga)),type="class")[,1]
factors_gbmlgg_cgga$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis=predict(object = cvfit_glm_lasso,newdata = as.matrix(t(cve_gene_gbmlgg_cgga)),type="class")[,1]
factors_gbmlgg_glass$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis=predict(object = cvfit_glm_lasso,newdata = as.matrix(t(cve_gene_gbmlgg_glass)),type="class")[,1]

factors_gbmlgg_tcga$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis=predict(object = cvfit_glm_ridge,newdata = as.matrix(t(cve_gene_gbmlgg_tcga)),type="class")[,1]
factors_gbmlgg_cgga$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis=predict(object = cvfit_glm_ridge,newdata = as.matrix(t(cve_gene_gbmlgg_cgga)),type="class")[,1]
factors_gbmlgg_glass$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis=predict(object = cvfit_glm_ridge,newdata = as.matrix(t(cve_gene_gbmlgg_glass)),type="class")[,1]

factors_gbmlgg_tcga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis=predict(object = cvfit_glm_elasticnet,newdata = as.matrix(t(cve_gene_gbmlgg_tcga)),type="class")[,1]
factors_gbmlgg_cgga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis=predict(object = cvfit_glm_elasticnet,newdata = as.matrix(t(cve_gene_gbmlgg_cgga)),type="class")[,1]
factors_gbmlgg_glass$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis=predict(object = cvfit_glm_elasticnet,newdata = as.matrix(t(cve_gene_gbmlgg_glass)),type="class")[,1]

factors_gbmlgg_tcga$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis=predict(object = cvfit_rf,newdata = as.matrix(t(cve_gene_gbmlgg_tcga)))
factors_gbmlgg_cgga$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis=predict(object = cvfit_rf,newdata = as.matrix(t(cve_gene_gbmlgg_cgga)))
factors_gbmlgg_glass$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis=predict(object = cvfit_rf,newdata = as.matrix(t(cve_gene_gbmlgg_glass)))

factors_gbmlgg_tcga$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis=predict(object = cvfit_xgb,newdata = as.matrix(t(cve_gene_gbmlgg_tcga)))
factors_gbmlgg_cgga$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis=predict(object = cvfit_xgb,newdata = as.matrix(t(cve_gene_gbmlgg_cgga)))
factors_gbmlgg_glass$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis=predict(object = cvfit_xgb,newdata = as.matrix(t(cve_gene_gbmlgg_glass)))


temp_roc_classify_glm_lasso=list(roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=1,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=1,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=1,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=1,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=1,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=1,"N/A"=0)))

temp_roc_classify_glm_ridge=list(roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=1,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=1,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=1,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=1,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=1,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=1,"N/A"=0)))

temp_roc_classify_glm_elasticnet=list(roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                      roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=1,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=1,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                      roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                      roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=1,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=1,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                      roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                      roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=1,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=1,"N/A"=0)))

temp_roc_classify_rf=list(roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=1,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=1,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=1,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=1,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=1,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=1,"N/A"=0)))

temp_roc_classify_xgb=list(roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=1,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=1,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=1,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=1,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=0,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=0,"N/A"=0)),
                                 roc(sapply(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=1,"N/A"=0),sapply(factors_gbmlgg_tcga$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=0,"Oligodendroglioma_G3"=1,"N/A"=0)))

figs2a=ggroc(temp_roc_classify_glm_elasticnet)+theme_classic()+ylab("Sensitivity")+xlab("Specificity")+labs(color = "Glioma Type")+scale_color_manual(labels = c(paste0("A (G2) (AUC=",signif(temp_roc_classify_glm_elasticnet[[1]]$auc[1],2),")"),paste0("A (G3) (AUC=",signif(temp_roc_classify_glm_elasticnet[[2]]$auc[1],2),")"),paste0("A (G4) (AUC=",signif(temp_roc_classify_glm_elasticnet[[3]]$auc[1],2),")"),paste0("GBM (G4) (AUC=",signif(temp_roc_classify_glm_elasticnet[[4]]$auc[1],2),")"),paste0("O (G2) (AUC=",signif(temp_roc_classify_glm_elasticnet[[5]]$auc[1],2),")"),paste0("O (G3) (AUC=",signif(temp_roc_classify_glm_elasticnet[[6]]$auc[1],2),")")), values = c("yellow2","orange2","firebrick2","purple2","lightblue2","blue2"))+theme(legend.position = c(0.725, 0.225))+geom_abline(slope=1, intercept = 1,linetype=3)+theme(text = element_text(size = 6),legend.key.height = unit(0.5, "lines"),legend.background = element_rect(fill = NA, colour = NA))


temp_acc_classify_glm_lasso=confusionMatrix(data=factor(factors_gbmlgg_tcga[factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis!="N/A",]$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis,levels=c("Glioblastoma_G4","Astrocytoma_G2","Astrocytoma_G3","Astrocytoma_G4","Oligodendroglioma_G2","Oligodendroglioma_G3")),reference=factor(factors_gbmlgg_tcga[factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis!="N/A",]$WHO_CNS5_2021_Diagnosis,levels=c("Glioblastoma_G4","Astrocytoma_G2","Astrocytoma_G3","Astrocytoma_G4","Oligodendroglioma_G2","Oligodendroglioma_G3")))
temp_acc_classify_glm_ridge=confusionMatrix(data=factor(factors_gbmlgg_tcga[factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis!="N/A",]$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis,levels=c("Glioblastoma_G4","Astrocytoma_G2","Astrocytoma_G3","Astrocytoma_G4","Oligodendroglioma_G2","Oligodendroglioma_G3")),reference=factor(factors_gbmlgg_tcga[factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis!="N/A",]$WHO_CNS5_2021_Diagnosis,levels=c("Glioblastoma_G4","Astrocytoma_G2","Astrocytoma_G3","Astrocytoma_G4","Oligodendroglioma_G2","Oligodendroglioma_G3")))
temp_acc_classify_glm_elasticnet=confusionMatrix(data=factor(factors_gbmlgg_tcga[factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis!="N/A",]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis,levels=c("Glioblastoma_G4","Astrocytoma_G2","Astrocytoma_G3","Astrocytoma_G4","Oligodendroglioma_G2","Oligodendroglioma_G3")),reference=factor(factors_gbmlgg_tcga[factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis!="N/A",]$WHO_CNS5_2021_Diagnosis,levels=c("Glioblastoma_G4","Astrocytoma_G2","Astrocytoma_G3","Astrocytoma_G4","Oligodendroglioma_G2","Oligodendroglioma_G3")))
temp_acc_classify_rf=confusionMatrix(data=factor(factors_gbmlgg_tcga[factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis!="N/A",]$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis,levels=c("Glioblastoma_G4","Astrocytoma_G2","Astrocytoma_G3","Astrocytoma_G4","Oligodendroglioma_G2","Oligodendroglioma_G3")),reference=factor(factors_gbmlgg_tcga[factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis!="N/A",]$WHO_CNS5_2021_Diagnosis,levels=c("Glioblastoma_G4","Astrocytoma_G2","Astrocytoma_G3","Astrocytoma_G4","Oligodendroglioma_G2","Oligodendroglioma_G3")))
temp_acc_classify_xgb=confusionMatrix(data=factor(factors_gbmlgg_tcga[factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis!="N/A",]$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis,levels=c("Glioblastoma_G4","Astrocytoma_G2","Astrocytoma_G3","Astrocytoma_G4","Oligodendroglioma_G2","Oligodendroglioma_G3")),reference=factor(factors_gbmlgg_tcga[factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis!="N/A",]$WHO_CNS5_2021_Diagnosis,levels=c("Glioblastoma_G4","Astrocytoma_G2","Astrocytoma_G3","Astrocytoma_G4","Oligodendroglioma_G2","Oligodendroglioma_G3")))


temp_roc_idh_glm_lasso=list(roc(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                            roc(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                            roc(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)))

temp_roc_idh_glm_ridge=list(roc(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                            roc(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                            roc(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)))

temp_roc_idh_glm_elasticnet=list(roc(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                            roc(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                            roc(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)))

temp_roc_idh_rf=list(roc(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                            roc(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                            roc(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)))

temp_roc_idh_xgb=list(roc(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                            roc(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                            roc(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)))

figs2b=ggroc(temp_roc_idh_glm_elasticnet)+theme_classic()+ylab("Sensitivity")+xlab("Specificity")+labs(color = "Dataset\n(IDH Status)")+scale_color_manual(labels = c(paste0("TCGA (AUC=",signif(temp_roc_idh_glm_elasticnet[[1]]$auc[1],2),")"),paste0("CGGA (AUC=",signif(temp_roc_idh_glm_elasticnet[[2]]$auc[1],2),")"),paste0("GLASS (AUC=",signif(temp_roc_idh_glm_elasticnet[[3]]$auc[1],2),")")), values = c("darkred","darkgreen","darkblue"))+theme(legend.position = c(0.7, 0.2))+geom_abline(slope=1, intercept = 1,linetype=3)+theme(text = element_text(size = 6),legend.key.height = unit(0.5, "lines"),legend.background = element_rect(fill = NA, colour = NA))


temp_acc_idh_glm_lasso=list(confusionMatrix(data=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),levels = c(0,1))),
                         confusionMatrix(data=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),levels = c(0,1))),
                         confusionMatrix(data=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),levels = c(0,1))))

temp_acc_idh_glm_ridge=list(confusionMatrix(data=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),levels = c(0,1))),
                        confusionMatrix(data=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),levels = c(0,1))),
                        confusionMatrix(data=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),levels = c(0,1))))

temp_acc_idh_glm_elasticnet=list(confusionMatrix(data=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),levels = c(0,1))),
                             confusionMatrix(data=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),levels = c(0,1))),
                             confusionMatrix(data=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),levels = c(0,1))))

temp_acc_idh_rf=list(confusionMatrix(data=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),levels = c(0,1))),
                 confusionMatrix(data=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),levels = c(0,1))),
                 confusionMatrix(data=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),levels = c(0,1))))

temp_acc_idh_xgb=list(confusionMatrix(data=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),levels = c(0,1))),
                  confusionMatrix(data=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),levels = c(0,1))),
                  confusionMatrix(data=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=1,"Astrocytoma_G3"=1,"Astrocytoma_G4"=1,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$IDH_Status!="N/A",]$IDH_Status, switch,"Mutant"=1,"Wildtype"=0),levels = c(0,1))))


temp_roc_codel_glm_lasso=list(roc(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                    roc(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                    roc(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)))

temp_roc_codel_glm_ridge=list(roc(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                              roc(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                              roc(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)))

temp_roc_codel_glm_elasticnet=list(roc(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                              roc(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                              roc(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)))

temp_roc_codel_rf=list(roc(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                              roc(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                              roc(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)))

temp_roc_codel_xgb=list(roc(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                              roc(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)),
                              roc(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1)))

figs2c=ggroc(temp_roc_codel_glm_elasticnet)+theme_classic()+ylab("Sensitivity")+xlab("Specificity")+labs(color = "Dataset (Chromosome\n1p/19q-codeletion)")+scale_color_manual(labels = c(paste0("TCGA (AUC=",signif(temp_roc_codel_glm_elasticnet[[1]]$auc[1],2),")"),paste0("CGGA (AUC=",signif(temp_roc_codel_glm_elasticnet[[2]]$auc[1],2),")"),paste0("GLASS (AUC=",signif(temp_roc_codel_glm_elasticnet[[3]]$auc[1],2),")")), values = c("darkred","darkgreen","darkblue"))+theme(legend.position = c(0.7, 0.2))+geom_abline(slope=1, intercept = 1,linetype=3)+theme(text = element_text(size = 6),legend.key.height = unit(0.5, "lines"),legend.background = element_rect(fill = NA, colour = NA))


temp_acc_codel_glm_lasso=list(confusionMatrix(data=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),levels = c(0,1))),
                         confusionMatrix(data=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),levels = c(0,1))),
                         confusionMatrix(data=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_Lasso_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),levels = c(0,1))))

temp_acc_codel_glm_ridge=list(confusionMatrix(data=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),levels = c(0,1))),
                         confusionMatrix(data=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),levels = c(0,1))),
                         confusionMatrix(data=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_Ridge_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),levels = c(0,1))))

temp_acc_codel_glm_elasticnet=list(confusionMatrix(data=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),levels = c(0,1))),
                              confusionMatrix(data=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),levels = c(0,1))),
                              confusionMatrix(data=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),levels = c(0,1))))

temp_acc_codel_rf=list(confusionMatrix(data=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),levels = c(0,1))),
                  confusionMatrix(data=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),levels = c(0,1))),
                  confusionMatrix(data=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_RandomForest_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),levels = c(0,1))))

temp_acc_codel_xgb=list(confusionMatrix(data=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_tcga[factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),levels = c(0,1))),
                   confusionMatrix(data=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_cgga[factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),levels = c(0,1))),
                   confusionMatrix(data=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Predicted_XGBoost_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"=0,"Astrocytoma_G3"=0,"Astrocytoma_G4"=0,"Glioblastoma_G4"=0,"Oligodendroglioma_G2"=1,"Oligodendroglioma_G3"=1),levels = c(0,1)),reference=factor(sapply(factors_gbmlgg_glass[factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion!="N/A",]$Chromosome_1p_19q_Codeletion, switch,"Codeleted"=1,"Non_Codeleted"=0),levels = c(0,1))))


classify_predicted_roc=data.frame(Machine_Learning_Algorithm=c(rep("GLM (Lasso)",6),rep("GLM (Ridge)",6),rep("GLM (Elastic Net)",6),rep("Random Forest",6),rep("XGBoost",6)),Glioma_Type=rep(c("A (G2)","A (G3)","A (G4)","GBM (G4)","O (G2)","O (G3)"),5),
                         AUC=c(temp_roc_classify_glm_lasso[[1]]$auc[1],temp_roc_classify_glm_lasso[[2]]$auc[1],temp_roc_classify_glm_lasso[[3]]$auc[1],temp_roc_classify_glm_lasso[[4]]$auc[1],temp_roc_classify_glm_lasso[[5]]$auc[1],temp_roc_classify_glm_lasso[[6]]$auc[1],temp_roc_classify_glm_ridge[[1]]$auc[1],temp_roc_classify_glm_ridge[[2]]$auc[1],temp_roc_classify_glm_ridge[[3]]$auc[1],temp_roc_classify_glm_ridge[[4]]$auc[1],temp_roc_classify_glm_ridge[[5]]$auc[1],temp_roc_classify_glm_ridge[[6]]$auc[1],temp_roc_classify_glm_elasticnet[[1]]$auc[1],temp_roc_classify_glm_elasticnet[[2]]$auc[1],temp_roc_classify_glm_elasticnet[[3]]$auc[1],temp_roc_classify_glm_elasticnet[[4]]$auc[1],temp_roc_classify_glm_elasticnet[[5]]$auc[1],temp_roc_classify_glm_elasticnet[[6]]$auc[1],temp_roc_classify_rf[[1]]$auc[1],temp_roc_classify_rf[[2]]$auc[1],temp_roc_classify_rf[[3]]$auc[1],temp_roc_classify_rf[[4]]$auc[1],temp_roc_classify_rf[[5]]$auc[1],temp_roc_classify_rf[[6]]$auc[1],temp_roc_classify_xgb[[1]]$auc[1],temp_roc_classify_xgb[[2]]$auc[1],temp_roc_classify_xgb[[3]]$auc[1],temp_roc_classify_xgb[[4]]$auc[1],temp_roc_classify_xgb[[5]]$auc[1],temp_roc_classify_xgb[[6]]$auc[1]))

idh_predicted_roc=data.frame(Machine_Learning_Algorithm=c(rep("GLM (Lasso)",3),rep("GLM (Ridge)",3),rep("GLM (Elastic Net)",3),rep("Random Forest",3),rep("XGBoost",3)),Dataset=rep(c("TCGA","CGGA","GLASS"),5),
                         AUC=c(temp_roc_idh_glm_lasso[[1]]$auc[1],temp_roc_idh_glm_lasso[[2]]$auc[1],temp_roc_idh_glm_lasso[[3]]$auc[1],temp_roc_idh_glm_ridge[[1]]$auc[1],temp_roc_idh_glm_ridge[[2]]$auc[1],temp_roc_idh_glm_ridge[[3]]$auc[1],temp_roc_idh_glm_elasticnet[[1]]$auc[1],temp_roc_idh_glm_elasticnet[[2]]$auc[1],temp_roc_idh_glm_elasticnet[[3]]$auc[1],temp_roc_idh_rf[[1]]$auc[1],temp_roc_idh_rf[[2]]$auc[1],temp_roc_idh_rf[[3]]$auc[1],temp_roc_idh_xgb[[1]]$auc[1],temp_roc_idh_xgb[[2]]$auc[1],temp_roc_idh_xgb[[3]]$auc[1]))

codel_predicted_roc=data.frame(Machine_Learning_Algorithm=c(rep("GLM (Lasso)",3),rep("GLM (Ridge)",3),rep("GLM (Elastic Net)",3),rep("Random Forest",3),rep("XGBoost",3)),Dataset=rep(c("TCGA","CGGA","GLASS"),5),
                         AUC=c(temp_roc_codel_glm_lasso[[1]]$auc[1],temp_roc_codel_glm_lasso[[2]]$auc[1],temp_roc_codel_glm_lasso[[3]]$auc[1],temp_roc_codel_glm_ridge[[1]]$auc[1],temp_roc_codel_glm_ridge[[2]]$auc[1],temp_roc_codel_glm_ridge[[3]]$auc[1],temp_roc_codel_glm_elasticnet[[1]]$auc[1],temp_roc_codel_glm_elasticnet[[2]]$auc[1],temp_roc_codel_glm_elasticnet[[3]]$auc[1],temp_roc_codel_rf[[1]]$auc[1],temp_roc_codel_rf[[2]]$auc[1],temp_roc_codel_rf[[3]]$auc[1],temp_roc_codel_xgb[[1]]$auc[1],temp_roc_codel_xgb[[2]]$auc[1],temp_roc_codel_xgb[[3]]$auc[1]))


classify_predicted_acc=data.frame(Machine_Learning_Algorithm=c(rep("GLM (Lasso)",6),rep("GLM (Ridge)",6),rep("GLM (Elastic Net)",6),rep("Random Forest",6),rep("XGBoost",6)),Glioma_Type=rep(c("A (G2)","A (G3)","A (G4)","GBM (G4)","O (G2)","O (G3)"),5),
                                  F1_Score=c(temp_acc_classify_glm_lasso$byClass[,7][1],temp_acc_classify_glm_lasso$byClass[,7][2],temp_acc_classify_glm_lasso$byClass[,7][3],temp_acc_classify_glm_lasso$byClass[,7][4],temp_acc_classify_glm_lasso$byClass[,7][5],temp_acc_classify_glm_lasso$byClass[,7][6],temp_acc_classify_glm_ridge$byClass[,7][1],temp_acc_classify_glm_ridge$byClass[,7][2],temp_acc_classify_glm_ridge$byClass[,7][3],temp_acc_classify_glm_ridge$byClass[,7][4],temp_acc_classify_glm_ridge$byClass[,7][5],temp_acc_classify_glm_ridge$byClass[,7][6],temp_acc_classify_glm_elasticnet$byClass[,7][1],temp_acc_classify_glm_elasticnet$byClass[,7][2],temp_acc_classify_glm_elasticnet$byClass[,7][3],temp_acc_classify_glm_elasticnet$byClass[,7][4],temp_acc_classify_glm_elasticnet$byClass[,7][5],temp_acc_classify_glm_elasticnet$byClass[,7][6],temp_acc_classify_rf$byClass[,7][1],temp_acc_classify_rf$byClass[,7][2],temp_acc_classify_rf$byClass[,7][3],temp_acc_classify_rf$byClass[,7][4],temp_acc_classify_rf$byClass[,7][5],temp_acc_classify_rf$byClass[,7][6],temp_acc_classify_xgb$byClass[,7][1],temp_acc_classify_xgb$byClass[,7][2],temp_acc_classify_xgb$byClass[,7][3],temp_acc_classify_xgb$byClass[,7][4],temp_acc_classify_xgb$byClass[,7][5],temp_acc_classify_xgb$byClass[,7][6]))

idh_predicted_acc=data.frame(Machine_Learning_Algorithm=c(rep("GLM (Lasso)",3),rep("GLM (Ridge)",3),rep("GLM (Elastic Net)",3),rep("Random Forest",3),rep("XGBoost",3)),Dataset=rep(c("TCGA","CGGA","GLASS"),5),
                             F1_Score=c(temp_acc_idh_glm_lasso[[1]]$byClass[7],temp_acc_idh_glm_lasso[[2]]$byClass[7],temp_acc_idh_glm_lasso[[3]]$byClass[7],temp_acc_idh_glm_ridge[[1]]$byClass[7],temp_acc_idh_glm_ridge[[2]]$byClass[7],temp_acc_idh_glm_ridge[[3]]$byClass[7],temp_acc_idh_glm_elasticnet[[1]]$byClass[7],temp_acc_idh_glm_elasticnet[[2]]$byClass[7],temp_acc_idh_glm_elasticnet[[3]]$byClass[7],temp_acc_idh_rf[[1]]$byClass[7],temp_acc_idh_rf[[2]]$byClass[7],temp_acc_idh_rf[[3]]$byClass[7],temp_acc_idh_xgb[[1]]$byClass[7],temp_acc_idh_xgb[[2]]$byClass[7],temp_acc_idh_xgb[[3]]$byClass[7]))

codel_predicted_acc=data.frame(Machine_Learning_Algorithm=c(rep("GLM (Lasso)",3),rep("GLM (Ridge)",3),rep("GLM (Elastic Net)",3),rep("Random Forest",3),rep("XGBoost",3)),Dataset=rep(c("TCGA","CGGA","GLASS"),5),
                               F1_Score=c(temp_acc_codel_glm_lasso[[1]]$byClass[7],temp_acc_codel_glm_lasso[[2]]$byClass[7],temp_acc_codel_glm_lasso[[3]]$byClass[7],temp_acc_codel_glm_ridge[[1]]$byClass[7],temp_acc_codel_glm_ridge[[2]]$byClass[7],temp_acc_codel_glm_ridge[[3]]$byClass[7],temp_acc_codel_glm_elasticnet[[1]]$byClass[7],temp_acc_codel_glm_elasticnet[[2]]$byClass[7],temp_acc_codel_glm_elasticnet[[3]]$byClass[7],temp_acc_codel_rf[[1]]$byClass[7],temp_acc_codel_rf[[2]]$byClass[7],temp_acc_codel_rf[[3]]$byClass[7],temp_acc_codel_xgb[[1]]$byClass[7],temp_acc_codel_xgb[[2]]$byClass[7],temp_acc_codel_xgb[[3]]$byClass[7]))

fig1a=ggplot(classify_predicted_acc, aes(x=Machine_Learning_Algorithm, y=F1_Score, fill=Glioma_Type))+geom_bar(stat = "identity",position = "dodge")+coord_cartesian(ylim=c(0.6,1))+theme_classic()+scale_fill_manual(values=c(`GBM (G4)`="purple2",`A (G2)`="yellow2",`A (G3)`="orange2",`A (G4)`="firebrick2",`O (G2)`="lightblue2",`O (G3)`="blue2"))+xlab("Machine Learning Algorithm")+ylab("F1 Score")+ labs(fill = "Glioma Type")+theme(text = element_text(size = 6),legend.key.height = unit(0.5, "lines"))
fig1b=ggplot(idh_predicted_acc, aes(x=Machine_Learning_Algorithm, y=F1_Score, fill=Dataset))+geom_bar(stat = "identity",position = "dodge")+coord_cartesian(ylim=c(0.6,1))+theme_classic()+scale_fill_manual(values = c(`TCGA`="darkred",`CGGA`="darkgreen",`GLASS`="darkblue"))+xlab("Machine Learning Algorithm")+ylab("F1 Score\n(IDH Status)")+labs(fill = "Dataset")+theme(text = element_text(size = 6),legend.key.height = unit(0.5, "lines"))
fig1c=ggplot(codel_predicted_acc, aes(x=Machine_Learning_Algorithm, y=F1_Score, fill=Dataset))+geom_bar(stat = "identity",position = "dodge")+coord_cartesian(ylim=c(0.6,1))+theme_classic()+scale_fill_manual(values = c(`TCGA`="darkred",`CGGA`="darkgreen",`GLASS`="darkblue"))+xlab("Machine Learning Algorithm")+ylab("F1 Score\n(Chromosome 1p/19q-codeletion)")+labs(fill = "Dataset")+theme(text = element_text(size = 6),legend.key.height = unit(0.5, "lines"))

figs1a=ggplot(classify_predicted_roc, aes(x=Machine_Learning_Algorithm, y=AUC, fill=Glioma_Type))+geom_bar(stat = "identity",position = "dodge")+coord_cartesian(ylim=c(0.6,1))+theme_classic()+scale_fill_manual(values=c(`GBM (G4)`="purple2",`A (G2)`="yellow2",`A (G3)`="orange2",`A (G4)`="firebrick2",`O (G2)`="lightblue2",`O (G3)`="blue2"))+xlab("Machine Learning Algorithm")+ylab("AUC")+ labs(fill = "Glioma Type")+theme(text = element_text(size = 6),legend.key.height = unit(0.5, "lines"))
figs1b=ggplot(idh_predicted_roc, aes(x=Machine_Learning_Algorithm, y=AUC, fill=Dataset))+geom_bar(stat = "identity",position = "dodge")+coord_cartesian(ylim=c(0.6,1))+theme_classic()+scale_fill_manual(values = c(`TCGA`="darkred",`CGGA`="darkgreen",`GLASS`="darkblue"))+xlab("Machine Learning Algorithm")+ylab("AUC\n(IDH Status)")+labs(fill = "Dataset")+theme(text = element_text(size = 6),legend.key.height = unit(0.5, "lines"))
figs1c=ggplot(codel_predicted_roc, aes(x=Machine_Learning_Algorithm, y=AUC, fill=Dataset))+geom_bar(stat = "identity",position = "dodge")+coord_cartesian(ylim=c(0.6,1))+theme_classic()+scale_fill_manual(values = c(`TCGA`="darkred",`CGGA`="darkgreen",`GLASS`="darkblue"))+xlab("Machine Learning Algorithm")+ylab("AUC\n(Chromosome 1p/19q-codeletion)")+labs(fill = "Dataset")+theme(text = element_text(size = 6),legend.key.height = unit(0.5, "lines"))

factors_gbmlgg_tcga$Dropped_Sample=F
factors_gbmlgg_tcga[(factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion=="Non_Codeleted"&!factors_gbmlgg_tcga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis%in%c("Astrocytoma_G2","Astrocytoma_G3","Astrocytoma_G4","Glioblastoma_G4")) | (factors_gbmlgg_tcga$Chromosome_1p_19q_Codeletion=="Codeleted"&factors_gbmlgg_tcga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis%in%c("Astrocytoma_G2","Astrocytoma_G3","Astrocytoma_G4","Glioblastoma_G4")) | (factors_gbmlgg_tcga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis=="Glioblastoma_G4"&factors_gbmlgg_tcga$IDH_Status=="Mutant") | (factors_gbmlgg_tcga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis!="Glioblastoma_G4"&factors_gbmlgg_tcga$IDH_Status=="Wildtype"),]$Dropped_Sample=T
factors_gbmlgg_cgga$Dropped_Sample=F
factors_gbmlgg_cgga[(factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion=="Non_Codeleted"&!factors_gbmlgg_cgga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis%in%c("Astrocytoma_G2","Astrocytoma_G3","Astrocytoma_G4","Glioblastoma_G4")) | (factors_gbmlgg_cgga$Chromosome_1p_19q_Codeletion=="Codeleted"&factors_gbmlgg_cgga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis%in%c("Astrocytoma_G2","Astrocytoma_G3","Astrocytoma_G4","Glioblastoma_G4")) | (factors_gbmlgg_cgga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis=="Glioblastoma_G4"&factors_gbmlgg_cgga$IDH_Status=="Mutant") | (factors_gbmlgg_cgga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis!="Glioblastoma_G4"&factors_gbmlgg_cgga$IDH_Status=="Wildtype"),]$Dropped_Sample=T
factors_gbmlgg_glass$Dropped_Sample=F
factors_gbmlgg_glass[(factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion=="Non_Codeleted"&!factors_gbmlgg_glass$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis%in%c("Astrocytoma_G2","Astrocytoma_G3","Astrocytoma_G4","Glioblastoma_G4")) | (factors_gbmlgg_glass$Chromosome_1p_19q_Codeletion=="Codeleted"&factors_gbmlgg_glass$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis%in%c("Astrocytoma_G2","Astrocytoma_G3","Astrocytoma_G4","Glioblastoma_G4")) | (factors_gbmlgg_glass$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis=="Glioblastoma_G4"&factors_gbmlgg_glass$IDH_Status=="Mutant") | (factors_gbmlgg_glass$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis!="Glioblastoma_G4"&factors_gbmlgg_glass$IDH_Status=="Wildtype"),]$Dropped_Sample=T


counts_tcga_og=table(factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis)
counts_tcga_pred=table(factors_gbmlgg_tcga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis,factors_gbmlgg_tcga$Dropped_Sample)
counts_cgga_pred=table(factors_gbmlgg_cgga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis,factors_gbmlgg_cgga$Dropped_Sample)
counts_glass_pred=table(factors_gbmlgg_glass$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis,factors_gbmlgg_glass$Dropped_Sample)

temp_table=data.frame(t(matrix(c(counts_tcga_og[1],counts_tcga_og[2],counts_tcga_og[3],counts_tcga_og[4],counts_tcga_og[6],counts_tcga_og[7],counts_tcga_og[5],0,
                      counts_tcga_pred[1,1],counts_tcga_pred[2,1],counts_tcga_pred[3,1],counts_tcga_pred[4,1],counts_tcga_pred[5,1],counts_tcga_pred[6,1],0,colSums(counts_tcga_pred)[2],
                      counts_cgga_pred[1,1],counts_cgga_pred[2,1],counts_cgga_pred[3,1],counts_cgga_pred[4,1],counts_cgga_pred[5,1],counts_cgga_pred[6,1],0,colSums(counts_cgga_pred)[2],
                      counts_glass_pred[1,1],counts_glass_pred[2,1],counts_glass_pred[3,1],counts_glass_pred[4,1],counts_glass_pred[5,1],counts_glass_pred[6,1],0,colSums(counts_glass_pred)[2]),ncol=4,nrow=8)))
temp_table=data.frame(c("Original TCGA Labels","Predicted TCGA Labels","Predicted CGGA Labels","Predicted GLASS Labels"),temp_table)
colnames(temp_table)=c("Dataset","Astrocytoma (Grade 2)","Astrocytoma (Grade 3)","Astrocytoma (Grade 4)","Glioblastoma (Grade 4)","Oligodendroglioma (Grade 2)","Oligodendroglioma (Grade 4)","Unlabelled (N/A)", "Misclassified (Dropped)")
temp_table[,c(2:ncol(temp_table))]=as.data.frame(mapply(function(x, y) paste0(x, " (", y, "%)"), temp_table[,c(2:ncol(temp_table))], round(100*temp_table[,c(2:ncol(temp_table))]/rowSums(temp_table[,c(2:ncol(temp_table))]),1)))
write.csv(temp_table,"Pub_Classify/Glioma_Label_Table.csv",row.names = F)


temp_coef=rbind(data.frame(Glioma_Type="Astrocytoma_G2",Ensembl_Gene_ID=names(coef(cvfit_glm_elasticnet)$Astrocytoma_G2),Estimate=coef(cvfit_glm_elasticnet)$Astrocytoma_G2),
      data.frame(Glioma_Type="Astrocytoma_G3",Ensembl_Gene_ID=names(coef(cvfit_glm_elasticnet)$Astrocytoma_G3),Estimate=coef(cvfit_glm_elasticnet)$Astrocytoma_G3),
      data.frame(Glioma_Type="Astrocytoma_G4",Ensembl_Gene_ID=names(coef(cvfit_glm_elasticnet)$Astrocytoma_G4),Estimate=coef(cvfit_glm_elasticnet)$Astrocytoma_G4),
      data.frame(Glioma_Type="Glioblastoma_G4",Ensembl_Gene_ID=names(coef(cvfit_glm_elasticnet)$Glioblastoma_G4),Estimate=coef(cvfit_glm_elasticnet)$Glioblastoma_G4),
      data.frame(Glioma_Type="Oligodendroglioma_G2",Ensembl_Gene_ID=names(coef(cvfit_glm_elasticnet)$Oligodendroglioma_G2),Estimate=coef(cvfit_glm_elasticnet)$Oligodendroglioma_G2),
      data.frame(Glioma_Type="Oligodendroglioma_G3",Ensembl_Gene_ID=names(coef(cvfit_glm_elasticnet)$Oligodendroglioma_G3),Estimate=coef(cvfit_glm_elasticnet)$Oligodendroglioma_G3))
temp_coef$id=1:nrow(temp_coef)
temp_coef=(merge(temp_coef,gene_labels,by="Ensembl_Gene_ID",all.x=T) %>% .[order(.$id),])[c(2,1,5,6,3)]
temp_coef[temp_coef$Ensembl_Gene_ID=="(Intercept)",]$HGNC_Symbol="(Intercept)"
temp_coef[temp_coef$Ensembl_Gene_ID=="(Intercept)",]$Entrez_Gene_ID="(Intercept)"
write.csv(temp_coef,"Pub_Classify/Glioma_Type_Coefficients.csv",row.names = F)

temp_data=rbind(data.frame(Sample_ID=rownames(factors_gbmlgg_tcga),Dataset="TCGA",Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis=factors_gbmlgg_tcga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis,Dropped_Sample=factors_gbmlgg_tcga$Dropped),
      data.frame(Sample_ID=rownames(factors_gbmlgg_cgga),Dataset="CGGA",Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis=factors_gbmlgg_cgga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis,Dropped_Sample=factors_gbmlgg_cgga$Dropped),
      data.frame(Sample_ID=rownames(factors_gbmlgg_glass),Dataset="GLASS",Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis=factors_gbmlgg_glass$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis,Dropped_Sample=factors_gbmlgg_glass$Dropped))
write.csv(temp_data,"Pub_Classify/Predicted_Glioma_Types.csv",row.names = F)

temp=factors_gbmlgg_tcga[,c("WHO_CNS5_2021_Diagnosis","Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis")]
colnames(temp)=c("WHO CNS5\nAnnotation","Predicted\nWHO CNS5 Annotation")
temp[temp=="Glioblastoma_G4"]="GBM (G4)"
temp[temp=="Astrocytoma_G2"]="A (G2)"
temp[temp=="Astrocytoma_G3"]="A (G3)"
temp[temp=="Astrocytoma_G4"]="A (G4)"
temp[temp=="Oligodendroglioma_G2"]="O (G2)"
temp[temp=="Oligodendroglioma_G3"]="O (G3)"
temp[temp=="N/A"]="Unlabelled"
temp$Dropped_Sample=factors_gbmlgg_tcga$Dropped_Sample
temp[temp$Dropped_Sample==T,]$`Predicted\nWHO CNS5 Annotation`="Misclassified"

tcga_sankey=temp %>% make_long(`WHO CNS5\nAnnotation`,`Predicted\nWHO CNS5 Annotation`)
tcga_sankey$node=factor(tcga_sankey$node,levels=rev(c("GBM (G4)","A (G2)","A (G3)","A (G4)","O (G2)","O (G3)","Unlabelled","Misclassified")))
fig2a=ggplot(tcga_sankey, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(flow.alpha = .6, node.color = "black") + geom_sankey_label(size = 1.5, color = "white", fill = "black") + scale_fill_viridis_d(drop = FALSE) + theme_sankey() + labs(x = NULL) + theme(legend.position = "none", plot.title = element_text(hjust = .5))+scale_fill_manual(values=c(`Unlabelled`="black",`Misclassified`="black",`GBM (G4)`="purple2",`A (G2)`="yellow2",`A (G3)`="orange2",`A (G4)`="firebrick2",`O (G2)`="lightblue2",`O (G3)`="blue2"))+theme(text = element_text(size = 6),plot.margin = margin(0.1, -1.5, 0.1, -1.5, "lines"))


temp=factors_gbmlgg_cgga[,c("CGGA_Diagnosis","Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis")]
colnames(temp)=c("CGGA\nAnnotation","Predicted\nWHO CNS5 Annotation")
temp[temp=="Glioblastoma_G4"]="GBM (G4)"
temp[temp=="Astrocytoma_G2"]="A (G2)"
temp[temp=="Astrocytoma_G3"]="A (G3)"
temp[temp=="Astrocytoma_G4"]="A (G4)"
temp[temp=="Oligodendroglioma_G2"]="O (G2)"
temp[temp=="Oligodendroglioma_G3"]="O (G3)"
temp[temp=="Glioblastoma_IV"]="GBM (IV)"
temp[temp=="Astrocytoma_II"]="A (II)"
temp[temp=="Astrocytoma_III"]="A (III)"
temp[temp=="Astrocytoma_IV"]="A (IV)"
temp[temp=="Oligodendroglioma_II"]="O (II)"
temp[temp=="Oligodendroglioma_III"]="O (III)"
temp[temp=="Anaplastic_Astrocytoma_III"]="AA (III)"
temp[temp=="Anaplastic_Oligoastrocytoma_III"]="AOA (III)"
temp[temp=="Anaplastic_Oligodendroglioma_III"]="AO (III)"
temp[temp=="Oligoastrocytoma_II"]="OA (II)"
temp[temp=="N/A"]="Unlabelled"
temp$Dropped_Sample=factors_gbmlgg_cgga$Dropped_Sample
temp[temp$Dropped_Sample==T,]$`Predicted\nWHO CNS5 Annotation`="Misclassified"

cgga_sankey=temp %>% make_long(`CGGA\nAnnotation`,`Predicted\nWHO CNS5 Annotation`)
cgga_sankey$node=factor(cgga_sankey$node,levels=rev(c("GBM (G4)","GBM (IV)","A (G2)","A (II)","A (G3)","A (III)","AA (III)","A (G4)","OA (II)","AOA (III)","O (G2)","O (II)","O (G3)","O (III)","AO (III)","Unlabelled","Misclassified")))
fig2d=ggplot(cgga_sankey, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(flow.alpha = .6, node.color = "black") + geom_sankey_label(size = 1.5, color = "white", fill = "black") + scale_fill_viridis_d(drop = FALSE) + theme_sankey() + labs(x = NULL) + theme(legend.position = "none", plot.title = element_text(hjust = .5))+scale_fill_manual(values=c(`Unlabelled`="black",`Misclassified`="black",`GBM (G4)`="purple2",`A (G2)`="yellow2",`A (G3)`="orange2",`A (G4)`="firebrick2",`O (G2)`="lightblue2",`O (G3)`="blue2",`GBM (IV)`="purple2",`A (II)`="yellow2",`A (III)`="orange2",`A (IV)`="firebrick2",`O (II)`="lightblue2",`O (III)`="blue2",`OA (II)`="palegreen2",`OA (III)`="seagreen2",`AA (III)`="orange1",`AA (II)`="yellow1",`AOA (II)`="palegreen1",`AOA (III)`="seagreen1",`AO (III)`="blue1"))+theme(text = element_text(size = 6),plot.margin = margin(0.1, -1.5, 0.1, -1.5, "lines"))

temp=factors_gbmlgg_glass[,c("GLASS_Diagnosis","Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis")]
colnames(temp)=c("GLASS\nAnnotation","Predicted\nWHO CNS5 Annotation")
temp[temp=="Glioblastoma_G4"]="GBM (G4)"
temp[temp=="Astrocytoma_G2"]="A (G2)"
temp[temp=="Astrocytoma_G3"]="A (G3)"
temp[temp=="Astrocytoma_G4"]="A (G4)"
temp[temp=="Oligodendroglioma_G2"]="O (G2)"
temp[temp=="Oligodendroglioma_G3"]="O (G3)"
temp[temp=="Glioblastoma_IV"]="GBM (IV)"
temp[temp=="Astrocytoma_II"]="A (II)"
temp[temp=="Astrocytoma_III"]="A (III)"
temp[temp=="Astrocytoma_IV"]="A (IV)"
temp[temp=="Oligodendroglioma_II"]="O (II)"
temp[temp=="Oligodendroglioma_III"]="O (III)"
temp[temp=="Anaplastic_Astrocytoma_II"]="AA (II)"
temp[temp=="Anaplastic_Astrocytoma_III"]="AA (III)"
temp[temp=="Oligoastrocytoma_II"]="OA (II)"
temp[temp=="Anaplastic_Oligoastrocytoma_II"]="AOA (II)"
temp[temp=="Anaplastic_Oligoastrocytoma_III"]="AOA (III)"
temp[temp=="Anaplastic_Oligodendroglioma_III"]="AO (III)"
temp[temp=="Diffuse_Astrocytoma_II"]="DA (II)"
temp[temp=="N/A"]="Unlabelled"
temp$Dropped_Sample=factors_gbmlgg_glass$Dropped_Sample
temp[temp$Dropped_Sample==T,]$`Predicted\nWHO CNS5 Annotation`="Misclassified"

glass_sankey=temp %>% make_long(`GLASS\nAnnotation`,`Predicted\nWHO CNS5 Annotation`)
glass_sankey$node=factor(glass_sankey$node,levels=rev(c("GBM (G4)","GBM (IV)","A (G2)","A (II)","AA (II)","DA (II)","A (G3)","A (III)","AA (III)","A (G4)","OA (II)","AOA (II)","AOA (III)","O (G2)","O (II)","O (G3)","O (III)","AO (III)","Misclassified")))
fig2g=ggplot(glass_sankey, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) + geom_sankey(flow.alpha = .6, node.color = "black") + geom_sankey_label(size = 1.5, color = "white", fill = "black") + scale_fill_viridis_d(drop = FALSE) + theme_sankey() + labs(x = NULL) + theme(legend.position = "none", plot.title = element_text(hjust = .5))+scale_fill_manual(values=c(`Misclassified`="black",`GBM (G4)`="purple2",`A (G2)`="yellow2",`A (G3)`="orange2",`A (G4)`="firebrick2",`O (G2)`="lightblue2",`O (G3)`="blue2",`GBM (IV)`="purple2",`A (II)`="yellow2",`A (III)`="orange2",`A (IV)`="firebrick2",`O (II)`="lightblue2",`O (III)`="blue2",`OA (II)`="palegreen2",`OA (III)`="seagreen2",`AA (III)`="orange1",`AA (II)`="yellow1",`AOA (II)`="palegreen1",`AOA (III)`="seagreen1",`AO (III)`="blue1",`DA (II)`="yellow3"))+theme(text = element_text(size = 6),plot.margin = margin(0.1, -1.5, 0.1, -1.5, "lines"))

##################################################################
#####UMAP (Bulk Glioma)#####
##################################################################
umap_tcga_gbmlgg=umap(t(cve_gene_gbmlgg_tcga[limma_filter(x,y,500),]))
df_umap <- data.frame(x = umap_tcga_gbmlgg$layout[,1], y = umap_tcga_gbmlgg$layout[,2], color1=factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis,Dropped=factors_gbmlgg_tcga$Dropped_Sample, color2=factors_gbmlgg_tcga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis)
df_umap[df_umap=="Glioblastoma_G4"]="GBM (G4)"
df_umap[df_umap=="Astrocytoma_G2"]="A (G2)"
df_umap[df_umap=="Astrocytoma_G3"]="A (G3)"
df_umap[df_umap=="Astrocytoma_G4"]="A (G4)"
df_umap[df_umap=="Oligodendroglioma_G2"]="O (G2)"
df_umap[df_umap=="Oligodendroglioma_G3"]="O (G3)"
df_umap[df_umap=="N/A"]="Unlabelled"
df_umap[df_umap$Dropped==T,]$color2="Misclassified"
df_umap$color1=factor(df_umap$color1,levels=c("GBM (G4)","A (G2)","A (G3)","A (G4)","O (G2)","O (G3)","Unlabelled","Misclassified"))
df_umap$color2=factor(df_umap$color2,levels=c("GBM (G4)","A (G2)","A (G3)","A (G4)","O (G2)","O (G3)","Unlabelled","Misclassified"))
fig2b=ggplot(df_umap, aes(x, y, colour = color1)) + geom_point(size=1)+theme_classic()+scale_color_discrete()+scale_color_manual(values=c(`Unlabelled`="black",`GBM (G4)`="purple2",`A (G2)`="yellow2",`A (G3)`="orange2",`A (G4)`="firebrick2",`O (G2)`="lightblue2",`O (G3)`="blue2"))+guides(color=guide_legend(title="WHO CNS5\nAnnotation", nrow = 3, byrow = TRUE))+xlab("UMAP1")+ylab("UMAP2")+theme(text = element_text(size = 6),legend.position = "top",legend.justification = "center",legend.key.height = unit(0.25, "lines"),legend.key.width = unit(0.1, "lines"),legend.margin=margin(c(0,25,0,0)))
fig2c=ggplot(df_umap, aes(x, y, colour = color2)) + geom_point(size=1)+theme_classic()+scale_color_discrete()+scale_color_manual(values=c(`Misclassified`="black",`GBM (G4)`="purple2",`A (G2)`="yellow2",`A (G3)`="orange2",`A (G4)`="firebrick2",`O (G2)`="lightblue2",`O (G3)`="blue2"))+guides(color=guide_legend(title="Predicted TCGA\nWHO CNS5\nAnnotation", nrow = 3, byrow = TRUE))+xlab("UMAP1")+ylab("UMAP2")+theme(text = element_text(size = 6),legend.position = "top",legend.justification = "center",legend.key.height = unit(0.25, "lines"),legend.key.width = unit(0.1, "lines"),legend.margin=margin(c(0,25,0,0)))

umap_cgga_gbmlgg=umap(t(cve_gene_gbmlgg_cgga[limma_filter(x,y,500),]))
df_umap <- data.frame(x = umap_cgga_gbmlgg$layout[,1], y = umap_cgga_gbmlgg$layout[,2], color1=factors_gbmlgg_cgga$CGGA_Diagnosis, color2=factors_gbmlgg_cgga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis,Dropped=factors_gbmlgg_cgga$Dropped_Sample)
df_umap[df_umap=="Glioblastoma_G4"]="GBM (G4)"
df_umap[df_umap=="Astrocytoma_G2"]="A (G2)"
df_umap[df_umap=="Astrocytoma_G3"]="A (G3)"
df_umap[df_umap=="Astrocytoma_G4"]="A (G4)"
df_umap[df_umap=="Oligodendroglioma_G2"]="O (G2)"
df_umap[df_umap=="Oligodendroglioma_G3"]="O (G3)"
df_umap[df_umap=="Glioblastoma_IV"]="GBM (IV)"
df_umap[df_umap=="Astrocytoma_II"]="A (II)"
df_umap[df_umap=="Astrocytoma_III"]="A (III)"
df_umap[df_umap=="Astrocytoma_IV"]="A (IV)"
df_umap[df_umap=="Oligodendroglioma_II"]="O (II)"
df_umap[df_umap=="Oligodendroglioma_III"]="O (III)"
df_umap[df_umap=="Anaplastic_Astrocytoma_III"]="AA (III)"
df_umap[df_umap=="Anaplastic_Oligoastrocytoma_III"]="AOA (III)"
df_umap[df_umap=="Anaplastic_Oligodendroglioma_III"]="AO (III)"
df_umap[df_umap=="Oligoastrocytoma_II"]="OA (II)"
df_umap[df_umap=="N/A"]="Unlabelled"
df_umap[df_umap$Dropped==T,]$color2="Misclassified"
df_umap$color1=factor(df_umap$color1,levels=c("GBM (G4)","GBM (IV)","A (G2)","A (II)","A (G3)","A (III)","AA (III)","A (G4)","OA (II)","AOA (III)","O (G2)","O (II)","O (G3)","O (III)","AO (III)","Unlabelled","Misclassified"))
df_umap$color2=factor(df_umap$color2,levels=c("GBM (G4)","GBM (IV)","A (G2)","A (II)","A (G3)","A (III)","AA (III)","A (G4)","OA (II)","AOA (III)","O (G2)","O (II)","O (G3)","O (III)","AO (III)","Unlabelled","Misclassified"))
fig2e=ggplot(df_umap, aes(x, y, colour = color1)) + geom_point(size=1)+theme_classic()+scale_color_discrete()+scale_color_manual(values=c(`Unlabelled`="black",`Misclassified`="black",`GBM (G4)`="purple2",`A (G2)`="yellow2",`A (G3)`="orange2",`A (G4)`="firebrick2",`O (G2)`="lightblue2",`O (G3)`="blue2",`GBM (IV)`="purple2",`A (II)`="yellow2",`A (III)`="orange2",`A (IV)`="firebrick2",`O (II)`="lightblue2",`O (III)`="blue2",`OA (II)`="palegreen2",`OA (III)`="seagreen2",`AA (III)`="orange1",`AA (II)`="yellow1",`AOA (II)`="palegreen1",`AOA (III)`="seagreen1",`AO (III)`="blue1"))+guides(color=guide_legend(title="CGGA\nAnnotation", nrow = 3, byrow = TRUE))+xlab("UMAP1")+ylab("UMAP2")+theme(text = element_text(size = 6),legend.position = "top",legend.justification = "center",legend.key.height = unit(0.25, "lines"),legend.key.width = unit(0.1, "lines"),legend.margin=margin(c(0,25,0,0)))
fig2f=ggplot(df_umap, aes(x, y, colour = color2)) + geom_point(size=1)+theme_classic()+scale_color_discrete()+scale_color_manual(values=c(`Unlabelled`="black",`Misclassified`="black",`GBM (G4)`="purple2",`A (G2)`="yellow2",`A (G3)`="orange2",`A (G4)`="firebrick2",`O (G2)`="lightblue2",`O (G3)`="blue2",`GBM (IV)`="purple2",`A (II)`="yellow2",`A (III)`="orange2",`A (IV)`="firebrick2",`O (II)`="lightblue2",`O (III)`="blue2",`OA (II)`="palegreen2",`OA (III)`="seagreen2",`AA (III)`="orange1",`AA (II)`="yellow1",`AOA (II)`="palegreen1",`AOA (III)`="seagreen1",`AO (III)`="blue1"))+guides(color=guide_legend(title="Predicted CGGA\nWHO CNS5\nAnnotation", nrow = 3, byrow = TRUE))+xlab("UMAP1")+ylab("UMAP2")+theme(text = element_text(size = 6),legend.position = "top",legend.justification = "center",legend.key.height = unit(0.25, "lines"),legend.key.width = unit(0.1, "lines"),legend.margin=margin(c(0,25,0,0)))

umap_glass_gbmlgg=umap(t(cve_gene_gbmlgg_glass[limma_filter(x,y,500),]))
df_umap <- data.frame(x = umap_glass_gbmlgg$layout[,1], y = umap_glass_gbmlgg$layout[,2], color1=factors_gbmlgg_glass$GLASS_Diagnosis, color2=factors_gbmlgg_glass$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis,Dropped=factors_gbmlgg_glass$Dropped_Sample)
df_umap[df_umap=="Glioblastoma_G4"]="GBM (G4)"
df_umap[df_umap=="Astrocytoma_G2"]="A (G2)"
df_umap[df_umap=="Astrocytoma_G3"]="A (G3)"
df_umap[df_umap=="Astrocytoma_G4"]="A (G4)"
df_umap[df_umap=="Oligodendroglioma_G2"]="O (G2)"
df_umap[df_umap=="Oligodendroglioma_G3"]="O (G3)"
df_umap[df_umap=="Glioblastoma_IV"]="GBM (IV)"
df_umap[df_umap=="Astrocytoma_II"]="A (II)"
df_umap[df_umap=="Astrocytoma_III"]="A (III)"
df_umap[df_umap=="Astrocytoma_IV"]="A (IV)"
df_umap[df_umap=="Oligodendroglioma_II"]="O (II)"
df_umap[df_umap=="Oligodendroglioma_III"]="O (III)"
df_umap[df_umap=="Anaplastic_Astrocytoma_II"]="AA (II)"
df_umap[df_umap=="Anaplastic_Astrocytoma_III"]="AA (III)"
df_umap[df_umap=="Oligoastrocytoma_II"]="OA (II)"
df_umap[df_umap=="Anaplastic_Oligoastrocytoma_II"]="AOA (II)"
df_umap[df_umap=="Anaplastic_Oligoastrocytoma_III"]="AOA (III)"
df_umap[df_umap=="Anaplastic_Oligodendroglioma_III"]="AO (III)"
df_umap[df_umap=="Diffuse_Astrocytoma_II"]="DA (II)"
df_umap[df_umap=="N/A"]="Unlabelled"
df_umap[df_umap$Dropped==T,]$color2="Misclassified"
df_umap$color1=factor(df_umap$color1,levels=c("GBM (G4)","GBM (IV)","A (G2)","A (II)","AA (II)","DA (II)","A (G3)","A (III)","AA (III)","A (G4)","OA (II)","AOA (II)","AOA (III)","O (G2)","O (II)","O (G3)","O (III)","AO (III)","Misclassified"))
df_umap$color2=factor(df_umap$color2,levels=c("GBM (G4)","GBM (IV)","A (G2)","A (II)","AA (II)","DA (II)","A (G3)","A (III)","AA (III)","A (G4)","OA (II)","AOA (II)","AOA (III)","O (G2)","O (II)","O (G3)","O (III)","AO (III)","Misclassified"))
fig2h=ggplot(df_umap, aes(x, y, colour = color1)) + geom_point(size=1)+theme_classic()+scale_color_discrete()+scale_color_manual(values=c(`Misclassified`="black",`GBM (G4)`="purple2",`A (G2)`="yellow2",`A (G3)`="orange2",`A (G4)`="firebrick2",`O (G2)`="lightblue2",`O (G3)`="blue2",`GBM (IV)`="purple2",`A (II)`="yellow2",`A (III)`="orange2",`A (IV)`="firebrick2",`O (II)`="lightblue2",`O (III)`="blue2",`OA (II)`="palegreen2",`OA (III)`="seagreen2",`AA (III)`="orange1",`AA (II)`="yellow1",`AOA (II)`="palegreen1",`AOA (III)`="seagreen1",`AO (III)`="blue1",`DA (II)`="yellow3"))+guides(color=guide_legend(title="GLASS\nAnnotation", nrow = 3, byrow = TRUE))+xlab("UMAP1")+ylab("UMAP2")+theme(text = element_text(size = 6),legend.position = "top",legend.justification = "center",legend.key.height = unit(0.25, "lines"),legend.key.width = unit(0.1, "lines"),legend.margin=margin(c(0,25,0,0)))
fig2i=ggplot(df_umap, aes(x, y, colour = color2)) + geom_point(size=1)+theme_classic()+scale_color_discrete()+scale_color_manual(values=c(`Misclassified`="black",`GBM (G4)`="purple2",`A (G2)`="yellow2",`A (G3)`="orange2",`A (G4)`="firebrick2",`O (G2)`="lightblue2",`O (G3)`="blue2",`GBM (IV)`="purple2",`A (II)`="yellow2",`A (III)`="orange2",`A (IV)`="firebrick2",`O (II)`="lightblue2",`O (III)`="blue2",`OA (II)`="palegreen2",`OA (III)`="seagreen2",`AA (III)`="orange1",`AA (II)`="yellow1",`AOA (II)`="palegreen1",`AOA (III)`="seagreen1",`AO (III)`="blue1",`DA (II)`="yellow3"))+guides(color=guide_legend(title="Predicted GLASS\nWHO CNS5\nAnnotation", nrow = 3, byrow = TRUE))+xlab("UMAP1")+ylab("UMAP2")+theme(text = element_text(size = 6),legend.position = "top",legend.justification = "center",legend.key.height = unit(0.25, "lines"),legend.key.width = unit(0.1, "lines"),legend.margin=margin(c(0,25,0,0)))

##################################################################
#####Survival (Bulk Glioma)#####
##################################################################
temp=rbind(data.frame(Dataset="TCGA",factors_gbmlgg_tcga[factors_gbmlgg_tcga$Recurrence_Status=="Primary",c("Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis","Dropped_Sample","OS","OS_Time")]),data.frame(Dataset="CGGA",factors_gbmlgg_cgga[factors_gbmlgg_cgga$Recurrence_Status=="Primary",c("Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis","Dropped_Sample","OS","OS_Time")]),data.frame(Dataset="GLASS",factors_gbmlgg_glass[factors_gbmlgg_glass$Recurrence_Status=="Primary",c("Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis","Dropped_Sample","OS")],OS_Time=30.417*factors_gbmlgg_glass[factors_gbmlgg_glass$Recurrence_Status=="Primary",c("OS_Time_Months")]))
temp[temp=="Glioblastoma_G4"]="GBM (G4)"
temp[temp=="Astrocytoma_G2"]="A (G2)"
temp[temp=="Astrocytoma_G3"]="A (G3)"
temp[temp=="Astrocytoma_G4"]="A (G4)"
temp[temp=="Oligodendroglioma_G2"]="O (G2)"
temp[temp=="Oligodendroglioma_G3"]="O (G3)"
temp[temp=="N/A"]="Unlabelled"
temp[temp$Dropped_Sample==T,]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis="Misclassified"
temp$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis=factor(temp$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis,levels=c("GBM (G4)","A (G2)","A (G3)","A (G4)","O (G2)","O (G3)","Unlabelled","Misclassified"))

fit1 <- survfit(Surv(OS_Time/30.417, OS) ~ Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, data = temp[temp$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis!="Misclassified",])
fig3a=(ggsurvplot(fit1,legend.title = "Glioma Type",surv.median.line = "hv",pval = TRUE, break.time.by = 12,risk.table = T, conf.int = F,palette = c("purple2","yellow2","orange2","firebrick2","lightblue2","blue2"),legend.labs=c("GBM (G4)", "A (G2)", "A (G3)", "A (G4)", "O (G2)", "O (G3)"),pval.coord = c(120, 0.95),test.for.trend = T,xlim=c(0,168),pval.size = 3,ggtheme = theme_classic(),size=0.5,censor.size=3))[[1]]+theme(text = element_text(size = 6),legend.key.height = unit(0.25, "lines"),legend.position = "top",legend.justification = "center")+xlab("Time (Months)")+ylab("Probability of Survival")


fit2 <- summary(coxph(Surv(OS_Time/30.417, OS) ~ Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis+strata(Dataset), data = temp[temp$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis!="Misclassified",]))
df=data.frame(Glioma_Type=c("GBM (G4)","A (G2)","A (G3)","A (G4)","O (G2)","O (G3)"),HR=c(1,fit2$coefficients[,2][1:5]),lower_CI=c(1,fit2$conf.int[,3][1:5]),upper_CI=c(1,fit2$conf.int[,4][1:5]),pval=c(1,fit2$coefficients[,5][1:5]))
df$Glioma_Type=paste0(df$Glioma_Type, " (HR=",round(df$HR,2),")\n(",round(df$lower_CI,2),"-",round(df$upper_CI,2),")")
df$Glioma_Type[1]="GBM (G4)\n(Reference)"
fig3b=ggplot(df, aes(x = factor(Glioma_Type,levels=rev(unique(Glioma_Type))), y = HR, ymin = lower_CI, ymax = upper_CI)) + geom_pointrange(size = 0.5) + geom_hline(yintercept = 1, linetype = "dashed", color = "red") + coord_flip() + scale_y_log10(breaks = c(0.05, 0.1, 0.2, 0.5, 1)) + labs( x = "Glioma Type", y = "Hazard Ratio (95% CI)") + theme_classic() + theme(panel.grid.minor = element_blank(), axis.text.y = element_text(hjust = 1),text = element_text(size = 6))


temp=data.frame(Dataset="TCGA",factors_gbmlgg_tcga[factors_gbmlgg_tcga$Recurrence_Status=="Primary",c("WHO_CNS5_2021_Diagnosis","Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis","Dropped_Sample","OS","OS_Time")])
temp[temp=="Glioblastoma_G4"]="GBM (G4)"
temp[temp=="Astrocytoma_G2"]="A (G2)"
temp[temp=="Astrocytoma_G3"]="A (G3)"
temp[temp=="Astrocytoma_G4"]="A (G4)"
temp[temp=="Oligodendroglioma_G2"]="O (G2)"
temp[temp=="Oligodendroglioma_G3"]="O (G3)"
temp[temp=="N/A"]="Unlabelled"
temp[temp$Dropped_Sample==T,]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis="Misclassified"
temp$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis=factor(temp$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis,levels=c("GBM (G4)","A (G2)","A (G3)","A (G4)","O (G2)","O (G3)"))
temp$WHO_CNS5_2021_Diagnosis=factor(temp$WHO_CNS5_2021_Diagnosis,levels=c("GBM (G4)","A (G2)","A (G3)","A (G4)","O (G2)","O (G3)"))

fit1 <- survfit(Surv(OS_Time/30.417, OS) ~ WHO_CNS5_2021_Diagnosis, data = temp[temp$WHO_CNS5_2021_Diagnosis!="Unlabelled",])
figs3a=(ggsurvplot(fit1,legend.title = "WHO CNS5\nAnnotation",surv.median.line = "hv",pval = TRUE, break.time.by = 12,risk.table = T, conf.int = F,palette = c("purple2","yellow2","orange2","firebrick2","lightblue2","blue2"),legend.labs=c("GBM (G4)", "A (G2)", "A (G3)", "A (G4)", "O (G2)", "O (G3)"),pval.coord = c(120, 0.95),test.for.trend = T,xlim=c(0,168),pval.size = 3,ggtheme = theme_classic(),size=0.5,censor.size=3))[[1]]+theme(text = element_text(size = 6),legend.key.height = unit(0.25, "lines"),legend.position = "top",legend.justification = "center")+xlab("Time (Months)")+ylab("Probability of Survival")

fit1 <- survfit(Surv(OS_Time/30.417, OS) ~ Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, data = temp[temp$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis!="Misclassified",])
figs3b=(ggsurvplot(fit1,legend.title = "Predicted\nWHO CNS5\nAnnotation",surv.median.line = "hv",pval = TRUE, break.time.by = 12,risk.table = T, conf.int = F,palette = c("purple2","yellow2","orange2","firebrick2","lightblue2","blue2"),legend.labs=c("GBM (G4)", "A (G2)", "A (G3)", "A (G4)", "O (G2)", "O (G3)"),pval.coord = c(120, 0.95),test.for.trend = T,xlim=c(0,168),pval.size = 3,ggtheme = theme_classic(),size=0.5,censor.size=3))[[1]]+theme(text = element_text(size = 6),legend.key.height = unit(0.25, "lines"),legend.position = "top",legend.justification = "center")+xlab("Time (Months)")+ylab("Probability of Survival")


temp=data.frame(Dataset="CGGA",factors_gbmlgg_cgga[factors_gbmlgg_cgga$Recurrence_Status=="Primary",c("CGGA_Diagnosis","Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis","Dropped_Sample","OS","OS_Time")])
temp[temp=="Glioblastoma_G4"]="GBM (G4)"
temp[temp=="Astrocytoma_G2"]="A (G2)"
temp[temp=="Astrocytoma_G3"]="A (G3)"
temp[temp=="Astrocytoma_G4"]="A (G4)"
temp[temp=="Oligodendroglioma_G2"]="O (G2)"
temp[temp=="Oligodendroglioma_G3"]="O (G3)"
temp[temp=="Glioblastoma_IV"]="GBM (IV)"
temp[temp=="Astrocytoma_II"]="A (II)"
temp[temp=="Astrocytoma_III"]="A (III)"
temp[temp=="Astrocytoma_IV"]="A (IV)"
temp[temp=="Oligodendroglioma_II"]="O (II)"
temp[temp=="Oligodendroglioma_III"]="O (III)"
temp[temp=="Anaplastic_Astrocytoma_III"]="AA (III)"
temp[temp=="Anaplastic_Oligoastrocytoma_III"]="AOA (III)"
temp[temp=="Anaplastic_Oligodendroglioma_III"]="AO (III)"
temp[temp=="Oligoastrocytoma_II"]="OA (II)"
temp[temp=="N/A"]="Unlabelled"
temp[temp$Dropped_Sample==T,]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis="Misclassified"
temp$CGGA_Diagnosis=factor(temp$CGGA_Diagnosis,levels=c("GBM (IV)","A (II)","AA (III)","OA (II)","AOA (III)","O (II)","AO (III)"))
temp$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis=factor(temp$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis,levels=c("GBM (G4)","A (G2)","A (G3)","A (G4)","O (G2)","O (G3)"))

fit1 <- survfit(Surv(OS_Time/30.417, OS) ~ CGGA_Diagnosis, data = temp)
figs3c=(ggsurvplot(fit1,legend.title = "CGGA\nAnnotation",surv.median.line = "hv",pval = TRUE, break.time.by = 12,risk.table = T, conf.int = F,palette = c("purple2","yellow2","orange1","palegreen2","seagreen1","lightblue2","blue1"),legend.labs=c("GBM (IV)", "A (II)", "AA (III)", "OA (II)", "AOA (II)", "O (II)", "AO (III)"),pval.coord = c(120, 0.95),test.for.trend = T,xlim=c(0,168),pval.size = 3,ggtheme = theme_classic(),size=0.5,censor.size=3))[[1]]+theme(text = element_text(size = 6),legend.key.height = unit(0.25, "lines"),legend.position = "top",legend.justification = "center")+xlab("Time (Months)")+ylab("Probability of Survival")

fit1 <- survfit(Surv(OS_Time/30.417, OS) ~ Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, data = temp[temp$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis!="Misclassified",])
figs3d=(ggsurvplot(fit1,legend.title = "Predicted\nWHO CNS5\nAnnotation",surv.median.line = "hv",pval = TRUE, break.time.by = 12,risk.table = T, conf.int = F,palette = c("purple2","yellow2","orange2","firebrick2","lightblue2","blue2"),legend.labs=c("GBM (G4)", "A (G2)", "A (G3)", "A (G4)", "O (G2)", "O (G3)"),pval.coord = c(120, 0.95),test.for.trend = T,xlim=c(0,168),pval.size = 3,ggtheme = theme_classic(),size=0.5,censor.size=3))[[1]]+theme(text = element_text(size = 6),legend.key.height = unit(0.25, "lines"),legend.position = "top",legend.justification = "center")+xlab("Time (Months)")+ylab("Probability of Survival")


temp=data.frame(Dataset="GLASS",factors_gbmlgg_glass[factors_gbmlgg_glass$Recurrence_Status=="Primary",c("GLASS_Diagnosis","Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis","Dropped_Sample","OS")],OS_Time=30.417*factors_gbmlgg_glass[factors_gbmlgg_glass$Recurrence_Status=="Primary",c("OS_Time_Months")])
temp[temp=="Glioblastoma_G4"]="GBM (G4)"
temp[temp=="Astrocytoma_G2"]="A (G2)"
temp[temp=="Astrocytoma_G3"]="A (G3)"
temp[temp=="Astrocytoma_G4"]="A (G4)"
temp[temp=="Oligodendroglioma_G2"]="O (G2)"
temp[temp=="Oligodendroglioma_G3"]="O (G3)"
temp[temp=="Glioblastoma_IV"]="GBM (IV)"
temp[temp=="Astrocytoma_II"]="A (II)"
temp[temp=="Astrocytoma_III"]="A (III)"
temp[temp=="Astrocytoma_IV"]="A (IV)"
temp[temp=="Oligodendroglioma_II"]="O (II)"
temp[temp=="Oligodendroglioma_III"]="O (III)"
temp[temp=="Anaplastic_Astrocytoma_II"]="AA (II)"
temp[temp=="Anaplastic_Astrocytoma_III"]="AA (III)"
temp[temp=="Oligoastrocytoma_II"]="OA (II)"
temp[temp=="Anaplastic_Oligoastrocytoma_II"]="AOA (II)"
temp[temp=="Anaplastic_Oligoastrocytoma_III"]="AOA (III)"
temp[temp=="Anaplastic_Oligodendroglioma_III"]="AO (III)"
temp[temp=="Diffuse_Astrocytoma_II"]="DA (II)"
temp[temp=="N/A"]="Unlabelled"
temp[temp$Dropped_Sample==T,]$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis="Misclassified"
temp$GLASS_Diagnosis=factor(temp$GLASS_Diagnosis,levels=c("GBM (IV)","A (III)","AA (II)","AA (III)","OA (II)","AOA (II)","AOA (III)","O (II)"))
temp$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis=factor(temp$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis,levels=c("GBM (G4)","A (G2)","A (G3)","A (G4)","O (G2)","O (G3)"))

fit1 <- survfit(Surv(OS_Time/30.417, OS) ~ GLASS_Diagnosis, data = temp)
figs3e=(ggsurvplot(fit1,legend.title = "GLASS\nAnnotation",surv.median.line = "hv",pval = TRUE, break.time.by = 12,risk.table = T, conf.int = F,palette = c("purple2","orange2","yellow1","orange1","palegreen2","palegreen1","seagreen1","lightblue2"),legend.labs=c("GBM (IV)", "A (III)", "AA (II)", "AA (III)", "OA (II)", "AOA (II)", "AOA (III)", "O (II)"),pval.coord = c(120, 0.95),test.for.trend = T,xlim=c(0,168),pval.size = 3,ggtheme = theme_classic(),size=0.5,censor.size=3))[[1]]+theme(text = element_text(size = 6),legend.key.height = unit(0.25, "lines"),legend.position = "top",legend.justification = "center")+xlab("Time (Months)")+ylab("Probability of Survival")

fit1 <- survfit(Surv(OS_Time/30.417, OS) ~ Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, data = temp[temp$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis!="Misclassified",])
figs3f=(ggsurvplot(fit1,legend.title = "Predicted\nWHO CNS5\nAnnotation",surv.median.line = "hv",pval = TRUE, break.time.by = 12,risk.table = T, conf.int = F,palette = c("purple2","yellow2","orange2","firebrick2","lightblue2","blue2"),legend.labs=c("GBM (G4)", "A (G2)", "A (G3)", "A (G4)", "O (G2)", "O (G3)"),pval.coord = c(120, 0.95),test.for.trend = T,xlim=c(0,168),pval.size = 3,ggtheme = theme_classic(),size=0.5,censor.size=3))[[1]]+theme(text = element_text(size = 6),legend.key.height = unit(0.25, "lines"),legend.position = "top",legend.justification = "center")+xlab("Time (Months)")+ylab("Probability of Survival")

##################################################################
#####DEGs (Bulk Glioma)#####
##################################################################
sub_factors_gbmlgg_tcga=factors_gbmlgg_tcga[!grepl("N/A",factors_gbmlgg_tcga$WHO_CNS5_2021_Diagnosis) & !is.na(factors_gbmlgg_tcga$Diagnosis_Age),]
sub_factors_gbmlgg_tcga$Scaled_Diagnosis_Age=scale(sub_factors_gbmlgg_tcga$Diagnosis_Age)
sub_cts_gene_gbmlgg_tcga=cts_gene_gbmlgg_tcga[,colnames(cts_gene_gbmlgg_tcga) %in% rownames(sub_factors_gbmlgg_tcga)]
sub_dge_gene_gbmlgg_tcga <- DGEList(sub_cts_gene_gbmlgg_tcga)
sub_dge_gene_gbmlgg_tcga <- calcNormFactors(sub_dge_gene_gbmlgg_tcga)

design <- ~0+WHO_CNS5_2021_Diagnosis+Recurrence_Status+Read_Length+Sex+Scaled_Diagnosis_Age
L <- makeContrastsDream(design, sub_factors_gbmlgg_tcga, contrasts = c(AstroII = "WHO_CNS5_2021_DiagnosisAstrocytoma_G2 - (WHO_CNS5_2021_DiagnosisAstrocytoma_G3+WHO_CNS5_2021_DiagnosisAstrocytoma_G4+WHO_CNS5_2021_DiagnosisGlioblastoma_G4+WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/5",
                                                                       AstroIII = "WHO_CNS5_2021_DiagnosisAstrocytoma_G3 - (WHO_CNS5_2021_DiagnosisAstrocytoma_G2+WHO_CNS5_2021_DiagnosisAstrocytoma_G4+WHO_CNS5_2021_DiagnosisGlioblastoma_G4+WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/5",
                                                                       AstroIV = "WHO_CNS5_2021_DiagnosisAstrocytoma_G4 - (WHO_CNS5_2021_DiagnosisAstrocytoma_G3+WHO_CNS5_2021_DiagnosisAstrocytoma_G2+WHO_CNS5_2021_DiagnosisGlioblastoma_G4+WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/5",
                                                                       GBMIV = "WHO_CNS5_2021_DiagnosisGlioblastoma_G4 - (WHO_CNS5_2021_DiagnosisAstrocytoma_G3+WHO_CNS5_2021_DiagnosisAstrocytoma_G4+WHO_CNS5_2021_DiagnosisAstrocytoma_G2+WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/5",
                                                                       OligoII = "WHO_CNS5_2021_DiagnosisOligodendroglioma_G2 - (WHO_CNS5_2021_DiagnosisAstrocytoma_G3+WHO_CNS5_2021_DiagnosisAstrocytoma_G4+WHO_CNS5_2021_DiagnosisGlioblastoma_G4+WHO_CNS5_2021_DiagnosisAstrocytoma_G2+WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/5",
                                                                       OligoIII = "WHO_CNS5_2021_DiagnosisOligodendroglioma_G3 - (WHO_CNS5_2021_DiagnosisAstrocytoma_G3+WHO_CNS5_2021_DiagnosisAstrocytoma_G4+WHO_CNS5_2021_DiagnosisGlioblastoma_G4+WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+WHO_CNS5_2021_DiagnosisAstrocytoma_G2)/5"))
vDream <- voomWithDreamWeights(sub_dge_gene_gbmlgg_tcga, design, sub_factors_gbmlgg_tcga, BPPARAM = SnowParam(workers = 10,exportglobals = F),plot = F)
mmfit <- dream(vDream, design, sub_factors_gbmlgg_tcga, L, BPPARAM = SnowParam(workers = 10,exportglobals = F))
efit <- eBayes(mmfit)
res_gene_tcga_gbmlgg_og_AstroII=data.frame(topTable(efit, coef="AstroII", sort.by = "none", n = Inf,lfc = 0))
res_gene_tcga_gbmlgg_og_AstroIII=data.frame(topTable(efit, coef="AstroIII", sort.by = "none", n = Inf,lfc = 0))
res_gene_tcga_gbmlgg_og_AstroIV=data.frame(topTable(efit, coef="AstroIV", sort.by = "none", n = Inf,lfc = 0))
res_gene_tcga_gbmlgg_og_GBMIV=data.frame(topTable(efit, coef="GBMIV", sort.by = "none", n = Inf,lfc = 0))
res_gene_tcga_gbmlgg_og_OligoII=data.frame(topTable(efit, coef="OligoII", sort.by = "none", n = Inf,lfc = 0))
res_gene_tcga_gbmlgg_og_OligoIII=data.frame(topTable(efit, coef="OligoIII", sort.by = "none", n = Inf,lfc = 0))
colnames(res_gene_tcga_gbmlgg_og_AstroII)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_tcga_gbmlgg_og_AstroIII)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_tcga_gbmlgg_og_AstroIV)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_tcga_gbmlgg_og_GBMIV)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_tcga_gbmlgg_og_OligoII)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_tcga_gbmlgg_og_OligoIII)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")


sub_factors_gbmlgg_tcga=factors_gbmlgg_tcga[factors_gbmlgg_tcga$Dropped!=T & !is.na(factors_gbmlgg_tcga$Diagnosis_Age) & !is.na(factors_gbmlgg_tcga$Sex),]
sub_factors_gbmlgg_tcga$Scaled_Diagnosis_Age=scale(sub_factors_gbmlgg_tcga$Diagnosis_Age)
sub_cts_gene_gbmlgg_tcga=cts_gene_gbmlgg_tcga[,colnames(cts_gene_gbmlgg_tcga) %in% rownames(sub_factors_gbmlgg_tcga)]
sub_dge_gene_gbmlgg_tcga <- DGEList(sub_cts_gene_gbmlgg_tcga)
sub_dge_gene_gbmlgg_tcga <- calcNormFactors(sub_dge_gene_gbmlgg_tcga)

design <- ~0+Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis+Recurrence_Status+Read_Length+Sex+Scaled_Diagnosis_Age
L <- makeContrastsDream(design, sub_factors_gbmlgg_tcga, contrasts = c(AstroII = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G2 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/5",
                                                                       AstroIII = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/5",
                                                                       AstroIV = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/5",
                                                                       GBMIV = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/5",
                                                                       OligoII = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/5",
                                                                       OligoIII = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G2)/5"))
vDream <- voomWithDreamWeights(sub_dge_gene_gbmlgg_tcga, design, sub_factors_gbmlgg_tcga, BPPARAM = SnowParam(workers = 10,exportglobals = F),plot = F)
mmfit <- dream(vDream, design, sub_factors_gbmlgg_tcga, L, BPPARAM = SnowParam(workers = 10,exportglobals = F))
efit <- eBayes(mmfit)
res_gene_tcga_gbmlgg_pred_AstroII=data.frame(topTable(efit, coef="AstroII", sort.by = "none", n = Inf,lfc = 0))
res_gene_tcga_gbmlgg_pred_AstroIII=data.frame(topTable(efit, coef="AstroIII", sort.by = "none", n = Inf,lfc = 0))
res_gene_tcga_gbmlgg_pred_AstroIV=data.frame(topTable(efit, coef="AstroIV", sort.by = "none", n = Inf,lfc = 0))
res_gene_tcga_gbmlgg_pred_GBMIV=data.frame(topTable(efit, coef="GBMIV", sort.by = "none", n = Inf,lfc = 0))
res_gene_tcga_gbmlgg_pred_OligoII=data.frame(topTable(efit, coef="OligoII", sort.by = "none", n = Inf,lfc = 0))
res_gene_tcga_gbmlgg_pred_OligoIII=data.frame(topTable(efit, coef="OligoIII", sort.by = "none", n = Inf,lfc = 0))
colnames(res_gene_tcga_gbmlgg_pred_AstroII)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_tcga_gbmlgg_pred_AstroIII)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_tcga_gbmlgg_pred_AstroIV)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_tcga_gbmlgg_pred_GBMIV)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_tcga_gbmlgg_pred_OligoII)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_tcga_gbmlgg_pred_OligoIII)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")


sub_factors_gbmlgg_cgga=factors_gbmlgg_cgga[factors_gbmlgg_cgga$Dropped!=T & !is.na(factors_gbmlgg_cgga$Diagnosis_Age) & !is.na(factors_gbmlgg_cgga$Sex),]
sub_factors_gbmlgg_cgga$Scaled_Diagnosis_Age=scale(sub_factors_gbmlgg_cgga$Diagnosis_Age)
sub_cts_gene_gbmlgg_cgga=cts_gene_gbmlgg_cgga[,colnames(cts_gene_gbmlgg_cgga) %in% rownames(sub_factors_gbmlgg_cgga)]
sub_dge_gene_gbmlgg_cgga <- DGEList(sub_cts_gene_gbmlgg_cgga)
sub_dge_gene_gbmlgg_cgga <- calcNormFactors(sub_dge_gene_gbmlgg_cgga)

design <- ~0+Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis+Recurrence_Status+CGGA_Dataset_Detailed+Sex+Scaled_Diagnosis_Age
L <- makeContrastsDream(design, sub_factors_gbmlgg_cgga, contrasts = c(AstroII = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G2 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/5",
                                                                       AstroIII = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/5",
                                                                       AstroIV = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/5",
                                                                       GBMIV = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/5",
                                                                       OligoII = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/5",
                                                                       OligoIII = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G2)/5"))
vDream <- voomWithDreamWeights(sub_dge_gene_gbmlgg_cgga, design, sub_factors_gbmlgg_cgga, BPPARAM = SnowParam(workers = 10,exportglobals = F),plot = F)
mmfit <- dream(vDream, design, sub_factors_gbmlgg_cgga, L, BPPARAM = SnowParam(workers = 10,exportglobals = F))
efit <- eBayes(mmfit)
res_gene_cgga_gbmlgg_pred_AstroII=data.frame(topTable(efit, coef="AstroII", sort.by = "none", n = Inf,lfc = 0))
res_gene_cgga_gbmlgg_pred_AstroIII=data.frame(topTable(efit, coef="AstroIII", sort.by = "none", n = Inf,lfc = 0))
res_gene_cgga_gbmlgg_pred_AstroIV=data.frame(topTable(efit, coef="AstroIV", sort.by = "none", n = Inf,lfc = 0))
res_gene_cgga_gbmlgg_pred_GBMIV=data.frame(topTable(efit, coef="GBMIV", sort.by = "none", n = Inf,lfc = 0))
res_gene_cgga_gbmlgg_pred_OligoII=data.frame(topTable(efit, coef="OligoII", sort.by = "none", n = Inf,lfc = 0))
res_gene_cgga_gbmlgg_pred_OligoIII=data.frame(topTable(efit, coef="OligoIII", sort.by = "none", n = Inf,lfc = 0))
colnames(res_gene_cgga_gbmlgg_pred_AstroII)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_cgga_gbmlgg_pred_AstroIII)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_cgga_gbmlgg_pred_AstroIV)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_cgga_gbmlgg_pred_GBMIV)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_cgga_gbmlgg_pred_OligoII)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_cgga_gbmlgg_pred_OligoIII)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")


sub_factors_gbmlgg_glass=factors_gbmlgg_glass[factors_gbmlgg_glass$Dropped!=T & !is.na(factors_gbmlgg_glass$Diagnosis_Age) & !is.na(factors_gbmlgg_glass$Sex),]
sub_factors_gbmlgg_glass$Scaled_Diagnosis_Age=scale(sub_factors_gbmlgg_glass$Diagnosis_Age)
sub_cts_gene_gbmlgg_glass=cts_gene_gbmlgg_glass[,colnames(cts_gene_gbmlgg_glass) %in% rownames(sub_factors_gbmlgg_glass)]
sub_dge_gene_gbmlgg_glass <- DGEList(sub_cts_gene_gbmlgg_glass)
sub_dge_gene_gbmlgg_glass <- calcNormFactors(sub_dge_gene_gbmlgg_glass)

design <- ~0+Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis+Recurrence_Status+Tissue_Source_Site_Detailed+Sex+Scaled_Diagnosis_Age
L <- makeContrastsDream(design, sub_factors_gbmlgg_glass, contrasts = c(AstroIII = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/4",
                                                                       AstroIV = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/4",
                                                                       GBMIV = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/4",
                                                                       OligoII = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3)/4",
                                                                       OligoIII = "Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G3 - (Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G3+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisAstrocytoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisGlioblastoma_G4+Predicted_GLM_ElasticNet_WHO_CNS5_2021_DiagnosisOligodendroglioma_G2)/4"))
vDream <- voomWithDreamWeights(sub_dge_gene_gbmlgg_glass, design, sub_factors_gbmlgg_glass, BPPARAM = SnowParam(workers = 10,exportglobals = F),plot = F)
mmfit <- dream(vDream, design, sub_factors_gbmlgg_glass, L, BPPARAM = SnowParam(workers = 10,exportglobals = F))
efit <- eBayes(mmfit)
res_gene_glass_gbmlgg_pred_AstroIII=data.frame(topTable(efit, coef="AstroIII", sort.by = "none", n = Inf,lfc = 0))
res_gene_glass_gbmlgg_pred_AstroIV=data.frame(topTable(efit, coef="AstroIV", sort.by = "none", n = Inf,lfc = 0))
res_gene_glass_gbmlgg_pred_GBMIV=data.frame(topTable(efit, coef="GBMIV", sort.by = "none", n = Inf,lfc = 0))
res_gene_glass_gbmlgg_pred_OligoII=data.frame(topTable(efit, coef="OligoII", sort.by = "none", n = Inf,lfc = 0))
res_gene_glass_gbmlgg_pred_OligoIII=data.frame(topTable(efit, coef="OligoIII", sort.by = "none", n = Inf,lfc = 0))
colnames(res_gene_glass_gbmlgg_pred_AstroIII)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_glass_gbmlgg_pred_AstroIV)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_glass_gbmlgg_pred_GBMIV)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_glass_gbmlgg_pred_OligoII)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")
colnames(res_gene_glass_gbmlgg_pred_OligoIII)=c("log2FC","log2Expression","t_Statistic","pval","padj","logOdds_B")

x_astroii=data.frame(t(data.frame(res_gene_tcga_gbmlgg_og_AstroII$padj<0.05,res_gene_tcga_gbmlgg_pred_AstroII$padj<0.05,res_gene_cgga_gbmlgg_pred_AstroII$padj<0.05)))
colnames(x_astroii)=rownames(sub_dge_gene_gbmlgg_cgga)
rownames(x_astroii)=c(paste0("Original TCGA: ",paste0(table(res_gene_tcga_gbmlgg_og_AstroII$padj<0.05)[2]), " DEGs"), paste0("Predicted TCGA: ",paste0(table(res_gene_tcga_gbmlgg_pred_AstroII$padj<0.05)[2]), " DEGs"), paste0("Predicted CGGA: ",paste0(table(res_gene_cgga_gbmlgg_pred_AstroII$padj<0.05)[2]), " DEGs"))
grp_astroii=x_astroii %>% as_tibble(rownames = "Annotation") %>% gather(Gene, Member, -Annotation) %>%filter(Member) %>% dplyr::select(-Member) %>% group_by(Gene) %>% dplyr::summarize(Annotation = list(Annotation))
figs4b=grp_astroii %>% ggplot(aes(x = Annotation))+scale_x_upset(n_intersections = 7)+geom_bar(fill="yellow2")+xlab("")+ylab("# of Intersecting DEGs")+theme_classic()+ theme(text = element_text(size = 6), axis.text.y=element_text(size = 6),legend.position="none")+theme_combmatrix(combmatrix.panel.point.color.fill = "yellow2",combmatrix.panel.line.color = "yellow2",combmatrix.label.height = unit(2,"lines"),combmatrix.label.make_space = T,combmatrix.panel.point.size = 1,combmatrix.panel.line.size = 0.3)

x_astroiii=data.frame(t(data.frame(res_gene_tcga_gbmlgg_og_AstroIII$padj<0.05,res_gene_tcga_gbmlgg_pred_AstroIII$padj<0.05,res_gene_cgga_gbmlgg_pred_AstroIII$padj<0.05,res_gene_glass_gbmlgg_pred_AstroIII$padj<0.05)))
colnames(x_astroiii)=rownames(sub_dge_gene_gbmlgg_cgga)
rownames(x_astroiii)=c(paste0("Original TCGA: ",paste0(table(res_gene_tcga_gbmlgg_og_AstroIII$padj<0.05)[2]), " DEGs"), paste0("Predicted TCGA: ",paste0(table(res_gene_tcga_gbmlgg_pred_AstroIII$padj<0.05)[2]), " DEGs"), paste0("Predicted CGGA: ",paste0(table(res_gene_cgga_gbmlgg_pred_AstroIII$padj<0.05)[2]), " DEGs"), paste0("Predicted GLASS: ",paste0(table(res_gene_glass_gbmlgg_pred_AstroIII$padj<0.05)[2]), " DEGs"))
grp_astroiii=x_astroiii %>% as_tibble(rownames = "Annotation") %>% gather(Gene, Member, -Annotation) %>%filter(Member) %>% dplyr::select(-Member) %>% group_by(Gene) %>% dplyr::summarize(Annotation = list(Annotation))
figs4c=grp_astroiii %>% ggplot(aes(x = Annotation))+scale_x_upset(n_intersections = 7)+geom_bar(fill="orange2")+xlab("")+ylab("# of Intersecting DEGs")+theme_classic()+ theme(text = element_text(size = 6), axis.text.y=element_text(size = 6),legend.position="none")+theme_combmatrix(combmatrix.panel.point.color.fill = "orange2",combmatrix.panel.line.color = "orange2",combmatrix.label.height = unit(2,"lines"),combmatrix.label.make_space = T,combmatrix.panel.point.size = 1,combmatrix.panel.line.size = 0.3)

x_astroiv=data.frame(t(data.frame(res_gene_tcga_gbmlgg_og_AstroIV$padj<0.05,res_gene_tcga_gbmlgg_pred_AstroIV$padj<0.05,res_gene_cgga_gbmlgg_pred_AstroIV$padj<0.05,res_gene_glass_gbmlgg_pred_AstroIV$padj<0.05)))
colnames(x_astroiv)=rownames(sub_dge_gene_gbmlgg_cgga)
rownames(x_astroiv)=c(paste0("Original TCGA: ",paste0(table(res_gene_tcga_gbmlgg_og_AstroIV$padj<0.05)[2]), " DEGs"), paste0("Predicted TCGA: ",paste0(table(res_gene_tcga_gbmlgg_pred_AstroIV$padj<0.05)[2]), " DEGs"), paste0("Predicted CGGA: ",paste0(table(res_gene_cgga_gbmlgg_pred_AstroIV$padj<0.05)[2]), " DEGs"), paste0("Predicted GLASS: ",paste0(table(res_gene_glass_gbmlgg_pred_AstroIV$padj<0.05)[2]), " DEGs"))
grp_astroiv=x_astroiv %>% as_tibble(rownames = "Annotation") %>% gather(Gene, Member, -Annotation) %>%filter(Member) %>% dplyr::select(-Member) %>% group_by(Gene) %>% dplyr::summarize(Annotation = list(Annotation))
figs4d=grp_astroiv %>% ggplot(aes(x = Annotation))+scale_x_upset(n_intersections = 7)+geom_bar(fill="firebrick2")+xlab("")+ylab("# of Intersecting DEGs")+theme_classic()+ theme(text = element_text(size = 6), axis.text.y=element_text(size = 6),legend.position="none")+theme_combmatrix(combmatrix.panel.point.color.fill = "firebrick2",combmatrix.panel.line.color = "firebrick2",combmatrix.label.height = unit(2,"lines"),combmatrix.label.make_space = T,combmatrix.panel.point.size = 1,combmatrix.panel.line.size = 0.3)

x_gbmiv=data.frame(t(data.frame(res_gene_tcga_gbmlgg_og_GBMIV$padj<0.05,res_gene_tcga_gbmlgg_pred_GBMIV$padj<0.05,res_gene_cgga_gbmlgg_pred_GBMIV$padj<0.05,res_gene_glass_gbmlgg_pred_GBMIV$padj<0.05)))
colnames(x_gbmiv)=rownames(sub_dge_gene_gbmlgg_cgga)
rownames(x_gbmiv)=c(paste0("Original TCGA: ",paste0(table(res_gene_tcga_gbmlgg_og_GBMIV$padj<0.05)[2]), " DEGs"), paste0("Predicted TCGA: ",paste0(table(res_gene_tcga_gbmlgg_pred_GBMIV$padj<0.05)[2]), " DEGs"), paste0("Predicted CGGA: ",paste0(table(res_gene_cgga_gbmlgg_pred_GBMIV$padj<0.05)[2]), " DEGs"), paste0("Predicted GLASS: ",paste0(table(res_gene_glass_gbmlgg_pred_GBMIV$padj<0.05)[2]), " DEGs"))
grp_gbmiv=x_gbmiv %>% as_tibble(rownames = "Annotation") %>% gather(Gene, Member, -Annotation) %>%filter(Member) %>% dplyr::select(-Member) %>% group_by(Gene) %>% dplyr::summarize(Annotation = list(Annotation))
figs4a=grp_gbmiv %>% ggplot(aes(x = Annotation))+scale_x_upset(n_intersections = 7)+geom_bar(fill="purple2")+xlab("")+ylab("# of Intersecting DEGs")+theme_classic()+ theme(text = element_text(size = 6), axis.text.y=element_text(size = 6),legend.position="none")+theme_combmatrix(combmatrix.panel.point.color.fill = "purple2",combmatrix.panel.line.color = "purple2",combmatrix.label.height = unit(2,"lines"),combmatrix.label.make_space = T,combmatrix.panel.point.size = 1,combmatrix.panel.line.size = 0.3)

x_oligoii=data.frame(t(data.frame(res_gene_tcga_gbmlgg_og_OligoII$padj<0.05,res_gene_tcga_gbmlgg_pred_OligoII$padj<0.05,res_gene_cgga_gbmlgg_pred_OligoII$padj<0.05,res_gene_glass_gbmlgg_pred_OligoII$padj<0.05)))
colnames(x_oligoii)=rownames(sub_dge_gene_gbmlgg_cgga)
rownames(x_oligoii)=c(paste0("Original TCGA: ",paste0(table(res_gene_tcga_gbmlgg_og_OligoII$padj<0.05)[2]), " DEGs"), paste0("Predicted TCGA: ",paste0(table(res_gene_tcga_gbmlgg_pred_OligoII$padj<0.05)[2]), " DEGs"), paste0("Predicted CGGA: ",paste0(table(res_gene_cgga_gbmlgg_pred_OligoII$padj<0.05)[2]), " DEGs"), paste0("Predicted GLASS: ",paste0(table(res_gene_glass_gbmlgg_pred_OligoII$padj<0.05)[2]), " DEGs"))
grp_oligoii=x_oligoii %>% as_tibble(rownames = "Annotation") %>% gather(Gene, Member, -Annotation) %>%filter(Member) %>% dplyr::select(-Member) %>% group_by(Gene) %>% dplyr::summarize(Annotation = list(Annotation))
figs4e=grp_oligoii %>% ggplot(aes(x = Annotation))+scale_x_upset(n_intersections = 7)+geom_bar(fill="lightblue2")+xlab("")+ylab("# of Intersecting DEGs")+theme_classic()+ theme(text = element_text(size = 6), axis.text.y=element_text(size = 6),legend.position="none")+theme_combmatrix(combmatrix.panel.point.color.fill = "lightblue2",combmatrix.panel.line.color = "lightblue2",combmatrix.label.height = unit(2,"lines"),combmatrix.label.make_space = T,combmatrix.panel.point.size = 1,combmatrix.panel.line.size = 0.3)

x_oligoiii=data.frame(t(data.frame(res_gene_tcga_gbmlgg_og_OligoIII$padj<0.05,res_gene_tcga_gbmlgg_pred_OligoIII$padj<0.05,res_gene_cgga_gbmlgg_pred_OligoIII$padj<0.05,res_gene_glass_gbmlgg_pred_OligoIII$padj<0.05)))
colnames(x_oligoiii)=rownames(sub_dge_gene_gbmlgg_cgga)
rownames(x_oligoiii)=c(paste0("Original TCGA: ",paste0(table(res_gene_tcga_gbmlgg_og_OligoIII$padj<0.05)[2]), " DEGs"), paste0("Predicted TCGA: ",paste0(table(res_gene_tcga_gbmlgg_pred_OligoIII$padj<0.05)[2]), " DEGs"), paste0("Predicted CGGA: ",paste0(table(res_gene_cgga_gbmlgg_pred_OligoIII$padj<0.05)[2]), " DEGs"), paste0("Predicted GLASS: ",paste0(table(res_gene_glass_gbmlgg_pred_OligoIII$padj<0.05)[2]), " DEGs"))
grp_oligoiii=x_oligoiii %>% as_tibble(rownames = "Annotation") %>% gather(Gene, Member, -Annotation) %>%filter(Member) %>% dplyr::select(-Member) %>% group_by(Gene) %>% dplyr::summarize(Annotation = list(Annotation))
figs4f=grp_oligoiii %>% ggplot(aes(x = Annotation))+scale_x_upset(n_intersections = 7)+geom_bar(fill="blue2")+xlab("")+ylab("# of Intersecting DEGs")+theme_classic()+ theme(text = element_text(size = 6), axis.text.y=element_text(size = 6),legend.position="none")+theme_combmatrix(combmatrix.panel.point.color.fill = "blue2",combmatrix.panel.line.color = "blue2",combmatrix.label.height = unit(2,"lines"),combmatrix.label.make_space = T,combmatrix.panel.point.size = 1,combmatrix.panel.line.size = 0.3)


figs5b=ggplot(melt(round(cor(data.frame("Original TCGA"=res_gene_tcga_gbmlgg_og_AstroII$t_Statistic,"Predicted TCGA"=res_gene_tcga_gbmlgg_pred_AstroII$t_Statistic,"Predicted CGGA"=res_gene_cgga_gbmlgg_pred_AstroII$t_Statistic)[limma_filter(x,y,500),],method = "spearman"),2)), aes(Var2, Var1, fill = value))+geom_tile(color = "white")+labs(x ="Moderated t-Statistic", y = "Moderated t-Statistic")+scale_fill_gradient2(low = "white", high = "yellow2", mid = "white", midpoint = 0, limit = c(0,1), space = "Lab", name="Spearman\nCorrelation") +  scale_y_discrete(labels=c("Original TCGA","Predicted TCGA","Predicted CGGA")) + scale_x_discrete(labels=c("Original TCGA","Predicted TCGA","Predicted CGGA")) +theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ coord_fixed() + geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) + theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.ticks = element_blank(),legend.position = "top",legend.justification = "left", legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.75, "lines"),legend.title = element_text(margin = margin(r = 25)),text = element_text(size = 6),legend.margin=margin(c(0,0,0,-25)))
figs5c=ggplot(melt(round(cor(data.frame("Original TCGA"=res_gene_tcga_gbmlgg_og_AstroIII$t_Statistic,"Predicted TCGA"=res_gene_tcga_gbmlgg_pred_AstroIII$t_Statistic,"Predicted CGGA"=res_gene_cgga_gbmlgg_pred_AstroIII$t_Statistic,"GLASS TCGA"=res_gene_glass_gbmlgg_pred_AstroIII$t_Statistic)[limma_filter(x,y,500),],method = "spearman"),2)), aes(Var2, Var1, fill = value))+geom_tile(color = "white")+labs(x ="Moderated t-Statistic", y = "Moderated t-Statistic")+scale_fill_gradient2(low = "white", high = "orange2", mid = "white", midpoint = 0, limit = c(0,1), space = "Lab", name="Spearman\nCorrelation") +  scale_y_discrete(labels=c("Original TCGA","Predicted TCGA","Predicted CGGA","Predicted GLASS")) + scale_x_discrete(labels=c("Original TCGA","Predicted TCGA","Predicted CGGA","Predicted GLASS")) +theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ coord_fixed() + geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) + theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.ticks = element_blank(),legend.position = "top",legend.justification = "left", legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.75, "lines"),legend.title = element_text(margin = margin(r = 25)),text = element_text(size = 6),legend.margin=margin(c(0,0,0,-25)))
figs5d=ggplot(melt(round(cor(data.frame("Original TCGA"=res_gene_tcga_gbmlgg_og_AstroIV$t_Statistic,"Predicted TCGA"=res_gene_tcga_gbmlgg_pred_AstroIV$t_Statistic,"Predicted CGGA"=res_gene_cgga_gbmlgg_pred_AstroIV$t_Statistic,"GLASS TCGA"=res_gene_glass_gbmlgg_pred_AstroIV$t_Statistic)[limma_filter(x,y,500),],method = "spearman"),2)), aes(Var2, Var1, fill = value))+geom_tile(color = "white")+labs(x ="Moderated t-Statistic", y = "Moderated t-Statistic")+scale_fill_gradient2(low = "white", high = "firebrick2", mid = "white", midpoint = 0, limit = c(0,1), space = "Lab", name="Spearman\nCorrelation") +  scale_y_discrete(labels=c("Original TCGA","Predicted TCGA","Predicted CGGA","Predicted GLASS")) + scale_x_discrete(labels=c("Original TCGA","Predicted TCGA","Predicted CGGA","Predicted GLASS")) +theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ coord_fixed() + geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) + theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.ticks = element_blank(),legend.position = "top",legend.justification = "left", legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.75, "lines"),legend.title = element_text(margin = margin(r = 25)),text = element_text(size = 6),legend.margin=margin(c(0,0,0,-25)))
figs5a=ggplot(melt(round(cor(data.frame("Original TCGA"=res_gene_tcga_gbmlgg_og_GBMIV$t_Statistic,"Predicted TCGA"=res_gene_tcga_gbmlgg_pred_GBMIV$t_Statistic,"Predicted CGGA"=res_gene_cgga_gbmlgg_pred_GBMIV$t_Statistic,"GLASS TCGA"=res_gene_glass_gbmlgg_pred_GBMIV$t_Statistic)[limma_filter(x,y,500),],method = "spearman"),2)), aes(Var2, Var1, fill = value))+geom_tile(color = "white")+labs(x ="Moderated t-Statistic", y = "Moderated t-Statistic")+scale_fill_gradient2(low = "white", high = "purple2", mid = "white", midpoint = 0, limit = c(0,1), space = "Lab", name="Spearman\nCorrelation") +  scale_y_discrete(labels=c("Original TCGA","Predicted TCGA","Predicted CGGA","Predicted GLASS")) + scale_x_discrete(labels=c("Original TCGA","Predicted TCGA","Predicted CGGA","Predicted GLASS")) +theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ coord_fixed() + geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) + theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.ticks = element_blank(),legend.position = "top",legend.justification = "left", legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.75, "lines"),legend.title = element_text(margin = margin(r = 25)),text = element_text(size = 6),legend.margin=margin(c(0,0,0,-25)))
figs5e=ggplot(melt(round(cor(data.frame("Original TCGA"=res_gene_tcga_gbmlgg_og_OligoII$t_Statistic,"Predicted TCGA"=res_gene_tcga_gbmlgg_pred_OligoII$t_Statistic,"Predicted CGGA"=res_gene_cgga_gbmlgg_pred_OligoII$t_Statistic,"GLASS TCGA"=res_gene_glass_gbmlgg_pred_OligoII$t_Statistic)[limma_filter(x,y,500),],method = "spearman"),2)), aes(Var2, Var1, fill = value))+geom_tile(color = "white")+labs(x ="Moderated t-Statistic", y = "Moderated t-Statistic")+scale_fill_gradient2(low = "white", high = "lightblue2", mid = "white", midpoint = 0, limit = c(0,1), space = "Lab", name="Spearman\nCorrelation") +  scale_y_discrete(labels=c("Original TCGA","Predicted TCGA","Predicted CGGA","Predicted GLASS")) + scale_x_discrete(labels=c("Original TCGA","Predicted TCGA","Predicted CGGA","Predicted GLASS")) +theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ coord_fixed() + geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) + theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.ticks = element_blank(),legend.position = "top",legend.justification = "left", legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.75, "lines"),legend.title = element_text(margin = margin(r = 25)),text = element_text(size = 6),legend.margin=margin(c(0,0,0,-25)))
figs5f=ggplot(melt(round(cor(data.frame("Original TCGA"=res_gene_tcga_gbmlgg_og_OligoIII$t_Statistic,"Predicted TCGA"=res_gene_tcga_gbmlgg_pred_OligoIII$t_Statistic,"Predicted CGGA"=res_gene_cgga_gbmlgg_pred_OligoIII$t_Statistic,"GLASS TCGA"=res_gene_glass_gbmlgg_pred_OligoIII$t_Statistic)[limma_filter(x,y,500),],method = "spearman"),2)), aes(Var2, Var1, fill = value))+geom_tile(color = "white")+labs(x ="Moderated t-Statistic", y = "Moderated t-Statistic")+scale_fill_gradient2(low = "white", high = "blue2", mid = "white", midpoint = 0, limit = c(0,1), space = "Lab", name="Spearman\nCorrelation") +  scale_y_discrete(labels=c("Original TCGA","Predicted TCGA","Predicted CGGA","Predicted GLASS")) + scale_x_discrete(labels=c("Original TCGA","Predicted TCGA","Predicted CGGA","Predicted GLASS")) +theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+ coord_fixed() + geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) + theme(panel.grid.major = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.ticks = element_blank(),legend.position = "top",legend.justification = "left", legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.75, "lines"),legend.title = element_text(margin = margin(r = 25)),text = element_text(size = 6),legend.margin=margin(c(0,0,0,-25)))

column_ha = HeatmapAnnotation(
  `Glioma Type`=as.character(c(rep(c("A (G2)","A (G3)","A (G4)","GBM (G4)","O (G2)","O (G3)"),2),c("A (G3)","A (G4)","GBM (G4)","O (G2)","O (G3)"))),
  `Dataset`=as.character(c(rep("TCGA",6),rep("CGGA",6),rep("GLASS",5))),
  col = list(`Glioma Type` = c(`A (G2)`="yellow2", `A (G3)`="orange2", `A (G4)`="firebrick2", `GBM (G4)`="purple2", `O (G2)`="lightblue2", `O (G3)`="blue2"),`Dataset` = c(`TCGA`="darkred", `CGGA`="darkgreen", `GLASS`="darkblue")),
  annotation_legend_param = list(`Glioma Type` = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5)),`Dataset` = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5))),
  annotation_name_gp = gpar(fontsize = 6)
  )
temp=data.frame(data.frame(Pred_TCGA_AstroII=res_gene_tcga_gbmlgg_pred_AstroII$t_Statistic,Pred_TCGA_AstroIII=res_gene_tcga_gbmlgg_pred_AstroIII$t_Statistic,Pred_TCGA_AstroIV=res_gene_tcga_gbmlgg_pred_AstroIV$t_Statistic,Pred_TCGA_GBMIV=res_gene_tcga_gbmlgg_pred_GBMIV$t_Statistic,Pred_TCGA_OligoII=res_gene_tcga_gbmlgg_pred_OligoII$t_Statistic,Pred_TCGA_OligoIII=res_gene_tcga_gbmlgg_pred_OligoIII$t_Statistic,Pred_CGGA_AstroII=res_gene_cgga_gbmlgg_pred_AstroII$t_Statistic,Pred_CGGA_AstroIII=res_gene_cgga_gbmlgg_pred_AstroIII$t_Statistic,Pred_CGGA_AstroIV=res_gene_cgga_gbmlgg_pred_AstroIV$t_Statistic,Pred_CGGA_GBMIV=res_gene_cgga_gbmlgg_pred_GBMIV$t_Statistic,Pred_CGGA_OligoII=res_gene_cgga_gbmlgg_pred_OligoII$t_Statistic,Pred_CGGA_OligoIII=res_gene_cgga_gbmlgg_pred_OligoIII$t_Statistic,Pred_GLASS_AstroIII=res_gene_glass_gbmlgg_pred_AstroIII$t_Statistic,Pred_GLASS_AstroIV=res_gene_glass_gbmlgg_pred_AstroIV$t_Statistic,Pred_GLASS_GBMIV=res_gene_glass_gbmlgg_pred_GBMIV$t_Statistic,Pred_GLASS_OligoII=res_gene_glass_gbmlgg_pred_OligoII$t_Statistic,Pred_GLASS_OligoIII=res_gene_glass_gbmlgg_pred_OligoIII$t_Statistic))
fig5a=ggplotify::as.grob(Heatmap(scale(temp)[limma_filter(x,y,500),],show_column_names=F,show_row_names = F,top_annotation=packLegend(column_ha),name="Scaled\nModerated\nt-Statistic",heatmap_legend_param = list(title_gp = gpar(fontsize = 6),labels_gp = gpar(fontsize = 5))))

##################################################################
#####Pathway Enrichment (Bulk Glioma)#####
##################################################################
dbs <- c("GO_Biological_Process_2025","CellMarker_2024")

enrich_gbmlgg_AstroII_up=enrichr(gene=gene_labels[gene_labels$Ensembl_Gene_ID %in% unique(colnames(x_astroii[,colSums(x_astroii)==3])) & gene_labels$Ensembl_Gene_ID %in% rownames(res_gene_tcga_gbmlgg_pred_AstroII[res_gene_tcga_gbmlgg_pred_AstroII$log2FC>0,]),]$HGNC_Symbol, dbs)
gc()
enrich_gbmlgg_AstroII_down=enrichr(gene=gene_labels[gene_labels$Ensembl_Gene_ID %in% unique(colnames(x_astroii[,colSums(x_astroii)==3])) & gene_labels$Ensembl_Gene_ID %in% rownames(res_gene_tcga_gbmlgg_pred_AstroII[res_gene_tcga_gbmlgg_pred_AstroII$log2FC<0,]),]$HGNC_Symbol, dbs)
gc()
enrich_gbmlgg_AstroII=rbind(data.frame(head(enrich_gbmlgg_AstroII_up[[1]][order(enrich_gbmlgg_AstroII_up[[1]]$P.value),],7),Sign="Upregulated"),data.frame(head(enrich_gbmlgg_AstroII_down[[1]][order(enrich_gbmlgg_AstroII_down[[1]]$P.value),],7),Sign="Downregulated"))
enrich_gbmlgg_AstroII$Term=sub('(.*)\\s\\(.*', '\\1',enrich_gbmlgg_AstroII$Term)
enrich_gbmlgg_AstroII$Adjusted.P.value=signif(enrich_gbmlgg_AstroII$Adjusted.P.value,2)
enrich_gbmlgg_AstroII$P.value=round(log10(enrich_gbmlgg_AstroII$Adjusted.P.value),1)
enrich_gbmlgg_AstroII$Term=str_wrap(enrich_gbmlgg_AstroII$Term, width = 40)
fig4b=plotEnrich(enrich_gbmlgg_AstroII, showTerms = 14, numChar = 100,xlab="Enriched Terms",ylab="Gene Count",orderBy = "P.value")+theme_classic()+ facet_grid(.~Sign)+scale_fill_continuous(high="yellow4",low="yellow1",labels = scales::label_number(accuracy = 1))+ theme(text = element_text(size = 6), plot.title = element_blank(),legend.position = "bottom",legend.justification = "left",legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.4, "lines"),legend.margin=margin(c(0,0,0,0)))+ guides(fill=guide_colourbar(title="Log10 Adjusted\np-value"))

enrich_gbmlgg_AstroIII_up=enrichr(gene=gene_labels[gene_labels$Ensembl_Gene_ID %in% unique(colnames(x_astroiii[,colSums(x_astroiii)>=3])) & gene_labels$Ensembl_Gene_ID %in% rownames(res_gene_tcga_gbmlgg_pred_AstroIII[res_gene_tcga_gbmlgg_pred_AstroIII$log2FC>0,]),]$HGNC_Symbol, dbs)
gc()
enrich_gbmlgg_AstroIII_down=enrichr(gene=gene_labels[gene_labels$Ensembl_Gene_ID %in% unique(colnames(x_astroiii[,colSums(x_astroiii)>=3])) & gene_labels$Ensembl_Gene_ID %in% rownames(res_gene_tcga_gbmlgg_pred_AstroIII[res_gene_tcga_gbmlgg_pred_AstroIII$log2FC<0,]),]$HGNC_Symbol, dbs)
gc()
enrich_gbmlgg_AstroIII=rbind(data.frame(head(enrich_gbmlgg_AstroIII_up[[1]][order(enrich_gbmlgg_AstroIII_up[[1]]$P.value),],7),Sign="Upregulated"),data.frame(head(enrich_gbmlgg_AstroIII_down[[1]][order(enrich_gbmlgg_AstroIII_down[[1]]$P.value),],7),Sign="Downregulated"))
enrich_gbmlgg_AstroIII$Term=sub('(.*)\\s\\(.*', '\\1',enrich_gbmlgg_AstroIII$Term)
enrich_gbmlgg_AstroIII$Adjusted.P.value=signif(enrich_gbmlgg_AstroIII$Adjusted.P.value,2)
enrich_gbmlgg_AstroIII$P.value=round(log10(enrich_gbmlgg_AstroIII$Adjusted.P.value),1)
enrich_gbmlgg_AstroIII$Term=str_wrap(enrich_gbmlgg_AstroIII$Term, width = 40)
fig4c=plotEnrich(enrich_gbmlgg_AstroIII, showTerms = 14, numChar = 100,xlab="Enriched Terms",ylab="Gene Count",orderBy = "P.value")+theme_classic()+ facet_grid(.~Sign)+scale_fill_continuous(high="orange4",low="orange1",labels = scales::label_number(accuracy = 1))+ theme(text = element_text(size = 6), plot.title = element_blank(),legend.position = "bottom",legend.justification = "left",legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.4, "lines"),legend.margin=margin(c(0,0,0,0)))+ guides(fill=guide_colourbar(title="Log10 Adjusted\np-value"))

enrich_gbmlgg_AstroIV_up=enrichr(gene=gene_labels[gene_labels$Ensembl_Gene_ID %in% unique(colnames(x_astroiv[,colSums(x_astroiv)>=3])) & gene_labels$Ensembl_Gene_ID %in% rownames(res_gene_tcga_gbmlgg_pred_AstroIV[res_gene_tcga_gbmlgg_pred_AstroIV$log2FC>0,]),]$HGNC_Symbol, dbs)
gc()
enrich_gbmlgg_AstroIV_down=enrichr(gene=gene_labels[gene_labels$Ensembl_Gene_ID %in% unique(colnames(x_astroiv[,colSums(x_astroiv)>=3])) & gene_labels$Ensembl_Gene_ID %in% rownames(res_gene_tcga_gbmlgg_pred_AstroIV[res_gene_tcga_gbmlgg_pred_AstroIV$log2FC<0,]),]$HGNC_Symbol, dbs)
gc()
enrich_gbmlgg_AstroIV=rbind(data.frame(head(enrich_gbmlgg_AstroIV_up[[1]][order(enrich_gbmlgg_AstroIV_up[[1]]$P.value),],7),Sign="Upregulated"),data.frame(head(enrich_gbmlgg_AstroIV_down[[1]][order(enrich_gbmlgg_AstroIV_down[[1]]$P.value),],7),Sign="Downregulated"))
enrich_gbmlgg_AstroIV$Term=sub('(.*)\\s\\(.*', '\\1',enrich_gbmlgg_AstroIV$Term)
enrich_gbmlgg_AstroIV$Adjusted.P.value=signif(enrich_gbmlgg_AstroIV$Adjusted.P.value,2)
enrich_gbmlgg_AstroIV$P.value=round(log10(enrich_gbmlgg_AstroIV$Adjusted.P.value),1)
enrich_gbmlgg_AstroIV$Term=str_wrap(enrich_gbmlgg_AstroIV$Term, width = 40)
fig4d=plotEnrich(enrich_gbmlgg_AstroIV, showTerms = 14, numChar = 100,xlab="Enriched Terms",ylab="Gene Count",orderBy = "P.value")+theme_classic()+ facet_grid(.~Sign)+scale_fill_continuous(high="firebrick4",low="firebrick1",labels = scales::label_number(accuracy = 1))+ theme(text = element_text(size = 6), plot.title = element_blank(),legend.position = "bottom",legend.justification = "left",legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.4, "lines"),legend.margin=margin(c(0,0,0,0)))+ guides(fill=guide_colourbar(title="Log10 Adjusted\np-value"))

enrich_gbmlgg_GBMIV_up=enrichr(gene=gene_labels[gene_labels$Ensembl_Gene_ID %in% unique(colnames(x_gbmiv[,colSums(x_gbmiv)>=3])) & gene_labels$Ensembl_Gene_ID %in% rownames(res_gene_tcga_gbmlgg_pred_GBMIV[res_gene_tcga_gbmlgg_pred_GBMIV$log2FC>0,]),]$HGNC_Symbol, dbs)
gc()
enrich_gbmlgg_GBMIV_down=enrichr(gene=gene_labels[gene_labels$Ensembl_Gene_ID %in% unique(colnames(x_gbmiv[,colSums(x_gbmiv)>=3])) & gene_labels$Ensembl_Gene_ID %in% rownames(res_gene_tcga_gbmlgg_pred_GBMIV[res_gene_tcga_gbmlgg_pred_GBMIV$log2FC<0,]),]$HGNC_Symbol, dbs)
gc()
enrich_gbmlgg_GBMIV=rbind(data.frame(head(enrich_gbmlgg_GBMIV_up[[1]][order(enrich_gbmlgg_GBMIV_up[[1]]$P.value),],7),Sign="Upregulated"),data.frame(head(enrich_gbmlgg_GBMIV_down[[1]][order(enrich_gbmlgg_GBMIV_down[[1]]$P.value),],7),Sign="Downregulated"))
enrich_gbmlgg_GBMIV$Term=sub('(.*)\\s\\(.*', '\\1',enrich_gbmlgg_GBMIV$Term)
enrich_gbmlgg_GBMIV$Adjusted.P.value=signif(enrich_gbmlgg_GBMIV$Adjusted.P.value,2)
enrich_gbmlgg_GBMIV$P.value=round(log10(enrich_gbmlgg_GBMIV$Adjusted.P.value),1)
enrich_gbmlgg_GBMIV$Term=str_wrap(enrich_gbmlgg_GBMIV$Term, width = 40)
fig4a=plotEnrich(enrich_gbmlgg_GBMIV, showTerms = 14, numChar = 100,xlab="Enriched Terms",ylab="Gene Count",orderBy = "P.value")+theme_classic()+ facet_grid(.~Sign)+scale_fill_continuous(high="purple4",low="purple1",labels = scales::label_number(accuracy = 1))+ theme(text = element_text(size = 6), plot.title = element_blank(),legend.position = "bottom",legend.justification = "left",legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.4, "lines"),legend.margin=margin(c(0,0,0,0)))+ guides(fill=guide_colourbar(title="Log10 Adjusted\np-value"))

enrich_gbmlgg_OligoII_up=enrichr(gene=gene_labels[gene_labels$Ensembl_Gene_ID %in% unique(colnames(x_oligoii[,colSums(x_oligoii)>=3])) & gene_labels$Ensembl_Gene_ID %in% rownames(res_gene_tcga_gbmlgg_pred_OligoII[res_gene_tcga_gbmlgg_pred_OligoII$log2FC>0,]),]$HGNC_Symbol, dbs)
gc()
enrich_gbmlgg_OligoII_down=enrichr(gene=gene_labels[gene_labels$Ensembl_Gene_ID %in% unique(colnames(x_oligoii[,colSums(x_oligoii)>=3])) & gene_labels$Ensembl_Gene_ID %in% rownames(res_gene_tcga_gbmlgg_pred_OligoII[res_gene_tcga_gbmlgg_pred_OligoII$log2FC<0,]),]$HGNC_Symbol, dbs)
gc()
enrich_gbmlgg_OligoII=rbind(data.frame(head(enrich_gbmlgg_OligoII_up[[1]][order(enrich_gbmlgg_OligoII_up[[1]]$P.value),],7),Sign="Upregulated"),data.frame(head(enrich_gbmlgg_OligoII_down[[1]][order(enrich_gbmlgg_OligoII_down[[1]]$P.value),],7),Sign="Downregulated"))
enrich_gbmlgg_OligoII$Term=sub('(.*)\\s\\(.*', '\\1',enrich_gbmlgg_OligoII$Term)
enrich_gbmlgg_OligoII$Adjusted.P.value=signif(enrich_gbmlgg_OligoII$Adjusted.P.value,2)
enrich_gbmlgg_OligoII$P.value=round(log10(enrich_gbmlgg_OligoII$Adjusted.P.value),1)
enrich_gbmlgg_OligoII$Term=str_wrap(enrich_gbmlgg_OligoII$Term, width = 40)
fig4e=plotEnrich(enrich_gbmlgg_OligoII, showTerms = 14, numChar = 100,xlab="Enriched Terms",ylab="Gene Count",orderBy = "P.value")+theme_classic()+ facet_grid(.~Sign)+scale_fill_continuous(high="lightblue4",low="lightblue1",labels = scales::label_number(accuracy = 1))+ theme(text = element_text(size = 6), plot.title = element_blank(),legend.position = "bottom",legend.justification = "left",legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.4, "lines"),legend.margin=margin(c(0,0,0,0)))+ guides(fill=guide_colourbar(title="Log10 Adjusted\np-value"))

enrich_gbmlgg_OligoIII_up=enrichr(gene=gene_labels[gene_labels$Ensembl_Gene_ID %in% unique(colnames(x_oligoiii[,colSums(x_oligoiii)>=3])) & gene_labels$Ensembl_Gene_ID %in% rownames(res_gene_tcga_gbmlgg_pred_OligoIII[res_gene_tcga_gbmlgg_pred_OligoIII$log2FC>0,]),]$HGNC_Symbol, dbs)
gc()
enrich_gbmlgg_OligoIII_down=enrichr(gene=gene_labels[gene_labels$Ensembl_Gene_ID %in% unique(colnames(x_oligoiii[,colSums(x_oligoiii)>=3])) & gene_labels$Ensembl_Gene_ID %in% rownames(res_gene_tcga_gbmlgg_pred_OligoIII[res_gene_tcga_gbmlgg_pred_OligoIII$log2FC<0,]),]$HGNC_Symbol, dbs)
gc()
enrich_gbmlgg_OligoIII=rbind(data.frame(head(enrich_gbmlgg_OligoIII_up[[1]][order(enrich_gbmlgg_OligoIII_up[[1]]$P.value),],7),Sign="Upregulated"),data.frame(head(enrich_gbmlgg_OligoIII_down[[1]][order(enrich_gbmlgg_OligoIII_down[[1]]$P.value),],7),Sign="Downregulated"))
enrich_gbmlgg_OligoIII$Term=sub('(.*)\\s\\(.*', '\\1',enrich_gbmlgg_OligoIII$Term)
enrich_gbmlgg_OligoIII$Adjusted.P.value=signif(enrich_gbmlgg_OligoIII$Adjusted.P.value,2)
enrich_gbmlgg_OligoIII$P.value=round(log10(enrich_gbmlgg_OligoIII$Adjusted.P.value),1)
enrich_gbmlgg_OligoIII$Term=str_wrap(enrich_gbmlgg_OligoIII$Term, width = 40)
fig4f=plotEnrich(enrich_gbmlgg_OligoIII, showTerms = 14, numChar = 100,xlab="Enriched Terms",ylab="Gene Count",orderBy = "P.value")+theme_classic()+ facet_grid(.~Sign)+scale_fill_continuous(high="blue4",low="blue1",labels = scales::label_number(accuracy = 1))+ theme(text = element_text(size = 6), plot.title = element_blank(),legend.position = "bottom",legend.justification = "left",legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.4, "lines"),legend.margin=margin(c(0,0,0,0)))+ guides(fill=guide_colourbar(title="Log10 Adjusted\np-value"))


enrich_gbmlgg_AstroII=data.frame(head(enrich_gbmlgg_AstroII_up[[2]][order(enrich_gbmlgg_AstroII_up[[2]]$P.value),] %>% .[grepl("Human",.$Term) & grepl("Cortex|Brain",.$Term),],7),Sign="Upregulated")
enrich_gbmlgg_AstroII$Term=sub('(.*)\\s\\(.*', '\\1',enrich_gbmlgg_AstroII$Term)
enrich_gbmlgg_AstroII$Adjusted.P.value=signif(enrich_gbmlgg_AstroII$Adjusted.P.value,2)
enrich_gbmlgg_AstroII$P.value=round(log10(enrich_gbmlgg_AstroII$Adjusted.P.value),1)
enrich_gbmlgg_AstroII$Term=str_wrap(gsub(" Human","",enrich_gbmlgg_AstroII$Term), width = 30)
fig5d=plotEnrich(enrich_gbmlgg_AstroII, showTerms = 7, numChar = 100,xlab="Enriched Terms",ylab="Gene Count",orderBy = "P.value")+theme_classic()+scale_fill_continuous(high="yellow4",low="yellow1",labels = scales::label_number(accuracy = 1))+ theme(text = element_text(size = 6), plot.title = element_blank(),legend.position = "bottom",legend.justification = "left",legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.4, "lines"),legend.margin=margin(c(0,0,0,-25)))+ guides(fill=guide_colourbar(title="Log10\nAdjusted\np-value"))

enrich_gbmlgg_AstroIII=data.frame(head(enrich_gbmlgg_AstroIII_up[[2]][order(enrich_gbmlgg_AstroIII_up[[2]]$P.value),] %>% .[grepl("Human",.$Term) & grepl("Cortex|Brain",.$Term),],7),Sign="Upregulated")
enrich_gbmlgg_AstroIII$Term=sub('(.*)\\s\\(.*', '\\1',enrich_gbmlgg_AstroIII$Term)
enrich_gbmlgg_AstroIII$Adjusted.P.value=signif(enrich_gbmlgg_AstroIII$Adjusted.P.value,2)
enrich_gbmlgg_AstroIII$P.value=round(log10(enrich_gbmlgg_AstroIII$Adjusted.P.value),1)
enrich_gbmlgg_AstroIII$Term=str_wrap(gsub(" Human","",enrich_gbmlgg_AstroIII$Term), width = 30)
fig5e=plotEnrich(enrich_gbmlgg_AstroIII, showTerms = 7, numChar = 100,xlab="Enriched Terms",ylab="Gene Count",orderBy = "P.value")+theme_classic()+scale_fill_continuous(high="orange4",low="orange1",labels = scales::label_number(accuracy = 1))+ theme(text = element_text(size = 6), plot.title = element_blank(),legend.position = "bottom",legend.justification = "left",legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.4, "lines"),legend.margin=margin(c(0,0,0,-25)))+ guides(fill=guide_colourbar(title="Log10\nAdjusted\np-value"))

enrich_gbmlgg_AstroIV=data.frame(head(enrich_gbmlgg_AstroIV_up[[2]][order(enrich_gbmlgg_AstroIV_up[[2]]$P.value),] %>% .[grepl("Human",.$Term) & grepl("Cortex|Brain",.$Term),],7),Sign="Upregulated")
enrich_gbmlgg_AstroIV$Term=sub('(.*)\\s\\(.*', '\\1',enrich_gbmlgg_AstroIV$Term)
enrich_gbmlgg_AstroIV$Adjusted.P.value=signif(enrich_gbmlgg_AstroIV$Adjusted.P.value,2)
enrich_gbmlgg_AstroIV$P.value=round(log10(enrich_gbmlgg_AstroIV$Adjusted.P.value),1)
enrich_gbmlgg_AstroIV$Term=str_wrap(gsub(" Human","",enrich_gbmlgg_AstroIV$Term), width = 30)
fig5f=plotEnrich(enrich_gbmlgg_AstroIV, showTerms = 7, numChar = 100,xlab="Enriched Terms",ylab="Gene Count",orderBy = "P.value")+theme_classic()+scale_fill_continuous(high="firebrick4",low="firebrick1",labels = scales::label_number(accuracy = 1))+ theme(text = element_text(size = 6), plot.title = element_blank(),legend.position = "bottom",legend.justification = "left",legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.4, "lines"),legend.margin=margin(c(0,0,0,-25)))+ guides(fill=guide_colourbar(title="Log10\nAdjusted\np-value"))

enrich_gbmlgg_GBMIV=data.frame(head(enrich_gbmlgg_GBMIV_up[[2]][order(enrich_gbmlgg_GBMIV_up[[2]]$P.value),] %>% .[grepl("Human",.$Term) & grepl("Cortex|Brain",.$Term),],7),Sign="Upregulated")
enrich_gbmlgg_GBMIV$Term=sub('(.*)\\s\\(.*', '\\1',enrich_gbmlgg_GBMIV$Term)
enrich_gbmlgg_GBMIV$Adjusted.P.value=signif(enrich_gbmlgg_GBMIV$Adjusted.P.value,2)
enrich_gbmlgg_GBMIV$P.value=round(log10(enrich_gbmlgg_GBMIV$Adjusted.P.value),1)
enrich_gbmlgg_GBMIV$Term=str_wrap(gsub(" Human","",enrich_gbmlgg_GBMIV$Term), width = 30)
fig5c=plotEnrich(enrich_gbmlgg_GBMIV, showTerms = 7, numChar = 100,xlab="Enriched Terms",ylab="Gene Count",orderBy = "P.value")+theme_classic()+scale_fill_continuous(high="purple4",low="purple1",labels = scales::label_number(accuracy = 1))+ theme(text = element_text(size = 6), plot.title = element_blank(),legend.position = "bottom",legend.justification = "left",legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.4, "lines"),legend.margin=margin(c(0,0,0,-25)))+ guides(fill=guide_colourbar(title="Log10\nAdjusted\np-value"))

enrich_gbmlgg_OligoII=data.frame(head(enrich_gbmlgg_OligoII_up[[2]][order(enrich_gbmlgg_OligoII_up[[2]]$P.value),] %>% .[grepl("Human",.$Term) & grepl("Cortex|Brain",.$Term),],7),Sign="Upregulated")
enrich_gbmlgg_OligoII$Term=sub('(.*)\\s\\(.*', '\\1',enrich_gbmlgg_OligoII$Term)
enrich_gbmlgg_OligoII$Adjusted.P.value=signif(enrich_gbmlgg_OligoII$Adjusted.P.value,2)
enrich_gbmlgg_OligoII$P.value=round(log10(enrich_gbmlgg_OligoII$Adjusted.P.value),1)
enrich_gbmlgg_OligoII$Term=str_wrap(gsub(" Human","",enrich_gbmlgg_OligoII$Term), width = 30)
fig5g=plotEnrich(enrich_gbmlgg_OligoII, showTerms = 7, numChar = 100,xlab="Enriched Terms",ylab="Gene Count",orderBy = "P.value")+theme_classic()+scale_fill_continuous(high="lightblue4",low="lightblue1",labels = scales::label_number(accuracy = 1))+ theme(text = element_text(size = 6), plot.title = element_blank(),legend.position = "bottom",legend.justification = "left",legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.4, "lines"),legend.margin=margin(c(0,0,0,-25)))+ guides(fill=guide_colourbar(title="Log10\nAdjusted\np-value"))

enrich_gbmlgg_OligoIII=data.frame(head(enrich_gbmlgg_OligoIII_up[[2]][order(enrich_gbmlgg_OligoIII_up[[2]]$P.value),] %>% .[grepl("Human",.$Term) & grepl("Cortex|Brain",.$Term),],7),Sign="Upregulated")
enrich_gbmlgg_OligoIII$Term=sub('(.*)\\s\\(.*', '\\1',enrich_gbmlgg_OligoIII$Term)
enrich_gbmlgg_OligoIII$Adjusted.P.value=signif(enrich_gbmlgg_OligoIII$Adjusted.P.value,2)
enrich_gbmlgg_OligoIII$P.value=round(log10(enrich_gbmlgg_OligoIII$Adjusted.P.value),1)
enrich_gbmlgg_OligoIII$Term=str_wrap(gsub(" Human","",enrich_gbmlgg_OligoIII$Term), width = 30)
fig5h=plotEnrich(enrich_gbmlgg_OligoIII, showTerms = 7, numChar = 100,xlab="Enriched Terms",ylab="Gene Count",orderBy = "P.value")+theme_classic()+scale_fill_continuous(high="blue4",low="blue1",labels = scales::label_number(accuracy = 1))+ theme(text = element_text(size = 6), plot.title = element_blank(),legend.position = "bottom",legend.justification = "left",legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.4, "lines"),legend.margin=margin(c(0,0,0,-25)))+ guides(fill=guide_colourbar(title="Log10\nAdjusted\np-value"))

##################################################################
#####Deconvolution (Bulk Glioma)#####
##################################################################
data_tpm_gbmlgg_tcga=data_gene_gbmlgg_tcga$abundance[rownames(data_gene_gbmlgg_tcga$abundance) %in% rownames(cts_gene_gbmlgg_tcga),]
rownames(data_tpm_gbmlgg_tcga)=gene_labels[match(rownames(data_tpm_gbmlgg_tcga), gene_labels$Ensembl_Gene_ID),]$HGNC_Symbol
data_tpm_gbmlgg_tcga=data_tpm_gbmlgg_tcga[rownames(data_tpm_gbmlgg_tcga)!="N/A" & !duplicated(rownames(data_tpm_gbmlgg_tcga)),]
corrected_data_tpm_gbmlgg_tcga <- 2^(removeBatchEffect(log2(data_tpm_gbmlgg_tcga+1),batch=factors_gbmlgg_tcga$Read_Length))-1

estimate_deconv=data.frame(t(data.frame(deconvolute(corrected_data_tpm_gbmlgg_tcga, "estimate"))))
estimate_deconv=estimate_deconv[-1,]
colnames(estimate_deconv)=c("ESTIMATE_Stromal_Score","ESTIMATE_Immune_Score","ESTIMATE_Estimate_Score","ESTIMATE_Tumor_Purity")
estimate_deconv[,c(1:4)]=signif(data.frame(sapply(data.frame(estimate_deconv[,c(1:4)]), as.numeric)),3)

factors_gbmlgg_tcga2=cbind(factors_gbmlgg_tcga,estimate_deconv)
factors_gbmlgg_tcga2$Predicted_Labels=sapply(factors_gbmlgg_tcga$Predicted_GLM_ElasticNet_WHO_CNS5_2021_Diagnosis, switch,"Astrocytoma_G2"="A (G2/3)","Astrocytoma_G3"="A (G2/3)","Astrocytoma_G4"="A (G4)","Glioblastoma_G4"="GBM (G4)","Oligodendroglioma_G2"="O (G2/3)","Oligodendroglioma_G3"="O (G2/3)")

write.csv(factors_gbmlgg_tcga2[c(1,32:33,35)],"Pub_Classify/TCGA_Glioma_ESTIMATE_Deconvolution.csv",row.names = F)

fig5b=ggplot(factors_gbmlgg_tcga2[factors_gbmlgg_tcga2$Dropped_Sample==F,], aes(x=Predicted_Labels, y=ESTIMATE_Tumor_Purity, fill=Predicted_Labels)) + geom_boxplot(outlier.shape=NA) + stat_compare_means(comparisons = list( c("A (G2/3)", "A (G4)"), c("A (G2/3)", "GBM (G4)"), c("A (G2/3)", "O (G2/3)"), c("A (G4)", "GBM (G4)"), c("A (G4)", "O (G2/3)"), c("GBM (G4)", "O (G2/3)") ), label = "p.signif",  method = "wilcox.test", p.adjust.method = "FDR",size = 1.5, vjust = 0.2, step.increase = 0.07, label.y=1)+theme_classic()+scale_fill_manual(values=c(`GBM (G4)`="purple2",`A (G2/3)`="chocolate2",`A (G4)`="firebrick2",`O (G2/3)`="dodgerblue2"))+guides(fill=guide_legend(title="Predicted\nWHO CNS5\nAnnotation", nrow = 2, byrow = TRUE))+xlab("Glioma Type")+ylab("Tumor Purity")+theme(text = element_text(size = 6),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.75, "lines"),legend.position = "top",legend.justification = "left")+ylim(0.65,1.25)

figs5g=ggplot(factors_gbmlgg_tcga2[factors_gbmlgg_tcga2$Dropped_Sample==F,], aes(x=Predicted_Labels, y=ESTIMATE_Stromal_Score, fill=Predicted_Labels)) + geom_boxplot(outlier.shape=NA) + stat_compare_means(comparisons = list( c("A (G2/3)", "A (G4)"), c("A (G2/3)", "GBM (G4)"), c("A (G2/3)", "O (G2/3)"), c("A (G4)", "GBM (G4)"), c("A (G4)", "O (G2/3)"), c("GBM (G4)", "O (G2/3)") ), label = "p.signif",  method = "wilcox.test", p.adjust.method = "FDR",size = 1.5, vjust = 0.2, step.increase = 0.07, label.y=750,hide.ns = F)+theme_classic()+scale_fill_manual(values=c(`GBM (G4)`="purple2",`A (G2/3)`="chocolate2",`A (G4)`="firebrick2",`O (G2/3)`="dodgerblue2"))+guides(fill=guide_legend(title="Predicted\nWHO CNS5\nAnnotation", nrow = 4, byrow = TRUE))+xlab("Glioma Type")+ylab("Stromal Score (ESTIMATE)")+theme(text = element_text(size = 6),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.75, "lines"))+ylim(-2000,2500)
figs5h=ggplot(factors_gbmlgg_tcga2[factors_gbmlgg_tcga2$Dropped_Sample==F,], aes(x=Predicted_Labels, y=ESTIMATE_Immune_Score, fill=Predicted_Labels)) + geom_boxplot(outlier.shape=NA) + stat_compare_means(comparisons = list( c("A (G2/3)", "A (G4)"), c("A (G2/3)", "GBM (G4)"), c("A (G2/3)", "O (G2/3)"), c("A (G4)", "GBM (G4)"), c("A (G4)", "O (G2/3)"), c("GBM (G4)", "O (G2/3)") ), label = "p.signif",  method = "wilcox.test", p.adjust.method = "FDR",size = 1.5, vjust = 0.2, step.increase = 0.07, label.y=1750)+theme_classic()+scale_fill_manual(values=c(`GBM (G4)`="purple2",`A (G2/3)`="chocolate2",`A (G4)`="firebrick2",`O (G2/3)`="dodgerblue2"))+guides(fill=guide_legend(title="Predicted\nWHO CNS5\nAnnotation", nrow = 4, byrow = TRUE))+xlab("Glioma Type")+ylab("Immune Score (ESTIMATE)")+theme(text = element_text(size = 6),axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),legend.key.width  = unit(0.75, "lines"),legend.key.height = unit(0.75, "lines"))+ylim(-1750,3500)

##################################################################
#####Final Figure Order#####
##################################################################
plot1_1=plot_grid(fig1a, fig1b, fig1c, labels = c('A', 'B', 'C'), label_size = 10, rel_widths = c(1, 1, 1),ncol = 1)

pdf("Pub_Classify/Fig1.pdf",height=4,width=6.5)
plot_grid("",plot1_1, ncol = 1,rel_heights=c(0.05,1))
dev.off()

plots1_1=plot_grid(figs1a, figs1b, figs1c, labels = c('A', 'B', 'C'), label_size = 10, rel_widths = c(1, 1, 1),ncol = 1)

pdf("Pub_Classify/SupFig1.pdf",height=4,width=6.5)
plot_grid("",plots1_1, ncol = 1,rel_heights=c(0.05,1))
dev.off()


plot2_1=plot_grid(fig2a, fig2b, fig2c, labels = c('A', 'B', 'C'), label_size = 10, rel_widths = c(0.75, 1, 1),ncol = 3)
plot2_2=plot_grid(fig2d, fig2e, fig2f, labels = c('D', 'E', 'F'), label_size = 10, rel_widths = c(0.75, 1, 1),ncol = 3)
plot2_3=plot_grid(fig2g, fig2h, fig2i, labels = c('G', 'H', 'I'), label_size = 10, rel_widths = c(0.75, 1, 1),ncol = 3)

pdf("Pub_Classify/Fig2.pdf",height=8.5,width=6.5)
plot_grid("",plot2_1,plot2_2,plot2_3, ncol = 1,rel_heights=c(0.05,1,1,1))
dev.off()

plots2_1=plot_grid(figs2a, figs2b, figs2c, labels = c('A', 'B', 'C'), label_size = 10, rel_widths = c(1, 1, 1),ncol = 3)

pdf("Pub_Classify/SupFig2.pdf",height=2.5,width=6.5)
plot_grid("",plots2_1, ncol = 1,rel_heights=c(0.05,1))
dev.off()


plot3_1=plot_grid(fig3a, fig3b, labels = c('A', 'B'), label_size = 10, rel_widths = c(1, 1),ncol = 2)

pdf("Pub_Classify/Fig3.pdf",height=3.5,width=6.5)
plot_grid("",plot3_1, ncol = 1,rel_heights=c(0.05,1))
dev.off()

plots3_1=plot_grid(figs3a, figs3b, labels = c('A', 'B'), label_size = 10, rel_widths = c(1.05, 1),ncol = 2)
plots3_2=plot_grid(figs3c, figs3d, labels = c('C', 'D'), label_size = 10, rel_widths = c(1.05, 1),ncol = 2)
plots3_3=plot_grid(figs3e, figs3f, labels = c('E', 'F'), label_size = 10, rel_widths = c(1.05, 1),ncol = 2)

pdf("Pub_Classify/SupFig3.pdf",height=7,width=6.5)
plot_grid("",plots3_1,plots3_2, plots3_3, ncol = 1,rel_heights=c(0.05,1,1,1))
dev.off()


plot4_1=plot_grid(fig4a, fig4b, labels = c('A', 'B'), label_size = 10, rel_widths = c(1, 1),ncol = 2)
plot4_2=plot_grid(fig4c, fig4d, labels = c('C', 'D'), label_size = 10, rel_widths = c(1, 1),ncol = 2)
plot4_3=plot_grid(fig4e, fig4f, labels = c('E', 'F'), label_size = 10, rel_widths = c(1, 1),ncol = 2)

pdf("Pub_Classify/Fig4.pdf",height=9,width=6.5)
plot_grid("",plot4_1,plot4_2,plot4_3, ncol = 1,rel_heights=c(0.05,1,1,1))
dev.off()

plots4_1=plot_grid(figs4a, figs4b, labels = c('A', 'B'), label_size = 10, rel_widths = c(1, 1),ncol = 2)
plots4_2=plot_grid(figs4c, figs4d, labels = c('C', 'D'), label_size = 10, rel_widths = c(1, 1),ncol = 2)
plots4_3=plot_grid(figs4e, figs4f, labels = c('E', 'F'), label_size = 10, rel_widths = c(1, 1),ncol = 2)

pdf("Pub_Classify/SupFig4.pdf",height=6,width=6.5)
plot_grid("",plots4_1,plots4_2,plots4_3, ncol = 1,rel_heights=c(0.05,1,1,1))
dev.off()


plot5_1=plot_grid(fig5a, fig5b, labels = c('A', 'B'), label_size = 10, rel_widths = c(1, 0.5),ncol = 2)
plot5_2=plot_grid(fig5c, fig5d, fig5e, labels = c('C' ,'D', 'E'), label_size = 10, rel_widths = c(1, 1, 1),ncol = 3)
plot5_3=plot_grid(fig5f, fig5g, fig5h, labels = c('F' ,'G', 'H'), label_size = 10, rel_widths = c(1, 1, 1),ncol = 3)

pdf("Pub_Classify/Fig5.pdf",height=8,width=6.5)
plot_grid("",plot5_1,plot5_2,plot5_3, ncol = 1,rel_heights=c(0.05,1.5,1,1))
dev.off()

plots5_1=plot_grid(figs5a, figs5b, figs5c, labels = c('A', 'B', 'C'), label_size = 10, rel_widths = c(1, 1, 1),ncol = 3)
plots5_2=plot_grid(figs5d, figs5e, figs5f, labels = c('D', 'E', 'F'), label_size = 10, rel_widths = c(1, 1, 1),ncol = 3)
plots5_3=plot_grid(figs5g, figs5h, labels = c('G', 'H'), label_size = 10, rel_widths = c(1, 1),ncol = 2)

pdf("Pub_Classify/SupFig5.pdf",height=8,width=6.5)
plot_grid("",plots5_1,plots5_2, plots5_3, ncol = 1,rel_heights=c(0.05,1,1,1))
dev.off()

tab_1=read.csv("Pub_Classify/Glioma_Label_Table.csv",check.names=F)

wb <- createWorkbook()
addWorksheet(wb, "Table 1")

writeData(wb, "Table 1", tab_1, startRow = 1, startCol = 1)

saveWorkbook(wb, file = "Pub_Classify/Tables.xlsx", overwrite = TRUE)


tabs_1=read.csv("Pub_Classify/TCGA_Glioma_Clinical_Info.csv")
tabs_1[is.na(tabs_1)]="N/A"
tabs_2=read.csv("Pub_Classify/CGGA_Glioma_Clinical_Info.csv")
tabs_2[is.na(tabs_2)]="N/A"
tabs_3=read.csv("Pub_Classify/GLASS_Glioma_Clinical_Info.csv")
tabs_3[is.na(tabs_3)]="N/A"
tabs_4=read.csv("Pub_Classify/Glioma_Type_Coefficients.csv")
tabs_5=read.csv("Pub_Classify/Predicted_Glioma_Types.csv")
tabs_6=read.csv("Pub_Classify/TCGA_Glioma_ESTIMATE_Deconvolution.csv")

wb <- createWorkbook()
addWorksheet(wb, "Supplemental Table 1")
addWorksheet(wb, "Supplemental Table 2")
addWorksheet(wb, "Supplemental Table 3")
addWorksheet(wb, "Supplemental Table 4")
addWorksheet(wb, "Supplemental Table 5")
addWorksheet(wb, "Supplemental Table 6")

writeData(wb, "Supplemental Table 1", tabs_1, startRow = 1, startCol = 1)
writeData(wb, "Supplemental Table 2", tabs_2, startRow = 1, startCol = 1)
writeData(wb, "Supplemental Table 3", tabs_3, startRow = 1, startCol = 1)
writeData(wb, "Supplemental Table 4", tabs_4, startRow = 1, startCol = 1)
writeData(wb, "Supplemental Table 5", tabs_5, startRow = 1, startCol = 1)
writeData(wb, "Supplemental Table 6", tabs_6, startRow = 1, startCol = 1)

saveWorkbook(wb, file = "Pub_Classify/Supplemental_Tables.xlsx", overwrite = TRUE)

