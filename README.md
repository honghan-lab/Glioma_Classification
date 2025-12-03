# Manuscript Title: [Gene Expression-Based WHO CNS5 Classification of Gliomas Using Machine Learning Approaches](TBD)
## Contents
### Code for generating data & associated figures (Glioma_Tumor_Classify_code.R; R script)
- Loading libraries
- Preparing covariate data
- Preparing gene counts
- Machine learning classifier training and validation for glioma type prediction 
- Generating UMAPs
- Differential expression of genes (DEG) analysis
- Pathway enrichment analysis
- Bulk deconvolution analysis
- Figure generation
### Trained models for glioma type prediction (RDS file format)
- Elastic net penalty GLM (Glioma_GLM_ElasticNet_Classifier.rds)
- Ridge penalty GLM (Glioma_GLM_Ridge_Classifier.rds)
- Lasso penalty GLM (Glioma_GLM_Lasso_Classifier.rds)
- Random forest (Glioma_RandomForest_Classifier.rds)
- XGBoost (Glioma_XGBoost_Classifier.rds)

