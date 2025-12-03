# Manuscript Title: [Gene Expression-Based WHO CNS5 Classification of Gliomas Using Machine Learning Approaches](TBD)
## Abstract
The 2021 WHO Classification of Tumors of the Central Nervous System (WHO CNS5) integrates molecular and histopathological features to refine glioma diagnosis and prognostication. However, most public transcriptomic datasets were annotated under earlier frameworks, impeding harmonized analyses and clinical translation. We developed a machine learning-based classifier using gene expression data to predict glioma types consistent with the WHO CNS5 criteria. Models were trained on the re-annotated TCGA cohort and validated in the independent CGGA and GLASS datasets. Among five algorithms tested, the elastic net generalized linear model (GLM) achieved optimal cross-cohort generalizability. Glioma type re-annotations were validated through survival, differential gene expression (DGE), pathway, cell deconvolution, and enrichment analyses. The classifier accurately predicted glioma types and grades with high sensitivity and specificity across all three cohorts, successfully re-annotating legacy samples lacking complete molecular data. Reclassified tumors showed consistent clustering and survival differences aligned with clinical expectations and WHO CNS5 categories. Transcriptomic profiling revealed type- and grade-specific enrichment of inflammatory, proliferative, and brain lineage programs, with GBM and lower-grade astrocytomas enriched for microglia- and macrophage-associated signatures, astrocytoma, grade 4 exhibiting proliferative neural progenitor-like signatures, and oligodendrogliomas showing strong oligodendrocyte lineage differentiation. This transcriptome-based framework provides a robust and scalable method for harmonizing glioma classifications across historical and contemporary datasets. By linking molecular diagnostics with biologically interpretable expression patterns, it supports retrospective analyses, multi-cohort integration, and the clinical translation of transcriptomic classifiers into modern neuro-oncologic practice. 
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

