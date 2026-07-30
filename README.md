# Manuscript Title: [Gene Expression-Based WHO CNS5 Classification of Gliomas Using Machine Learning Approaches](TBD)
## Abstract
The 2021 WHO Classification of Tumors of the Central Nervous System (WHO CNS5) integrates molecular alterations with histopathology to refine glioma classification. However, many widely used public transcriptomic datasets were generated and annotated under earlier frameworks, limiting cross-cohort integration and retrospective analyses using current diagnostic criteria. We developed a transcriptome-based machine learning framework to assign adult-type diffuse gliomas to WHO CNS5-consistent types and grades using bulk RNA sequencing (RNA-seq) data alone. Models were trained using WHO CNS5-reannotated TCGA gliomas and externally evaluated in the CGGA and GLASS cohorts, excluding rare glioma entities. Among five evaluated algorithms, an elastic net generalized linear model demonstrated robust cross-cohort performance and interpretability. Predicted annotations were assessed using survival analyses, differential gene expression, pathway enrichment, and tumor microenvironment deconvolution. Notably, the classifier enabled systematic re-annotation of legacy glioma samples, including cases with incomplete molecular metadata. Predicted glioma classes exhibited coherent transcriptomic organization and survival stratification across datasets. Comparative analyses identified relative enrichment of immune-associated and extracellular matrix pathways in glioblastoma, proliferative programs in astrocytoma, grade 4, and lineage-restricted oligodendrocyte-associated programs in oligodendrogliomas. Furthermore, re-annotation systematically uncovered previously underappreciated divergence in molecular programs that distinguish astrocytoma, grades 2-4 and oligodendroglioma, grades 2-3. This study provides a scalable transcriptome-based approach for harmonizing adult-type diffuse glioma annotations across historical and contemporary RNA-seq datasets. By enabling WHO CNS5-consistent re-annotation, biologically coherent stratification, and molecular and cellular characterization of glioma types and grades, this framework facilitates cross-cohort integration and retrospective discovery in neuro-oncology research.
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

