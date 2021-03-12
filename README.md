# PPSV
This repository consists of source codes and data sets of the Protein-Protein Similarity Vector (PPSV) to predict drug-disease associations. The folder of ‘Code’ consists of the code for measuring PPSV scores and the code for predicting the score of drug-disease associations. Additionally, the folder of ‘Datasets’ consists of the original data sets such as drug-target, target-disease, and drug-disease.

# Requirements
The code is written in R (R-4.0.3) and Python (version 3.8) languages. The requirement packages compose of ‘readxl’, ‘readr’, ‘igraph’, ‘protr’, ‘Biostrings’, ‘testit’ in R language and ‘os’, ‘numpy, ‘pandas’, including ‘sklearn’ for Python language.

# Code
This folder consists of the code for measuring PPSV scores and the code for predicting the score of drug-disease associations.
## PPSV_scores.R
This code is to account the Protein-Protein Similarity Vector (PPSV) which contains seven scores of biological meanings between proteins (Uniprot_id) as follows:
1.	The neighboring similarity score
2.	The closeness score
3.	Local alignment score
4.	Global alignment score
5.	Pathway score
6.	Functional score
7.	Druggable property score
## PPSV_model.py
This code is for predicting drug-diseases associations based on separate drug-disease matrix.

# Datasets
This folder consists of the original data sets such as drug-target, target-disease, and drug-disease that used in our experiments. It got 14,264 known drug-disease associations with 1,317 approved drugs and 478 diseases. Moreover, we provided the example of drug-disease matrix for users who would like to replicate the method. The descriptions of all data are described as follows.
1.	‘drug_disease.txt’ : the drug-disease associations with 1,317 approved drugs (“DrugBankIDs”) and 478 diseases ("DiseaseID” in MESH term).
2.	‘drug_target.txt’ : the drug-target associations with 1,317 approved drugs.
3.	‘target_disease.txt’: the target-disease associations with 478 diseases.
4.	‘sequence_unp.txt’ : the information of proteins ‘ sequences.
5.	‘ppi_data.rar’ : the protein-protein interaction data.
6.	‘pathway_data.txt’ : the information of metabolic pathways.
7.	‘Drug-disease matrix.txt’ : the example of drug-disease matrix for carcinoma disease.
