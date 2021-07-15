# minerva-adjuvant-nsclc
Reference scripts for study "Genomic Signatures Define Three Subtypes of EGFR-Mutant Stage II-III NSCLC With Distinct Adjuvant Therapy Outcomes"

## run 
All analysis are running under R version 3.5.0

0.load_required_package.R -- To install packages used in this analysis.

1.minerva_score_establishment.R -- To perform gene-by-treatment interaction test, select genetic features and assign MINERVA score for each patient.

2.cross_validation.R -- To construct and evaluate MINERVA score by cross validation.

3.leaveoneout.R -- To construct and evaluate MINERVA score by leave one out cross validation.


## data
multigene-171_combined.csv -- data used for survival analysis and MINERVA construction.
