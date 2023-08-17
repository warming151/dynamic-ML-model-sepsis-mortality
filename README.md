# dynamic-ML-model-sepsis-mortality
Coding for paper "A machine learning model derived from analysis of time-course gene-expression datasets reveal temporally stable gene markers predictive of sepsis mortality"

## Running step
Download the public dataset from GEO Database.

"Limma" R package for obtaining bulk DEGs.

"masig.R" for obtaining temperal DEGs.

Once get the DEGs, create the matrix (n\*m) (n is all the samples and m is the DEGs we selected)

Inputting the matrix and labels into "dynamic_svm.ipynb" can obtain the prediction results.

