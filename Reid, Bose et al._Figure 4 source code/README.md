# Predicting Response to Treatment from Pre-treatment Metabolic Profile of Bone Marrow Samples

## Synopsis

The file [Reid, Bose et al_Figure 4 data set_Bone_Marrow.csv](Reid,%20Bose%20et%20al_Figure%204%20data%20set_Bone_Marrow.csv) contains the metabolic profile of bone marrow samples taken from patients before treatment. Columns represent patients and indicate whether or not they responded to treatment.

The R script [Reid, Bose et al_Figure 4CEF.R](Reid,%20Bose%20et%20al_Figure%204CEF.R) computes the Mutual Information Coefficient (MIC) for each metabolite profile across patients and the treatment response vector. It then takes a set of metabolites to train and test a logistic regression model predicting response to treatment. 


## Requirements
The `minerva` package is used in the calculations of the MIC. The modeling step makes use of the following packages: `glmnet, sets, pROC, fscaret, fbroc`.

## Usages

The code can be run as is or in individual chunks, and will save the results into files named `MIC_BaselineMarrowRevised_Tx.csv` and `ROC_Curves.csv`.

## Contributors

**Peter G Mikhael**

+ [http://github.com/pgmikhael](http://github.com/pgmikhael)

## License

This software is released under the [MIT License](LICENSE-MIT).
