# Bivariate-pseudo-observations

This repository contains the code for the simulations and data analysis from the paper "Pseudo-Observations for Bivariate Survival Data".

## Description
The file "psuedo_observations_functions.R" contains the functions that compute the pseudo-observations based on either the Dabrowska estimator or the Lin and Ying estimator of the bivariate survival function.
Additionally, this repository contains two folders, one folder for the data analysis example, and another folder for the simulations. 
### data analysis
This folder contains two files: (1) the diabetic retinopathy data as a CSV file, and (2) the code for calculating the covariate-adjusted 5-year joint survival probability based on the proposed approach (Section 4 of our paper).
### simulations
This folder contains two files: (1) functions to generate the simulated data (bivariate logistic and bivariate lognormal), and (2) an example of running the simulations for the bivariate logistic data and the fixed time point selection approach (Section 3.1 of our paper).

