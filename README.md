# Bivariate-pseudo-observations

This repository contains the code for the simulations and data analysis from the paper "Pseudo-Observations for Bivariate Survival Data".

## Description
The file "psuedo_observations_functions.R" contains the functions that compute the pseudo-observations based on either the Dabrowska estimator or the Lin and Ying estimator of the bivariate survival function.
Additionally, this repository contains two folders, one folder for the data analysis example, and another folder for the simulations. 
### data analysis
This folder contains four files: (1) the diabetic retinopathy data as a CSV file, (2) the code for calculating the covariate-adjusted 5-year joint survival probability based on the proposed approach (Section 4 of our paper), and (3)-(4) a tutorial explaining how to use the approach to estimate covariate-adjusted conditional probabilities (RMD file + HTML file containing interactive 3D figures).
### simulations
This folder contains two files: (1) functions to generate the simulated data (bivariate logistic and bivariate lognormal), and (2) an example of running the simulations for the bivariate logistic data and the fixed time point selection approach (Section 3.1 of our paper). Additional simulation settings can be derived similarly.

## Example - diabetic retinopathy study
1. Load the data and the pseudo observations functions:
```{r}
library(geepack)
obs <- read.csv("retinopathy_data.csv")
source("pseudo_observations_functions.R")
```

2. Compute the pseudo-observations for the 5-year (60-month) joint survival probability based on the Dabrowska estimator:
```{r}
t0 <- c(60,60)
obs <- PO_func_dabrowska(obs,t0)
```

3. Fit a regression model using a GEE with a logit link:
```{r}
fit <- geese(PO~age+mean_risk+type,
             id=id, data=obs,scale.fix=TRUE,family=gaussian,
             jack=TRUE, mean.link="logit",corstr="independence")
```
