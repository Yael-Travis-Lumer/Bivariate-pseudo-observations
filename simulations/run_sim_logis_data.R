library(MASS)
library(mvtnorm)
library(survival)
library(tidyr)
library(zoo)

source("generate data.R")
source("pseudo_observations_functions.R")

###-----Set parameter values and sample the univariate covariate-----
a1 <- 1
a2 <- 3
lambda <- 1
beta <- 2
n <- 200
set.seed(123)
Z <- runif(n,0.5,1.5) 

###----Select time points and calculate the corresponding true conditional bivariate survival at these time points------
a <- c(0.5,1)
b <- c(0.7,1.2,1.5)
t0 <- expand.grid(a,b)
t0_merged <- paste0(t0[,1],",",t0[,2], sep = "")
S_mat <- matrix(0,nrow = n,ncol=nrow(t0))
c_vec <- exp(beta*Z)
for (k in 1:nrow(t0)){
  S_mat[,k] <- S(t1 = t0[k,1],t2=t0[k,2],a1 = a1,a2 = a2, lambda = lambda, c = c_vec)
}
S_DF <- as.data.frame(S_mat)
colnames(S_DF) <- paste0("S_",1:nrow(t0),"")
summary(S_DF)


###----Simulations---- 
m <- 500
betas_df_LY <- data.frame(matrix(NA,m,1+nrow(t0)))
colnames(betas_df_LY) <- c(paste0("beta0_",1:nrow(t0),""),"Z")
betas_df_Dab <- betas_df_LY
var_betas_df_Dab <- var_betas_df_LY <- list(m)
est_S_logit_list_LY <- list(m)
est_S_logit_list_Dab <- est_S_logit_list_LY
est_S_logit_LY <- matrix(0,n,nrow(t0))
colnames(est_S_logit_LY) <- paste0("S_hat_",1:nrow(t0),"")
est_S_logit_Dab <- est_S_logit_LY

for (i in 1:m){
  #generate data
  obs <- gen_logit_data(a1, a2, lambda ,beta ,Z, seed=i)
  
  #calculate the POs for the two estimators
  obs_LY <- joint_surv_PO_mult(obs,t0)
  obs_Dab <- PO_func_dabrowska_mult_new(obs,t0)
  
  #reshape POs
  obs_long_LY <- pivot_longer(obs_LY,cols = starts_with("PO_t"))
  obs_long_Dab <- pivot_longer(obs_Dab,cols = starts_with("PO_t"))
  obs_long_LY$times <- rep(t0_merged,n)
  obs_long_Dab$times <- rep(t0_merged,n)
  
  #fit two regression models (one per estimator) using a GEE with a logit link
  set.seed(m+i)
  fit_LY <- geese(value~factor(times, levels=t0_merged)+Z,
                    id=id, data=obs_long_LY,scale.fix=TRUE,family=gaussian,
                    jack=TRUE, mean.link="logit",corstr="independence",)
  betas_LY <- fit_LY$beta
  var_betas_LY <- fit_LY$vbeta
  
    
  set.seed(2*m+i)
  fit_Dab <- geese(value~factor(times, levels=t0_merged)+Z,
                     id=id, data=obs_long_Dab,scale.fix=TRUE,family=gaussian,
                     jack=TRUE, mean.link="logit",corstr="independence",)
  betas_Dab <- fit_Dab$beta
  var_betas_Dab <- fit_Dab$vbeta
  
  #store results  
  betas_df_LY[i,] <- betas_LY
  betas_df_Dab[i,] <- betas_Dab
  var_betas_df_LY[[i]] <- data.frame(var_betas_LY)
  var_betas_df_Dab[[i]] <- data.frame(var_betas_Dab)
  
  
  #calculate the estimated conditional joint survival for each time point and for each type of estimator
    for (j in 1:nrow(t0)){
      if (j>1){
        est_S_logit_LY[,j] <- exp(betas_LY[1]+betas_LY[j]+betas_LY[nrow(t0)+1]*Z)/(1+exp(betas_LY[1]+betas_LY[j]+betas_LY[nrow(t0)+1]*Z))
        est_S_logit_Dab[,j] <- exp(betas_Dab[1]+betas_Dab[j]+betas_Dab[nrow(t0)+1]*Z)/(1+exp(betas_Dab[1]+betas_Dab[j]+betas_Dab[nrow(t0)+1]*Z))
      } else {
        est_S_logit_LY[,j] <- exp(betas_LY[1]+betas_LY[nrow(t0)+1]*Z)/(1+exp(betas_LY[1]+betas_LY[nrow(t0)+1]*Z))
        est_S_logit_Dab[,j] <- exp(betas_Dab[1]+betas_Dab[nrow(t0)+1]*Z)/(1+exp(betas_Dab[1]+betas_Dab[nrow(t0)+1]*Z))
      }
    }
    est_S_logit_list_LY[[i]] <- est_S_logit_LY
    est_S_logit_list_Dab[[i]] <- est_S_logit_Dab
  
}

#save results
save(betas_df_LY, var_betas_df_LY, est_S_logit_list_LY, file="logistic_model_grid_points_LY_500.Rdata")
save(betas_df_Dab, var_betas_df_Dab ,est_S_logit_list_Dab, file="logistic_model_grid_points_Dab_500.Rdata")
