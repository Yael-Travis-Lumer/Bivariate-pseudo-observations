library(pseudo)
library(mhazard)
library(Rfast)
library(survival)

###----POs for the Dabrowska estimator at a single time point----------
# obs: a data frame containing the observed times (obs1,obs2) and their corresponding indicators (delta1,delta2)
# t0: a single bivariate time point
PO_func_dabrowska <- function(obs,t0){
  n <- nrow(obs)
  #compute the estimate for the bi-variate survival
  S_dabrowska <- npSurv2(Y1 = obs$obs1, Delta1 = obs$delta1,Y2 = obs$obs2, Delta2 = obs$delta2,newT1 = t0[1],newT2 = t0[2])
  S_hat <- S_dabrowska$Fhat_est
  #compute leave-one-out estimate - using for loop 
  for (i in 1:n){
    obs_i <- obs[-i,]
    S_dabrowska_i <- npSurv2(Y1 = obs_i$obs1, Delta1 = obs_i$delta1,Y2 = obs_i$obs2, Delta2 = obs_i$delta2,newT1 = t0[1],newT2 = t0[2])
    obs[i,"S_hat_i"] <- S_dabrowska_i$Fhat_est
  }
  #compute POs
  obs$PO <- n*S_hat-(n-1)*obs$S_hat_i
  obs$id <- 1:n
  return(obs)
}


###----POs for the Dabrowska estimator at multiple time points----------
# obs: a data frame containing the observed times (obs1,obs2) and their corresponding indicators (delta1,delta2)
# t0: a data frame of size (K*2) with K bivariate time point. The two columns correspond to the two coordinates (t1,t2), and each of the K rows corresponds to a different time point.
PO_func_dabrowska_mult_new <- function(obs,t0){
  n <- nrow(obs)
  old_name <- colnames(obs)
  m <- ncol(obs)
  l <- nrow(t0)
  obs_temp <- obs
  obs_temp[,(m+1):(m+l)] <- NA
  S_dabrowska <- npSurv2(Y1 = obs$obs1, Delta1 = obs$delta1,Y2 = obs$obs2, Delta2 = obs$delta2,newT1 = t0[,1],newT2 = t0[,2])
  S_hat <- diag(S_dabrowska$Fhat_est)
  S_hat_mat  <- do.call("rbind", replicate(n,S_hat,simplify = FALSE))
  #compute leave-one-out estimate - using for loop 
  for (i in 1:n){
    obs_i <- obs[-i,]
    S_dabrowska_i <- npSurv2(Y1 = obs_i$obs1, Delta1 = obs_i$delta1,Y2 = obs_i$obs2, Delta2 = obs_i$delta2,newT1 = t0[,1],newT2 = t0[,2])
    obs_temp[i,(m+1):(m+l)] <- as.list(diag(S_dabrowska_i$Fhat_est))
  }
  #compute POs
  obs[,(m+1):(m+l)] <- n*S_hat_mat-(n-1)*obs_temp[,(m+1):(m+l)]
  colnames(obs) <- c(old_name,paste0("PO_t",1:l,""))
  obs$id <- 1:n
  return(obs)
}

###-----Helper functions from the pseudo R package----
#These functions are used for computing the LOO Kaplan-Meier estimate used in the denominator of the Lin and Ying estimator
"surv.omit" <-
  function(pseudo,tmax){
    
    # calculate Kaplan - Meier, leave one out 
    howmany <- nrow(pseudo)
    
    td <- pseudo$time[pseudo$event==1]
    lt.temp <- c(td[-1],td[length(td)]+1)
    lt <- which(td!=lt.temp)
    
    #km - i
    Y1 <- matrix(howmany:1,byrow=TRUE,ncol=howmany,nrow=howmany)
    Y2 <- matrix((howmany-1):0,byrow=TRUE,ncol=howmany,nrow=howmany)
    Y <- upper.tri(Y1,diag=FALSE)*Y1+lower.tri(Y2,diag=TRUE)*Y2
    N <- matrix(pseudo$event,byrow=TRUE,ncol=howmany,nrow=howmany)
    Ndiag <- diag(diag(N))
    N <- N - Ndiag
    
    kmji <- (Y-N)/Y
    
    km <- t(apply(kmji,1,cumprod))
    
    if(!missing(tmax)){
      tt <- matrix(pseudo$time,byrow=TRUE,nrow=nrow(pseudo),ncol=nrow(pseudo))
      #diag(tt) <- c(diag(tt[-nrow(pseudo),-1]),tmax)
      diag(tt) <- c(0,diag(tt[-1,-nrow(pseudo)]))
      tt <- tt[,pseudo$event==1,drop=FALSE]
      tt <- tt[,lt,drop=FALSE]
      tt <- cbind(rep(0,nrow(pseudo)),tt,rep(tmax,nrow(pseudo)))
      tt <- t(apply(tt,1,diff))
    }
    
    
    #corrected value for the last time - last value carried forward 
    aje <- which(is.na(km[howmany,]))
    if(length(aje)>0){
      kir <- min(aje)
      km[howmany,kir:ncol(km)] <- km[howmany,kir-1] 
    }
    
    #only for deaths, one value per tie
    km <- km[,pseudo$event==1,drop=FALSE]
    km <- km[,lt,drop=FALSE]
    if(!missing(tmax)){
      km <- apply(cbind(rep(1,nrow(pseudo)),km)*tt,1,sum)
    }
    km	
  }

"surv.tot" <-
  function(pseudo,tmax){
    
    # calculate Kaplan - Meier, all cases
    
    howmany <- nrow(pseudo)
    
    td <- pseudo$time[pseudo$event==1]
    lt.temp <- c(td[-1],td[length(td)]+1)
    lt <- which(td!=lt.temp)
    
    #km - i
    Y <- howmany:1
    N <- pseudo$event
    
    kmji <- (Y-N)/Y
    
    km <- cumprod(kmji)
    
    if(!missing(tmax)){
      tt <- pseudo$time[pseudo$event==1]
      tt <- tt[lt]
      tt <- c(0,tt,tmax)
      tt <- diff(tt)
    }
    
    #only for deaths, one value per tie
    km <- km[pseudo$event==1]
    km <- km[lt]
    if(!missing(tmax)){
      km <- sum(c(1,km)*tt)
    }
    km	
  }


###--------POs for the Lin and Ying estimator----------
# obs: a data frame containing the observed times (obs1,obs2) and their corresponding indicators (delta1,delta2)
# t0: a data frame of size (K*2) with K bivariate time point. The two columns correspond to the two coordinates (t1,t2), and each of the K rows corresponds to a different time point.
joint_surv_PO_mult <- function(obs, t0){
  old_name <- colnames(obs)
  m <- ncol(obs)
  time <- pmax(obs$obs1,obs$obs2)
  event <- 1-obs$delta1*obs$delta2
  tmax <- pmax(t0[,1],t0[,2])
  tmax <- sort(tmax)
  ltmax <- length(tmax)
  howmany <- length(time)
  pseudo <- data.frame(id = 1:howmany, time = time, event = event)
  pseudo <- pseudo[order(pseudo$time, -pseudo$event), ]
  tu <- unique(pseudo$time[pseudo$event == 1])
  ltu <- length(tu)
  tu <- matrix(tu, byrow = TRUE, ncol = ltu, nrow = ltmax)
  tgiven <- matrix(tmax, byrow = FALSE, ncol = ltu, nrow = ltmax)
  inx <- apply(tgiven >= tu, 1, sum)
  KM.omit <- surv.omit(pseudo)
  KM.omit <- KM.omit[, inx, drop = FALSE]
  KM.omit <- data.frame(KM.omit[order(pseudo$id),])
  KM.tot <- surv.tot(pseudo)[inx]
  KM.tot <- matrix(KM.tot, byrow = TRUE, nrow = howmany, ncol = length(tmax))
  KM.tot <- data.frame(KM.tot[order(pseudo$id),])
  for (k in 1:ltmax){
    ind <- obs$obs1>t0[k,1]&obs$obs2>t0[k,2]
    LOO_sum <- sum(ind)-ind
    POs <- (KM.omit[,k]*ind +(KM.omit[,k]-KM.tot[,k])*LOO_sum)/(KM.tot[,k]*KM.omit[,k])
    obs[,m+k] <- POs
  }
  colnames(obs) <- c(old_name,paste0("PO_t",1:ltmax,""))
  return(obs)
}


###--------POs for the Lin and Ying estimator adjusted for bivariate censoring----------
# obs: a data frame containing the observed times (obs1,obs2) and their corresponding indicators (delta1,delta2)
# t0: a data frame of size (K*2) with K bivariate time point. The two columns correspond to the two coordinates (t1,t2), and each of the K rows corresponds to a different time point.
PO_func_bivar_cens_mult <- function(obs,t0){
  old_names <- colnames(obs)
  p <- ncol(obs)
  k <- nrow(t0)
  km_fit1 <- survfit(Surv(obs$obs1, 1-obs$delta1) ~ 1, data=obs)
  km_fit2 <- survfit(Surv(obs$obs2, 1-obs$delta2) ~ 1, data=obs)
  S_hat <- data.frame(matrix(NA,1,k))
  for (i in 1:k){
    t1 <- t0[i,1]
    t2 <- t0[i,2]
    
    #compute the estimate for the bi-variate survival
    obs$ind <- obs$obs1>t1&obs$obs2>t2
    numerator <- mean(obs$ind)
    
    survs1 <- km_fit1$surv[km_fit1$time<=t1]
    survs2 <- km_fit2$surv[km_fit2$time<=t2]
    denominator <- survs1[length(survs1)]*survs2[length(survs2)]
    
    S_hat[i] <- numerator/denominator
  }
  S_hat_mat <- do.call("rbind", replicate(n, S_hat, simplify = FALSE))
  #compute leave-one-out estimate - using for loop
  LOO <- data.frame(matrix(NA,n,k))
  for (i in 1:n){
    obs_i <- obs[-i,]
    km_fit1 <- survfit(Surv(obs_i$obs1, 1-obs_i$delta1) ~ 1, data=obs_i)
    km_fit2 <- survfit(Surv(obs_i$obs2, 1-obs_i$delta2) ~ 1, data=obs_i)
    for (j in 1:k){
      t1 <- t0[j,1]
      t2 <- t0[j,2]
      obs_i$ind_i <- obs_i$obs1>t1&obs_i$obs2>t2
      numerator <- mean(obs_i$ind_i)
      survs1 <- km_fit1$surv[km_fit1$time<=t1]
      survs2 <- km_fit2$surv[km_fit2$time<=t2]
      denominator <- survs1[length(survs1)]*survs2[length(survs2)]
      LOO[i,j] <- numerator/denominator
    }
  }
  #compute POs
  obs[,(p+1):(p+k)] <- n*S_hat_mat-(n-1)*LOO
  colnames(obs) <- c(old_names, paste0("PO",1:k,""))
  obs$id <- 1:n
  return(obs)
}

