
###-------Bivariate logistic failure times-----

#Bivariate survival
S <- function(t1,t2,a1,a2,lambda,c){ #c=exp(beta*Z)
  return(exp(-lambda*log(a1*t1+a2*t2))*c/(1+exp(-lambda*log(a1*t1+a2*t2))*c))
}

#Univariate CDF
F1 <- function(t,a1,lambda,c){
  return(1-(c*(a1*t)^(-lambda))/(1+c*(a1*t)^(-lambda)))
}
#inverse of univariate CDF
F1_inv <- function(u,a1,lambda,c){
  return((1/a1)*((1-u)/(c*u))^(-1/lambda))
}

#Conditional CDF
F21 <- function(t2,t1,a1,a2,lambda,c){
  return(1-((a1*t1)/(a1*t1+a2*t2))*S(t1,t2,a1,a2,lambda,c)*(1-S(t1,t2,a1,a2,lambda,c))/(S(t1,0,a1,a2,lambda,c)*(1-S(t1,0,a1,a2,lambda,c))))
}

#Find the inverse conditional CDF by setting to zero the following difference
f <- function(x,v,t1i,ci){
  return(F21(t2=x,t1=t1i,a1=a1,a2=a2,lambda=lambda,c=ci)-v)
}

#generate bivariate logistic data
gen_logit_data <- function(a1,a2,lambda,beta, Z,seed){
  n <- length(Z) #Z is a numeric vector of size n (a univariate covariate)
  c_vec <- exp(beta*Z)
  Failure_times <- data.frame(T1=NA,T2=NA,covar=NA)
  set.seed(seed)
  unifs <- runif(2*n)
  for (j in 1:n) {
    u <- unifs[j]
    v <- unifs[n+j]
    t1i <- F1_inv(u,a1=a1,lambda=lambda,c=c_vec[j])
    root <- uniroot(f, c(0,10^7),v=v,t1=t1i,c=c_vec[j])
    t2i <- root$root
    Failure_times[j,] <- c(t1i,t2i,Z[j])
  }
  #add univariate censoring
  set.seed(1000+seed)
  C <- rexp(n, rate=0.3)
  obs <- data.frame(obs1=pmin(Failure_times[,1],C),delta1=as.numeric(Failure_times[,1]<=C),obs2=pmin(Failure_times[,2],C), delta2=as.numeric(Failure_times[,2]<=C), C=C, Z=Failure_times[,3])
  obs$id <- 1:n
  return(obs)
}

#generate bivariate logistic data with bivariate censoring
gen_logit_data_bivar_cens <- function(a1,a2,lambda,beta, Z,seed,low_cens=FALSE){
  n <- length(Z) #Z is a numeric vector of size n (a univariate covariate)
  c_vec <- exp(beta*Z)
  Failure_times <- data.frame(T1=NA,T2=NA,covar=NA)
  set.seed(seed)
  unifs <- runif(2*n)
  for (j in 1:n) {
    u <- unifs[j]
    v <- unifs[n+j]
    t1i <- F1_inv(u,a1=a1,lambda=lambda,c=c_vec[j])
    root <- uniroot(f, c(0,10^7),v=v,t1=t1i,c=c_vec[j])
    t2i <- root$root
    Failure_times[j,] <- c(t1i,t2i,Z[j])
  }
  #add univariate censoring
  set.seed(1000+seed)
  C1 <- rexp(n, rate=0.3)
  set.seed(10000+seed)
  C2 <- rexp(n, rate=0.2)
  if (isTRUE(low_cens)){
    set.seed(1000+seed)
    C1 <- rexp(n, rate=0.1)
    set.seed(10000+seed)
    C2 <- rexp(n, rate=0.05)
  }
  obs <- data.frame(obs1=pmin(Failure_times[,1],C1),delta1=as.numeric(Failure_times[,1]<=C1),obs2=pmin(Failure_times[,2],C2), delta2=as.numeric(Failure_times[,2]<=C2), Z=Failure_times[,3])
  obs$id <- 1:n
  return(obs)
}


###-------Bivariate lognormal failure times-----

#generate bivariate lognormal data
gen_LN_data <- function(Sigma ,gamma1=0.2,gamma2=1, Z,seed){
  n <- length(Z) #Z is a numeric vector of size n (a univariate covariate)
  mu <- data.frame(mu1=gamma1*Z,mu2=gamma2*Z)
  set.seed(seed)
  Y2 <- t(apply(X=mu,MARGIN = 1,FUN = mvrnorm,n=1,Sigma=Sigma))
  Failure_times <- exp(Y2)
  #add univariate censoring
  set.seed(1000+seed)
  C <- rexp(n, rate=0.3)
  obs <- data.frame(obs1=pmin(Failure_times[,1],C),delta1=as.numeric(Failure_times[,1]<=C),obs2=pmin(Failure_times[,2],C), delta2=as.numeric(Failure_times[,2]<=C), C=C, Z=Z)
  obs$id <- 1:n
  return(obs)
}