args <- commandArgs(trailingOnly = TRUE)
seed=as.numeric(args[1])
cat(paste(" Seed =", seed,"\n") ) # Print the random seed 
set.seed(seed)
library(fda)
library(mvnfast)  # to use rmvn and dmvn function to draw samples
library(gtools)   # to use rdirichlet function

#-------------------------------------# 
#-------------------------------------# 
#-------------------------------------# 
#Part 1: Simulation# 
#-------------------------------------# 
#-------------------------------------# 
#-------------------------------------# 
# I groups, each with ni patients, each with T measurements. There are K clusters.
# 12 groups, each with 50 patients, each with 20 measurements. There are 2 clusters.

ni=50
I=12
T=20
K=2

# Group ID and patient ID
Group=c()
for (i in 1:I){Group=c(Group,rep(i,ni*T))}
Subject=c()
for (i in 1:(I*ni)){Subject=c(Subject,rep(i,T))}

# spline
tobs = seq(0,1,length = T+1)[2:(T+1)]
Time=c(rep(tobs,I*ni))
nobs = length(tobs)
knots = c(0.05,0.33,0.67,1);
norder = 4;
nbasis = length(knots) + norder - 2;
basis = create.bspline.basis(c(min(tobs),max(tobs)),nbasis,norder,knots);
X_spline   = eval.basis(Time, basis);

# spline part parameters
tau_u<-1
tau_w<-1
tau_e<-1
pi<-c(0.5, 0.5)
v=rnorm(I,0,tau_u^(-1/2))
v_long=c()
for (i in 1:I){v_long=c(v_long,rep(v[i],ni*T))}
w=rnorm(I*ni,0,tau_w^(-1/2))
w_long=c()
for (i in 1:(i*ni)){w_long=c(w_long,rep(w[i],T))}
epsilon=rnorm(I*ni*T,0,tau_e^(-1/2))
alpha1=c(1,1,1,1,1,1)
alpha2=c(1,2,3,4,5,6)
ci.true=t(matrix(c(rep(c(1,0),6), rep(c(0,1),6)),ncol=2,byrow = TRUE))
muktl=rep(0,I*ni*T)
for (i in 1:I){if(ci.true[1,i]==1){muktl[((i-1)*ni*T+1):(i*ni*T)]=X_spline[((i-1)*ni*T+1):(i*ni*T),]%*%alpha1} else {muktl[((i-1)*ni*T+1):(i*ni*T)]=X_spline[((i-1)*ni*T+1):(i*ni*T),]%*%alpha2}}

# biomarker measurements
z=muktl + v_long + w_long + epsilon

# survival time
lambda<-1.96
r<-1.5
theta_1<-0
theta_2<-0

# Weibull latent event times
v <- runif(ni*I)
t=rep(0,ni*I)
for (j in 1:I){
  for (i in ((j-1)*ni+1):(j*ni)){if (ci.true[1,j]==1) {t[i] <- (- log(v[i]) / (lambda * exp(theta_1)))^(1 / r)} else {t[i] <- (- log(v[i]) / (lambda * exp(theta_2)))^(1 / r)}}
}

# censoring times
C <- rexp(n=ni*I, rate=0.01)

# follow-up times and event indicators
observed_time_short <- pmin(t, C)
observed_time=c()
for (i in 1:(ni*I)){observed_time=c(observed_time,rep(observed_time_short[i],T))}
omega_short <- as.numeric(t <= C)
omega=c()
for (i in 1:(ni*I)){omega=c(omega,rep(omega_short[i],T))}

# create data frame
longdata_full <- cbind(Group,Subject,Time,z,omega,observed_time)










#-------------------------------------# 
#-------------------------------------# 
#-------------------------------------# 
#Part 2: Proposed Model Interim 1# 
#-------------------------------------# 
#-------------------------------------# 
#-------------------------------------# 

nknot=2
mcsim=8000
nburnin=4000
ncluster=2
ci.true=t(ci.true)
#-------------------------------------#


#-------------------------------------#
# A. Longitudinal part data organize  #
#-------------------------------------#

# interim analysis 1
longdata <- longdata_full[Subject<=ni*Group-20,]
data<-longdata

#---------------------------------------------#
# I. Spline part creation                     #
#---------------------------------------------#

#-------------------------------------------------------------#
# II. extract dimension parameters from the data:             #
# number of observations of each patient (N_obs);             #
# number of total observations(N_tot)= N_subj*length(t);      #
# number of subjects in the data (N_subj);                    #
# number of groups in the data (N_grp);                       #
# number of patients per group (N_pt);                        #
#-------------------------------------------------------------#
# rep(number of time points, number of individuals)
N_obs<-as.numeric(as.vector(table(data[,2])))
# total number of observations
N_tot<-sum(N_obs)
# number of patients in the data (N_subj)
N_subj<-as.numeric(length(unique(data[,2])))
# number of groups in the data (N_grp)
N_grp<-as.numeric(length(unique(data[,1])))
N_grp_y<-N_grp
# number of patients per group as a vector
N_pt<-rep(N_subj/N_grp,N_grp)
# number of patients per group
N_pts<-N_subj/N_grp
# subject ID indicator (c(1,1,1,1,1,2,2,2,2,2,.....))
sub<-as.numeric(match(data[,2], unique(data[,2])))
# group ID indicator
grp<-as.numeric(match(data[,1], unique(data[,1])))
grp_short<-c()
for (i in 1:N_subj){
  grp_short=c(grp_short,grp[1+(i-1)*(N_tot/N_subj)])
}

#------------------------------------------------#
# III. Vectors to store the posterior samples    #
# of parameters and initialize them:             #
# alpha: all the spline coefficients             #
# beta: all the r.e., v_i and w_ij               #
# tau_u: precision of v_i, group-specific r.e.   #
# tau_w: precision of w_ij, subject-specfic r.e. #
# tau_e: precision of e_ijt, measurement error   #
#------------------------------------------------#

alpha<-array(NA, dim=c(mcsim, 6, ncluster))
beta<-matrix(NA, nrow=mcsim, ncol=N_grp+N_subj)
tau_u<-tau_w<-tau_e<-rep(NA,mcsim)
lambda<-rep(NA,mcsim)
r<-rep(NA,mcsim)
theta_2<-rep(NA,mcsim)
# Initialize
alpha[1, ,]<- rep(1.5, 6)
beta[1,]<-rep(0,N_grp+N_subj)
tau_u[1]<-tau_w[1]<-tau_e[1]<- 1
lambda[1]<- 1.5
r[1]<- 1.05
theta_1 <- 0
theta_2[1]<- -0.75

# Design matrix for longitudinal part
# Dimension of X_spline: (N_tot, 3+nknot)
tobs = seq(0,1,length = T+1)[2:(T+1)]
Time=c(rep(tobs,N_grp*N_pts))
nobs = length(tobs)
knots = c(0.05,0.33,0.67,1);
norder = 4;
nbasis = length(knots) + norder - 2;
basis = create.bspline.basis(c(min(tobs),max(tobs)),nbasis,norder,knots);
X_spline   = eval.basis(Time, basis);
# Dimension of X_u: (N_tot, N_grp)
X_u <-NULL
for (i in 1: N_grp){
  grp.temp<-grp
  grp.temp[which(grp==i)]<-1
  grp.temp[which(grp!=i)]<-0
  X_u<-cbind(X_u,grp.temp)
}
# Dimension of X_w: (N_tot, N_subj)
X_w<-matrix(0, nrow=N_tot,ncol=N_subj)
for (i in 1:N_subj){
  X_w[which(sub==i),i]<-1
}
# Dimension of X: (N_tot, 6+N_grp+N_subj)
X<-cbind(X_spline, X_u, X_w)
X_random<-cbind(X_u, X_w)

#-------------------------------------#
# B. response part data organize        #
#-------------------------------------#  
# biomarker measurements
z<-as.numeric(data[,4]) # dependent variable, vector format
# mean biomarker measurement for each patient
z_mean<-c()
for (i in 1:N_subj){
  z_mean=c(z_mean,mean(z[which(sub==i)]))
}

# censoring indicator variable
omega<-as.numeric(data[,5])
omega_short<-c()
for (i in 1:N_subj){
  omega_short=c(omega_short,omega[1+(i-1)*(N_tot/N_subj)])
}

# observed time variable
t<-as.numeric(data[,6])
# a vector of t which does not contain replicates for each patient
t_short<-c()
for (i in 1:N_subj){
  t_short=c(t_short,t[1+(i-1)*(N_tot/N_subj)])
}


#------------------------------------------------------------------#  
#----------------------------------------------#
# C. latent indicator part                     #
# ci: latent indicator for each subgroup       #
# pi: Pr(c_ik=1), group i belongs to cluster k #
# v.prob: allocation probability (tau_ik)      #
#----------------------------------------------#    
pi.prior<-rep(2,ncluster) # prior specification, use 2 instead of 1/ncluster to avoid sampled p close to 0
pi<-array(dim=c(mcsim,ncluster,N_grp))
v.prob<-array(dim=c(mcsim,N_grp,ncluster))
ci<-array(dim=c(mcsim,ncluster,N_grp))

for (i in 1:N_grp){
  pi[1, ,i]<-c(0.5, 0.5)
  ci[1, ,i]<-ci.true[i,]
  v.prob[1,i,]<-rep(1/ncluster, ncluster)
}


#----------------------------------------------#
# D. MCMC part                                 #
#                                              #
#----------------------------------------------#    
a<-0.001; b<-0.001

alpha.subgrp<-array(dim=c((mcsim), 6, N_grp))
pi.subgrp<-matrix(NA,nrow=(mcsim), ncol=N_grp)

set.seed(seed)
for (iter in 2:mcsim){
  #--------------#  
  # update ci    #
  #--------------#  
  for (i in 1:N_grp){
    
    idx.temp<-c(i, (sum(N_pt[1:i])+N_grp-(N_pt[i]-1)):(N_grp+sum(N_pt[1:i])))
    sig.long <- (tau_e[iter-1]^(-1)*diag(nrow=N_pt[i]*N_obs[1], ncol=N_pt[i]*N_obs[1]))
    dens.long<-numeric(ncluster)
    
    mean.long_grp1 <- X_spline[grp==i,]%*%alpha[iter-1, ,1]+X_random[grp==i,idx.temp]%*%beta[iter-1,idx.temp]
    mean.long_grp2 <- X_spline[grp==i,]%*%alpha[iter-1, ,2]+X_random[grp==i,idx.temp]%*%beta[iter-1,idx.temp]
    dens.long[1]<-dmvn(z[grp==i], mean.long_grp1, sig.long, log=TRUE) + sum(omega_short[grp_short==i]*log(exp(-lambda[iter-1]*t_short[grp_short==i]^r[iter-1]*exp(theta_1))*(lambda[iter-1]*exp(theta_1)*r[iter-1]*t_short[grp_short==i]^(r[iter-1]-1))) + (1-omega_short[grp_short==i])*(-lambda[iter-1]*t_short[grp_short==i]^r[iter-1]*exp(theta_1))) #take log
    dens.long[2]<-dmvn(z[grp==i], mean.long_grp2, sig.long, log=TRUE) + sum(omega_short[grp_short==i]*log(exp(-lambda[iter-1]*t_short[grp_short==i]^r[iter-1]*exp(theta_2[iter-1]))*(lambda[iter-1]*exp(theta_2[iter-1])*r[iter-1]*t_short[grp_short==i]^(r[iter-1]-1))) + (1-omega_short[grp_short==i])*(-lambda[iter-1]*t_short[grp_short==i]^r[iter-1]*exp(theta_2[iter-1]))) #take log
    
    if ((dens.long[1]-dens.long[2])>10^2){
      v.prob[iter,i,]<-c(1,0)
    } else if ((dens.long[1]-dens.long[2])< -10^2) {
      v.prob[iter,i,]<-c(0,1)
    } else {
      num<-(pi[iter-1,1,i]*exp(dens.long[1]-dens.long[2]))
      v.prob[iter,i,]<-c(num/(num+pi[iter-1,2,i]),1-num/(num+pi[iter-1,2,i]))
    }
    
    ci[iter,,i]<-as.vector(rmultinom(1, 1, prob=v.prob[iter,i,]))
    
    #-------------# 
    # update pi   #
    #-------------#
    
    # there is a specific pi for each group i
    pi[iter, ,i] <-rdirichlet(1, (ci[iter, ,i]+pi.prior))
  }
  
  #---------------------------#  
  # update alpha              #
  #---------------------------#   
  for (k in 1:ncluster){
    clust.temp<-which(ci[iter,k ,]==1) #subgroup indicator of those belonging to cluster k
    sigma.spline<-diag(c(rep(10^4, 6)), nrow=6, ncol=6)
    if (length(clust.temp)==0){
      mean.alpha<-matrix(rep(0,6),nrow=6,ncol=1)
      sigma.alpha<-solve(tau_e[iter-1]*(diag(6))+solve(sigma.spline))
      
    } else {
      
      myxspline<-NULL
      myxrandom<-NULL
      myz<-NULL
      for (j in 1:length(clust.temp)){
        myxspline<-rbind(myxspline, X_spline[grp==clust.temp[j],])    #longitudinal part
        myxrandom<-rbind(myxrandom, X_random[grp==clust.temp[j],])    #longitudinal part
        myz<-c(myz, z[grp==clust.temp[j]])                            #longitudinal part
      }
      
      sigma.alpha<-solve(tau_e[iter-1]*crossprod(myxspline)+solve(sigma.spline))
      
      mean.alpha<-tau_e[iter-1]*sigma.alpha%*%crossprod(myxspline,(myz-myxrandom%*%matrix(beta[iter-1,],ncol=1)))
      
    } #end else
    
    # update alpha
    alpha[iter, ,k] <-rmvn(1, mean.alpha, sigma.alpha)
    
  }# end k
  
  #---------------------------#  
  # update beta and tau_e     #
  #---------------------------#
  res_sum<-0
  for (i in 1: N_grp){
    sigma.random<-diag(c(1/tau_u[iter-1], rep(1/tau_w[iter-1], N_pt[i])), nrow=1+N_pt[i], ncol=1+N_pt[i])
    idx.temp<-c(i, (sum(N_pt[1:i])+N_grp-(N_pt[i]-1)):(N_grp+sum(N_pt[1:i])))
    sigma.beta<-solve(solve(sigma.random)+crossprod(X_random[grp==i,idx.temp])*tau_e[iter-1])
    
    k.temp<-which(ci[iter, ,i]==1)
    mu.beta<-tau_e[iter-1]*sigma.beta%*%crossprod((X_random[grp==i,idx.temp]),(z[grp==i]-X_spline[grp==i,]%*%alpha[iter, ,k.temp]))
    
    beta[iter,idx.temp]<-rmvn(1, mu.beta, sigma.beta)
    
    res_temp<-z[grp==i]-X_spline[grp==i,]%*%alpha[iter, ,k.temp]-X_random[grp==i,idx.temp]%*%beta[iter,idx.temp]
    res_sum<-res_sum+sum(res_temp^2)
    
    # re-organize results
    pi.subgrp[(iter),i]<-pi[iter,k.temp,i]
    alpha.subgrp[(iter),,i]<-alpha[iter, , k.temp]
  }
  
  tau_e[iter]<-rgamma(1, a+0.5*N_tot, scale=(b+0.5*res_sum)^(-1))
  #---------------------------#    
  # update tau_u#
  #---------------------------#
  tau_u[iter]<-rgamma(1, a+0.5*N_grp, scale=(b+0.5*sum(beta[iter, (1:N_grp)]^2))^(-1))
  
  #---------------------------#
  # update tau_w #
  #---------------------------#
  tau_w[iter]<-rgamma(1, a+0.5*N_subj, scale=(b+0.5*sum(beta[iter, ((N_grp+1):(N_grp+N_subj))]^2))^(-1))
  
  #----------------------------------#
  # update theta_2, lambda, r#
  #----------------------------------#
  cluster.ind=rep(0,N_subj)
  for (i in 1:N_grp){if (ci[iter,1,i]==1){cluster.ind[((i-1)*N_pt[1]+1):(i*N_pt[1])]=rep(1,N_pt[1])} else {cluster.ind[((i-1)*N_pt[1]+1):(i*N_pt[1])]=rep(2,N_pt[1])}}
  
  # log density of theta_2
  logden.theta_2 <- function(theta_2, lambda, r) {
    return(sum(omega_short[cluster.ind==2]*log(exp(-lambda*t_short[cluster.ind==2]^r*exp(theta_2))*(lambda*exp(theta_2)*r*t_short[cluster.ind==2]^(r-1))) + (1-omega_short[cluster.ind==2])*(-lambda*t_short[cluster.ind==2]^r*exp(theta_2))) + dunif(theta_2, -1, 1, log=TRUE))
  }
  # log density of lambda
  logden.lambda <- function(lambda, theta_2, r) {
    return(sum(omega_short[cluster.ind==1]*log(exp(-lambda*t_short[cluster.ind==1]^r*exp(theta_1))*(lambda*exp(theta_1)*r*t_short[cluster.ind==1]^(r-1))) + (1-omega_short[cluster.ind==1])*(-lambda*t_short[cluster.ind==1]^r*exp(theta_1)))
           + sum(omega_short[cluster.ind==2]*log(exp(-lambda*t_short[cluster.ind==2]^r*exp(theta_2))*(lambda*exp(theta_2)*r*t_short[cluster.ind==2]^(r-1))) + (1-omega_short[cluster.ind==2])*(-lambda*t_short[cluster.ind==2]^r*exp(theta_2))) + dgamma(lambda, shape=0.1, rate=0.1, log=TRUE))
  }
  # log density of r
  logden.r <- function(r, theta_2, lambda) {
    return(sum(omega_short[cluster.ind==1]*log(exp(-lambda*t_short[cluster.ind==1]^r*exp(theta_1))*(lambda*exp(theta_1)*r*t_short[cluster.ind==1]^(r-1))) + (1-omega_short[cluster.ind==1])*(-lambda*t_short[cluster.ind==1]^r*exp(theta_1)))
           + sum(omega_short[cluster.ind==2]*log(exp(-lambda*t_short[cluster.ind==2]^r*exp(theta_2))*(lambda*exp(theta_2)*r*t_short[cluster.ind==2]^(r-1))) + (1-omega_short[cluster.ind==2])*(-lambda*t_short[cluster.ind==2]^r*exp(theta_2))) + dgamma(r, shape=0.1, rate=0.1, log=TRUE))
  }
  
  # update theta_2
  theta_2.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_2.proposed = theta_2.proposal(theta_2[iter-1])
  prob.theta_2 = exp(logden.theta_2(theta_2.proposed, lambda[iter-1], r[iter-1])-logden.theta_2(theta_2[iter-1], lambda[iter-1], r[iter-1]))
  if (runif(1) < prob.theta_2){
    theta_2[iter] = theta_2.proposed
  }else{
    theta_2[iter] = theta_2[iter-1]
  }
  # update lambda
  lambda.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  lambda.proposed = lambda.proposal(lambda[iter-1])
  prob.lambda = exp(logden.lambda(lambda.proposed, theta_2[iter], r[iter-1])-logden.lambda(lambda[iter-1], theta_2[iter], r[iter-1]))
  if (runif(1) < prob.lambda){
    lambda[iter] = lambda.proposed
  }else{
    lambda[iter] = lambda[iter-1]
  }
  # update r
  r.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r.proposed = r.proposal(r[iter-1])
  prob.r = exp(logden.r(r.proposed, theta_2[iter], lambda[iter])-logden.r(r[iter-1], theta_2[iter], lambda[iter]))
  if (runif(1) < prob.r){
    r[iter] = r.proposed
  }else{
    r[iter] = r[iter-1]
  }
}

# Qf = 0.05, Q = 0.8
# check
print("Proposed Model Interim 1")
for (i in 1:I){print(mean(ci[(nburnin+1):mcsim,1,i]))}
sum(((1/lambda)^(1/r)*(log(2))^(1/r))[(nburnin+1):mcsim]>0.6)/(mcsim-nburnin)
sum(((1/lambda/exp(theta_2))^(1/r)*(log(2))^(1/r))[(nburnin+1):mcsim]>0.6)/(mcsim-nburnin)

longdata <- longdata_full
if (sum(((1/lambda)^(1/r)*(log(2))^(1/r))[(nburnin+1):mcsim]>0.6)/(mcsim-nburnin)<0.05){
  for (i in 1:I){if (mean(ci[(nburnin+1):mcsim,1,i])>0.5){
    longdata <- longdata[longdata[,1] !=i,]
  }}}
if (sum(((1/lambda/exp(theta_2))^(1/r)*(log(2))^(1/r))[(nburnin+1):mcsim]>0.6)/(mcsim-nburnin)<0.05){
  for (i in 1:I){if (mean(ci[(nburnin+1):mcsim,1,i])<=0.5){
    longdata <- longdata[(longdata[,1] !=i),]
  }}}
ci.true <- t(matrix(rep(c(0.5,0.5),length(unique(longdata[,1]))),ncol=2,byrow = TRUE))










#-------------------------------------# 
#-------------------------------------# 
#-------------------------------------# 
#Part 3: Proposed Model Interim 2# 
#-------------------------------------# 
#-------------------------------------# 
#-------------------------------------# 
if (length(unique(longdata[,1]))!=0){
  nknot=2
  mcsim=8000
  nburnin=4000
  ncluster=2
  ci.true=t(ci.true)
  #-------------------------------------#
  
  
  #-------------------------------------#
  # A. Longitudinal part data organize  #
  #-------------------------------------#
  data<-longdata
  
  #---------------------------------------------#
  # I. Spline part creation                     #
  #---------------------------------------------#
  
  #-------------------------------------------------------------#
  # II. extract dimension parameters from the data:             #
  # number of observations of each patient (N_obs);             #
  # number of total observations(N_tot)= N_subj*length(t);      #
  # number of subjects in the data (N_subj);                    #
  # number of groups in the data (N_grp);                       #
  # number of patients per group (N_pt);                        #
  #-------------------------------------------------------------#
  # rep(number of time points, number of individuals)
  N_obs<-as.numeric(as.vector(table(data[,2])))
  # total number of observations
  N_tot<-sum(N_obs)
  # number of patients in the data (N_subj)
  N_subj<-as.numeric(length(unique(data[,2])))
  # number of groups in the data (N_grp)
  N_grp<-as.numeric(length(unique(data[,1])))
  N_grp_y<-N_grp
  # number of patients per group as a vector
  N_pt<-rep(N_subj/N_grp,N_grp)
  # number of patients per group
  N_pts<-N_subj/N_grp
  # subject ID indicator (c(1,1,1,1,1,2,2,2,2,2,.....))
  sub<-as.numeric(match(data[,2], unique(data[,2])))
  # group ID indicator
  grp<-as.numeric(match(data[,1], unique(data[,1])))
  grp_short<-c()
  for (i in 1:N_subj){
    grp_short=c(grp_short,grp[1+(i-1)*(N_tot/N_subj)])
  }
  
  #------------------------------------------------#
  # III. Vectors to store the posterior samples    #
  # of parameters and initialize them:             #
  # alpha: all the spline coefficients             #
  # beta: all the r.e., v_i and w_ij               #
  # tau_u: precision of v_i, group-specific r.e.   #
  # tau_w: precision of w_ij, subject-specfic r.e. #
  # tau_e: precision of e_ijt, measurement error   #
  #------------------------------------------------#
  
  alpha<-array(NA, dim=c(mcsim, 6, ncluster))
  beta<-matrix(NA, nrow=mcsim, ncol=N_grp+N_subj)
  tau_u<-tau_w<-tau_e<-rep(NA,mcsim)
  lambda<-rep(NA,mcsim)
  r<-rep(NA,mcsim)
  theta_2<-rep(NA,mcsim)
  # Initialize
  alpha[1, ,]<- rep(1.5, 6)
  beta[1,]<-rep(0,N_grp+N_subj)
  tau_u[1]<-tau_w[1]<-tau_e[1]<- 1
  lambda[1]<- 1.5
  r[1]<- 1.05
  theta_1 <- 0
  theta_2[1]<- -0.75
  
  # Design matrix for longitudinal part
  # Dimension of X_spline: (N_tot, nbasis)
  tobs = seq(0,1,length = T+1)[2:(T+1)]
  Time=c(rep(tobs,N_grp*N_pts))
  nobs = length(tobs)
  knots = c(0.05,0.33,0.67,1);
  norder = 4;
  nbasis = length(knots) + norder - 2;
  basis = create.bspline.basis(c(min(tobs),max(tobs)),nbasis,norder,knots);
  X_spline   = eval.basis(Time, basis);
  # Dimension of X_u: (N_tot, N_grp)
  X_u <-NULL
  for (i in 1: N_grp){
    grp.temp<-grp
    grp.temp[which(grp==i)]<-1
    grp.temp[which(grp!=i)]<-0
    X_u<-cbind(X_u,grp.temp)
  }
  # Dimension of X_w: (N_tot, N_subj)
  X_w<-matrix(0, nrow=N_tot,ncol=N_subj)
  for (i in 1:N_subj){
    X_w[which(sub==i),i]<-1
  }
  # Dimension of X: (N_tot, 6+N_grp+N_subj)
  X<-cbind(X_spline, X_u, X_w)
  X_random<-cbind(X_u, X_w)
  
  #-------------------------------------#
  # B. response part data organize        #
  #-------------------------------------#  
  # biomarker measurements
  z<-as.numeric(data[,4]) # dependent variable, vector format
  # mean biomarker measurement for each patient
  z_mean<-c()
  for (i in 1:N_subj){
    z_mean=c(z_mean,mean(z[which(sub==i)]))
  }
  
  # censoring indicator variable
  omega<-as.numeric(data[,5])
  omega_short<-c()
  for (i in 1:N_subj){
    omega_short=c(omega_short,omega[1+(i-1)*(N_tot/N_subj)])
  }
  
  # observed time variable
  t<-as.numeric(data[,6])
  # a vector of t which does not contain replicates for each patient
  t_short<-c()
  for (i in 1:N_subj){
    t_short=c(t_short,t[1+(i-1)*(N_tot/N_subj)])
  }
  
  
  #------------------------------------------------------------------#  
  #----------------------------------------------#
  # C. latent indicator part                     #
  # ci: latent indicator for each subgroup       #
  # pi: Pr(c_ik=1), group i belongs to cluster k #
  # v.prob: allocation probability (tau_ik)      #
  #----------------------------------------------#    
  pi.prior<-rep(2,ncluster) # prior specification, use 2 instead of 1/ncluster to avoid sampled p close to 0
  pi<-array(dim=c(mcsim,ncluster,N_grp))
  v.prob<-array(dim=c(mcsim,N_grp,ncluster))
  ci<-array(dim=c(mcsim,ncluster,N_grp))
  
  for (i in 1:N_grp){
    pi[1, ,i]<-c(0.5, 0.5)
    ci[1, ,i]<-ci.true[i,]
    v.prob[1,i,]<-rep(1/ncluster, ncluster)
  }
  
  
  #----------------------------------------------#
  # D. MCMC part                                 #
  #                                              #
  #----------------------------------------------#    
  a<-0.001; b<-0.001
  
  alpha.subgrp<-array(dim=c((mcsim), 6, N_grp))
  pi.subgrp<-matrix(NA,nrow=(mcsim), ncol=N_grp)
  
  set.seed(seed)
  for (iter in 2:mcsim){
    #--------------#  
    # update ci    #
    #--------------#  
    for (i in 1:N_grp){
      
      idx.temp<-c(i, (sum(N_pt[1:i])+N_grp-(N_pt[i]-1)):(N_grp+sum(N_pt[1:i])))
      sig.long <- (tau_e[iter-1]^(-1)*diag(nrow=N_pt[i]*N_obs[1], ncol=N_pt[i]*N_obs[1]))
      dens.long<-numeric(ncluster)
      
      mean.long_grp1 <- X_spline[grp==i,]%*%alpha[iter-1, ,1]+X_random[grp==i,idx.temp]%*%beta[iter-1,idx.temp]
      mean.long_grp2 <- X_spline[grp==i,]%*%alpha[iter-1, ,2]+X_random[grp==i,idx.temp]%*%beta[iter-1,idx.temp]
      dens.long[1]<-dmvn(z[grp==i], mean.long_grp1, sig.long, log=TRUE) + sum(omega_short[grp_short==i]*log(exp(-lambda[iter-1]*t_short[grp_short==i]^r[iter-1]*exp(theta_1))*(lambda[iter-1]*exp(theta_1)*r[iter-1]*t_short[grp_short==i]^(r[iter-1]-1))) + (1-omega_short[grp_short==i])*(-lambda[iter-1]*t_short[grp_short==i]^r[iter-1]*exp(theta_1))) #take log
      dens.long[2]<-dmvn(z[grp==i], mean.long_grp2, sig.long, log=TRUE) + sum(omega_short[grp_short==i]*log(exp(-lambda[iter-1]*t_short[grp_short==i]^r[iter-1]*exp(theta_2[iter-1]))*(lambda[iter-1]*exp(theta_2[iter-1])*r[iter-1]*t_short[grp_short==i]^(r[iter-1]-1))) + (1-omega_short[grp_short==i])*(-lambda[iter-1]*t_short[grp_short==i]^r[iter-1]*exp(theta_2[iter-1]))) #take log
      
      if ((dens.long[1]-dens.long[2])>10^2){
        v.prob[iter,i,]<-c(1,0)
      } else if ((dens.long[1]-dens.long[2])< -10^2) {
        v.prob[iter,i,]<-c(0,1)
      } else {
        num<-(pi[iter-1,1,i]*exp(dens.long[1]-dens.long[2]))
        v.prob[iter,i,]<-c(num/(num+pi[iter-1,2,i]),1-num/(num+pi[iter-1,2,i]))
      }
      
      ci[iter,,i]<-as.vector(rmultinom(1, 1, prob=v.prob[iter,i,]))
      
      #-------------# 
      # update pi   #
      #-------------#
      
      # there is a specific pi for each group i
      pi[iter, ,i] <-rdirichlet(1, (ci[iter, ,i]+pi.prior))
    }
    
    #---------------------------#  
    # update alpha              #
    #---------------------------#   
    for (k in 1:ncluster){
      clust.temp<-which(ci[iter,k ,]==1) #subgroup indicator of those belonging to cluster k
      sigma.spline<-diag(c(rep(10^4, 6)), nrow=6, ncol=6)
      if (length(clust.temp)==0){
        mean.alpha<-matrix(rep(0,6),nrow=6,ncol=1)
        sigma.alpha<-solve(tau_e[iter-1]*(diag(6))+solve(sigma.spline))
        
      } else {
        
        myxspline<-NULL
        myxrandom<-NULL
        myz<-NULL
        for (j in 1:length(clust.temp)){
          myxspline<-rbind(myxspline, X_spline[grp==clust.temp[j],])    #longitudinal part
          myxrandom<-rbind(myxrandom, X_random[grp==clust.temp[j],])    #longitudinal part
          myz<-c(myz, z[grp==clust.temp[j]])                            #longitudinal part
        }
        
        sigma.alpha<-solve(tau_e[iter-1]*crossprod(myxspline)+solve(sigma.spline))
        
        mean.alpha<-tau_e[iter-1]*sigma.alpha%*%crossprod(myxspline,(myz-myxrandom%*%matrix(beta[iter-1,],ncol=1)))
        
      } #end else
      
      # update alpha
      alpha[iter, ,k] <-rmvn(1, mean.alpha, sigma.alpha)
      
    }# end k
    
    #---------------------------#  
    # update beta and tau_e     #
    #---------------------------#
    res_sum<-0
    for (i in 1: N_grp){
      sigma.random<-diag(c(1/tau_u[iter-1], rep(1/tau_w[iter-1], N_pt[i])), nrow=1+N_pt[i], ncol=1+N_pt[i])
      idx.temp<-c(i, (sum(N_pt[1:i])+N_grp-(N_pt[i]-1)):(N_grp+sum(N_pt[1:i])))
      sigma.beta<-solve(solve(sigma.random)+crossprod(X_random[grp==i,idx.temp])*tau_e[iter-1])
      
      k.temp<-which(ci[iter, ,i]==1)
      mu.beta<-tau_e[iter-1]*sigma.beta%*%crossprod((X_random[grp==i,idx.temp]),(z[grp==i]-X_spline[grp==i,]%*%alpha[iter, ,k.temp]))
      
      beta[iter,idx.temp]<-rmvn(1, mu.beta, sigma.beta)
      
      res_temp<-z[grp==i]-X_spline[grp==i,]%*%alpha[iter, ,k.temp]-X_random[grp==i,idx.temp]%*%beta[iter,idx.temp]
      res_sum<-res_sum+sum(res_temp^2)
      
      # re-organize results
      pi.subgrp[(iter),i]<-pi[iter,k.temp,i]
      alpha.subgrp[(iter),,i]<-alpha[iter, , k.temp]
    }
    
    tau_e[iter]<-rgamma(1, a+0.5*N_tot, scale=(b+0.5*res_sum)^(-1))
    #---------------------------#    
    # update tau_u#
    #---------------------------#
    tau_u[iter]<-rgamma(1, a+0.5*N_grp, scale=(b+0.5*sum(beta[iter, (1:N_grp)]^2))^(-1))
    
    #---------------------------#
    # update tau_w #
    #---------------------------#
    tau_w[iter]<-rgamma(1, a+0.5*N_subj, scale=(b+0.5*sum(beta[iter, ((N_grp+1):(N_grp+N_subj))]^2))^(-1))
    
    #----------------------------------#
    # update theta_2, lambda, r#
    #----------------------------------#
    cluster.ind=rep(0,N_subj)
    for (i in 1:N_grp){if (ci[iter,1,i]==1){cluster.ind[((i-1)*N_pt[1]+1):(i*N_pt[1])]=rep(1,N_pt[1])} else {cluster.ind[((i-1)*N_pt[1]+1):(i*N_pt[1])]=rep(2,N_pt[1])}}
    
    # log density of theta_2
    logden.theta_2 <- function(theta_2, lambda, r) {
      return(sum(omega_short[cluster.ind==2]*log(exp(-lambda*t_short[cluster.ind==2]^r*exp(theta_2))*(lambda*exp(theta_2)*r*t_short[cluster.ind==2]^(r-1))) + (1-omega_short[cluster.ind==2])*(-lambda*t_short[cluster.ind==2]^r*exp(theta_2))) + dunif(theta_2, -1, 1, log=TRUE))
    }
    # log density of lambda
    logden.lambda <- function(lambda, theta_2, r) {
      return(sum(omega_short[cluster.ind==1]*log(exp(-lambda*t_short[cluster.ind==1]^r*exp(theta_1))*(lambda*exp(theta_1)*r*t_short[cluster.ind==1]^(r-1))) + (1-omega_short[cluster.ind==1])*(-lambda*t_short[cluster.ind==1]^r*exp(theta_1)))
             + sum(omega_short[cluster.ind==2]*log(exp(-lambda*t_short[cluster.ind==2]^r*exp(theta_2))*(lambda*exp(theta_2)*r*t_short[cluster.ind==2]^(r-1))) + (1-omega_short[cluster.ind==2])*(-lambda*t_short[cluster.ind==2]^r*exp(theta_2))) + dgamma(lambda, shape=0.1, rate=0.1, log=TRUE))
    }
    # log density of r
    logden.r <- function(r, theta_2, lambda) {
      return(sum(omega_short[cluster.ind==1]*log(exp(-lambda*t_short[cluster.ind==1]^r*exp(theta_1))*(lambda*exp(theta_1)*r*t_short[cluster.ind==1]^(r-1))) + (1-omega_short[cluster.ind==1])*(-lambda*t_short[cluster.ind==1]^r*exp(theta_1)))
             + sum(omega_short[cluster.ind==2]*log(exp(-lambda*t_short[cluster.ind==2]^r*exp(theta_2))*(lambda*exp(theta_2)*r*t_short[cluster.ind==2]^(r-1))) + (1-omega_short[cluster.ind==2])*(-lambda*t_short[cluster.ind==2]^r*exp(theta_2))) + dgamma(r, shape=0.1, rate=0.1, log=TRUE))
    }
    
    # update theta_2
    theta_2.proposal = function(x){return(rnorm(1,x,0.01))}
    theta_2.proposed = theta_2.proposal(theta_2[iter-1])
    prob.theta_2 = exp(logden.theta_2(theta_2.proposed, lambda[iter-1], r[iter-1])-logden.theta_2(theta_2[iter-1], lambda[iter-1], r[iter-1]))
    if (runif(1) < prob.theta_2){
      theta_2[iter] = theta_2.proposed
    }else{
      theta_2[iter] = theta_2[iter-1]
    }
    # update lambda
    lambda.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
    lambda.proposed = lambda.proposal(lambda[iter-1])
    prob.lambda = exp(logden.lambda(lambda.proposed, theta_2[iter], r[iter-1])-logden.lambda(lambda[iter-1], theta_2[iter], r[iter-1]))
    if (runif(1) < prob.lambda){
      lambda[iter] = lambda.proposed
    }else{
      lambda[iter] = lambda[iter-1]
    }
    # update r
    r.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
    r.proposed = r.proposal(r[iter-1])
    prob.r = exp(logden.r(r.proposed, theta_2[iter], lambda[iter])-logden.r(r[iter-1], theta_2[iter], lambda[iter]))
    if (runif(1) < prob.r){
      r[iter] = r.proposed
    }else{
      r[iter] = r[iter-1]
    }
  }
  
  # Qf = 0.05, Q = 0.8
  # check
  print("Proposed Model Interim 2")
  for (i in 1:length(unique(longdata[,1]))){print(mean(ci[(nburnin+1):mcsim,1,i]))}
  print(sum(((1/lambda)^(1/r)*(log(2))^(1/r))[(nburnin+1):mcsim]>0.6)/(mcsim-nburnin))
  print(sum(((1/lambda/exp(theta_2))^(1/r)*(log(2))^(1/r))[(nburnin+1):mcsim]>0.6)/(mcsim-nburnin))
}










#-------------------------------------# 
#-------------------------------------# 
#-------------------------------------# 
#Part 4: Independent Approach Interim 1# 
#-------------------------------------# 
#-------------------------------------# 
#-------------------------------------# 
longdata <- longdata_full[Subject<=ni*Group-20,]

nknot=2
mcsim=8000
nburnin=4000
#-------------------------------------#


#-------------------------------------#
# A. Longitudinal part data organize  #
#-------------------------------------#
data<-longdata

#---------------------------------------------#
# I. Spline part creation                     #
#---------------------------------------------#

#-------------------------------------------------------------#
# II. extract dimension parameters from the data:             #
# number of observations of each patient (N_obs);             #
# number of total observations(N_tot)= N_subj*length(t);      #
# number of subjects in the data (N_subj);                    #
# number of groups in the data (N_grp);                       #
# number of patients per group (N_pt);                        #
#-------------------------------------------------------------#
# rep(number of time points, number of individuals)
N_obs<-as.numeric(as.vector(table(data[,2])))
# total number of observations
N_tot<-sum(N_obs)
# number of patients in the data (N_subj)
N_subj<-as.numeric(length(unique(data[,2])))
# number of groups in the data (N_grp)
N_grp<-as.numeric(length(unique(data[,1])))
N_grp_y<-N_grp
# number of patients per group as a vector
N_pt<-rep(N_subj/N_grp,N_grp)
# number of patients per group
N_pts<-N_subj/N_grp
# subject ID indicator (c(1,1,1,1,1,2,2,2,2,2,.....))
sub<-as.numeric(match(data[,2], unique(data[,2])))
# group ID indicator
grp<-as.numeric(match(data[,1], unique(data[,1])))
grp_short<-c()
for (i in 1:N_subj){
  grp_short=c(grp_short,grp[1+(i-1)*(N_tot/N_subj)])
}

#------------------------------------------------#
# III. Vectors to store the posterior samples    #
# of parameters and initialize them:             #
# alpha: all the spline coefficients             #
# beta: all the r.e., v_i and w_ij               #
# tau_u: precision of v_i, group-specific r.e.   #
# tau_w: precision of w_ij, subject-specfic r.e. #
# tau_e: precision of e_ijt, measurement error   #
#------------------------------------------------#

tau_u<-tau_w<-tau_e<-rep(NA,mcsim)
r_1<-r_2<-r_3<-r_4<-r_5<-r_6<-r_7<-r_8<-r_9<-r_10<-r_11<-r_12<-rep(NA,mcsim)
theta_1<-theta_2<-theta_3<-theta_4<-theta_5<-theta_6<-theta_7<-theta_8<-theta_9<-theta_10<-theta_11<-theta_12<-rep(NA,mcsim)
# Initialize
tau_u[1]<-tau_w[1]<-tau_e[1]<- 1
lambda<- 1.96
r_1[1]<-r_2[1]<-r_3[1]<-r_4[1]<-r_5[1]<-r_6[1]<-r_7[1]<-r_8[1]<-r_9[1]<-r_10[1]<-r_11[1]<-r_12[1]<- 1.05
theta_1[1]<-theta_2[1]<-theta_3[1]<-theta_4[1]<-theta_5[1]<-theta_6[1]<-theta_7[1]<-theta_8[1]<-theta_9[1]<-theta_10[1]<-theta_11[1]<-theta_12[1]<- -0.75

#-------------------------------------#
# B. response part data organize        #
#-------------------------------------#  
# biomarker measurements
z<-as.numeric(data[,4]) # dependent variable, vector format
# mean biomarker measurement for each patient
z_mean<-c()
for (i in 1:N_subj){
  z_mean=c(z_mean,mean(z[which(sub==i)]))
}

# censoring indicator variable
omega<-as.numeric(data[,5])
omega_short<-c()
for (i in 1:N_subj){
  omega_short=c(omega_short,omega[1+(i-1)*(N_tot/N_subj)])
}

# observed time variable
t<-as.numeric(data[,6])
# a vector of t which does not contain replicates for each patient
t_short<-c()
for (i in 1:N_subj){
  t_short=c(t_short,t[1+(i-1)*(N_tot/N_subj)])
}


#------------------------------------------------------------------#  


#----------------------------------------------#
# D. MCMC part                                 #
#                                              #
#----------------------------------------------#  
set.seed(seed)
for (iter in 2:mcsim){
  
  #----------------------------------#
  # update theta and r#
  #----------------------------------#
  # log density of theta
  logden.theta <- function(theta, lambda, r, i) {
    return(sum(omega_short[grp_short==i]*log(exp(-lambda*t_short[grp_short==i]^r*exp(theta))*(lambda*exp(theta)*r*t_short[grp_short==i]^(r-1))) + (1-omega_short[grp_short==i])*(-lambda*t_short[grp_short==i]^r*exp(theta)))
           + dunif(theta, -2, 0, log=TRUE))
  }
  
  # log density of r
  logden.r <- function(r, theta, lambda, i) {
    return(sum(omega_short[grp_short==i]*log(exp(-lambda*t_short[grp_short==i]^r*exp(theta))*(lambda*exp(theta)*r*t_short[grp_short==i]^(r-1))) + (1-omega_short[grp_short==i])*(-lambda*t_short[grp_short==i]^r*exp(theta)))
           + dgamma(r, shape=0.1, rate=0.1, log=TRUE))
  }
  # update theta_1
  theta_1.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_1.proposed = theta_1.proposal(theta_1[iter-1])
  prob.theta_1 = exp(logden.theta(theta_1.proposed, lambda, r_1[iter-1], 1)-logden.theta(theta_1[iter-1], lambda, r_1[iter-1], 1))
  if (runif(1) < prob.theta_1){
    theta_1[iter] = theta_1.proposed
  }else{
    theta_1[iter] = theta_1[iter-1]
  }
  # update theta_2
  theta_2.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_2.proposed = theta_2.proposal(theta_2[iter-1])
  prob.theta_2 = exp(logden.theta(theta_2.proposed, lambda, r_2[iter-1], 2)-logden.theta(theta_2[iter-1], lambda, r_2[iter-1], 2))
  if (runif(1) < prob.theta_2){
    theta_2[iter] = theta_2.proposed
  }else{
    theta_2[iter] = theta_2[iter-1]
  }
  # update theta_3
  theta_3.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_3.proposed = theta_3.proposal(theta_3[iter-1])
  prob.theta_3 = exp(logden.theta(theta_3.proposed, lambda, r_3[iter-1], 3)-logden.theta(theta_3[iter-1], lambda, r_3[iter-1], 3))
  if (runif(1) < prob.theta_3){
    theta_3[iter] = theta_3.proposed
  }else{
    theta_3[iter] = theta_3[iter-1]
  }
  # update theta_4
  theta_4.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_4.proposed = theta_4.proposal(theta_4[iter-1])
  prob.theta_4 = exp(logden.theta(theta_4.proposed, lambda, r_4[iter-1], 4)-logden.theta(theta_4[iter-1], lambda, r_4[iter-1], 4))
  if (runif(1) < prob.theta_4){
    theta_4[iter] = theta_4.proposed
  }else{
    theta_4[iter] = theta_4[iter-1]
  }
  # update theta_5
  theta_5.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_5.proposed = theta_5.proposal(theta_5[iter-1])
  prob.theta_5 = exp(logden.theta(theta_5.proposed, lambda, r_5[iter-1], 5)-logden.theta(theta_5[iter-1], lambda, r_5[iter-1], 5))
  if (runif(1) < prob.theta_5){
    theta_5[iter] = theta_5.proposed
  }else{
    theta_5[iter] = theta_5[iter-1]
  }
  # update theta_6
  theta_6.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_6.proposed = theta_6.proposal(theta_6[iter-1])
  prob.theta_6 = exp(logden.theta(theta_6.proposed, lambda, r_6[iter-1], 6)-logden.theta(theta_6[iter-1], lambda, r_6[iter-1], 6))
  if (runif(1) < prob.theta_6){
    theta_6[iter] = theta_6.proposed
  }else{
    theta_6[iter] = theta_6[iter-1]
  }
  # update theta_7
  theta_7.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_7.proposed = theta_7.proposal(theta_7[iter-1])
  prob.theta_7 = exp(logden.theta(theta_7.proposed, lambda, r_7[iter-1], 7)-logden.theta(theta_7[iter-1], lambda, r_7[iter-1], 7))
  if (runif(1) < prob.theta_7){
    theta_7[iter] = theta_7.proposed
  }else{
    theta_7[iter] = theta_7[iter-1]
  }
  # update theta_8
  theta_8.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_8.proposed = theta_8.proposal(theta_8[iter-1])
  prob.theta_8 = exp(logden.theta(theta_8.proposed, lambda, r_8[iter-1], 8)-logden.theta(theta_8[iter-1], lambda, r_8[iter-1], 8))
  if (runif(1) < prob.theta_8){
    theta_8[iter] = theta_8.proposed
  }else{
    theta_8[iter] = theta_8[iter-1]
  }
  # update theta_9
  theta_9.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_9.proposed = theta_9.proposal(theta_9[iter-1])
  prob.theta_9 = exp(logden.theta(theta_9.proposed, lambda, r_9[iter-1], 9)-logden.theta(theta_9[iter-1], lambda, r_9[iter-1], 9))
  if (runif(1) < prob.theta_9){
    theta_9[iter] = theta_9.proposed
  }else{
    theta_9[iter] = theta_9[iter-1]
  }
  # update theta_10
  theta_10.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_10.proposed = theta_10.proposal(theta_10[iter-1])
  prob.theta_10 = exp(logden.theta(theta_10.proposed, lambda, r_10[iter-1], 10)-logden.theta(theta_10[iter-1], lambda, r_10[iter-1], 10))
  if (runif(1) < prob.theta_10){
    theta_10[iter] = theta_10.proposed
  }else{
    theta_10[iter] = theta_10[iter-1]
  }
  # update theta_11
  theta_11.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_11.proposed = theta_11.proposal(theta_11[iter-1])
  prob.theta_11 = exp(logden.theta(theta_11.proposed, lambda, r_11[iter-1], 11)-logden.theta(theta_11[iter-1], lambda, r_11[iter-1], 11))
  if (runif(1) < prob.theta_11){
    theta_11[iter] = theta_11.proposed
  }else{
    theta_11[iter] = theta_11[iter-1]
  }
  # update theta_12
  theta_12.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_12.proposed = theta_12.proposal(theta_12[iter-1])
  prob.theta_12 = exp(logden.theta(theta_12.proposed, lambda, r_12[iter-1], 12)-logden.theta(theta_12[iter-1], lambda, r_12[iter-1], 12))
  if (runif(1) < prob.theta_12){
    theta_12[iter] = theta_12.proposed
  }else{
    theta_12[iter] = theta_12[iter-1]
  }
  # update r_1
  r_1.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_1.proposed = r_1.proposal(r_1[iter-1])
  prob.r_1 = exp(logden.r(r_1.proposed, theta_1[iter], lambda, 1)-logden.r(r_1[iter-1], theta_1[iter], lambda, 1))
  if (runif(1) < prob.r_1){
    r_1[iter] = r_1.proposed
  }else{
    r_1[iter] = r_1[iter-1]
  }
  # update r_2
  r_2.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_2.proposed = r_2.proposal(r_2[iter-1])
  prob.r_2 = exp(logden.r(r_2.proposed, theta_2[iter], lambda, 2)-logden.r(r_2[iter-1], theta_2[iter], lambda, 2))
  if (runif(1) < prob.r_2){
    r_2[iter] = r_2.proposed
  }else{
    r_2[iter] = r_2[iter-1]
  }
  # update r_3
  r_3.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_3.proposed = r_3.proposal(r_3[iter-1])
  prob.r_3 = exp(logden.r(r_3.proposed, theta_3[iter], lambda, 3)-logden.r(r_3[iter-1], theta_3[iter], lambda, 3))
  if (runif(1) < prob.r_3){
    r_3[iter] = r_3.proposed
  }else{
    r_3[iter] = r_3[iter-1]
  }
  # update r_4
  r_4.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_4.proposed = r_4.proposal(r_4[iter-1])
  prob.r_4 = exp(logden.r(r_4.proposed, theta_4[iter], lambda, 4)-logden.r(r_4[iter-1], theta_4[iter], lambda, 4))
  if (runif(1) < prob.r_4){
    r_4[iter] = r_4.proposed
  }else{
    r_4[iter] = r_4[iter-1]
  }
  # update r_5
  r_5.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_5.proposed = r_5.proposal(r_5[iter-1])
  prob.r_5 = exp(logden.r(r_5.proposed, theta_5[iter], lambda, 5)-logden.r(r_5[iter-1], theta_5[iter], lambda, 5))
  if (runif(1) < prob.r_5){
    r_5[iter] = r_5.proposed
  }else{
    r_5[iter] = r_5[iter-1]
  }
  # update r_6
  r_6.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_6.proposed = r_6.proposal(r_6[iter-1])
  prob.r_6 = exp(logden.r(r_6.proposed, theta_6[iter], lambda, 6)-logden.r(r_6[iter-1], theta_6[iter], lambda, 6))
  if (runif(1) < prob.r_6){
    r_6[iter] = r_6.proposed
  }else{
    r_6[iter] = r_6[iter-1]
  }
  # update r_7
  r_7.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_7.proposed = r_7.proposal(r_7[iter-1])
  prob.r_7 = exp(logden.r(r_7.proposed, theta_7[iter], lambda, 7)-logden.r(r_7[iter-1], theta_7[iter], lambda, 7))
  if (runif(1) < prob.r_7){
    r_7[iter] = r_7.proposed
  }else{
    r_7[iter] = r_7[iter-1]
  }
  # update r_8
  r_8.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_8.proposed = r_8.proposal(r_8[iter-1])
  prob.r_8 = exp(logden.r(r_8.proposed, theta_8[iter], lambda, 8)-logden.r(r_8[iter-1], theta_8[iter], lambda, 8))
  if (runif(1) < prob.r_8){
    r_8[iter] = r_8.proposed
  }else{
    r_8[iter] = r_8[iter-1]
  }
  # update r_9
  r_9.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_9.proposed = r_9.proposal(r_9[iter-1])
  prob.r_9 = exp(logden.r(r_9.proposed, theta_9[iter], lambda, 9)-logden.r(r_9[iter-1], theta_9[iter], lambda, 9))
  if (runif(1) < prob.r_9){
    r_9[iter] = r_9.proposed
  }else{
    r_9[iter] = r_9[iter-1]
  }
  # update r_10
  r_10.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_10.proposed = r_10.proposal(r_10[iter-1])
  prob.r_10 = exp(logden.r(r_10.proposed, theta_10[iter], lambda, 10)-logden.r(r_10[iter-1], theta_10[iter], lambda, 10))
  if (runif(1) < prob.r_10){
    r_10[iter] = r_10.proposed
  }else{
    r_10[iter] = r_10[iter-1]
  }
  # update r_11
  r_11.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_11.proposed = r_11.proposal(r_11[iter-1])
  prob.r_11 = exp(logden.r(r_11.proposed, theta_11[iter], lambda, 11)-logden.r(r_11[iter-1], theta_11[iter], lambda, 11))
  if (runif(1) < prob.r_11){
    r_11[iter] = r_11.proposed
  }else{
    r_11[iter] = r_11[iter-1]
  }
  # update r_12
  r_12.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_12.proposed = r_12.proposal(r_12[iter-1])
  prob.r_12 = exp(logden.r(r_12.proposed, theta_12[iter], lambda, 12)-logden.r(r_12[iter-1], theta_12[iter], lambda, 12))
  if (runif(1) < prob.r_12){
    r_12[iter] = r_12.proposed
  }else{
    r_12[iter] = r_12[iter-1]
  }
}

# check
print("Independent Approach Interim 1")
mean(((1/lambda/exp(theta_1))^(1/r_1)*(log(2))^(1/r_1))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_2))^(1/r_2)*(log(2))^(1/r_2))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_3))^(1/r_3)*(log(2))^(1/r_3))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_4))^(1/r_4)*(log(2))^(1/r_4))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_5))^(1/r_5)*(log(2))^(1/r_5))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_6))^(1/r_6)*(log(2))^(1/r_6))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_7))^(1/r_7)*(log(2))^(1/r_7))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_8))^(1/r_8)*(log(2))^(1/r_8))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_9))^(1/r_9)*(log(2))^(1/r_9))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_10))^(1/r_10)*(log(2))^(1/r_10))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_11))^(1/r_11)*(log(2))^(1/r_11))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_12))^(1/r_12)*(log(2))^(1/r_12))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
longdata <- longdata_full
if (mean(((1/lambda/exp(theta_1))^(1/r_1)*(log(2))^(1/r_1))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)<0.05){longdata <- longdata[(longdata[,1] !=1),]}
if (mean(((1/lambda/exp(theta_2))^(1/r_2)*(log(2))^(1/r_2))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)<0.05){longdata <- longdata[(longdata[,1] !=2),]}
if (mean(((1/lambda/exp(theta_3))^(1/r_3)*(log(2))^(1/r_3))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)<0.05){longdata <- longdata[(longdata[,1] !=3),]}
if (mean(((1/lambda/exp(theta_4))^(1/r_4)*(log(2))^(1/r_4))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)<0.05){longdata <- longdata[(longdata[,1] !=4),]}
if (mean(((1/lambda/exp(theta_5))^(1/r_5)*(log(2))^(1/r_5))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)<0.05){longdata <- longdata[(longdata[,1] !=5),]}
if (mean(((1/lambda/exp(theta_6))^(1/r_6)*(log(2))^(1/r_6))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)<0.05){longdata <- longdata[(longdata[,1] !=6),]}
if (mean(((1/lambda/exp(theta_7))^(1/r_7)*(log(2))^(1/r_7))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)<0.05){longdata <- longdata[(longdata[,1] !=7),]}
if (mean(((1/lambda/exp(theta_8))^(1/r_8)*(log(2))^(1/r_8))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)<0.05){longdata <- longdata[(longdata[,1] !=8),]}
if (mean(((1/lambda/exp(theta_9))^(1/r_9)*(log(2))^(1/r_9))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)<0.05){longdata <- longdata[(longdata[,1] !=9),]}
if (mean(((1/lambda/exp(theta_10))^(1/r_10)*(log(2))^(1/r_10))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)<0.05){longdata <- longdata[(longdata[,1] !=10),]}
if (mean(((1/lambda/exp(theta_11))^(1/r_11)*(log(2))^(1/r_11))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)<0.05){longdata <- longdata[(longdata[,1] !=11),]}
if (mean(((1/lambda/exp(theta_12))^(1/r_12)*(log(2))^(1/r_12))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)<0.05){longdata <- longdata[(longdata[,1] !=12),]}










#-------------------------------------# 
#-------------------------------------# 
#-------------------------------------# 
#Part 5: Independent Approach Interim 2# 
#-------------------------------------# 
#-------------------------------------# 
#-------------------------------------# 
#-------------------------------------# 
nknot=2
mcsim=8000
nburnin=4000
#-------------------------------------#


#-------------------------------------#
# A. Longitudinal part data organize  #
#-------------------------------------#
data<-longdata

#---------------------------------------------#
# I. Spline part creation                     #
#---------------------------------------------#

#-------------------------------------------------------------#
# II. extract dimension parameters from the data:             #
# number of observations of each patient (N_obs);             #
# number of total observations(N_tot)= N_subj*length(t);      #
# number of subjects in the data (N_subj);                    #
# number of groups in the data (N_grp);                       #
# number of patients per group (N_pt);                        #
#-------------------------------------------------------------#
# rep(number of time points, number of individuals)
N_obs<-as.numeric(as.vector(table(data[,2])))
# total number of observations
N_tot<-sum(N_obs)
# number of patients in the data (N_subj)
N_subj<-as.numeric(length(unique(data[,2])))
# number of groups in the data (N_grp)
N_grp<-as.numeric(length(unique(data[,1])))
N_grp_y<-N_grp
# number of patients per group as a vector
N_pt<-rep(N_subj/N_grp,N_grp)
# number of patients per group
N_pts<-N_subj/N_grp
# subject ID indicator (c(1,1,1,1,1,2,2,2,2,2,.....))
sub<-as.numeric(match(data[,2], unique(data[,2])))
# group ID indicator
grp<-as.numeric(match(data[,1], unique(data[,1])))
grp_short<-c()
for (i in 1:N_subj){
  grp_short=c(grp_short,grp[1+(i-1)*(N_tot/N_subj)])
}

#------------------------------------------------#
# III. Vectors to store the posterior samples    #
# of parameters and initialize them:             #
# alpha: all the spline coefficients             #
# beta: all the r.e., v_i and w_ij               #
# tau_u: precision of v_i, group-specific r.e.   #
# tau_w: precision of w_ij, subject-specfic r.e. #
# tau_e: precision of e_ijt, measurement error   #
#------------------------------------------------#

tau_u<-tau_w<-tau_e<-rep(NA,mcsim)
r_1<-r_2<-r_3<-r_4<-r_5<-r_6<-r_7<-r_8<-r_9<-r_10<-r_11<-r_12<-rep(NA,mcsim)
theta_1<-theta_2<-theta_3<-theta_4<-theta_5<-theta_6<-theta_7<-theta_8<-theta_9<-theta_10<-theta_11<-theta_12<-rep(NA,mcsim)
# Initialize
tau_u[1]<-tau_w[1]<-tau_e[1]<- 1
lambda<- 1.96
r_1[1]<-r_2[1]<-r_3[1]<-r_4[1]<-r_5[1]<-r_6[1]<-r_7[1]<-r_8[1]<-r_9[1]<-r_10[1]<-r_11[1]<-r_12[1]<- 1.05
theta_1[1]<-theta_2[1]<-theta_3[1]<-theta_4[1]<-theta_5[1]<-theta_6[1]<-theta_7[1]<-theta_8[1]<-theta_9[1]<-theta_10[1]<-theta_11[1]<-theta_12[1]<- -0.75

#-------------------------------------#
# B. response part data organize        #
#-------------------------------------#  
# biomarker measurements
z<-as.numeric(data[,4]) # dependent variable, vector format
# mean biomarker measurement for each patient
z_mean<-c()
for (i in 1:N_subj){
  z_mean=c(z_mean,mean(z[which(sub==i)]))
}

# censoring indicator variable
omega<-as.numeric(data[,5])
omega_short<-c()
for (i in 1:N_subj){
  omega_short=c(omega_short,omega[1+(i-1)*(N_tot/N_subj)])
}

# observed time variable
t<-as.numeric(data[,6])
# a vector of t which does not contain replicates for each patient
t_short<-c()
for (i in 1:N_subj){
  t_short=c(t_short,t[1+(i-1)*(N_tot/N_subj)])
}


#------------------------------------------------------------------#  


#----------------------------------------------#
# D. MCMC part                                 #
#                                              #
#----------------------------------------------#  
set.seed(seed)
for (iter in 2:mcsim){
  
  #----------------------------------#
  # update theta and r#
  #----------------------------------#
  # log density of theta
  logden.theta <- function(theta, lambda, r, i) {
    return(sum(omega_short[grp_short==i]*log(exp(-lambda*t_short[grp_short==i]^r*exp(theta))*(lambda*exp(theta)*r*t_short[grp_short==i]^(r-1))) + (1-omega_short[grp_short==i])*(-lambda*t_short[grp_short==i]^r*exp(theta)))
           + dunif(theta, -2, 0, log=TRUE))
  }
  
  # log density of r
  logden.r <- function(r, theta, lambda, i) {
    return(sum(omega_short[grp_short==i]*log(exp(-lambda*t_short[grp_short==i]^r*exp(theta))*(lambda*exp(theta)*r*t_short[grp_short==i]^(r-1))) + (1-omega_short[grp_short==i])*(-lambda*t_short[grp_short==i]^r*exp(theta)))
           + dgamma(r, shape=0.1, rate=0.1, log=TRUE))
  }
  # update theta_1
  theta_1.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_1.proposed = theta_1.proposal(theta_1[iter-1])
  prob.theta_1 = exp(logden.theta(theta_1.proposed, lambda, r_1[iter-1], 1)-logden.theta(theta_1[iter-1], lambda, r_1[iter-1], 1))
  if (runif(1) < prob.theta_1){
    theta_1[iter] = theta_1.proposed
  }else{
    theta_1[iter] = theta_1[iter-1]
  }
  # update theta_2
  theta_2.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_2.proposed = theta_2.proposal(theta_2[iter-1])
  prob.theta_2 = exp(logden.theta(theta_2.proposed, lambda, r_2[iter-1], 2)-logden.theta(theta_2[iter-1], lambda, r_2[iter-1], 2))
  if (runif(1) < prob.theta_2){
    theta_2[iter] = theta_2.proposed
  }else{
    theta_2[iter] = theta_2[iter-1]
  }
  # update theta_3
  theta_3.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_3.proposed = theta_3.proposal(theta_3[iter-1])
  prob.theta_3 = exp(logden.theta(theta_3.proposed, lambda, r_3[iter-1], 3)-logden.theta(theta_3[iter-1], lambda, r_3[iter-1], 3))
  if (runif(1) < prob.theta_3){
    theta_3[iter] = theta_3.proposed
  }else{
    theta_3[iter] = theta_3[iter-1]
  }
  # update theta_4
  theta_4.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_4.proposed = theta_4.proposal(theta_4[iter-1])
  prob.theta_4 = exp(logden.theta(theta_4.proposed, lambda, r_4[iter-1], 4)-logden.theta(theta_4[iter-1], lambda, r_4[iter-1], 4))
  if (runif(1) < prob.theta_4){
    theta_4[iter] = theta_4.proposed
  }else{
    theta_4[iter] = theta_4[iter-1]
  }
  # update theta_5
  theta_5.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_5.proposed = theta_5.proposal(theta_5[iter-1])
  prob.theta_5 = exp(logden.theta(theta_5.proposed, lambda, r_5[iter-1], 5)-logden.theta(theta_5[iter-1], lambda, r_5[iter-1], 5))
  if (runif(1) < prob.theta_5){
    theta_5[iter] = theta_5.proposed
  }else{
    theta_5[iter] = theta_5[iter-1]
  }
  # update theta_6
  theta_6.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_6.proposed = theta_6.proposal(theta_6[iter-1])
  prob.theta_6 = exp(logden.theta(theta_6.proposed, lambda, r_6[iter-1], 6)-logden.theta(theta_6[iter-1], lambda, r_6[iter-1], 6))
  if (runif(1) < prob.theta_6){
    theta_6[iter] = theta_6.proposed
  }else{
    theta_6[iter] = theta_6[iter-1]
  }
  # update theta_7
  theta_7.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_7.proposed = theta_7.proposal(theta_7[iter-1])
  prob.theta_7 = exp(logden.theta(theta_7.proposed, lambda, r_7[iter-1], 7)-logden.theta(theta_7[iter-1], lambda, r_7[iter-1], 7))
  if (runif(1) < prob.theta_7){
    theta_7[iter] = theta_7.proposed
  }else{
    theta_7[iter] = theta_7[iter-1]
  }
  # update theta_8
  theta_8.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_8.proposed = theta_8.proposal(theta_8[iter-1])
  prob.theta_8 = exp(logden.theta(theta_8.proposed, lambda, r_8[iter-1], 8)-logden.theta(theta_8[iter-1], lambda, r_8[iter-1], 8))
  if (runif(1) < prob.theta_8){
    theta_8[iter] = theta_8.proposed
  }else{
    theta_8[iter] = theta_8[iter-1]
  }
  # update theta_9
  theta_9.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_9.proposed = theta_9.proposal(theta_9[iter-1])
  prob.theta_9 = exp(logden.theta(theta_9.proposed, lambda, r_9[iter-1], 9)-logden.theta(theta_9[iter-1], lambda, r_9[iter-1], 9))
  if (runif(1) < prob.theta_9){
    theta_9[iter] = theta_9.proposed
  }else{
    theta_9[iter] = theta_9[iter-1]
  }
  # update theta_10
  theta_10.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_10.proposed = theta_10.proposal(theta_10[iter-1])
  prob.theta_10 = exp(logden.theta(theta_10.proposed, lambda, r_10[iter-1], 10)-logden.theta(theta_10[iter-1], lambda, r_10[iter-1], 10))
  if (runif(1) < prob.theta_10){
    theta_10[iter] = theta_10.proposed
  }else{
    theta_10[iter] = theta_10[iter-1]
  }
  # update theta_11
  theta_11.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_11.proposed = theta_11.proposal(theta_11[iter-1])
  prob.theta_11 = exp(logden.theta(theta_11.proposed, lambda, r_11[iter-1], 11)-logden.theta(theta_11[iter-1], lambda, r_11[iter-1], 11))
  if (runif(1) < prob.theta_11){
    theta_11[iter] = theta_11.proposed
  }else{
    theta_11[iter] = theta_11[iter-1]
  }
  # update theta_12
  theta_12.proposal = function(x){return(rnorm(1,x,0.01))}
  theta_12.proposed = theta_12.proposal(theta_12[iter-1])
  prob.theta_12 = exp(logden.theta(theta_12.proposed, lambda, r_12[iter-1], 12)-logden.theta(theta_12[iter-1], lambda, r_12[iter-1], 12))
  if (runif(1) < prob.theta_12){
    theta_12[iter] = theta_12.proposed
  }else{
    theta_12[iter] = theta_12[iter-1]
  }
  # update r_1
  r_1.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_1.proposed = r_1.proposal(r_1[iter-1])
  prob.r_1 = exp(logden.r(r_1.proposed, theta_1[iter], lambda, 1)-logden.r(r_1[iter-1], theta_1[iter], lambda, 1))
  if (runif(1) < prob.r_1){
    r_1[iter] = r_1.proposed
  }else{
    r_1[iter] = r_1[iter-1]
  }
  # update r_2
  r_2.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_2.proposed = r_2.proposal(r_2[iter-1])
  prob.r_2 = exp(logden.r(r_2.proposed, theta_2[iter], lambda, 2)-logden.r(r_2[iter-1], theta_2[iter], lambda, 2))
  if (runif(1) < prob.r_2){
    r_2[iter] = r_2.proposed
  }else{
    r_2[iter] = r_2[iter-1]
  }
  # update r_3
  r_3.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_3.proposed = r_3.proposal(r_3[iter-1])
  prob.r_3 = exp(logden.r(r_3.proposed, theta_3[iter], lambda, 3)-logden.r(r_3[iter-1], theta_3[iter], lambda, 3))
  if (runif(1) < prob.r_3){
    r_3[iter] = r_3.proposed
  }else{
    r_3[iter] = r_3[iter-1]
  }
  # update r_4
  r_4.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_4.proposed = r_4.proposal(r_4[iter-1])
  prob.r_4 = exp(logden.r(r_4.proposed, theta_4[iter], lambda, 4)-logden.r(r_4[iter-1], theta_4[iter], lambda, 4))
  if (runif(1) < prob.r_4){
    r_4[iter] = r_4.proposed
  }else{
    r_4[iter] = r_4[iter-1]
  }
  # update r_5
  r_5.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_5.proposed = r_5.proposal(r_5[iter-1])
  prob.r_5 = exp(logden.r(r_5.proposed, theta_5[iter], lambda, 5)-logden.r(r_5[iter-1], theta_5[iter], lambda, 5))
  if (runif(1) < prob.r_5){
    r_5[iter] = r_5.proposed
  }else{
    r_5[iter] = r_5[iter-1]
  }
  # update r_6
  r_6.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_6.proposed = r_6.proposal(r_6[iter-1])
  prob.r_6 = exp(logden.r(r_6.proposed, theta_6[iter], lambda, 6)-logden.r(r_6[iter-1], theta_6[iter], lambda, 6))
  if (runif(1) < prob.r_6){
    r_6[iter] = r_6.proposed
  }else{
    r_6[iter] = r_6[iter-1]
  }
  # update r_7
  r_7.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_7.proposed = r_7.proposal(r_7[iter-1])
  prob.r_7 = exp(logden.r(r_7.proposed, theta_7[iter], lambda, 7)-logden.r(r_7[iter-1], theta_7[iter], lambda, 7))
  if (runif(1) < prob.r_7){
    r_7[iter] = r_7.proposed
  }else{
    r_7[iter] = r_7[iter-1]
  }
  # update r_8
  r_8.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_8.proposed = r_8.proposal(r_8[iter-1])
  prob.r_8 = exp(logden.r(r_8.proposed, theta_8[iter], lambda, 8)-logden.r(r_8[iter-1], theta_8[iter], lambda, 8))
  if (runif(1) < prob.r_8){
    r_8[iter] = r_8.proposed
  }else{
    r_8[iter] = r_8[iter-1]
  }
  # update r_9
  r_9.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_9.proposed = r_9.proposal(r_9[iter-1])
  prob.r_9 = exp(logden.r(r_9.proposed, theta_9[iter], lambda, 9)-logden.r(r_9[iter-1], theta_9[iter], lambda, 9))
  if (runif(1) < prob.r_9){
    r_9[iter] = r_9.proposed
  }else{
    r_9[iter] = r_9[iter-1]
  }
  # update r_10
  r_10.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_10.proposed = r_10.proposal(r_10[iter-1])
  prob.r_10 = exp(logden.r(r_10.proposed, theta_10[iter], lambda, 10)-logden.r(r_10[iter-1], theta_10[iter], lambda, 10))
  if (runif(1) < prob.r_10){
    r_10[iter] = r_10.proposed
  }else{
    r_10[iter] = r_10[iter-1]
  }
  # update r_11
  r_11.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_11.proposed = r_11.proposal(r_11[iter-1])
  prob.r_11 = exp(logden.r(r_11.proposed, theta_11[iter], lambda, 11)-logden.r(r_11[iter-1], theta_11[iter], lambda, 11))
  if (runif(1) < prob.r_11){
    r_11[iter] = r_11.proposed
  }else{
    r_11[iter] = r_11[iter-1]
  }
  # update r_12
  r_12.proposal = function(x){return(runif(1,max(0,x-0.1),x+0.1))}
  r_12.proposed = r_12.proposal(r_12[iter-1])
  prob.r_12 = exp(logden.r(r_12.proposed, theta_12[iter], lambda, 12)-logden.r(r_12[iter-1], theta_12[iter], lambda, 12))
  if (runif(1) < prob.r_12){
    r_12[iter] = r_12.proposed
  }else{
    r_12[iter] = r_12[iter-1]
  }
}

# check
print("Independent Approach Interim 2")
mean(((1/lambda/exp(theta_1))^(1/r_1)*(log(2))^(1/r_1))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_2))^(1/r_2)*(log(2))^(1/r_2))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_3))^(1/r_3)*(log(2))^(1/r_3))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_4))^(1/r_4)*(log(2))^(1/r_4))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_5))^(1/r_5)*(log(2))^(1/r_5))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_6))^(1/r_6)*(log(2))^(1/r_6))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_7))^(1/r_7)*(log(2))^(1/r_7))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_8))^(1/r_8)*(log(2))^(1/r_8))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_9))^(1/r_9)*(log(2))^(1/r_9))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_10))^(1/r_10)*(log(2))^(1/r_10))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_11))^(1/r_11)*(log(2))^(1/r_11))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)
mean(((1/lambda/exp(theta_12))^(1/r_12)*(log(2))^(1/r_12))[(nburnin+1):mcsim]>0.6, na.rm=TRUE)