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
K=3

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
ci.true=t(matrix(c(rep(c(1,0,0),8), rep(c(0,1,0),2), rep(c(0,0,1),2)),ncol=3,byrow = TRUE))
muktl=rep(0,I*ni*T)
for (i in 1:I){if(ci.true[1,i]==1){muktl[((i-1)*ni*T+1):(i*ni*T)]=X_spline[((i-1)*ni*T+1):(i*ni*T),]%*%alpha1} else {muktl[((i-1)*ni*T+1):(i*ni*T)]=X_spline[((i-1)*ni*T+1):(i*ni*T),]%*%alpha2}}

# biomarker measurements
z=muktl + v_long + w_long + epsilon

# survival time
lambda<-1.96
r<-1.5
theta_1<-0
theta_2<- -0.505
theta_1.5<- -0.275

# Weibull latent event times
v <- runif(ni*I)
t=rep(0,ni*I)
for (j in 1:I){
  for (i in ((j-1)*ni+1):(j*ni)){
    if (ci.true[1,j]==1) {t[i] <- (- log(v[i]) / (lambda * exp(theta_1)))^(1 / r)}
    else if (ci.true[2,j]==1) {t[i] <- (- log(v[i]) / (lambda * exp(theta_1.5)))^(1 / r)}
    else {t[i] <- (- log(v[i]) / (lambda * exp(theta_2)))^(1 / r)}}
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
#Part 2: BHM Interim 1# 
#-------------------------------------# 
#-------------------------------------# 
#-------------------------------------# 

nknot=2
mcsim=8000
nburnin=4000
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

lambda<-matrix(NA, nrow=mcsim, ncol=I)
r<-matrix(NA, nrow=mcsim, ncol=I)
theta<-matrix(NA, nrow=mcsim, ncol=I)
hyperparameter_1<-rep(NA,mcsim)
hyperparameter_2<-rep(NA,mcsim)
hyperparameter_3<-rep(NA,mcsim)
hyperparameter_4<-rep(NA,mcsim)
hyperparameter_5<-rep(NA,mcsim)
hyperparameter_6<-rep(NA,mcsim)
# Initialize
lambda[1,]<- rep(1.5, I)
r[1,]<- rep(1.05, I)
theta[1,]<- rep(-0.75, I)
hyperparameter_1[1]<- -1
hyperparameter_2[1]<- 1
hyperparameter_3[1]<- 0.1
hyperparameter_4[1]<- 0.1
hyperparameter_5[1]<- 0.1
hyperparameter_6[1]<- 0.1


#-------------------------------------#
# B. response part data organize      #
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


#----------------------------------------------#
# D. MCMC part                                 #
#                                              #
#----------------------------------------------#    

set.seed(seed)
for (iter in 2:mcsim){
  #----------------------------------#
  # update theta, lambda, r#
  #----------------------------------#
  
  # log density
  logden <- function(theta, lambda, r, hyperparameter_1, hyperparameter_2, hyperparameter_3, hyperparameter_4, hyperparameter_5, hyperparameter_6) {
    logdensity=rep(NA,I)
    for (j in 1:I){logdensity[j] = sum(omega_short[grp_short==j]*log(exp(-lambda[j]*t_short[grp_short==j]^r[j]*exp(theta[j]))*(lambda[j]*exp(theta[j])*r[j]*t_short[grp_short==j]^(r[j]-1))) + (1-omega_short[grp_short==j])*(-lambda[j]*t_short[grp_short==j]^r[j]*exp(theta[j]))) + dnorm(theta[j], hyperparameter_1, hyperparameter_2, log=TRUE) + dgamma(lambda[j], shape=hyperparameter_3, rate=hyperparameter_4, log=TRUE) + dgamma(r[j], shape=hyperparameter_5, rate=hyperparameter_6, log=TRUE)}
    return(sum(logdensity)+dnorm(hyperparameter_1, 0, 10, log=TRUE) + dunif(hyperparameter_2, 0, 10, log=TRUE) + dunif(hyperparameter_3, 0, 1, log=TRUE) + dunif(hyperparameter_4, 0, 1, log=TRUE) + dunif(hyperparameter_5, 0, 1, log=TRUE) + dunif(hyperparameter_6, 0, 1, log=TRUE))
  }
  # update theta
  theta.proposal = function(x){return(rnorm(1,x,0.01))}
  i=1
  {
    theta.proposed = theta.proposal(theta[iter-1,i])
    prob.theta = exp(logden(c(theta.proposed, theta[iter-1,(i+1):I]), lambda[iter-1,], r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                     -logden(c(theta[iter-1,i:I]), lambda[iter-1,], r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
    )
    if (runif(1) < prob.theta){
      theta[iter,i] = theta.proposed
    }else{
      theta[iter,i] = theta[iter-1,i]
    }
  }
  for (i in 2:(I-1)){
    theta.proposed = theta.proposal(theta[iter-1,i])
    prob.theta = exp(logden(c(theta[iter,1:(i-1)], theta.proposed, theta[iter-1,(i+1):I]), lambda[iter-1,], r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                     -logden(c(theta[iter,1:(i-1)], theta[iter-1,i:I]), lambda[iter-1,], r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
    )
    if (runif(1) < prob.theta){
      theta[iter,i] = theta.proposed
    }else{
      theta[iter,i] = theta[iter-1,i]
    }
  }
  i=I
  {
    theta.proposed = theta.proposal(theta[iter-1,i])
    prob.theta = exp(logden(c(theta[iter,1:(i-1)], theta.proposed), lambda[iter-1,], r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                     -logden(c(theta[iter,1:(i-1)], theta[iter-1,i:I]), lambda[iter-1,], r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
    )
    if (runif(1) < prob.theta){
      theta[iter,i] = theta.proposed
    }else{
      theta[iter,i] = theta[iter-1,i]
    }
  }
  # update lambda
  lambda.proposal = function(x){return(max(0.01,rnorm(1,x,0.01)))}
  i=1
  {
    lambda.proposed = lambda.proposal(lambda[iter-1,i])
    prob.lambda = exp(logden(theta[iter,], c(lambda.proposed, lambda[iter-1,(i+1):I]), r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                      -logden(theta[iter,], c(lambda[iter-1,i:I]), r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
    )
    if (runif(1) < prob.lambda){
      lambda[iter,i] = lambda.proposed
    }else{
      lambda[iter,i] = lambda[iter-1,i]
    }
  }
  for (i in 2:(I-1)){
    lambda.proposed = lambda.proposal(lambda[iter-1,i])
    prob.lambda = exp(logden(theta[iter,], c(lambda[iter,1:(i-1)], lambda.proposed, lambda[iter-1,(i+1):I]), r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                      -logden(theta[iter,], c(lambda[iter,1:(i-1)], lambda[iter-1,i:I]), r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
    )
    if (runif(1) < prob.lambda){
      lambda[iter,i] = lambda.proposed
    }else{
      lambda[iter,i] = lambda[iter-1,i]
    }
  }
  i=I
  {
    lambda.proposed = lambda.proposal(lambda[iter-1,i])
    prob.lambda = exp(logden(theta[iter,], c(lambda[iter,1:(i-1)], lambda.proposed), r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                      -logden(theta[iter,], c(lambda[iter,1:(i-1)], lambda[iter-1,i:I]), r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
    )
    if (runif(1) < prob.lambda){
      lambda[iter,i] = lambda.proposed
    }else{
      lambda[iter,i] = lambda[iter-1,i]
    }
  }
  # update r
  r.proposal = function(x){return(max(0.01,rnorm(1,x,0.01)))}
  i=1
  {
    r.proposed = r.proposal(r[iter-1,i])
    prob.r = exp(logden(theta[iter,], lambda[iter,], c(r.proposed, r[iter-1,(i+1):I]), hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                 -logden(theta[iter,], lambda[iter,], c(r[iter-1,i:I]), hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
    )
    if (runif(1) < prob.r){
      r[iter,i] = r.proposed
    }else{
      r[iter,i] = r[iter-1,i]
    }
  }
  for (i in 2:(I-1)){
    r.proposed = r.proposal(r[iter-1,i])
    prob.r = exp(logden(theta[iter,], lambda[iter,], c(r[iter,1:(i-1)], r.proposed, r[iter-1,(i+1):I]), hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                 -logden(theta[iter,], lambda[iter,], c(r[iter,1:(i-1)], r[iter-1,i:I]), hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
    )
    if (runif(1) < prob.r){
      r[iter,i] = r.proposed
    }else{
      r[iter,i] = r[iter-1,i]
    }
  }
  i=I
  {
    r.proposed = r.proposal(r[iter-1,i])
    prob.r = exp(logden(theta[iter,], lambda[iter,], c(r[iter,1:(i-1)], r.proposed),  hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                 -logden(theta[iter,], lambda[iter,], c(r[iter,1:(i-1)], r[iter-1,i:I]), hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
    )
    if (runif(1) < prob.r){
      r[iter,i] = r.proposed
    }else{
      r[iter,i] = r[iter-1,i]
    }
  }
  # update hyperparameter_1
  hyperparameter_1.proposal = function(x){return(rnorm(1,x,0.01))}
  hyperparameter_1.proposed = hyperparameter_1.proposal(hyperparameter_1[iter-1])
  prob.hyperparameter_1 = exp(logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1.proposed, hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                              -logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
  )
  if (runif(1) < prob.hyperparameter_1){
    hyperparameter_1[iter] = hyperparameter_1.proposed
  }else{
    hyperparameter_1[iter] = hyperparameter_1[iter-1]
  }
  # update hyperparameter_2
  hyperparameter_2.proposal = function(x){return(max(0.0001,rnorm(1,x,0.01)))}
  hyperparameter_2.proposed = hyperparameter_2.proposal(hyperparameter_2[iter-1])
  prob.hyperparameter_2 = exp(logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2.proposed, hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                              -logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
  )
  if (runif(1) < prob.hyperparameter_2){
    hyperparameter_2[iter] = hyperparameter_2.proposed
  }else{
    hyperparameter_2[iter] = hyperparameter_2[iter-1]
  }
  # update hyperparameter_3
  hyperparameter_3.proposal = function(x){return(max(0.0001,rnorm(1,x,0.01)))}
  hyperparameter_3.proposed = hyperparameter_3.proposal(hyperparameter_3[iter-1])
  prob.hyperparameter_3 = exp(logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter], hyperparameter_3.proposed, hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                              -logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
  )
  if (runif(1) < prob.hyperparameter_3){
    hyperparameter_3[iter] = hyperparameter_3.proposed
  }else{
    hyperparameter_3[iter] = hyperparameter_3[iter-1]
  }
  # update hyperparameter_4
  hyperparameter_4.proposal = function(x){return(max(0.0001,rnorm(1,x,0.01)))}
  hyperparameter_4.proposed = hyperparameter_4.proposal(hyperparameter_4[iter-1])
  prob.hyperparameter_4 = exp(logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter], hyperparameter_3[iter], hyperparameter_4.proposed, hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                              -logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter], hyperparameter_3[iter], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
  )
  if (runif(1) < prob.hyperparameter_4){
    hyperparameter_4[iter] = hyperparameter_4.proposed
  }else{
    hyperparameter_4[iter] = hyperparameter_4[iter-1]
  }
  # update hyperparameter_5
  hyperparameter_5.proposal = function(x){return(max(0.0001,rnorm(1,x,0.01)))}
  hyperparameter_5.proposed = hyperparameter_5.proposal(hyperparameter_5[iter-1])
  prob.hyperparameter_5 = exp(logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter], hyperparameter_3[iter], hyperparameter_4[iter], hyperparameter_5.proposed, hyperparameter_6[iter-1])
                              -logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter], hyperparameter_3[iter], hyperparameter_4[iter], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
  )
  if (runif(1) < prob.hyperparameter_5){
    hyperparameter_5[iter] = hyperparameter_5.proposed
  }else{
    hyperparameter_5[iter] = hyperparameter_5[iter-1]
  }
  # update hyperparameter_6
  hyperparameter_6.proposal = function(x){return(max(0.0001,rnorm(1,x,0.01)))}
  hyperparameter_6.proposed = hyperparameter_6.proposal(hyperparameter_6[iter-1])
  prob.hyperparameter_6 = exp(logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter], hyperparameter_3[iter], hyperparameter_4[iter], hyperparameter_5[iter], hyperparameter_6.proposed)
                              -logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter], hyperparameter_3[iter], hyperparameter_4[iter], hyperparameter_5[iter], hyperparameter_6[iter-1])
  )
  if (runif(1) < prob.hyperparameter_6){
    hyperparameter_6[iter] = hyperparameter_6.proposed
  }else{
    hyperparameter_6[iter] = hyperparameter_6[iter-1]
  }
}

# Qf = 0.05, Q = 0.8
# check
print("BHM Interim 1")
for (i in 1:I){print(mean(((1/lambda[(nburnin+1):mcsim,i]/exp(theta[(nburnin+1):mcsim,i]))^(1/r[(nburnin+1):mcsim,i])*(log(2))^(1/r[(nburnin+1):mcsim,i]))>0.6, na.rm=TRUE))}
longdata <- longdata_full
if (mean(((1/lambda[(nburnin+1):mcsim,i]/exp(theta[(nburnin+1):mcsim,i]))^(1/r[(nburnin+1):mcsim,i])*(log(2))^(1/r[(nburnin+1):mcsim,i]))>0.6, na.rm=TRUE)<0.05){longdata <- longdata[(longdata[,1] !=i),]}

# print("mean of theta")
# for (i in 1:I){print(mean(theta[(nburnin+1):mcsim,i]))}
# print("mean of lambda")
# for (i in 1:I){print(mean(lambda[(nburnin+1):mcsim,i]))}
# print("mean of r")
# for (i in 1:I){print(mean(r[(nburnin+1):mcsim,i]))}
# print("mean of hyperparameters")
# mean(hyperparameter_1[(nburnin+1):mcsim])
# mean(hyperparameter_2[(nburnin+1):mcsim])
# mean(hyperparameter_3[(nburnin+1):mcsim])
# mean(hyperparameter_4[(nburnin+1):mcsim])
# mean(hyperparameter_5[(nburnin+1):mcsim])
# mean(hyperparameter_6[(nburnin+1):mcsim])
# print("mean of median survival time")
# for (i in 1:I){print(mean((1/lambda[(nburnin+1):mcsim,i]/exp(theta[(nburnin+1):mcsim,i]))^(1/r[(nburnin+1):mcsim,i])*(log(2))^(1/r[(nburnin+1):mcsim,i])))}










#-------------------------------------# 
#-------------------------------------# 
#-------------------------------------# 
#Part 3: BHM Interim 2# 
#-------------------------------------# 
#-------------------------------------# 
#-------------------------------------# 
if (length(unique(longdata[,1]))!=0){
  nknot=2
  mcsim=8000
  nburnin=4000
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
  
  lambda<-matrix(NA, nrow=mcsim, ncol=I)
  r<-matrix(NA, nrow=mcsim, ncol=I)
  theta<-matrix(NA, nrow=mcsim, ncol=I)
  hyperparameter_1<-rep(NA,mcsim)
  hyperparameter_2<-rep(NA,mcsim)
  hyperparameter_3<-rep(NA,mcsim)
  hyperparameter_4<-rep(NA,mcsim)
  hyperparameter_5<-rep(NA,mcsim)
  hyperparameter_6<-rep(NA,mcsim)
  # Initialize
  lambda[1,]<- rep(1.5, I)
  r[1,]<- rep(1.05, I)
  theta[1,]<- rep(-0.75, I)
  hyperparameter_1[1]<- -1
  hyperparameter_2[1]<- 1
  hyperparameter_3[1]<- 0.1
  hyperparameter_4[1]<- 0.1
  hyperparameter_5[1]<- 0.1
  hyperparameter_6[1]<- 0.1
  
  
  #-------------------------------------#
  # B. response part data organize      #
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
  
  
  #----------------------------------------------#
  # D. MCMC part                                 #
  #                                              #
  #----------------------------------------------#    
  
  set.seed(seed)
  for (iter in 2:mcsim){
    #----------------------------------#
    # update theta, lambda, r#
    #----------------------------------#
    # log density
    logden <- function(theta, lambda, r, hyperparameter_1, hyperparameter_2, hyperparameter_3, hyperparameter_4, hyperparameter_5, hyperparameter_6) {
      logdensity=rep(NA,I)
      for (j in 1:I){logdensity[j] = sum(omega_short[grp_short==j]*log(exp(-lambda[j]*t_short[grp_short==j]^r[j]*exp(theta[j]))*(lambda[j]*exp(theta[j])*r[j]*t_short[grp_short==j]^(r[j]-1))) + (1-omega_short[grp_short==j])*(-lambda[j]*t_short[grp_short==j]^r[j]*exp(theta[j]))) + dnorm(theta[j], hyperparameter_1, hyperparameter_2, log=TRUE) + dgamma(lambda[j], shape=hyperparameter_3, rate=hyperparameter_4, log=TRUE) + dgamma(r[j], shape=hyperparameter_5, rate=hyperparameter_6, log=TRUE)}
      return(sum(logdensity)+dnorm(hyperparameter_1, 0, 10, log=TRUE) + dunif(hyperparameter_2, 0, 10, log=TRUE) + dunif(hyperparameter_3, 0, 1, log=TRUE) + dunif(hyperparameter_4, 0, 1, log=TRUE) + dunif(hyperparameter_5, 0, 1, log=TRUE) + dunif(hyperparameter_6, 0, 1, log=TRUE))
    }
    # update theta
    theta.proposal = function(x){return(rnorm(1,x,0.01))}
    i=1
    {
      theta.proposed = theta.proposal(theta[iter-1,i])
      prob.theta = exp(logden(c(theta.proposed, theta[iter-1,(i+1):I]), lambda[iter-1,], r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                       -logden(c(theta[iter-1,i:I]), lambda[iter-1,], r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
      )
      if (runif(1) < prob.theta){
        theta[iter,i] = theta.proposed
      }else{
        theta[iter,i] = theta[iter-1,i]
      }
    }
    for (i in 2:(I-1)){
      theta.proposed = theta.proposal(theta[iter-1,i])
      prob.theta = exp(logden(c(theta[iter,1:(i-1)], theta.proposed, theta[iter-1,(i+1):I]), lambda[iter-1,], r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                       -logden(c(theta[iter,1:(i-1)], theta[iter-1,i:I]), lambda[iter-1,], r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
      )
      if (runif(1) < prob.theta){
        theta[iter,i] = theta.proposed
      }else{
        theta[iter,i] = theta[iter-1,i]
      }
    }
    i=I
    {
      theta.proposed = theta.proposal(theta[iter-1,i])
      prob.theta = exp(logden(c(theta[iter,1:(i-1)], theta.proposed), lambda[iter-1,], r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                       -logden(c(theta[iter,1:(i-1)], theta[iter-1,i:I]), lambda[iter-1,], r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
      )
      if (runif(1) < prob.theta){
        theta[iter,i] = theta.proposed
      }else{
        theta[iter,i] = theta[iter-1,i]
      }
    }
    # update lambda
    lambda.proposal = function(x){return(max(0.01,rnorm(1,x,0.01)))}
    i=1
    {
      lambda.proposed = lambda.proposal(lambda[iter-1,i])
      prob.lambda = exp(logden(theta[iter,], c(lambda.proposed, lambda[iter-1,(i+1):I]), r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                        -logden(theta[iter,], c(lambda[iter-1,i:I]), r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
      )
      if (runif(1) < prob.lambda){
        lambda[iter,i] = lambda.proposed
      }else{
        lambda[iter,i] = lambda[iter-1,i]
      }
    }
    for (i in 2:(I-1)){
      lambda.proposed = lambda.proposal(lambda[iter-1,i])
      prob.lambda = exp(logden(theta[iter,], c(lambda[iter,1:(i-1)], lambda.proposed, lambda[iter-1,(i+1):I]), r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                        -logden(theta[iter,], c(lambda[iter,1:(i-1)], lambda[iter-1,i:I]), r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
      )
      if (runif(1) < prob.lambda){
        lambda[iter,i] = lambda.proposed
      }else{
        lambda[iter,i] = lambda[iter-1,i]
      }
    }
    i=I
    {
      lambda.proposed = lambda.proposal(lambda[iter-1,i])
      prob.lambda = exp(logden(theta[iter,], c(lambda[iter,1:(i-1)], lambda.proposed), r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                        -logden(theta[iter,], c(lambda[iter,1:(i-1)], lambda[iter-1,i:I]), r[iter-1,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
      )
      if (runif(1) < prob.lambda){
        lambda[iter,i] = lambda.proposed
      }else{
        lambda[iter,i] = lambda[iter-1,i]
      }
    }
    # update r
    r.proposal = function(x){return(max(0.01,rnorm(1,x,0.01)))}
    i=1
    {
      r.proposed = r.proposal(r[iter-1,i])
      prob.r = exp(logden(theta[iter,], lambda[iter,], c(r.proposed, r[iter-1,(i+1):I]), hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                   -logden(theta[iter,], lambda[iter,], c(r[iter-1,i:I]), hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
      )
      if (runif(1) < prob.r){
        r[iter,i] = r.proposed
      }else{
        r[iter,i] = r[iter-1,i]
      }
    }
    for (i in 2:(I-1)){
      r.proposed = r.proposal(r[iter-1,i])
      prob.r = exp(logden(theta[iter,], lambda[iter,], c(r[iter,1:(i-1)], r.proposed, r[iter-1,(i+1):I]), hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                   -logden(theta[iter,], lambda[iter,], c(r[iter,1:(i-1)], r[iter-1,i:I]), hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
      )
      if (runif(1) < prob.r){
        r[iter,i] = r.proposed
      }else{
        r[iter,i] = r[iter-1,i]
      }
    }
    i=I
    {
      r.proposed = r.proposal(r[iter-1,i])
      prob.r = exp(logden(theta[iter,], lambda[iter,], c(r[iter,1:(i-1)], r.proposed),  hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                   -logden(theta[iter,], lambda[iter,], c(r[iter,1:(i-1)], r[iter-1,i:I]), hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
      )
      if (runif(1) < prob.r){
        r[iter,i] = r.proposed
      }else{
        r[iter,i] = r[iter-1,i]
      }
    }
    # update hyperparameter_1
    hyperparameter_1.proposal = function(x){return(rnorm(1,x,0.01))}
    hyperparameter_1.proposed = hyperparameter_1.proposal(hyperparameter_1[iter-1])
    prob.hyperparameter_1 = exp(logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1.proposed, hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                                -logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter-1], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
    )
    if (runif(1) < prob.hyperparameter_1){
      hyperparameter_1[iter] = hyperparameter_1.proposed
    }else{
      hyperparameter_1[iter] = hyperparameter_1[iter-1]
    }
    # update hyperparameter_2
    hyperparameter_2.proposal = function(x){return(max(0.0001,rnorm(1,x,0.01)))}
    hyperparameter_2.proposed = hyperparameter_2.proposal(hyperparameter_2[iter-1])
    prob.hyperparameter_2 = exp(logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2.proposed, hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                                -logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter-1], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
    )
    if (runif(1) < prob.hyperparameter_2){
      hyperparameter_2[iter] = hyperparameter_2.proposed
    }else{
      hyperparameter_2[iter] = hyperparameter_2[iter-1]
    }
    # update hyperparameter_3
    hyperparameter_3.proposal = function(x){return(max(0.0001,rnorm(1,x,0.01)))}
    hyperparameter_3.proposed = hyperparameter_3.proposal(hyperparameter_3[iter-1])
    prob.hyperparameter_3 = exp(logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter], hyperparameter_3.proposed, hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                                -logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter], hyperparameter_3[iter-1], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
    )
    if (runif(1) < prob.hyperparameter_3){
      hyperparameter_3[iter] = hyperparameter_3.proposed
    }else{
      hyperparameter_3[iter] = hyperparameter_3[iter-1]
    }
    # update hyperparameter_4
    hyperparameter_4.proposal = function(x){return(max(0.0001,rnorm(1,x,0.01)))}
    hyperparameter_4.proposed = hyperparameter_4.proposal(hyperparameter_4[iter-1])
    prob.hyperparameter_4 = exp(logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter], hyperparameter_3[iter], hyperparameter_4.proposed, hyperparameter_5[iter-1], hyperparameter_6[iter-1])
                                -logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter], hyperparameter_3[iter], hyperparameter_4[iter-1], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
    )
    if (runif(1) < prob.hyperparameter_4){
      hyperparameter_4[iter] = hyperparameter_4.proposed
    }else{
      hyperparameter_4[iter] = hyperparameter_4[iter-1]
    }
    # update hyperparameter_5
    hyperparameter_5.proposal = function(x){return(max(0.0001,rnorm(1,x,0.01)))}
    hyperparameter_5.proposed = hyperparameter_5.proposal(hyperparameter_5[iter-1])
    prob.hyperparameter_5 = exp(logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter], hyperparameter_3[iter], hyperparameter_4[iter], hyperparameter_5.proposed, hyperparameter_6[iter-1])
                                -logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter], hyperparameter_3[iter], hyperparameter_4[iter], hyperparameter_5[iter-1], hyperparameter_6[iter-1])
    )
    if (runif(1) < prob.hyperparameter_5){
      hyperparameter_5[iter] = hyperparameter_5.proposed
    }else{
      hyperparameter_5[iter] = hyperparameter_5[iter-1]
    }
    # update hyperparameter_6
    hyperparameter_6.proposal = function(x){return(max(0.0001,rnorm(1,x,0.01)))}
    hyperparameter_6.proposed = hyperparameter_6.proposal(hyperparameter_6[iter-1])
    prob.hyperparameter_6 = exp(logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter], hyperparameter_3[iter], hyperparameter_4[iter], hyperparameter_5[iter], hyperparameter_6.proposed)
                                -logden(theta[iter,], lambda[iter,], r[iter,], hyperparameter_1[iter], hyperparameter_2[iter], hyperparameter_3[iter], hyperparameter_4[iter], hyperparameter_5[iter], hyperparameter_6[iter-1])
    )
    if (runif(1) < prob.hyperparameter_6){
      hyperparameter_6[iter] = hyperparameter_6.proposed
    }else{
      hyperparameter_6[iter] = hyperparameter_6[iter-1]
    }
  }
  
  # Qf = 0.05, Q = 0.8
  # check
  print("BHM Interim 2")
  for (i in 1:I){print(mean(((1/lambda[(nburnin+1):mcsim,i]/exp(theta[(nburnin+1):mcsim,i]))^(1/r[(nburnin+1):mcsim,i])*(log(2))^(1/r[(nburnin+1):mcsim,i]))>0.6, na.rm=TRUE))}
}




