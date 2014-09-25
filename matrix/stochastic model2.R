## Modification of code from Jesse & Heesterbeek J. Theor Biol (2011) 275 12-20 # #
# Initially has S I R in patch and matrix
rm(list = ls())

## Define some parameters for patch host
b_p=0.5 #birth rate of patch host
d_p=0.01 #death rate of patch host
s_p=0.1  #recovery rate of patch host
nu_p=0.05 #loss of immunity of patch host
alpha_p=0.3 #death rate of infected patch hosts
k_p=0.01 #density dependence in the absence of infection for patch host

## And the same parameters for the matrix host
b_m=0.5 #birth rate of matrix host
d_m=0.01 #death rate of matrix host
s_m=0.1  #recovery rate of matrix host
nu_m=0.05 #loss of immunity of matrix host
alpha_m=0.2 #death rate of infected matrix hosts
k_m=0.0005 #density dependence in the absence of infection for matrix host
## now some infection rates
beta_pp=0.001 #infection rate in patches
beta_mm=0.001  #infection rate in matrix
beta_mp=0.0001  #infection rate from matrix to patch
beta_pm=0.0001  #infection rate from patch to matrix


##########################################################
# set up matrices to take outputs (SIR) for each host type
hosts_p=matrix(NA,nrow=1000,ncol=3)
colnames(hosts_p)<-c("S","I","R")
hosts_m=matrix(NA,nrow=1000,ncol=3)
colnames(hosts_m)<-c("S","I","R")
#set up initial values
hosts_p[1,1]=10
hosts_p[1,2]=0
hosts_p[1,3]=0
hosts_m[1,1]=200
hosts_m[1,2]=10
hosts_m[1,3]=0
fadeout_patch=0
#now run the guts of the model
for(t in 1:999) {
  #for the patch host first
  #newborns produced from poisson process
  hosts_p[t+1,1]=hosts_p[t,1]+rpois(1,(hosts_p[t,1]+hosts_p[t,2]+hosts_p[t,3])*b_p)
  #probability of survival  in absence of disease
  surv_ps=exp(-(d_p+k_p*(hosts_p[t,1]+hosts_p[t,2]+hosts_p[t,3])))
  #and now apply binomial prob to produce survivors
  hosts_p[t+1,1]=rbinom(1,hosts_p[t+1,1],surv_ps)
  #probability of becoming infected
  infprob=(1-exp(-hosts_p[t+1,1]*(beta_pp*hosts_p[t,2]+beta_mp*hosts_m[t,2])))
  #number infected from binomial
  inf_p=rbinom(1,hosts_p[t+1,1],infprob)
  #so decrement susceptibles
  hosts_p[t+1,1]=hosts_p[t+1,1]-inf_p
  #and increment infecteds
  hosts_p[t+1,2]=hosts_p[t,2]+inf_p
  #prob of survival if infected
  surv_pi=exp(-(d_p+alpha_p+k_p*(hosts_p[t,1]+hosts_p[t,2]+hosts_p[t,3])))
  #and calculate survivors from binomial
  hosts_p[t+1,2]=rbinom(1,hosts_p[t+1,2],surv_pi)
  #calculate number of recovereds
  rec_p=rbinom(1,hosts_p[t+1,2],1-exp(-s_p))
  #so decrement infecteds
  hosts_p[t+1,2]=hosts_p[t+1,2]-rec_p
  #and increment recovered class
  hosts_p[t+1,3]=hosts_p[t,3]+rec_p
  #these suffer same death rate as susceptibles
  hosts_p[t+1,3]=rbinom(1,hosts_p[t+1,3],surv_ps)
  #and some of them lose immunity
  lossimm_p=rbinom(1,hosts_p[t+1,3],1-exp(-nu_p))
  #so decrement the recovered class
  hosts_p[t+1,3]=hosts_p[t+1,3]-lossimm_p
  #and put them back into susceptibles
  hosts_p[t+1,1]=hosts_p[t+1]+lossimm_p
  ############all done, now repeat for matrix species
  
  #newborns produced from poisson process
  hosts_m[t+1,1]=hosts_m[t,1]+rpois(1,(hosts_p[t,1]+hosts_p[t,2]+hosts_p[t,3])*b_m)
  #probability of survival  in absence of disease
  surv_ms=exp(-(d_m+k_m*(hosts_m[t,1]+hosts_m[t,2]+hosts_m[t,3])))
  #and now apply binomial prob to produce survivors
  hosts_m[t+1,1]=rbinom(1,hosts_m[t+1,1],surv_ms)
  #probability of becoming infected
  infprob=(1-exp(-hosts_m[t+1,1]*(beta_mm*hosts_m[t,2]+beta_pm*hosts_p[t+1,2])))
  #number infected from binomial
  inf_m=rbinom(1,hosts_m[t+1,1],infprob)
  #so decrement susceptibles
  hosts_m[t+1,1]=hosts_m[t+1,1]-inf_m
  #and increment infecteds
  hosts_m[t+1,2]=hosts_m[t,2]+inf_m
  #prob of survival if infected
  surv_mi=exp(-(d_m+alpha_m+k_m*(hosts_m[t,1]+hosts_m[t,2]+hosts_m[t,3])))
  #and calculate survivors from binomial
  hosts_m[t+1,2]=rbinom(1,hosts_m[t+1,2],surv_mi)
  #calculate number of recovereds
  rec_m=rbinom(1,hosts_m[t+1,2],1-exp(-s_m))
  #so decrement infecteds
  hosts_m[t+1,2]=hosts_m[t+1,2]-rec_m
  #and increment recovered class
  hosts_m[t+1,3]=hosts_m[t,3]+rec_m
  #these suffer same death rate as susceptibles
  hosts_m[t+1,3]=rbinom(1,hosts_m[t+1,3],surv_ms)
  #and some of them lose immunity
  lossimm_m=rbinom(1,hosts_m[t+1,3],1-exp(-nu_m))
  #so decrement the recovered class
  hosts_m[t+1,3]=hosts_m[t+1,3]-lossimm_m
  #and put them back into susceptibles
  hosts_m[t+1,1]=hosts_m[t+1]+lossimm_m
   
  
}
  
plot(1:1000,(hosts_m[,1]+hosts_m[,2]+hosts_m[,3]),type="l",ylab="population",xlab="time",ylim=c(0,300))
lines(1:1000,(hosts_p[,1]+hosts_p[,2]+hosts_p[,3]),type="l",col="red")

plot(1:1000,hosts_m[,2]/(hosts_m[,1]+hosts_m[,2]+hosts_m[,3]),type="l",ylab="prevalence",xlab="time")
lines(1:1000,hosts_p[,2]/(hosts_p[,1]+hosts_p[,2]+hosts_p[,3]),type="l",col="red")
#count fadeouts
fade=0
for (t in 1:999){
if(hosts_p[t,2]==0) fade<-fade+1
}
