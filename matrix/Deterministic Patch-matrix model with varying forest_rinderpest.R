######### 
# Rinderpest Patch-Matrix model with slider for patch to matrix transmission
# Christina Faust Version October 8 2014
########
rm(list = ls())


library(deSolve)
library(manipulate)
library(graphics)
#see word document "pathogen spillover for fragmented landscapes" for equation and parameter definitions
#note that kappa describes position on density-dependent (0) or  frequency (1) dependent transmission continuum
#beta.pp is pathogen transmission within patch
#beta.mm is pathogen transmission within matrix
#beta.pm is transmission from patch to matrix (ie spillover)

# f is the proportion forested (ranges from 1 to 0)
# 1. scaled by edge; edge= 1+cos(f*(pi*3/2)-2.5); ranges from 0 to 2 across values of f; symmetrical)
# 2. decreases k.p (reduce carrying capacity in patch habitat, because of habitat loss) 
# 3. relative increase in k.m (increase carrying capacity in matrix habitat, as more matrix created)

patch.matrix.model.f <- function(Time, State, Parameters) {
  with(as.list(c(State, Parameters)), {
    N.p<-S.p+I.p+R.p
    N.m<-S.m+I.m+R.m
    
    dS.p<-N.p*b.p*(1-N.p*k.p/f) - ((S.p*beta.pp*I.p)/N.p^kappa+(epsilon*beta.mp*I.m*S.p)/(N.p+epsilon*N.m)^kappa) - d.p*S.p + gamma.p*R.p
    dI.p<-(beta.pp*S.p*I.p)/N.p^kappa+(epsilon*beta.mp*S.p*I.m)/(N.p+epsilon*N.m)^kappa - I.p*(alpha.p+d.p+sigma.p)
    dR.p<-I.p*sigma.p-R.p*(d.p+gamma.p)
    
    dS.m<-N.m*b.m*(1-N.m*k.m/(1.5-f)) - ((beta.mm*S.m*I.m)/N.m^kappa+(epsilon*beta.pm*S.p*I.p)/(epsilon*N.p+N.m)^kappa) - d.m*S.m + gamma.m*R.m
    dI.m<-((beta.mm*S.m*I.m)/N.m^kappa+(epsilon*beta.pm*I.p*S.p)/(epsilon*N.p+N.m)^kappa) - I.m*(alpha.m+d.m+sigma.m)
    dR.m<-I.m*sigma.m-R.m*(d.m+gamma.m)
    
   return(list(c(dS.p, dI.p, dR.p,dS.m, dI.m, dR.m)))
  })
}



#Between species contact rate varying by deforestation; can change functional response
f=seq(0,1,by=0.01)
epsilon=(1+cos(f*(pi*3/2)-1.5))
plot(f, epsilon, bty='n', type='l', col='darkgreen', lwd=4,
     xlim=c(0,1), xlab= 'proportion forested', 
     ylim=c(0,2), ylab= expression(paste("forested edge (",epsilon,")")),
     main="cosine function")
#
f=seq(0,1,by=0.01)
epsilon=2- ((3*f^3)/(1+f^(4)))
plot(f, epsilon, bty='n', type='l', col='darkgreen', lwd=4,
     xlim=c(0,1), xlab= 'proportion forested', 
     #ylim=c(0,2), 
     ylab= expression(paste("forested edge (",epsilon,")")),
     main='increasing edge')

epsilon= numeric(length(f))

for (i in 1:length(f)){
  if(f[i]<0.2){
    epsilon[i]=10*f[i]
  }else{
    epsilon[i]=1/(2.5*f[i])
  } 
}
plot(f, epsilon, bty='n', type='l', col='darkgreen', lwd=4,
     xlim=c(0,1), xlab= 'proportion forested', 
     #ylim=c(0,2), 
     ylab= expression(paste("forested edge (",epsilon,")")),
     main='max at 2')

f=seq(0,1,by=0.01)
epsilon=(1.5+cos(3-f*(pi*8/2)))/1.5
plot(f, epsilon, bty='n', type='l', col='darkgreen', lwd=4,
     xlim=c(0,1), xlab= 'proportion forested', 
     ylim=c(0,2), ylab= expression(paste("forested edge (",epsilon,")")),
     main="cosine function")
#


### output of times ###
times <- seq(0, 100, by = 0.01)

plot.fun<-function(f){

### set some parameters
  #parameters from 
params<-c(b.p=.5,d.p=0.1,k.p=0.001,beta.pp=0.05,beta.pm=0.001,gamma.p=0,alpha.p=0.001,
          sigma.p=0.05,epsilon=2- ((3*f^3)/(1+f^(4))), f=f,
          b.m=.1,d.m=0.02,k.m=0.0001,beta.mm=0.001,beta.mp=0.0,gamma.m=0,alpha.m=0.01,
          sigma.m=0.05,kappa=0)
initial.values<-c(S.p=20,I.p=1,R.p=0,S.m=100,I.m=0,R.m=0)
print(system.time(
  out<-ode(func=patch.matrix.model.f,y=initial.values,parms=params,times=times))) # if you are having LSODA error messages you may need to switch to ode45 solver; insert, method="ode45" into the  ode function)
head(out,n=2)
par(mfrow=c(2,1))

#plot the patch hosts
matplot(out[,"time"],out[,2:4],type="l",xlab="time",ylab="number", bty='n',
        main="Patch Hosts",lwd=2, ylim=c(0,500))
legend("topright",c("susc","inf","rec"),col=1:3,lty=1:3)

#plot the matrix hosts
matplot(out[,"time"],out[,5:7],type="l",xlab="time",ylab="number", bty='n',
        main="Matrix Hosts",lwd=2, ylim=c(0,5000))
legend("topright",c("susc","inf","rec"),col=1:3,lty=1:3)

}

#to plot graphs with varying 'forested' areas
manipulate(plot.fun(f),f=slider(0.0000001,1.0)) #can't divide by zero
f<-slider(0.0000001,1.0)



#### R0; secondary infections in a completely naive population

#functions for calculating matrix parts for within and between species transmissions and transitions
ngm.within.dd<-function(beta, N1, alpha,sigma,d){
  z = (beta*N1)/(alpha+sigma+d)
  return(z)
}
ngm.between.dd<-function(epsilon, beta, N1, N2, alpha, sigma, d){
  z = (beta*(epsilon*N1^2/(epsilon*N1+N2))/(alpha+sigma+d))
  return(z)
}


#R0 at different f; just by calculating a pathogen that comes into a completely naive population at carrying capacity for that proportion forested
f=seq(0.01,1,by=0.01) #vector for different forested proportions
params<-c(d.p=0.1,k.p=0.007,beta.pp=0.01,beta.pm=0.001,gamma.p=0.03,alpha.p=0.001,sigma.p=0.05,#patch parameters that are relevant
          d.m=0.02,k.m=0.006,beta.mm=0.001,beta.mp=0.0,gamma.m=0.05,alpha.m=0.01,sigma.m=0.05) #matrix
epsilon=(1+cos(f*(pi*3/2)-2)) #vector of edges at different forested levels
S.m = (1.1-f)/(params[["k.m"]]) #calculating carrying capacity at each level of 
S.p = (f)/(params[["k.p"]])

#empty vectors for DD R0
R0.p.dd= numeric(length(f)) #assuming no matrix species present (or not transmissible)
R0.m.dd= numeric(length(f)) #assuming no patch species present (or not transmissible)
R0.combin.dd= numeric(length(f))
ngm.dd<-matrix(NA,nrow=2,ncol=2) #next generation matrix 2x2 for 2 species

#loop to calculate each R0 and NGM at varying forest levels
for (i in 1:length(f)){
  ngm1.11 <- ngm.within.dd(params[["beta.mm"]],S.m[i],params[["alpha.m"]],params[["sigma.m"]], params[["d.m"]])
  ngm1.22 <- ngm.within.dd(params[["beta.pp"]],S.p[i],params[["alpha.p"]],params[["sigma.p"]], params[["d.p"]])
  ngm1.12 <- ngm.between.dd(epsilon,params[["beta.mp"]], S.m[i],S.p[i], params[["alpha.m"]],params[["sigma.m"]], params[["d.m"]])
  ngm1.21 <- ngm.between.dd(epsilon,params[["beta.pm"]], S.p[i], S.m[i], params[["alpha.p"]],params[["sigma.p"]], params[["d.p"]])
  ngm.dd<-matrix(c(ngm1.11,ngm1.21,ngm1.12,ngm1.22), nrow=2,ncol=2)
  eigen.ngm.dd<-eigen(ngm.dd)
  R0.combin.dd[i] = eigen.ngm.dd$values[[1]]
  R0.p.dd[i] = (params[["beta.pp"]]*S.p[i]/(params[["alpha.p"]]+params[["gamma.p"]]+params[["sigma.p"]]))
  R0.m.dd[i] = (params[["beta.mm"]]*S.m[i]/(params[["alpha.m"]]+params[["gamma.m"]]+params[["sigma.m"]]))
}


plot(f, R0.combin.dd, type="l",lwd=4, 
     ylim=c(0,20),bty='n', 
     ylab="R0",xlab="proportion forested")
lines(f, R0.m.dd, lwd=4, col='darkgoldenrod')
lines(f, R0.p.dd, lwd=4, col='forestgreen')
abline(h=1, lty=2, lwd=4, col='darkred')
lines(f, 7.5*epsilon, lty=3, lwd=1, col='grey44')
legend("topleft", c("combined" , "within matrix" ,"within patch","relative edge"),
       lty=c(1,1,1,3), lwd=3, col=c("black", "darkgoldenrod", "forestgreen", "grey44"), bty='n' )

