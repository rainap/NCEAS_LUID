######### 
# Deterministic Patch-Matrix model with slider for patch to matrix transmission
#Hamish McCallum Version September 24 2014
# Adjusted by CFaust (documentation on git hub)
########
rm(list = ls())


library(deSolve)
library(manipulate)

patch.matrix.model <- function(Time, State, Parameters) {
  with(as.list(c(State, Parameters)), {
    N.p<-S.p+I.p+R.p
    dS.p<-N.p*(b.p-N.p*k.p) - S.p*(beta.pp*I.p+beta.mp*I.m)/N.p^kappa - S.p*d.p +gamma.p*I.p
    dI.p<-S.p*(beta.pp*I.p+beta.mp*I.m)/N.p^kappa-I.p*(alpha.p+d.p+sigma.p+N.p*k.p)
    dR.p<-I.p*sigma.p-R.p*(d.p+gamma.p+N.p*k.p)
    
    N.m<-S.m+I.m+R.m
    dS.m<-S.m*(b.m-d.m-N.m*k.m)-S.m*(beta.mm*I.m+beta.pm*I.p)/N.m^kappa+gamma.m*I.m
    dI.m<-S.m*(beta.mm*I.m+beta.pm*I.p)/N.m^kappa-I.m*(alpha.m+d.m+sigma.m+N.m*k.m)
    dR.m<-I.m*sigma.m-R.m*(d.m+gamma.m+N.m*k.m)
    
   return(list(c(dS.p, dI.p, dR.p,dS.m, dI.m, dR.m)))
  })
}

### output of times ###
times <- seq(0, 100, by = 1)

plot.fun<-function(beta.pm.sl){

### set some parameters
params<-c(b.p=.5, #birth rate for p (patch) species
          d.p=0.1,
          k.p=0.001, #1/carrying capacity
          beta.pp=0.01, #within patch species transmission rate
          beta.pm=beta.pm.sl, #variable patch to matrix transmission
          gamma.p=0.05, #loss of immunity R->S
          alpha.p=0.2, #disease-induced mortality
          sigma.p=0.05, #recovery rate I->R
          b.m=.1, #Matrix parameters
          d.m=0.02,
          k.m=0.0001,
          beta.mm=0.0001,
          beta.mp=0.0,
          gamma.m=0.05,
          alpha.m=0.0001,
          sigma.m=0.05,
          kappa=0) # can adjust for frequency vs. density dependent

initial.values<-c(S.p=20,I.p=1,R.p=0,S.m=100,I.m=0,R.m=0)

print(system.time(
  out<-ode(func=patch.matrix.model,y=initial.values,parms=params,times=times)))
head(out,n=3)
par(mfrow=c(2,1))
#plot the patch hosts
matplot(out[,"time"],out[,2:4],type="l",xlab="time",ylab="number",
        main="Patch Hosts",lwd=2)
legend("topright",c("susc","inf","rec"),col=1:3,lty=1:3)

#plot the matrix hosts
matplot(out[,"time"],out[,5:7],type="l",xlab="time",ylab="number",
        main="Matrix Hosts",lwd=2)
legend("topright",c("susc","inf","rec"),col=1:3,lty=1:3)


}

manipulate(plot.fun(beta.pm.sl),beta.pm.sl=slider(0,0.01))


},

beta.pm.sl<-slider(0,0.1)
)

