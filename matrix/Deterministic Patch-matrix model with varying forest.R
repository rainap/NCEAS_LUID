######### 
# Deterministic Patch-Matrix model with slider for patch to matrix transmission
# Hamish McCallum Version September 24 2014
########
rm(list = ls())


library(deSolve)
library(manipulate)
#see word document "pathogen spillover for fragmented landscapes" for equation and parameter definitions
# most of this is pretty obvious
#note that kappa describes position on density-dependent frequency dependent transmission continuum
#beta.pp is pathogen transmission within patch
#beta.mm is pathogen transmission within matrix
#beta.pm is transmission from patch to matrix (ie spillover)

# f is the proportion forested (ranges from 1 to 0)
# 1  increases and then decreases beta.pm and beta.mp (scaled by edge; edge= 1+cos(f*(pi*3/2)-2.5); ranges from 0 to 2 across values of f; symmetrical)
# 2  decreases k.p (reduce carrying capacity in patch habitat, because of habitat loss) 
# 3 increases k.m (increase carrying capacity in matrix habitat, as more matrix created)

patch.matrix.model.f <- function(Time, State, Parameters) {
  with(as.list(c(State, Parameters)), {
    N.p<-S.p+I.p+R.p
    dS.p<-S.p*(b.p-d.p-N.p*k.p/f)-S.p*(beta.pp*I.p+beta.mp*I.m*edge)/N.p^kappa+gamma.p*I.p
    dI.p<-S.p*(beta.pp*I.p+beta.mp*I.m*edge)/N.p^kappa-I.p*(alpha.p+d.p+sigma.p+N.p*k.p/f)
    dR.p<-I.p*sigma.p-R.p*(d.p+gamma.p+N.p*k.p/f)
    
    N.m<-S.m+I.m+R.m
    dS.m<-S.m*(b.m-d.m-N.m*k.m/(2-f))-S.m*(beta.mm*I.m+beta.pm*I.p*edge)/N.m^kappa+gamma.m*I.m
    dI.m<-S.m*(beta.mm*I.m+beta.pm*I.p*edge)/N.m^kappa-I.m*(alpha.m+d.m+sigma.m+N.m*k.m/(2-f))
    dR.m<-I.m*sigma.m-R.m*(d.m+gamma.m+N.m*k.m/(2-f))
    
   return(list(c(dS.p, dI.p, dR.p,dS.m, dI.m, dR.m)))
  })
}

### output of times ###
times <- seq(0, 100, by = 1)

plot.fun<-function(f){

### set some parameters
params<-c(b.p=.5,d.p=0.1,k.p=0.001,beta.pp=0.01,beta.pm=0.001,gamma.p=0.05,alpha.p=0.001,
          sigma.p=0.05,edge=(1+cos(f*(pi*3/2)-2.5)), f=f,
          b.m=.1,d.m=0.02,k.m=0.0001,beta.mm=0.0001,beta.mp=0.0,gamma.m=0.05,alpha.m=0.01,
          sigma.m=0.05,kappa=0)
initial.values<-c(S.p=20,I.p=2,R.p=0,S.m=100,I.m=0,R.m=0)
print(system.time(
  out<-ode(func=patch.matrix.model.f,y=initial.values,parms=params,times=times)))
head(out,n=3)
par(mfrow=c(2,1))
#plot the patch hosts
matplot(out[,"time"],out[,2:4],type="l",xlab="time",ylab="number",
        main="Patch Hosts",lwd=2, ylim=c(0,150))
legend("topright",c("susc","inf","rec"),col=1:3,lty=1:3)

#plot the matrix hosts
matplot(out[,"time"],out[,5:7],type="l",xlab="time",ylab="number",
        main="Matrix Hosts",lwd=2, ylim=c(0,600))
legend("topright",c("susc","inf","rec"),col=1:3,lty=1:3)


}

manipulate(plot.fun(f),f=slider(0.0000001,1.0))


},

f<-slider(0,1.0)
)

