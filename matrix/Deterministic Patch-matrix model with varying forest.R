######### 
# Deterministic Patch-Matrix model with slider for patch to matrix transmission
# Christina Faust Version October 8 2014
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
    N.m<-S.m+I.m+R.m
    
    dS.p<-N.p*b.p*(1-N.p*k.p/f) - ((S.p*beta.pp*I.p)/N.p^kappa+(epsilon*beta.mp*I.m)/(N.p+N.m)^kappa) - d.p*S.p + gamma.p*R.p
    dI.p<-(S.p*beta.pp*I.p)/N.p^kappa+(epsilon*beta.mp*I.m)/(N.p+N.m)^kappa - I.p*(alpha.p+d.p+sigma.p)
    dR.p<-I.p*sigma.p-R.p*(d.p+gamma.p)
    
    dS.m<-N.m*(b.m-N.m*k.m/(2-f)) - ((S.m*beta.mm*I.m)/N.m^kappa+(epsilon*beta.pm*I.p)/(N.p+N.m)^kappa) - d.m*S.m + gamma.m*R.m
    dI.m<-((S.m*beta.mm*I.m)/N.m^kappa+(epsilon*beta.pm*I.p)/(N.p+N.m)^kappa) - I.m*(alpha.m+d.m+sigma.m)
    dR.m<-I.m*sigma.m-R.m*(d.m+gamma.m)
    
   return(list(c(dS.p, dI.p, dR.p,dS.m, dI.m, dR.m)))
  })
}

### output of times ###
times <- seq(0, 100, by = 0.01)

plot.fun<-function(f){

### set some parameters
params<-c(b.p=.5,d.p=0.1,k.p=0.001,beta.pp=0.01,beta.pm=0.001,gamma.p=0.05,alpha.p=0.001,
          sigma.p=0.05,epsilon=(1+cos(f*(pi*3/2)-2.5)), f=f,
          b.m=.1,d.m=0.02,k.m=0.0001,beta.mm=0.0001,beta.mp=0.0,gamma.m=0.05,alpha.m=0.01,
          sigma.m=0.05,kappa=0)
initial.values<-c(S.p=20,I.p=1,R.p=0,S.m=100,I.m=0,R.m=0)
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

#plot the R0
plot(type, R0, xlim=c(0,5))
abline(h=1, lty=3)


}

manipulate(plot.fun(f),f=slider(0.0000001,1.0))


},

f<-slider(0,1.0)
)

ngm1.1<-function(beta, alpha, gamma){
  z = beta/(alpha+beta+gamma)
  return(z)
}
ngm1.2<-function(beta, edge){
  z = beta*edge
  return(z)
}


f=0.5
params<-c(b.p=.5,d.p=0.1,k.p=0.001,beta.pp=0.01,beta.pm=0.001,gamma.p=0.05,alpha.p=0.001,
          sigma.p=0.05,edge=(1+cos(f*(pi*3/2)-2)), f=f,
          b.m=.1,d.m=0.02,k.m=0.0001,beta.mm=0.0001,beta.mp=0.0,gamma.m=0.05,alpha.m=0.01,
          sigma.m=0.05,kappa=0)
ngm1.11 <- ngm1.1(params[["beta.mm"]], params[["alpha.m"]], params[["gamma.m"]])
ngm1.22 <- ngm1.1(params[["beta.pp"]], params[["alpha.p"]], params[["gamma.p"]])
ngm1.12 <- ngm1.2(params[["beta.pm"]],params[["edge"]])
ngm1.21 <- ngm1.2(params[["beta.mp"]],params[["edge"]])

ngm_f1<-matrix(c(ngm1.11,ngm1.21,ngm1.12,ngm1.22), nrow=2,ncol=2)
eigen(ngm_f1)

initial.values<-c(S.p=20,I.p=1,R.p=0,S.m=100,I.m=0,R.m=0)


R0= (params[["beta.pp"]]*S.p*I.p/(params[["alpha.p"]]+params[["gamma.p"]]))

f=seq(0,1,by=0.01)
edge=(1+cos(f*(pi*3/2)-2))
plot(f, edge, bty='n', type='l', col='darkgreen', lwd=4,
     xlim=c(0,1), xlab= 'proportion forested', 
     ylim=c(0,2), ylab= expression(paste("forested edge ( ", epsilon, " )")))

