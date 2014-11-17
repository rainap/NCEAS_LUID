######### 
# Plasmodium knowlesi Patch-Matrix model with slider for patch to matrix transmission
# Christina Faust Version October 8 2014
########
rm(list = ls())


library(deSolve)
library(manipulate)
library(graphics)
# f is the proportion forested (ranges from 1 to 0)
# 1. scaled by edge; edge= 1+cos(f*(pi*3/2)-2.5); ranges from 0 to 2 across values of f; symmetrical)
#Forest edge varying by deforestation; can change functional response
f=seq(0,1,by=0.01)
epsilon=(1+cos(f*(pi*3/2)-1.5))
epsilon=2- ((3*f^3)/(1+f^(4)))
epsilon= numeric(length(f))
for (i in 1:length(f)){
  if(f[i]<0.2){
    epsilon[i]=10*f[i]
  }else{
    epsilon[i]=1/(2.5*f[i])
  } 
}
epsilon=(1.5+cos(3-f*(pi*8/2)))/1.5

#R0 for equilibrium
#R0 at different f; just by calculating a pathogen that comes into a completely naive population at carrying capacity for that proportion forested
f=seq(0.01,1,by=0.01) #vector for different forested proportions
params<-c(d.p=0.1,k.p=0.007,beta.pp=0.15,beta.pm=0.15,gamma.p=0.03,alpha.p=0.001,sigma.p=0.05,
          d.m=0.02,k.m=0.006,beta.mm=0.05,beta.mp=0.1,gamma.m=0.05,alpha.m=0.01,sigma.m=0.05) 

epsilon=(1+cos(f*(pi*3/2)-2)) #vector of edges at different forested levels

ngm.within.fd<-function(beta, alpha,sigma,d){
  z = (beta)/(alpha+sigma+d)
  return(z)
}
ngm.between.fd<-function(epsilon, beta, alpha, sigma, d){
  z = (epsilon*beta)/(alpha+sigma+d)
  return(z)
}
f=seq(0.01,1,by=0.01) #vector for different forested proportions
params<-c(b.p=0.3,d.p=0.1,k.p=0.007,beta.pp=0.2,beta.pm=0.1,gamma.p=0.03,alpha.p=0.001,sigma.p=0.05,
          b.m=0.05, d.m=0.02,k.m=0.006,beta.mm=0.002,beta.mp=0.1,gamma.m=0.05,alpha.m=0.01,sigma.m=0.05)

epsilon=(1+cos(f*(pi*3/2)-2)) #vector of edges at different forested levels
#R0; plotting within and between host R0
R0.p.fd= numeric(length(f)) #assuming no matrix species present (or not transmissible)
R0.m.fd= numeric(length(f)) #assuming no patch species present (or not transmissible)
R0.combin.fd= numeric(length(f))
ngm.fd<-matrix(NA,nrow=2,ncol=2) #next generation matrix 2x2 for 2 species

for (i in 1:length(f)){
  ngm1.11 <- ngm.within.fd(params[["beta.mm"]],params[["alpha.m"]],params[["sigma.m"]], params[["d.m"]])
  ngm1.22 <- ngm.within.fd(params[["beta.pp"]],params[["alpha.p"]],params[["sigma.p"]], params[["d.p"]])
  ngm1.12 <- ngm.between.fd(epsilon[[i]],params[["beta.mp"]], params[["alpha.m"]],params[["sigma.m"]], params[["d.m"]])
  ngm1.21 <- ngm.between.fd(epsilon[[i]],params[["beta.pm"]], params[["alpha.p"]],params[["sigma.p"]], params[["d.p"]])
  ngm.fd<-matrix(c(ngm1.11,ngm1.21,ngm1.12,ngm1.22), nrow=2,ncol=2)
  eigen.ngm.fd<-eigen(ngm.fd)
  R0.combin.fd[i] = eigen.ngm.fd$values[[1]]
  R0.p.fd[i] = (params[["beta.pp"]]/(params[["alpha.p"]]+params[["gamma.p"]]+params[["sigma.p"]]))
  R0.m.fd[i] = (params[["beta.mm"]]/(params[["alpha.m"]]+params[["gamma.m"]]+params[["sigma.m"]]))
}
par(bg='white') 
plot(f, R0.combin.fd, type="l",lwd=4, 
     ylim=c(0,5),bty='n', 
     ylab="R0",xlab="proportion forested")
lines(f, R0.m.fd, lwd=4, col='darkgoldenrod')
lines(f, R0.p.fd, lwd=4, col='forestgreen')
abline(h=1, lty=2, lwd=4, col='darkred')
#lines(f, 7.5*epsilon, lty=3, lwd=1, col='grey44')
legend("topleft", c("combined" , "within matrix" ,"within patch"),
       lty=c(1,1,1,3), lwd=3, col=c("black", "darkgoldenrod", "forestgreen"), bty='n' )


#ODES for dynamical model
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
params<-c(d.p=0.1,k.p=0.007,beta.pp=0.2,beta.pm=0.1,gamma.p=0.03,alpha.p=0.001,sigma.p=0.05,
          d.m=0.02,k.m=0.006,beta.mm=0.02,beta.mp=0.1,gamma.m=0.05,alpha.m=0.01,sigma.m=0.05) 
#iterations
times <- seq(0, 200, by = 0.01)
plot.fun<-function(f){
    params<-c(b.p=0.3,d.p=0.1,k.p=0.007,beta.pp=0.2,beta.pm=0.1,gamma.p=0.03,alpha.p=0.001,sigma.p=0.05,
              b.m=0.05, d.m=0.02,k.m=0.006,beta.mm=0.002,beta.mp=0.1,gamma.m=0.05,alpha.m=0.01,sigma.m=0.05,
              f=f,epsilon=(1+cos(f*(pi*3/2)-2)),kappa=1) 
    initial.values<-c(S.p=20,I.p=2,R.p=0,S.m=100,I.m=0,R.m=0)
    
    print(system.time(
      out<-ode(func=patch.matrix.model.f,y=initial.values,parms=params,times=times))) 
    head(out,n=4)
    par(mfrow=c(2,2))
    matplot(out[,"time"],out[,2:4],type="l",xlab="time",ylab="number", bty='n',cex=0.8,
            main="Patch Hosts",lwd=2, ylim=c(0,100),col=c("black","darkred","forestgreen"))
    legend("topright",c("susc","inf","rec"),col=c("black","darkred","forestgreen"),pch=20,bty='n',cex=0.7)
    matplot(out[,"time"],(out[,3]/(out[,2]+out[,3]+out[,4])),type="l",xlab="time",ylab="number", bty='n',
            main="Prop. Infected Patch Hosts",lwd=2, ylim=c(0,1),cex=0.8) 
    matplot(out[,"time"],out[,5:7],type="l",xlab="time",ylab="number", bty='n',
            main="Matrix Hosts",lwd=2, ylim=c(0,200),cex=0.8,col=c("black","darkred","forestgreen"))
    legend("topright",c("susc","inf","rec"),col=c("black","darkred","forestgreen"),pch=20, bty='n',cex=0.7)
    matplot(out[,"time"],(out[,6]/(out[,5]+out[,6]+out[,7])),type="l",xlab="time",ylab="number", bty='n',
            main="Prop. Infected Matrix Hosts",lwd=2, ylim=c(0,1),cex=0.8)    
}

manipulate(plot.fun(f),f=slider(0.0000001,1.0)) #can't divide by zero
f<-slider(0.0000001,1.0)

