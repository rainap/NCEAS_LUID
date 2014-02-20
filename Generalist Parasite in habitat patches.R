######### 
# BASIC CODE FOR A GENERALIST PATHOGEN; 3 HOSTS
########
rm(list = ls())

##  Set range of pathogen Ros
Amax = 100 # maximum number of unique patches
Area=numeric(Amax)
SpecAcc=numeric(Amax)

#  Set a range of Area values
Area[1]=1.0
for(i in 1: (Amax/10) ) { 
  Area[i+1]=Area[i]+1
}
for(i in (Amax/10):(Amax/4)) { 
  Area[i+1]=Area[i]+2
}
for(i in (Amax/4):(Amax/2)) { 
  Area[i+1]=Area[i]+5
}    
for(i in (Amax/2):(Amax*3/4)) { 
  Area[i+1]=Area[i]+10
}
for(i in (Amax*3/4):(Amax-1)) { 
  Area[i+1]=Area[i]+20
}
Area=Area/40

#Hosts
Hosts=3  #number of host species

Mass=numeric(Hosts)
# Set a range of body masses in Kg
Mass[1:3]=c(0.1, 10.0,100.0)

#Empty vectors for parameters
Dens = numeric(Hosts) #assume each host has a unique density that is fixed
betaD = numeric(Hosts) #density dependent transmission; each host has a unique rate that is fixed
betaF = numeric(Hosts) #frequency dependent transmission; each host has a unique rate that is fixed
r = numeric(Hosts) #reproduction rates for each host are dependent of mass (births - deaths)
d = numeric(Hosts) #death rates for each host; dependent on mass
b = numeric(Hosts)
alphaF = numeric(Hosts) #disease induced mortality for Frequency dependent pathogens
alphaD = numeric(Hosts) #disease induced mortality for density dependent pathogens
delta = numeric(Hosts)
gamma = numeric(Hosts)

#  Density-dependent transmission
md = 26  #26x the background mortality rate
for (i in 1:Hosts) {
  Dens[i]=16.2*Mass[i]^-0.70
  betaD[i]=0.045*md*Mass[i]^0.44 
  r[i]=0.6*Mass[i]^-0.26  # births - deaths
  d[i]=0.4*Mass[i]^-0.27 #  
  b[i]=r[i]+d[i]
  delta[i] = (r[i]/Dens[i])  #real way to represent b-d/N*
  alphaD[i]= d[i]*md
  gamma[i]=0.1*Mass[i]^-0.27
}

#  Frequency-dependent transmission
mf = 1.12 #1.12 
for (i in 1:Hosts) {
  betaF[i]=1.3*mf*Mass[i]^-0.26
  alphaF[i]= d[i]*mf
}

#Calculate carrying capacity for a patch
PatchK=matrix(NA, ncol = Amax, nrow = Hosts) #Number of Host individuals in a given species+patch; carrycapacity for patch
for(i in 1:Amax)   {
  SpecAcc[i]=0
  for(j in 1:Hosts) {
    PatchK[j,i]=Dens[j]*Area[i] #fix scales
    if (PatchK[j,i]>1) PatchK[j,i] else PatchK[j,i]=0
    if (PatchK[j,i]>1) SpecAcc[i]=SpecAcc[i]+1 else SpecAcc[i]
  }
}  

#ODEs
library(deSolve)
state<- c(S1=PatchK[1,50],
          I1=1,
          S2=PatchK[2,50],
          I2=1,
          S3=PatchK[3,50],
          I3=1 
            )
parameters <- c(betaF,
                b,
                d,
                alphaF,
                gamma)


threehostpatch <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS1<- b[1]*(S1+I1)-(betaF[1]*S1*I1+betaF[2]*S1*I2+betaF[3]*S1*I3)-d[1]*S1+ gamma[1]*I1
    dI1<- (betaF[1]*S1*I1+betaF[2]*S1*I2+betaF[3]*S1*I3) - (d[1]+alphaF[1])*I1- gamma[1]*I1
    dS2<- b[2]*(S2+I2)-(betaF[1]*S2*I1+betaF[2]*S2*I2+betaF[3]*S2*I3)-d[1]*S1+ gamma[2]*I2
    dI2<- (betaF[1]*S2*I1+betaF[2]*S2*I2+betaF[3]*S2*I3) - (d[2]+alphaF[2])*I2- gamma[2]*I2
    dS3<- b[3]*(S3+I3)-(betaF[1]*S3*I1+betaF[2]*S3*I2+betaF[3]*S3*I3)-d[3]*S1+ gamma[3]*I3
    dI3<- (betaF[1]*S3*I1+betaF[2]*S3*I2+betaF[3]*S3*I3) - (d[3]+alphaF[3])*I3- gamma[3]*I3
    list(c(dS1, dI1, dS2, dI2, dS3, dI3))
  }) 
}

### output of times ###
times <- seq(0, 365, by = 0.01)

out <- as.data.frame(lsoda(y= state, times = times, func = threehostpatch, parms = parameters))
out$time = NULL

#GRAPHS

par(mar=c(4,4,2,2))
layout(matrix(1:6,2,3))
layout.show(6)

plot(times, out$S1, type= 'l', col="peru", 
      xlab = 'days',main= 'Species 1', bty="n", xlim=c(0,10))
lines(times, out$I1, type= 'l', col="peru", lty=2)
plot(times, out$S2, type= 'l', col="darkgreen", xlim=c(0,10), 
     xlab = 'days',main= 'Species 2', bty="n")
lines(times, out$I2, type= 'l', col="darkgreen", lty=2)
plot(times, out$S3, type= 'l', col="navy", 
     xlab = 'days',main= 'Species 3', bty="n")
lines(times, out$I3, type= 'l', col="navy", lty=2)
lines(times, out$S3+out$I3, pch=20, col="navy")
