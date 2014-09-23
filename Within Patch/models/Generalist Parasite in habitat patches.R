######### 
# BASIC CODE FOR A GENERALIST PATHOGEN; 3 HOSTS
########
rm(list = ls())
library(deSolve)
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
Mass[1:3]=c(0.1, 1.0,10.0)

#Empty vectors for parameters
Dens = numeric(Hosts) #assume each host has a unique density that is fixed
betaD = numeric(Hosts) #density dependent transmission; each host has a unique rate that is fixed
betaF = numeric(Hosts) #frequency dependent transmission; each host has a unique rate that is fixed
r = numeric(Hosts) #reproduction rates for each host are dependent of mass (births - deaths)
d = numeric(Hosts) #death rates for each host; dependent on mass
b = numeric(Hosts) #birth rates r+d
alphaF = numeric(Hosts) #disease induced mortality for Frequency dependent pathogens
alphaD = numeric(Hosts) #disease induced mortality for density dependent pathogens
delta = numeric(Hosts)
gamma = numeric(Hosts)

#  Density-dependent transmission
md = 26  #26x the background mortality rate
for (i in 1:Hosts) {
  Dens[i]=16.2*Mass[i]^-0.70
  betaD[i]=0.0045*md*Mass[i]^0.44 
  r[i]=0.6*Mass[i]^-0.26  # births - deaths
  d[i]=0.4*Mass[i]^-0.27 # deaths
  b[i]=r[i]+d[i]
  delta[i] = (r[i]/Dens[i])  #real way to represent b-d/N*
  alphaD[i]= d[i]*md # disease induced mortality
  gamma[i]=0.05*Mass[i]^-0.27 #
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
# library(deSolve)
# state<- c(S1=PatchK[1,50],
#           I1=1,
#           S2=PatchK[2,50],
#           I2=1,
#           S3=PatchK[3,50],
#           I3=1 
#             )
# parameters <- c(betaF,
#                 b,
#                 d,
#                 alphaF,
#                 gamma)
# 
# threehostpatchF <- function(t, state, parameters) {
#   with(as.list(c(state, parameters)), {
#     dS1 <- b[1]*(S1+I1)-(betaF[1]*S1*I1)/(S1*I1)+(betaF[2]*S1*I2)/(S1+I2) + (betaF[3]*S1*I3)/(S1+I3))-d[1]*S1+gamma[1]*I1
#     dI1<- ((betaF[1]*S1*I1)/(S1+I1+)+(betaF[2]*S1*I2)/(S1+I2) + (betaF[3]*S1*I3)/(S1+I3)) - (d[1]+alphaF[1]+gamma[1])*I1
#     dS2 <- b[2]*(S2+I2)-((betaF[1]*S2*I1)/(S2*I1)+(betaF[2]*S2*I2)/(S2+I2)+(betaF[3]*S2*I3)/(S2+I3))-d[2]*S2+gamma[2]*I2
#     dI2<- ((betaF[1]*S2*I1)/(S2+I1)+(betaF[2]*S2*I2)/(S2+I2)+ (betaF[3]*S2*I3))/(S2+I3) - (d[2]+alphaF[2]+gamma[2])*I2
#     dS3 <- b[3]*(S3+I3)-((betaF[1]*S3*I1)/(S3+I1)+(betaF[2]*S3*I2)/(S3+I2)+(betaF[3]*S3*I3)/(S3+I3))-d[3]*S3+gamma[3]*I3
#     dI3<- ((betaF[1]*S3*I1)/(S3+I1)+(betaF[2]*S3*I2)/(S3+I2)+ (betaF[3]*S3*I3)/(S3+I3)) - (d[3]+alphaF[3]+gamma[3])*I3
#     list(c(dS1, dI1, dS2, dI2, dS3, dI3))
#   })
# }
# 
# ### output of times ###
# times <- seq(0, 10, by = 0.1)
# 
# out <- as.data.frame(lsoda(y= state, times = times, func = threehostpatchF, parms = parameters))
# out$time = NULL
# 
# 
# par(mar=c(4,4,2,2))
# layout(matrix(1:3,1, 3))
# layout.show(6)
# 
# plot(times, out$S1, type='l',  xlab = 'years', main= 'Species 1', bty='n', col="darkgreen")
# lines(times, out$I1, type='l', lty= 3, col="darkgreen")
# 
# plot(times, out$S2, type= 'l', lty=1, col="peru", xlab = 'days', ylim= c(0,30), main= 'Species2', bty='n')
# lines(times, out$I2, type= 'l', lty=3, col="darkorange3")
# 
# plot(times, out$S3, type= 'l', col="navy", xlab = 'days', ylim= c(0,30), main= 'Species3', bty='n')
# lines(times, out$I3, type= 'l', lty=3, col="navy")
# 
# plot(times, out$I1/(out$S1+out$I1), type='l',col="darkgreen", main="Prevalence", bty='n', ylim=c(0,1))
# lines(times, out$I2/(out$S2+out$I2), type= 'l', lty=3, col="darkorange3")
# lines(times,out$I3/(out$S3+out$I3), type='l', lty=3, col= "navy")
# 


Area=1
####
library(deSolve)
state<- c(S1=if (Dens[1]*Area>1) Dens[1]*Area else 0,
          I1=1,
          S2=if (Dens[2]*Area>1) Dens[2]*Area else 0,
          I2=0,
          S3=if (Dens[3]*Area>1) Dens[1]*Area else 0,
          I3=0 
)
parameters <- c(betaD,
                b,
                d,
                alphaD,
                gamma,
                Dens,
                Area)

threehostpatchD <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
  dS1 <- b[1]*(S1+I1)*(1-((S1+I1)/(Dens[1]*Area))) - ((betaD[1]*S1*I1) + (betaD[2]*S1*I2) + (betaD[3]*S1*I3)) - d[1]*S1 + gamma[1]*I1
  dI1 <- ((betaD[1]*S1*I1) + (betaD[2]*S1*I2) + (betaD[3]*S1*I3)) - (d[1]+alphaD[1]+gamma[1])*I1
  dS2 <- b[2]*(S2+I2)*(1-((S2+I2)/(Dens[2]*Area))) - ((betaD[1]*S2*I1) + (betaD[2]*S2*I2) + (betaD[3]*S2*I3)) - d[2]*S2 + gamma[2]*I2
  dI2 <- ((betaD[1]*S2*I1) + (betaD[2]*S2*I2) + (betaD[3]*S2*I3)) - (d[2]+alphaD[2]+gamma[2])*I2
  dS3 <- b[3]*(S3+I3)*(1-((S3+I3)/(Dens[1]*Area))) - ((betaD[1]*S3*I1) +  (betaD[2]*S3*I2) + (betaD[3]*S3*I3)) - d[3]*S3 + gamma[3]*I3
  dI3 <- ((betaD[1]*S3*I1) +  (betaD[2]*S3*I2) + (betaD[3]*S3*I3)) - (d[3]+alphaD[3]+gamma[3])*I3
  list(c(dS1, dI1, dS2, dI2, dS3, dI3))
  })
}

### output of times ###
times <- seq(0, 10, by = 0.5)

out <- as.data.frame(lsoda(y= state, times = times, func = threehostpatchD, parms = parameters))
out$time = NULL

# par(mar=c(4,4,2,2))
# layout(matrix(1:1,1, 1))
# layout.show(6)

plot(times, out$S1, type='l',  xlab = 'years', main= 'Species 1', bty='n', col="darkgreen", 
     ylim= c(0,Dens[1]*Area), ylab= 'individuals')
lines(times, out$I1, type='l', lty= 3, col="darkgreen")
#lines(times, out$S1+ out$I1, col= 'black')
legend("topright", c("susceptibles", "infected"), col="darkgreen", lty= c(1,3), bty='n')

plot(times, out$S2, type= 'l', lty=1, col="peru", xlab = 'days', 
     ylim= c(0,15), main= 'Species2', bty='n',
     ylab= 'individuals')
lines(times, out$I2, type= 'l', lty=3, col="peru")
legend("topright", c("susceptibles", "infected"), col="peru", lty= c(1,3), bty='n')

plot(times, out$S3, type= 'l', col="navy", xlab = 'days', ylim= c(0,4), main= 'Species3', bty='n')
lines(times, out$I3, type= 'l', lty=3, col="navy")
abline(h=PatchK[3,50])
legend("topright", c("susceptibles", "infected"), col="navy", lty= c(1,3), bty='n')

plot(times, out$I1/(out$S1+out$I1), type='l',col="darkgreen", 
     ylab="prevalence", bty='n', ylim=c(0,1))
lines(times, out$I2/(out$S2+out$I2), type= 'l', col="darkorange3")
lines(times,out$I3/(out$S3+out$I3), type='l', col= "navy")
legend("topright", c("sp1", "sp2", "sp3"), bty='n', lty=1, col=c("darkgreen","darkorange3", "navy"))


#########
#code to find equilibrium in each patch size
### Calculating state variables at time = 40 ####
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
Area=Area
times <- seq(0, 10, by = 0.1)
# Species1 #
storeS1 = {};
storeI1 = {};
storeS2 = {};
storeI2 = {};
storeS3 = {};
storeI3 = {};
for (i in Area) {
  parameters["Area"] = i;
  out <- as.data.frame(ode(y= state, times = times, func = threehostpatchD, parms = parameters));
  storeS1 = c(storeS1, out$S1[100]);
  storeI1 = c(storeI1, out$I1[100]);
  storeS2 = c(storeS2, out$S2[100]);
  storeI2 = c(storeI2, out$I2[100]);
  storeS3 = c(storeS3, out$S3[100]);
  storeI3 = c(storeI3, out$I3[100]);
}

plot(Area, (storeS1+storeI1)/(Area*Dens[1]), type = "l", 
     col= "gray40", xlab= "patch size", ylab= "proportion of carrying capacity",
     ylim=c(0, 1), bty='n', cex=0.4,pch=19)
lines(Area, (storeS2+storeI2)/(Area*Dens[2]), type = "l", 
     col= "navy")
lines(Area, (storeS3+storeI3)/(Area*Dens[3]), type = "l", 
      col= "forestgreen")
legend("topright", c( "sp1", "sp2","sp3"),
       col=c('darkgrey','navy','forestgreen'), 
       lty=1, bty='n') 

plot(Area, storeS1, type = "l", col='grey40', bty='n', ylim=c(0,max(storeS1)))
lines(Area,storeS2, col='navy')    
lines(Area, storeS3, col='forestgreen')
lines(Area,storeI1, col="darkred")

plot(Area,storeI1/ (storeS1 + storeI1), type= 'l', 
     col="gray40", ylim= c(0, 1.0), ylab='proportion infected', 
     xlab = 'area', bty='n')
lines(Area,storeI2/(storeS2 + storeI2),  col="navy")
lines(Area, storeI3/(storeS3 + storeI3), col="forestgreen")
