######### 
# BASIC CODE FOR A GENERALIST PATHOGEN; 
# varying within and between host transmission
# c faust created sept 24 for sesync
rm(list = ls())

#Libraries
library(deSolve)
library(lattice)
library(rgl)

#Hosts
Hosts=3  #number of host species

Mass=numeric(Hosts)
Mass[1:3]=c(0.1, 1.0,10.0) # Set a range of body masses in Kg

Dens = numeric(Hosts) #assume each host has a unique density that is fixed
betaD = numeric(Hosts) #density dependent transmission; each host has a unique rate that is fixed
betaF = numeric(Hosts) #frequency dependent transmission; each host has a unique rate that is fixed
r = numeric(Hosts) #reproduction rates for each host are dependent of mass (births - deaths)
d = numeric(Hosts) #death rates for each host; dependent on mass
b = numeric(Hosts) #birth rates r+d
alphaF = numeric(Hosts) #disease induced mortality for Frequency dependent pathogens
alphaD = numeric(Hosts) #disease induced mortality for density dependent pathogens
delta = numeric(Hosts)
gamma = numeric(Hosts) #recovery rate

#slider for within and between host 
w = 0.1

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
 

####
Area=1
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
                w,
                alphaD,
                gamma,
                Dens,
                Area)

threehostpatchD <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
  dS1 <- b[1]*(S1+I1)*(1-((S1+I1)/(Dens[1]*Area))) - ((betaD[1]*S1*I1) + (w*betaD[2]*S1*I2) + (w*betaD[3]*S1*I3)) - d[1]*S1 + gamma[1]*I1
  dI1 <- ((betaD[1]*S1*I1) + (w*betaD[2]*S1*I2) + (w*betaD[3]*S1*I3)) - (d[1]+alphaD[1]+gamma[1])*I1
  dS2 <- b[2]*(S2+I2)*(1-((S2+I2)/(Dens[2]*Area))) - ((w*betaD[1]*S2*I1) + (betaD[2]*S2*I2) + (w*betaD[3]*S2*I3)) - d[2]*S2 + gamma[2]*I2
  dI2 <- ((w*betaD[1]*S2*I1) + (betaD[2]*S2*I2) + (w*betaD[3]*S2*I3)) - (d[2]+alphaD[2]+gamma[2])*I2
  dS3 <- b[3]*(S3+I3)*(1-((S3+I3)/(Dens[1]*Area))) - ((w*betaD[1]*S3*I1) +  (w*betaD[2]*S3*I2) + (betaD[3]*S3*I3)) - d[3]*S3 + gamma[3]*I3
  dI3 <- ((w*betaD[1]*S3*I1) +  (w*betaD[2]*S3*I2) + (betaD[3]*S3*I3)) - (d[3]+alphaD[3]+gamma[3])*I3
  list(c(dS1, dI1, dS2, dI2, dS3, dI3))
  })
}

### output of times ###
times <- seq(0, 100, by = 1)

out <- as.data.frame(ode(y= state, times = times, func = threehostpatchD, parms = parameters,method="ode45" ))
out$time = NULL

plot(times, out$S1, type='l',  xlab = 'years', main= 'Species 1', bty='n', col="darkgreen", 
     ylim= c(0,Dens[1]*Area), ylab= 'individuals')
lines(times, out$I1, type='l', lty= 3, col="darkgreen")
#lines(times, out$S1+ out$I1, col= 'black')
legend("topright", c("susceptibles", "infected"), col="darkgreen", lty= c(1,3), bty='n')
lines(times, out$S2, type= 'l', lty=1, col="peru", xlab = 'days', 
      main= 'Species2', bty='n',
     ylab= 'individuals')
lines(times, out$I2, type= 'l', lty=3, col="peru")
#legend("topright", c("susceptibles", "infected"), col="peru", lty= c(1,3), bty='n')
lines(times, out$S3, type= 'l', col="navy", xlab = 'days',  main= 'Species3', bty='n')
lines(times, out$I3, type= 'l', lty=3, col="navy")
#legend("topright", c("susceptibles", "infected"), col="navy", lty= c(1,3), bty='n')

plot(times, out$I1/(out$S1+out$I1), type='l',col="darkgreen", 
     ylab="prevalence", bty='n', ylim=c(0,1))
lines(times, out$I2/(out$S2+out$I2), type= 'l', col="darkorange3")
lines(times,out$I3/(out$S3+out$I3), type='l', col= "navy")
legend("topright", c("sp1", "sp2", "sp3"), bty='n', lty=1, col=c("darkgreen","darkorange3", "navy"))


#########
#code to find equilibrium in each patch size
### Calculating state variables at time = 40 ####
Amax = 100 # maximum number of unique patches
Area=numeric(Amax) #size of each patch
SpecAcc=numeric(Amax) # number of species in each patch
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
w=seq(0,1, by=0.1)
times <- seq(0, 50, by = 1)
# Species1 #
storeS1 = as.data.frame(matrix(NA,length(w),length(Area)));
names(storeS1)=Area
storeI1 = as.data.frame(matrix(NA,length(w),length(Area)));
names(storeI1)=Area
storeS2 = as.data.frame(matrix(NA,length(w),length(Area)));
names(storeS2)=Area
storeI2 = as.data.frame(matrix(NA,length(w),length(Area)));
names(storeI2)=Area
storeS3 = as.data.frame(matrix(NA,length(w),length(Area)));
names(storeS3)=Area
storeI3 = as.data.frame(matrix(NA,length(w),length(Area)));
names(storeI3)=Area

for (i in Area) {
  parameters["Area"] = i;
  for (j in w){
    parameters["w"] = j;
    out <- as.data.frame(ode(y= state, times = times, func = threehostpatchD, parms = parameters, method="ode45"));
    storeS1[j,i] = out$S1[length(times)];
    storeI1[j,i] = out$I1[length(times)];
    storeS2[j,i] = out$S2[length(times)];
    storeI2[j,i] = out$I2[length(times)];
    storeS3[j,i] = out$S3[length(times)];
    storeI3[j,i] = out$I3[length(times)];
  }
}

plot(w, storeS1+storeI1, type='p',  col= "gray40", pch=20, cex=0.5,
     bty='n', ylim=c(0, max(storeS1+storeI1)), 
     ylab='population', xlab="between host probability of transmission",
     main='host pop. dependent on between species transmission')
points(w, storeS2+storeI2, col= "navy", pch=20, cex=0.5)
points(w, storeS3+storeI3, col= "forestgreen", pch=20, cex=0.5)

plot(w, storeI1/(storeS1+storeI1), type='p',  col= "gray40", pch=20, cex=0.5,
     bty='n', ylim=c(0, 1), 
     ylab='population', xlab="between host probability of transmission",
     main='prevalence of infection')
points(w, storeI2/(storeS2+storeI2), col= "navy", pch=20, cex=0.5)
points(w, storeI3/(storeS3+storeI3), col= "forestgreen", pch=20, cex=0.5)


plot(Area, (storeS1+storeI1)/(Area*Dens[1]), type = "l", 
     col= "gray40", xlab= "patch size", ylab= "proportion of carrying capacity",
     ylim=c(0, 1), bty='n', cex=0.4,pch=19, log='x', main="populations of hosts across patch size")
lines(Area, (storeS2+storeI2)/(Area*Dens[2]), type = "l", 
     col= "navy")
lines(Area, (storeS3+storeI3)/(Area*Dens[3]), type = "l", 
      col= "forestgreen")
legend("topright", c( "sp1", "sp2","sp3"),
       col=c('darkgrey','navy','forestgreen'), 
       lty=1, bty='n') 

plot(Area, storeS1, type = "l", 
     col= "gray40", xlab= "patch size", ylab= "individuals",
     bty='n', cex=0.4,pch=19, log='x', ylim=c(0,max(storeS1)))
lines(Area, storeI1, lty=3, col= "gray40")
lines(Area, storeS2,  col= "navy")
lines(Area, storeI2, lty=3, col= "navy")
lines(Area, storeS3,  col= "forestgreen")
lines(Area, storeI3, lty=3, col= "forestgreen")
legend("topleft", c( "sp1", "sp2","sp3"),
       col=c('darkgrey','navy','forestgreen'), 
       lty=1, bty='n') 


plot(Area,storeI1/ (storeS1 + storeI1), type= 'l', 
     col="gray40", ylim= c(0, 1.0), ylab='proportion infected', 
     xlab = 'area', bty='n', log='x')
lines(Area,storeI2/(storeS2 + storeI2),  col="navy")
lines(Area, storeI3/(storeS3 + storeI3), col="forestgreen")

