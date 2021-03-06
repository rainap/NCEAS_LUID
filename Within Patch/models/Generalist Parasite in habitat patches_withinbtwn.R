######### 
# CODE FOR A GENERALIST PATHOGEN; 
# varying within and between host transmission
# c faust created sept 24 for sesync land use, economics, & infectious disease workshop
# based on allometrically scaled model (see others in github repo)
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


#### Input for DD model
Area=1
state<- c(S1=if (Dens[1]*Area>1) Dens[1]*Area else 0,
          I1=1,
          S2=if (Dens[2]*Area>1) Dens[2]*Area else 0,
          I2=0,
          S3=if (Dens[3]*Area>1) Dens[1]*Area else 0,
          I3=0 
)
parametersD <- c(betaD,
                b,
                d,
                w,
                alphaD,
                gamma,
                Dens,
                Area)

threehostpatchD <- function(t, state, parameters) {
  with(as.list(c(state, parametersD)), {
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

out <- as.data.frame(ode(y= state, times = times, func = threehostpatchD, parms = parametersD,method="ode45" ))
out$time = NULL

plot(times, out$S1, type='l',  xlab = 'years', main= 'Species 1', bty='n', col="darkgreen", 
      ylim= c(0,Dens[1]*Area), ylab= 'individuals')
# lines(times, out$I1, type='l', lty= 3, col="darkgreen")
# #lines(times, out$S1+ out$I1, col= 'black')
# legend("topright", c("susceptibles", "infected"), col="darkgreen", lty= c(1,3), bty='n')
# lines(times, out$S2, type= 'l', lty=1, col="peru", xlab = 'days', 
#       main= 'Species2', bty='n',
#      ylab= 'individuals')
# lines(times, out$I2, type= 'l', lty=3, col="peru")
# #legend("topright", c("susceptibles", "infected"), col="peru", lty= c(1,3), bty='n')
# lines(times, out$S3, type= 'l', col="navy", xlab = 'days',  main= 'Species3', bty='n')
# lines(times, out$I3, type= 'l', lty=3, col="navy")
# #legend("topright", c("susceptibles", "infected"), col="navy", lty= c(1,3), bty='n')
# 
# plot(times, out$I1/(out$S1+out$I1), type='l',col="darkgreen", 
#      ylab="prevalence", bty='n', ylim=c(0,1))
# lines(times, out$I2/(out$S2+out$I2), type= 'l', col="darkorange3")
# lines(times,out$I3/(out$S3+out$I3), type='l', col= "navy")
# legend("topright", c("sp1", "sp2", "sp3"), bty='n', lty=1, col=c("darkgreen","darkorange3", "navy"))


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

storeS1 = {};
storeI1 = {};
storeS2 = {};
storeI2 = {};
storeS3 = {};
storeI3 = {};
plotArea = {};
plotw= {};

for (i in Area) {
  parametersD["Area"] = i;
  for (j in w){
    parametersD["w"] = j;
    out <- as.data.frame(ode(y= state, times = times, func = threehostpatchD, parms = parametersD, method="ode45"));
    storeS1 = c(storeS1, out$S1[length(times)]);
    storeI1 = c(storeI1, out$I1[length(times)]);
    storeS2 = c(storeS2, out$S2[length(times)]);
    storeI2 = c(storeI2, out$I2[length(times)]);
    storeS3 = c(storeS3, out$S3[length(times)]);
    storeI3 = c(storeI3, out$I3[length(times)]);
    plotArea = c(plotArea, i)
    plotw= c(plotw, j)
  }
}
# 
# col.l <- colorRampPalette(c('blue', 'cyan','green','yellow','orange','red'))
# wireframe(storeS1 ~ plotArea * plotw, main="Host1-smallsize-.5kg",
#           drape = TRUE,
#           #perspective = FALSE, 
#           #aspect = c(10,10), 
#           #colorkey = TRUE, pretty=T,
#           #col.regions= c(   'blue',  'cyan')
#           #form = "full", col.regions=col.l
#             )
# wireframe(storeS2 ~ plotArea * plotw, main="Host2-midsize-1kg",
#           drape= T , zlim=c(-1000,max(storeS2,na.rm=T)))
# wireframe(storeS3 ~ plotArea * plotw, main="Host3-large-5kg",
#           drape= T)
# wireframe((storeI1+storeI2+storeI3) ~ plotArea * plotw, 
#           drape= T, main="total infecteds",
#           zlab='number of infecteds', xlab= "patch size", ylab="between species trans")
# theseCol = colorRampPalette(c("red", "white","blue"))
# wireframe((storeI1/(storeI1+storeS1)) ~ plotArea * plotw, 
#           drape= T, main="prevalence in sp1",
#           #col.regions=theseCol,
#           xlab= "patch size", ylab="between species trans",
#           zlab='prevalence')
# library(graphics)
# 
# filled.contour(plotArea,plotw, storeS2)
#                #xlim = range(x, finite = TRUE),
#                #ylim = range(y, finite = TRUE),
#                #zlim = range(z, finite = TRUE),
#                levels = pretty(storeS2, nlevels), nlevels = 20,
#                color.palette = cm.colors,
#                col = color.palette(length(levels) - 1),
#                plot.title, plot.axes, key.title, key.axes,
#                asp = NA, xaxs = "i", yaxs = "i", las = 1,
#                axes = TRUE, frame.plot = axes, ...)

DD3hostgenw<-as.data.frame(cbind(plotArea,plotw,storeS1,storeI1,storeS2,storeI2,storeS3,storeI3))
write.csv(DD3hostgenw, file = "DD3hostgenw.csv",row.names=FALSE)

wireframe(storeS2 ~ plotArea * plotw, data=DD3hostgenw, 
          main="DDHost S2-midsize-1kg", drape = TRUE)
wireframe(storeI2 ~ plotArea * plotw, data=DD3hostgenw, 
          main="DDHost I2-midsize-1kg", drape = TRUE)
wireframe(storeS3 ~ plotArea * plotw, data=DD3hostgenw,
          main="DDHost S3-large-10kg", drape = TRUE, zlim=c(0,10))
wireframe(storeI3 ~ plotArea * plotw, data=DD3hostgenw,
          main="DDHost I3-large-10kg", drape = TRUE)
wireframe(storeI1 ~ plotArea * plotw,data=DD3hostgenw, 
          main="DDHost I1-small-0.1kg", drape = TRUE)
wireframe(storeS1 ~ plotArea * plotw,data=DD3hostgenw, 
          main="DDHost S1-small-0.1kg",drape = TRUE)

wireframe((storeI1+storeI2+storeI3)/(storeI1+storeS1+storeI2+storeI3+storeS2+storeS3)  ~ plotArea * plotw,
          data=DD3hostgenw, main="prevalence in all", drape = TRUE)
wireframe((storeI1+storeI2+storeI3)~ plotArea * plotw,data=DD3hostgenw, 
          main="number of infected individuals", drape = TRUE)

########
# Frequency dependent transmission

#same state variables as DD

mf = 3 #1.12 
for (i in 1:Hosts) {
  betaF[i]=1.3*mf*Mass[i]^-0.26
  alphaF[i]= d[i]*mf
}
Area=100
parametersF <- c(betaF,
                 b,
                 d,
                 w,
                 alphaF,
                 gamma,
                 Dens,
                 Area)

#ODE for FD transmission

threehostpatchF <- function(t, state, parameters) {
  with(as.list(c(state, parametersF)), {
    dS1 <- b[1]*(S1+I1)*(1-((S1+I1)/(Dens[1]*Area))) - ((betaF[1]*S1*I1/(S1+I1+S2+I2+S3+I3)) + (w*betaF[2]*S1*I2/(S1+I1+S2+I2+S3+I3)) + (w*betaF[3]*S1*I3/(S1+I1+S2+I2+S3+I3))) - d[1]*S1 + gamma[1]*I1
    dI1 <- ((betaF[1]*S1*I1/(S1+I1+S2+I2+S3+I3)) + (w*betaF[2]*S1*I2/(S1+I1+S2+I2+S3+I3)) + (w*betaF[3]*S1*I3/(S1+I1+S2+I2+S3+I3))) - (d[1]+alphaF[1]+gamma[1])*I1
    dS2 <- b[2]*(S2+I2)*(1-((S2+I2)/(Dens[2]*Area))) - ((w*betaF[1]*S2*I1/(S1+I1+S2+I2+S3+I3)) + (betaF[2]*S2*I2/(S1+I1+S2+I2+S3+I3)) + (w*betaF[3]*S2*I3/(S1+I1+S2+I2+S3+I3))) - d[2]*S2 + gamma[2]*I2
    dI2 <- ((w*betaF[1]*S2*I1/(S1+I1+S2+I2+S3+I3)) + (betaF[2]*S2*I2/(S1+I1+S2+I2+S3+I3)) + (w*betaF[3]*S2*I3/(S1+I1+S2+I2+S3+I3))) - (d[2]+alphaF[2]+gamma[2])*I2
    dS3 <- b[3]*(S3+I3)*(1-((S3+I3)/(Dens[1]*Area))) - ((w*betaF[1]*S3*I1/(S1+I1+S2+I2+S3+I3)) +  (w*betaF[2]*S3*I2/(S1+I1+S2+I2+S3+I3)) + (betaF[3]*S3*I3/(S1+I1+S2+I2+S3+I3))) - d[3]*S3 + gamma[3]*I3
    dI3 <- ((w*betaF[1]*S3*I1/(S1+I1+S2+I2+S3+I3)) +  (w*betaF[2]*S3*I2/(S1+I1+S2+I2+S3+I3)) + (betaF[3]*S3*I3/(S1+I1+S2+I2+S3+I3))) - (d[3]+alphaF[3]+gamma[3])*I3
    list(c(dS1, dI1, dS2, dI2, dS3, dI3))
  })
}

### output of times ###
times <- seq(0, 100, by = 1)

out <- as.data.frame(ode(y= state, times = times, func = threehostpatchF, parms = parametersF, method="ode45" ))
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
     
### Looking a frequency dependence across the landscape
#########
#code to find equilibrium in each patch size
### Calculating state variables at time = 40 ####
Amax=100
Area=numeric(Amax)
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

times <- seq(0, 40, by = 1)

storeS1 = {};
storeI1 = {};
storeS2 = {};
storeI2 = {};
storeS3 = {};
storeI3 = {};
plotArea = {};
plotw= {};

for (i in Area) {
  parametersF["Area"] = i;
  for (j in w){
    parametersF["w"] = j;
    out <- as.data.frame(ode(y= state, times = times, func = threehostpatchF, parms = parametersF, method="ode45"));
    storeS1 = c(storeS1, out$S1[length(times)]);
    storeI1 = c(storeI1, out$I1[length(times)]);
    storeS2 = c(storeS2, out$S2[length(times)]);
    storeI2 = c(storeI2, out$I2[length(times)]);
    storeS3 = c(storeS3, out$S3[length(times)]);
    storeI3 = c(storeI3, out$I3[length(times)]);
    plotArea = c(plotArea, i)
    plotw= c(plotw, j)
  }
}


FD3hostgenw<-as.data.frame(cbind(plotArea,plotw,storeS1,storeI1,storeS2,storeI2,storeS3,storeI3))
write.csv(FD3hostgenw, file = "FD3hostgenw.csv",row.names=FALSE)
FD3hostgenw$storeS2

wireframe(storeS2 ~ plotArea * plotw, data=FD3hostgenw, main="Host S2-midsize-1kg")
wireframe(storeI2 ~ plotArea * plotw, data=FD3hostgenw, main="Host I2-midsize-1kg")
wireframe(storeS3 ~ plotArea * plotw, data=FD3hostgenw,main="Host S3-large-10kg")
wireframe(storeI3 ~ plotArea * plotw, data=FD3hostgenw,main="Host I3-large-10kg")
wireframe(storeI1 ~ plotArea * plotw,data=FD3hostgenw, main="Host I1-small-0.1kg")
wireframe(storeS1 ~ plotArea * plotw,data=FD3hostgenw, main="Host S1-small-0.1kg")

wireframe((storeI1+storeI2+storeI3)/(storeI1+storeS1+storeI2+storeI3+storeS2+storeS3)  ~ plotArea * plotw,data=FD3hostgenw, main="prevalence in all")
wireframe((storeI1+storeI2+storeI3)~ plotArea * plotw,data=FD3hostgenw, main="number of infected individuals")
