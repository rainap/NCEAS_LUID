##        Andy Dobson - EEB, Princeton University
##
#############################################################################
##                                                                          #
##   Code for model to examine shrinking habitat and pathogen burdens       #
##                                                                          #
#############################################################################
rm(list = ls())
#  install.package("ggplot2")
library("ggplot2")

##  Set range of pathogen Ros
Amax = 100 # maximum number of unique patches
Area=numeric(Amax)
SpecRich=numeric(Amax)

#  Set a range of Area values
Area[1]=0.1
for(i in 1: (Amax/10) ) { 
  Area[i+1]=Area[i]+0.2
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
print(Area)

#Hosts
Hosts=12  #number of host species

Mass=numeric(Hosts)
# Set a range of body masses in Kg
Mass[1:12]=c(0.05,0.10,0.20,0.333,0.5, 1.0, 2.0, 3.333, 5.0, 10.0, 20.0, 50.0)

Dens = numeric(Hosts) #assume each host has a unique density that is fixed
BetaD = numeric(Hosts) #density dependent transmission; each host has a unique rate that is fixed
BetaF = numeric(Hosts) #frequency dependent transmission; each host has a unique rate that is fixed
r = numeric(Hosts) #reproduction rates for each host are dependent of mass (births - deaths)
d = numeric(Hosts) #death rates for each host; dependent on mass
alphaF = numeric(Hosts) #disease induced mortality for Frequency dependent pathogens
alphaD = numeric(Hosts) #disease induced mortality for density dependent pathogens
delta = numeric(Hosts)
R0F=numeric(Hosts)
R0D= numeric(Hosts)

# Set Carrying capacity for each host species (#/km^2)
#  and also r and d for host transmission rate for pathogen
#  m is increase in host mortality rate

#  Density-dependent transmission
md = 26  #26x the background mortality rate
for (i in 1:Hosts) {
    Dens[i]=16.2*Mass[i]^-0.70
    BetaD[i]=0.045*md*Mass[i]^0.44 
     r[i]=0.6*Mass[i]^-0.26  # births - deaths
     d[i]=0.4*Mass[i]^-0.27 #  
     delta[i] = (r[i]/Dens[i])  #real way to represent b-d/N*
     alphaD[i]= d[i]*md
     R0D[i] = BetaD[i]/(d[i]+alphaD[i])
    }



#  Frequency-dependent transmission

mf = 1.12 #1.12 
for (i in 1:Hosts) {
    BetaF[i]=1.3*mf*Mass[i]^-0.26
    alphaF[i]= d[i]*mf
    R0F[i] = BetaF[i]/(d[i]+alphaF[i])
    }


plot(Mass, R0F, log='x',  bty="n", ylim=c(0,3), ylab='R0')
points(Mass, R0D*Dens, pch=20)
legend('topleft', pch=c(1,20), bty='n',
       c("frequency dependent", "density dependent"))

plot(Mass, r, main="birth, death, delta by body size",
     log='x',  bty="n", ylim=c(0,1.2), xlab='',
     las=1, pch=20, xlim=c(0.1,50), col="red")
points(Mass, d, pch=20)
points(Mass, delta,log='x', pch=20, col="navy")

plot(Mass, Dens, log='x',  bty="n", 
     las=1, pch=20, xlim=c(0.1,50), ylab='density (individuals/km2)')

plot(Mass, BetaD, log='x', 
     main="Beta for hosts, same R0",
     bty="n", las=1, pch=20, ylim=c(0,7), xlim=c(0.1,50))
points(Mass, BetaF)
legend('topleft', pch=c(1,20), bty='n',
       c("frequency dependent", "density dependent"))

##########################################################
#       Determine abundance in each patch
#            and then calculate number infected hosts for frequency dependent pathogens.
#
PatchK=matrix(NA, ncol = Amax, nrow = Hosts) #Number of Host individuals in a given species+patch; carrycapacity for patch
DisPatchKF=matrix(NA, ncol = Amax, nrow = Hosts)
IAbunF=matrix(NA, ncol = Amax, nrow = Hosts) #Number of Infected Individuals
PrevF=matrix(NA, ncol = Amax, nrow = Hosts) # Prevalence 

for(i in 1:Amax)   {
     SpecRich[i]=0
     for(j in 1:Hosts) {
     PatchK[j,i]=Dens[j]*Area[i] #fix scales
     if (PatchK[j,i]>1) PatchK[j,i] else PatchK[j,i]=0
     DisPatchKF[j,i]=((r[j]+d[j]+alphaF[j]-BetaF[j])/delta[j])*Area[i]
     if (DisPatchKF[j,i]>1) DisPatchKF[j,i] else DisPatchKF[j,i]=0
     if (DisPatchKF[j,i]>1) SpecRich[i]=SpecRich[i]+1 else SpecRich[i]
     PrevF[j,i]=(r[j]-(delta[j]*(DisPatchKF[j,i]/Area[i])))/BetaF[j]
     IAbunF[j,i]=PrevF[j,i]*DisPatchKF[j,i]
     }
   }  

SpecRichComp=numeric(Amax)
CompTest=numeric(Hosts)
CompPatchDens=matrix(NA, ncol = Amax, nrow = Hosts)  # carrying capacity in patch with competition
CompAbund=matrix(NA, ncol = Amax, nrow = Hosts)
for(i in 1:Amax)   {
  SpecRichComp[i]=0
  for(j in 2:Hosts) {
    #PatchK[j,i]=Dens[j]*Area[i] #fix scales
    #if (PatchK[j,i]>1) PatchK[j,i] else PatchK[j,i]=0
    CompTest[j]= 1-(Mass[j]-Mass[j-1])^2/(Mass[j]*Mass[j-1])
    if (CompTest[j]>1) CompTest[j] else CompTest[j]=0
    if (CompTest[j]>1) SpecRichComp[i]=SpecRichComp[i]+1 else SpecRichComp[i]
    CompAbund[j,i]=CompAbund[j,i]*Area[i]
  }
}

test=matrix(NA, ncol = 2, nrow = Hosts)  
for (i in 1:Hosts){
  if ((i-1)>0) test[i,1]=(Dens[i]*(1-(Mass[i]-Mass[i-1]))) else test[i,1]=Dens[i]
  if ((i+1)<13) test[i,2]=(Dens[i]*((1-(Mass[i]-Mass[i-1]))+(1-(Mass[i]-Mass[i+1])))) else test[i,2]=Dens[i]
}
test[,1]

i=12
if ((i+1)<13) test[i,2]=(Dens[i]*((1-(Mass[i]-Mass[i-1]))+(1-(Mass[i]-Mass[i+1])))) else test[i,2]=Dens[i]
test[i,2]

1-(abs(Mass[1]+Mass[(1+1)])/(Mass[1]*Mass[(1+1)]))
test[,1]
i=3
if (Mass[i-1]>0) test[i,2]=Dens[i]*(1-(Mass[i]-Mass[i-1])) else test[i,2]=Dens[i]
test[i,2]
# Assuming there 
##   Boring plot of the number of hosts in each patch

par(mfrow=c(1,1))  
plot(Area, PatchK[1,],  ylab='Hosts', xlab='Patch Area(km^2)',
  col ="red", log = c("x", 'y'), #ylim = c(0.1, 300), 
  main = 'Number of individuals in patches of different size',
  pch=20)
points(Area, DisPatchKF[1,], col="red", cex=.5)
points(Area, DisPatchKD[1,], col="red", cex=.5,pch=2)

points(Area, PatchK[2,], col = "blue",pch=20)
points(Area, PatchK[5,], col = "green", pch=20) 
points(Area, PatchK[6,], col= "black", pch=20)
points(Area, PatchK[9,], col= "orange", pch=20) 
leg.text<-c("5gram", "100gram ", "500gram ", "1kg ", "5kg")
legend("topleft",leg.text,lty=rep(1,4),col=c( "red", "blue", "green", "black", "orange"),bty="n")
 

##  Do a quick plot of the  number of infected hosts        
 
plot(Area, IAbunF[1,], ylab='Infected Hosts', xlab='Patch Area(km^2)',
  col ="red", log = "x", #ylim = c(0.1, 100), 
  main = 'Number in patches of different size',
  pch=20)
points(Area, IAbunF[2,], col = "blue", pch=20)
points(Area, IAbunF[5,], col = "green", pch=20) 
points(Area, IAbunF[6,], col= "black", pch=20)
points(Area, IAbunF[9,], col= "orange", pch=20) 
legend("topleft",leg.text,lty=rep(1,4),
       col=c("red", "blue", "green", "black", "orange"),bty="n")
 

## Boring plot of the prevalence data
plot(Area, PrevF[1,], ylab='Prevalence', xlab='Patch Area(km^2)',
  col ="red", log = "x", #ylim = c(0.01, 1.0), 
  main = 'Prevalence in patches of different size',
  pch=20)
points(Area, PrevF[2,], col = "blue", pch=20)
points(Area, PrevF[5,], col = "green", pch=20) 
points(Area, PrevF[6,], col= "black", pch=20)
points(Area, PrevF[9,], col= "orange", pch=20) 
leg.text<-c("5gram", "100gram ", "500gram ", "1kg ", "5kg")
legend("topleft",leg.text,lty=rep(1,4),
       col=c("red", "blue", "green", "black", "orange"),bty="n")



##########################################################
#       DENSITY dependent pathogens
#
DisPatchKD=matrix(NA, ncol = Amax, nrow = Hosts)
IAbunD=matrix(NA, ncol = Amax, nrow = Hosts) #Number of Infected Individuals
PrevD=matrix(NA, ncol = Amax, nrow = Hosts) # Prevalence 

for(i in 1:Amax)   {
  SpecRich[i]=0
  for(j in 1:Hosts) {
    PatchK[j,i]=Dens[j]*Area[i] 
    if (PatchK[j,i]>1) PatchK[j,i] else PatchK[j,i]=0
    DisPatchKD[j,i]=((r[j]+d[j]+alphaD[j])/(BetaD[j]+delta[j]))*Area[i]
    if (DisPatchKD[j,i]>1) DisPatchKD[j,i] else DisPatchKD[j,i]=0
    if (DisPatchKD[j,i]>1) SpecRich[i]=SpecRich[i]+1 else SpecRich[i]
    PrevD[j,i]=((r[j]-(delta[j]*((r[j]+d[j]+alphaD[j])/(BetaD[j]+delta[j]))))/BetaD[j])/
      ((r[j]+d[j]+alphaD[j])/(BetaD[j]+delta[j]))
    IAbunD[j,i]=PrevD[j,i]*DisPatchKD[j,i]
  }
}  

##   Boring plot of the number of hosts in each patch

par(mfrow=c(1,1))  
plot(Area, PatchK[1,], ylab='Hosts', xlab='Patch Area(km^2)',
     col ="red", log = "x", #ylim = c(0.1, 300), 
     main = 'Number of individuals in patches of different size',
     pch=20)
plot(Area, DisPatchKD[1,]/Area, col='red')
points(Area, PatchK[2,], col = "blue",pch=20)
points(Area, PatchK[5,], col = "green", pch=20) 
points(Area, PatchK[6,], col= "black", pch=20)
points(Area, PatchK[9,], col= "orange", pch=20) 
leg.text<-c("5gram", "100gram ", "500gram ", "1kg ", "5kg")
legend("topleft",leg.text,lty=rep(1,4),col=c( "red", "blue", "green", "black", "orange"),bty="n")


##  Do a quick plot of the  number of infected hosts        

plot(Area, IAbunD[1,], ylab='Infected Hosts', xlab='Patch Area(km^2)',
     col ="red", log = "x",# ylim = c(0.1, 100), 
     main = 'Number in patches of different size',
     pch=20)
points(Area, IAbunD[2,], col = "blue", pch=20)
points(Area, IAbunD[5,], col = "green", pch=20) 
points(Area, IAbunD[6,], col= "black", pch=20)
points(Area, IAbunD[9,], col= "orange", pch=20) 
legend("topleft",leg.text,lty=rep(1,4),
       col=c("red", "blue", "green", "black", "orange"),bty="n")


## Boring plot of the prevalence data
plot(Area, PrevD[1,], ylab='Prevalence', xlab='Patch Area(km^2)',
     col ="red", log = "x", ylim = c(0,0.1), 
     main = 'Prevalence in patches of different size',
     pch=20)
points(Area, PrevD[2,], col = "blue", pch=20)
points(Area, PrevD[5,], col = "green", pch=20) 
points(Area, PrevD[6,], col= "black", pch=20)
points(Area, PrevD[9,], col= "orange", pch=20) 
leg.text<-c("5gram", "100gram ", "500gram ", "1kg ", "5kg")
legend("topleft",leg.text,lty=rep(1,4),
       col=c("red", "blue", "green", "black", "orange"),bty="n")



###########################################################################
##
##      Species Area curve

par(mfrow=c(1,1))
plot(Area, SpecRich, ylab="Species", xlab="Area of Patch", log="x",
      main= "Species Area Curve", pch=20)
      
## Now add the number of infected together to get at net risk 

SumDen=numeric(Amax)
ISumDen=numeric(Amax)
AvgPrev=numeric(Amax)

for(i in 1:Amax)   {
     SumDen[i]=0.0
     ISumDen[i]=0.0
     for(j in 1:Hosts) {
     SumDen[i]=PatchK[j,i]+SumDen[i]
#
#    Pathogen prevalence
#     
     ISumDen[i]=IAbunF[j,i]+ISumDen
     }
   AvgPrev[i]=ISumDen[i]/SumDen[i]
   }  

##  Total number of animals and number infected
# 
par(mfrow=c(1,1))  
plot(Area, SumDen, ylab='Numbers', xlab=' Patch Area(km^2)',
  col ="blue", log = "xy", ylim = c(1, 100000), main = 'Numbers in patches of different size' )
points(Area, ISumDen, col = "red")

leg.text<-c("Total Hosts", "Number Infectious")
legend("topleft",leg.text,lty=rep(1,4),col=c("blue", "red"),bty="n")

##  Proportion of all animals that are infected
#
par(mfrow=c(1,1)) 
plot(Area, AvgPrev, ylab='Prevalence', xlab=' Patch Area(km^2)',
  col ="red", log = "x", ylim = c(0.0005, 0.05), main = 'Overall Prevalence ' ) 

###
##   Add in 3d surface if number hosts and infected
##

require(grDevices) # for trans3d
##                   -----------

op <- par(bg = "white")
persp(Mass, Area, Prev, theta = 30, phi = 30, expand = 0.5, col = "lightblue")
persp(Mass, Area, Prev, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
       ltheta = 120, shade = 0.75, ticktype = "detailed",
       xlab = "Mass of Host", ylab = "Area of Patch", zlab = "Infecteds"
) -> res
round(res, 3)

####

###########################################################
out <- qplot(Area,PatchK[1,],color=qsec)
 + qplot(Area, PatchK[2,], color =qsec)
 
 out
 

## Superimpose plots on top of each other

qplot(Area,Prev[1,],color="red")
qplot(Area,Prev[2,],color="blue")

# p.tmp
   

######################################################################
##
##   Explore alternative formulations
##
########################################################################



