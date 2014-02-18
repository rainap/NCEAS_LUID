##        Andy Dobson - EEB, Princeton University
##
#############################################################################
##                                                                          #
##   Code for model to examine shrinking habitat and pathogen burdens       #
##                                                                          #
#############################################################################

#  install.package("ggplot2")
library("ggplot2")

##  Set range of Areas
Amax = 100
SpecArea=numeric(Amax)

# Set a range of Area values
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

#Set up 12 Hosts
Hosts=12
# Set a range of body masses  in Kg
Mass=numeric(Hosts)
Mass[1]=0.05
Mass[2]=0.10
Mass[3]=0.20
Mass[4]=0.333
Mass[5]=0.5
Mass[6]=1.0
Mass[7]=2.
Mass[8]=3.333
Mass[9]=5.0
Mass[10]=10.0
Mass[11]=20.0
Mass[12]=50.0

#Parameters for each host species
Dens = numeric(Hosts) #assume each host has a unique density that is fixed
Beta = numeric(Hosts) #each host has a unique transmission rate
r = numeric(Hosts) #reproduction rates for each host
d = numeric(Hosts) #death rates

Prev=matrix(NA, ncol = Amax, nrow = Hosts)
HDen=matrix(NA, ncol = Amax, nrow = Hosts)
IDen=matrix(NA, ncol = Amax, nrow = Hosts)


# Set Carrying capacity for each host species (#/km^2)
#  and also r and d for host transmission rate for pathogen
#  m is increase in host mortality rate

#  Density-dependent transmission
m = 26
for (i in 1:Hosts) {
    Dens[i]=16.2*Mass[i]^-0.70
    Beta[i]=0.0247*m*Mass[i]^0.44
     r[i]=0.6*Mass[i]^-0.26
     d[i]=0.4*Mass[i]^-0.27
    }

#  Frequency-dependent transmission
m = 1.12
for (i in 1:Hosts) {
    Dens[i]=16.2*Mass[i]^-0.70
    Beta[i]=0.4*m*Mass[i]^-0.26
     r[i]=0.6*Mass[i]^-0.26
     d[i]=0.4*Mass[i]^-0.27
    }


##########################################################
#
#       Set abundance in each patch
#            and then calculate number infected hosts
#

for(i in 1:Amax)   {
     SpecArea[i]=0
     for(j in 1:Hosts) {
     HDen[j,i]=Dens[j]*Area[i]/500
     if (HDen[j,i]>1) HDen[j,i] else HDen[j,i]=0
     if  (HDen[j,i]>1) SpecArea[i]=SpecArea[i]+1 else SpecArea[i]
#
#    Pathogen prevalence
#     
     Prev[j,i]=(r[j]/(m+d[j]))/(1+(r[j]/(m+d[j])))
     IDen[j,i]=Prev[j,i]*HDen[j,i]
     
     }
   }  

##   Boring plot of the number of hosts in each patch

par(mfrow=c(1,1))  
plot(Area, HDen[1,], ylab='Hosts', xlab='Patch Area(km^2)',
  col ="red", log = "x", ylim = c(0.1, 1000), main = 'Number in patches of different size' )
points(Area, HDen[2,], col = "blue")
points(Area, HDen[5,], col = "green") 
points(Area, HDen[6,], col= "black")
points(Area, HDen[9,], col= "orange") 
leg.text<-c("5gram", "100gram ", "500gram ", "1kg ", "5kg")
legend("topleft",leg.text,lty=rep(1,4),col=c( "red", "blue", "green", "black", "orange"),bty="n")


##  Do a quick plot of the  number of infected hosts        

par(mfrow=c(1,1))  
plot(Area, IDen[1,], ylab='Infected Hosts', xlab='Patch Area(km^2)',
  col ="red", log = "x", ylim = c(0.1, 100), main = 'Number in patches of different size' )
points(Area, IDen[2,], col = "blue")
points(Area, IDen[5,], col = "green") 
points(Area, IDen[6,], col= "black")
points(Area, IDen[9,], col= "orange") 
leg.text<-c("5gram", "100gram ", "500gram ", "1kg ", "5kg")
legend("topleft",leg.text,lty=rep(1,4),col=c("red", "blue", "green", "black", "orange"),bty="n")
 

##  ##   Boring plot of the prevalence data

par(mfrow=c(1,1))  
plot(Area, Prev[1,], ylab='Prevalence', xlab='Patch Area(km^2)',
  col ="red", log = "x", ylim = c(0.01, 0.5), main = 'Prevalence in patches of different size' )
points(Area, Prev[2,], col = "blue")
points(Area, Prev[5,], col = "green") 
points(Area, Prev[6,], col= "black")
points(Area, Prev[9,], col= "orange") 
leg.text<-c("5gram", "100gram ", "500gram ", "1kg ", "5kg")
legend("topleft",leg.text,lty=rep(1,4),col=c("red", "blue", "green", "black", "orange"),bty="n")


###########################################################################
##
##      Species Area curve

par(mfrow=c(1,1))
plot(Area, SpecArea, ylab="Species", xlab="Area of Patch", log="x",
      main= "Species Area Curve")
      
## Now add the number of infected together to get at net risk 

SumDen=numeric(Amax)
ISumDen=numeric(Amax)
AvgPrev=numeric(Amax)

for(i in 1:Amax)   {
     SumDen[i]=0.0
     ISumDen[i]=0.0
     for(j in 1:Hosts) {
     SumDen[i]=HDen[j,i]+SumDen[i]
#
#    Pathogen prevalence
#     
     ISumDen[i]=IDen[j,i]+ISumDen
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
out <- qplot(Area,HDen[1,],color=qsec)
 + qplot(Area, HDen[2,], color =qsec)
 
 out
 

## Superimpose plots on top of each other

qplot(Area,Prev[1,],color=red)
qplot(Area,Prev[2,],color=blue)

# p.tmp
   

######################################################################
##
##   Explore alternative formulations
##
########################################################################



