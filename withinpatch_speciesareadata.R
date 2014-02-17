rm(list = ls())
setwd("/Users/cfaust/Dropbox/NCEAS/Within Patch Host Dynamics")
library(Hmisc)
library(plyr)
library(scales)
#import data
bodysize <- read.csv("bodysizepatch_17Feb2014.csv", header=T)

bodysize$bodyclasses<- as.numeric(cut2(bodysize$body_size_grams, 
                                       cuts=c(0, 10, 10^2, 10^3, 10^4, 10^5, 10^6)))
plot(log(bodysize$patch_area+1), log(bodysize$density),
     col=c("black","darkred","purple","darkgoldenrod","navy", "lightseagreen")[bodysize$type],
     ylab="density, individuals/hectare", xlab="log of patch area")
legend("topleft", col=c("black","darkred","purple","darkgoldenrod","navy", "lightseagreen")[bodysize$type],bodysize$type)

plot(log(bodysize$patch_area+1), log(bodysize$body_size_grams),
     col=c("black","darkred","purple","darkgoldenrod","navy", "lightseagreen")[bodysize$type],
     ylab="dbody size log", xlab="log of patch area")

plot(log(bodysize$patch_area), bodysize$bodyclasses,
     col=c("black","darkred","purple","darkgoldenrod","navy", "lightseagreen")[bodysize$type],
     ylab="dbody size log", xlab="log of patch area")
symbols(log(bodysize$patch_area), bodysize$density, 
        circles= bodysize$bodyclasses,
        pch=19,
        ylab="density", xlab="log of patch area")
colors=c("firebrick3", "cadetblue","steelblue4", 
                       "coral", "darkgoldenrod1", "darkolivegreen4")
symbols(x=log(bodysize$patch_area), y=log(bodysize$density+0.1), 
        circles=bodysize$bodyclasses, inches=1/5,  
        bg= alpha(colors,0.5)[bodysize$type], 
        pch=16,
        bty="n", las=1, cex.lab=1)
symbols(x=log(bodysize$density+0.1), y=bodysize$bodyclass, 
        circles=bodysize$bodyclasses, inches=1/5,  
        bg= alpha(colors,0.5)[bodysize$type], 
        pch=16,
        bty="n", las=1, cex.lab=1)
legend("topright", bty="n",names(table(bodysize$type)), 
       pch=20, col=colors)
symbols(x=log10(bodysize$patch_area), y=(bodysize$density+1), 
        circles=bodysize$bodyclasses, inches=1/5,  
        bg= alpha(colors,0.5)[bodysize$type], 
        xlim= c(2, 5),pch=16,
        bty="n", las=1, cex.lab=1)
legend("topleft", bty="n",names(table(bodysize$type)), 
       pch=20, col=colors)

