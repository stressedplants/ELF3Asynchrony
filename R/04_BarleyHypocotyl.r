#########################################
#
#. Step 4: ELF3 equivalent in barley

library(fda)

cols=c(rgb(1, 0.5,0.5,1),rgb(1, 0.7,0.2,1), rgb(0.8, 0.6,1,1))
barleyGrowth=read.csv("data/barleyGrowthCurves.csv", header=T)
times=c(1,3,5,8,10,12,15,17)
maxMass=max(barleyGrowth[,paste("W", times, sep="")])


sapply(unique(barleyGrowth$Genotype), function(g){
  ids=which(barleyGrowth$Genotype==g)
  
  plot(c(), xlim=c(0, 20), ylim=c(0, maxMass), main=g)
  
  apply(barleyGrowth[ids,], 1, function(i){
    lines(times[1:8],i[paste("W", times[1:8], sep="")])
  })
  
  
})


plot(c(), xlim=c(0, 17), ylim=c(0, maxMass), xlab="Day", ylab="Weight (g)")

colsForGenos=cols
names(colsForGenos)=c("Bowman", "B289", "B290")

sapply(names(colsForGenos), function(g){
  ids=which(barleyGrowth$Genotype==g)
  
  apply(barleyGrowth[ids,], 1, function(i){
    lines(times[1:8],i[paste("W", times[1:8], sep="")], col=colsForGenos[g])
  })
  
  
})
legend(0, 2.5, names(colsForGenos), col=colsForGenos, lwd=2, bty='n')

#tidy for ggplots
dfBarleyGrowth=do.call(rbind, apply(barleyGrowth, 1, function(i){
  data.frame("Genotype"=rep(i["Genotype"],length(times)), 
             "Day"=times,
             "Weight"=i[paste("W", times, sep="")])
}))
dfBarleyGrowth$Genotype=as.factor(dfBarleyGrowth$Genotype)
dfBarleyGrowth$Day=as.factor(dfBarleyGrowth$Day)
dfBarleyGrowth$Weight=as.numeric(dfBarleyGrowth$Weight)


#TODO: fix x-axis
ggplot(dfBarleyGrowth[which(dfBarleyGrowth$Genotype %in% names(colsForGenos)),], aes(x=Day, y=Weight, fill=Genotype))+geom_boxplot() + theme_minimal() 

renames=c("Bowman", "eam10.m (Hvlux)", "eam5.x (HvPHYC)", "eam8.k (Hvelf3)", "eam8.w (Hvelf3)", "Antonella")
names(renames)=c("Bowman", "B284", "B285", "B289", "B290", "Antonella")
dfBarleyGrowth$Genotype=factor(renames[dfBarleyGrowth$Genotype], levels=renames)
ggplot(dfBarleyGrowth, aes(x=Day, y=Weight, fill=Genotype))+geom_boxplot() + theme_minimal() + ylab("Weight (g)")
ggsave(file="plots/04_01_genotypeGrowthCurveAllBarley.pdf", width=7.5, height=5)

########################################
#do parallel analysis as with Arabidopsis


yt=lapply(unique(barleyGrowth[,"Genotype"]), function(g){
  barleyGrowth[which(barleyGrowth[,"Genotype"]==g),paste("W", times[1:8], sep="")]
})
#get range of means
range=sapply(yt, function(i){
  temp=colMeans(i)
  plot(temp)
  c(min(temp), max(temp))
})

#sig diff stdevs?
centredFinal=lapply(yt, function(i){
  i[,"W17"]-mean(i[,"W17"])
})

#I only care about ELF3
unadjusted=c(var.test(centredFinal[[4]], centredFinal[[6]], alternative="greater")$p.value,
             var.test(centredFinal[[5]], centredFinal[[6]], alternative="greater")$p.value)

adjusted=p.adjust(unadjusted)

#variance over time
sdOverTime=sapply(unique(barleyGrowth$Genotype), function(g){
  apply(barleyGrowth[which(barleyGrowth$Genotype==g),paste("W", times[1:8], sep="")], 2, function(i){sd(i)})
})

meanOverTime=sapply(unique(barleyGrowth$Genotype), function(g){
  apply(barleyGrowth[which(barleyGrowth$Genotype==g),paste("W", times[1:8], sep="")], 2, function(i){mean(i)})
})

plot(meanOverTime["W5",], sdOverTime["W5",])

plot(c(), xlim=c(1, 17), ylim=c(0,0.3))

apply(sdOverTime, 2, function(i){
  lines(times, i)
})

nbasis <- 8
Lfdobj <- 2
lambda <- c(10^(-4))
rng <- c(min(range), max(range))
hgtbasis <- create.bspline.basis(rng, nbasis)
#growfdPar <- fdPar(hgtbasis, Lfdobj, lambda)

##ii) What is g^(-1)(y(t))?
#install.packages("fdapace")
library(fdapace)
invGtofYt=lapply(1:length(yt), function(id){
  i=yt[[id]]
  print(id)
  
  #temp=colMeans(i)
  
  # Ly=lapply(c(1:(dim(i)[1])), function(j){
  #   as.numeric(i[j,])
  # })
  # Lt=lapply(c(1:(dim(i)[1])), function(j){
  #   times
  # })
  temp=colMeans(i) #GetMeanCurve(Ly, Lt)$mu
  #print(temp)
  
  rng <- c(min(i), max(i))
  hgtbasis <- create.bspline.basis(rng, nbasis)
  growfdPar <- fdPar(hgtbasis, Lfdobj, lambda[1])
  
  temp2=smooth.basis(temp, times, growfdPar)$fd
  plot(temp2, main=id)
  points(temp, c(1:(dim(yt[[1]])[2])))
  apply(i, 1, function(j){
    predict.fd(temp2, j)
  })
  
})

##iii) plot g^(-1)(y(t))
#Recall: yt=list(Ws2HeightLD, Ws2HeightSD, elf3HeightsLD, elf3HeightsSD)

sapply(invGtofYt, function(i){
  plot(c(), xlim=c(0, 17), ylim=c(min(i), max(i)))
  apply(i, 2, function(j){
    lines(times, j)
  })
})


####include
names(invGtofYt)=unique(barleyGrowth$Genotype)
plot(c(), xlim=c(0, 17), ylim=c(min(unlist(invGtofYt)), max(unlist(invGtofYt))), xlab="Chronological time (day)", ylab="Biomass scaled by mean growth curve")
sapply(c("Bowman", "B289", "B290"), function(id){
  i=invGtofYt[[id]]
  apply(i, 2, function(j){
    lines(times, j, col=colsForGenos[id])
  })
})


#align? 0=mx+b, -b/m
xints=lapply(invGtofYt, function(i){
  plot(c(), xlim=c(0, 20), ylim=c(min(i), max(i)))
  apply(i, 2, function(j){
    cfs=lm(j~times)$coef
    xint=-cfs[1]/cfs[2]
    lines(times-xint, j)
    xint
  })
})
#supp
names(xints)=names(invGtofYt)
names(xints)=renames[names(xints)]
pdf(file="plots/04_02_xInterceptsBarley.pdf", width=7, height=5)
par(mar = c(10, 5, 0.5, 0.5))
boxplot(xints, ylab="X-intercept (predicted germination day)", xlab="", las=3)
abline(h=c(-1, 1), lty=2)
dev.off()
###Can we use the initial fresh weight as an indicator of starting point variation?
pdf(file="plots/04_03_xInterceptsVsIFW.pdf", width=7, height=5)
plot(barleyGrowth[,"W1"], unlist(xints), xlab="Initial Fresh Weight (g)", ylab="X-intercept (day)") ####Definitely, YES!!!
dev.off()
###What do the growth curves look like once aligned?

png(file="plots/04_04_barleyGrowthCurvesBioAge.png", width=7, height=7, res=500, unit='in')
par(mfrow = c(2, 2))
#not aligned
plot(c(), xlim=c(0, 20), ylim=c(0, maxMass), xlab="Day", ylab="Weight (g)")

colsForGenos=cols
names(colsForGenos)=c("Bowman", "B289", "B290")
names(xints)=unique(barleyGrowth[,"Genotype"])
sapply(names(colsForGenos), function(g){
  print(g)
  ids=which(barleyGrowth$Genotype==g)
  
  sapply(c(1:length(ids)), function(is){
    i=barleyGrowth[ids[is],]
    lines(times[1:8],i[paste("W", times[1:8], sep="")], col=colsForGenos[g])
    
    # lines(times[1:8]-xints[[g]][is],i[paste("W", times[1:8], sep="")], col=colsForGenos[g])
  })
  
  
})
legend(0, 2.5, c("Bowman", "eam8.k (HvELF3)", "eam8.w (HvELF3)"), col=colsForGenos, lwd=2, bty='n')

#aligned
plot(c(), xlim=c(0, 20), ylim=c(0, maxMass), xlab="Day", ylab="Weight (g)")

colsForGenos=cols
names(colsForGenos)=c("Bowman", "B289", "B290")
names(xints)=unique(barleyGrowth[,"Genotype"])
sapply(names(colsForGenos), function(g){
  print(g)
  ids=which(barleyGrowth$Genotype==g)
  
  sapply(c(1:length(ids)), function(is){
    i=barleyGrowth[ids[is],]
    #lines(times[1:8],i[paste("W", times[1:8], sep="")], col=colsForGenos[g])
    
    lines(times[1:8]-xints[[g]][is],i[paste("W", times[1:8], sep="")], col=colsForGenos[g])
  })
  
  
})
###Unaligned biological age
names(invGtofYt)=unique(barleyGrowth$Genotype)
plot(c(), xlim=c(0, 17), ylim=c(min(unlist(invGtofYt)), max(unlist(invGtofYt))), xlab="Chronological time (day)", ylab="Biomass scaled by mean growth curve")
sapply(c("Bowman", "B289", "B290"), function(id){
  i=invGtofYt[[id]]
  apply(i, 2, function(j){
    lines(times, j, col=colsForGenos[id])
  })
})
#aligned biological age
plot(c(), xlim=c(0, 17), ylim=c(min(unlist(invGtofYt)), max(unlist(invGtofYt))), xlab="Chronological time (day)", ylab="Biomass scaled by mean growth curve")
sapply(c("Bowman", "B289", "B290"), function(id){
  i=invGtofYt[[id]]
  #find smooth function of yt
  sapply(1:length(xints[[id]]), function(jd){
    j=i[,jd]
    lines(times-xints[[id]][jd], j, col=colsForGenos[id])
  })
})

dev.off()


pdf(file="plots/04_05_barleyGrowthCurvesBioAge.pdf", width=3.8, height=7)
par(mfrow = c(3, 1))
par(mar = c(5, 5, 1.3, 1))
####Draw histograms at day 10, aligned and not aligned
#not aligned:
names(yt)=names(invGtofYt)
temp=lapply(c("Bowman", "B289", "B290"), function(id){
  i=yt[[id]][,"W10"]
})
names(temp)=c("Bowman", "eam8.k (Hvelf3)", "eam8.w (Hvelf3)")
boxplot(temp, col=colsForGenos, xlab="Genotypes", ylab="Weight (g)", main="Day 10 (unaligned)")

dfTemp=data.frame(size=c(temp[[1]], 
                         temp[[2]],
                         temp[[3]]),
                  genotype=c(rep("Bowman", length(temp[[1]])),
                             rep("eam8.k", length(temp[[2]])),
                             rep("eam8.w", length(temp[[3]]))))
TukeyHSD(aov(size~genotype,data=dfTemp))

pvals=c(var.test(temp[[1]], temp[[2]])$p.value,
        var.test(temp[[1]], temp[[3]])$p.value)
p.adjust(pvals)



#aligned: 
names(xints)=names(yt)
tempAligned=lapply(c("Bowman", "B289", "B290"), function(id){
  i=yt[[id]]
  rng <- c(-10, 30)
  hgtbasis <- create.bspline.basis(rng, nbasis)
  growfdPar <- fdPar(hgtbasis, Lfdobj, lambda[1])
  idRev=id
  print(idRev)
  sapply(1:length(xints[[idRev]]), function(jd){
    j=as.numeric(i[jd,])
    print(length(j))
    temp2=smooth.basis(times-xints[[idRev]][jd],j, growfdPar )$fd
    predict.fd(temp2, c(10))
  })
  
})
names(tempAligned)=c("Bowman", "eam8.k (Hvelf3)", "eam8.w (Hvelf3)")
boxplot(tempAligned, col=colsForGenos, xlab="Genotypes", ylab="Weight (g)", main="Day 10 (aligned)")

dfTemp=data.frame(size=c(tempAligned[[1]], 
                         tempAligned[[2]],
                         tempAligned[[3]]),
                  genotype=c(rep("Bowman", length(tempAligned[[1]])),
                             rep("eam8.k", length(tempAligned[[2]])),
                             rep("eam8.w", length(tempAligned[[3]]))))
TukeyHSD(aov(size~genotype,data=dfTemp))

pvals=c(var.test(tempAligned[[1]], tempAligned[[2]])$p.value,
        var.test(tempAligned[[1]], tempAligned[[3]])$p.value)
p.adjust(pvals)


xintsSub=list("Bowman"=xints[["Bowman"]],"eam8.k (Hvelf3)"=xints[["B289"]], "eam8.w (Hvelf3)"=xints[["B290"]])
boxplot(xintsSub, col=colsForGenos, xlab="Genotypes", ylab="Germination day")
abline(h=c(-1, 1), lty=2)


dev.off()


###What do biological times look like for plants that are within one day of each other in terms of germination time?
names(xints)=names(yt)
#aligned biological age
plot(c(), xlim=c(0, 17), ylim=c(min(unlist(invGtofYt)), max(unlist(invGtofYt))), xlab="Chronological time (day)", ylab="Biomass scaled by mean growth curve")
sapply(c("Bowman", "B289", "B290"), function(id){
  i=invGtofYt[[id]][,which(xints[[id]]>(-0.5) & xints[[id]]<(0.5))]
  #find smooth function of yt
  sapply(1:length(which(xints[[id]]>(-0.5) & xints[[id]]<(0.5))), function(jd){
    j=i[,jd]
    lines(times-xints[[id]][jd], j, col=colsForGenos[id])
  })
})


###Are daily growth rates following a gamma distribution?
dailyGrowths=sapply(invGtofYt, function(i){
  apply(i, 2, function(j){
    (j[2:8]-j[1:7])/(times[2:8]-times[1:7])
  })
})

slopesPos=lapply(dailyGrowths, function(i){
  temp=as.numeric(i)
  temp=temp[which(temp>0)]
})





#get value of "maturation rate" for each genotype
#Let's get the value of the maturation 



pdf(file=paste("plots/04_06_biologicalAgeBarley.pdf", sep=""), height=9.2, width=4.5)
par(mfcol = c(2,1), oma = c(0.1, 0.1, 0.1, 0.1))

cols_by_cond=c("lightgreen", "forestgreen", "lightblue", "cornflowerblue")
plot(c(), xlim=c(74,(74+143)), ylim=c(74, 300), xlab="Chronological time (h)", ylab="Biomass scaled by mean growth curve")
sapply(c(1:4), function(id){
  i=invGtofYt[[id]]
  apply(i, 2, function(j){
    lines(c(74:(74+143))[seq(1, 143, 24)]+12, j[seq(1, 143, 24)]+74, col=cols_by_cond[id])
  })
})

legend(80, 300, c("Ws-2 SD", "Ws-2 LD", "elf3-4 SD", "elf3-4 LD"), lwd=1, col=c("forestgreen", "lightgreen", "cornflowerblue", "lightblue"), bty='n')



###colour-code and 



dev.off()


