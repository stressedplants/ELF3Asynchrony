##############################
#
# Step 2: Analyse hypocotyl data in Ws-2 and elf3-4

#0. Load data
data=read.csv("data/expt1_rawdata_phyBelf3growth.csv", stringsAsFactors = F, header=F)

###process data:
condition=paste(data[1,])
genotypes=paste(data[2,])
replicate=paste(data[3,])

colnames(data)=paste(condition, genotypes, replicate, sep="_")
colnames(data)[1]="ImageN"
colnames(data)[2]="Time(h)"
colnames(data)[3]="ZT"

dataNum=data[-c(1,2,3),]

dataNum=apply(dataNum, c(1,2), function(i){as.numeric(i)})

ids=which(genotypes=="Ws2")
data_WS2=dataNum[,c(1:3, ids)]



###### Plot the growth curves for each genotype


pdf(file=paste("plots/02_01_hypocotylHeights_Ws-2_elf3only.pdf", sep=""), height=9.2, width=4.5)
par(mfcol = c(2,1), oma = c(0.1, 0.1, 0.1, 0.1))


genos=c("Ws-2", "elf3-4")
names(genos)=c("Ws2", "elf3")

library(fda)
nbasis=200
basisobj <- create.bspline.basis(c(74, 74+143), nbasis)
lambda   <- 0.00001
fdParobj <- fdPar(basisobj, 2, lambda)

heights=lapply(c("Ws2", "elf3"), function(geno){
  ids=which(genotypes==geno & condition %in% c("SD", "LD"))
  condition_sub=condition[c(1:3, ids)]
  data_WS2=dataNum[,c(1:3, ids)]
  
  #growth rate
  # plot(data_WS2[,2], data_WS2[,4], type="l")
  
  #height
  height=sapply(c(2:dim(data_WS2)[1]), function(i){
    colSums(data_WS2[1:i, ])
  })
  
  # plot(data_WS2[2:dim(data_WS2)[1],2], height[4,], type="l")
  font=1
  if(geno!="Ws2"){
    plot(c(), lwd=2, cex.lab=font, cex.axis=font, cex.main=font, cex.sub=font,
         xlim=c(min(data_WS2[2:dim(data_WS2)[1],2]), 
                max(data_WS2[2:dim(data_WS2)[1],2])), 
         ylim=c(0, 12), main=substitute(italic(g), list(g = genos[geno])), xlab="Time (h)", ylab="Hypocotyl length (mm)")
  }else{
    plot(c(), lwd=2, cex.lab=font, cex.axis=font, cex.main=font, cex.sub=font,
         xlim=c(min(data_WS2[2:dim(data_WS2)[1],2]), 
                max(data_WS2[2:dim(data_WS2)[1],2])), 
         ylim=c(0, 12), main=substitute(g, list(g = genos[geno])), xlab="Time (h)", ylab="Hypocotyl length (mm)")
  }
  if(geno=="Ws2"){
    cols_by_cond=c("forestgreen", "lightgreen") #, "darkmagenta", "orange")
  }else{
    cols_by_cond=c("cornflowerblue", "lightblue")
  }
  names(cols_by_cond)=c("SD", "LD") #, "LDtoSD", "SDtoLD")
  sapply(4:length(height[,143]), function(i){
    lines(data_WS2[2:dim(data_WS2)[1],2], height[i,], col=cols_by_cond[condition_sub[i]])
  })
  if(geno=="Ws2"){
    legend(65, 17, legend= c("SD", "LD"), col = cols_by_cond[c("SD", "LD")], lty=1, bty='n')
  }
  
  list("data_WS2"=data_WS2, "height"=height, "condition_sub"=condition_sub)
})


dev.off()

######
#A: what is the "biological age curve"? g(f(t))=y(t), where g(t) is the average curve, f(t) is the biological age curve, and y(t) is the observed growth curve
# Therefore, to find f(t), I need g^(-1))y(t)

##i) What is y(t)?
Ws2HeightLD=heights[[1]][[2]][which(heights[[1]][[3]]=="LD"),]
Ws2HeightSD=heights[[1]][[2]][which(heights[[1]][[3]]=="SD"),]
elf3HeightsLD=heights[[2]][[2]][which(heights[[2]][[3]]=="LD"),]
elf3HeightsSD=heights[[2]][[2]][which(heights[[2]][[3]]=="SD"),]

yt=list(Ws2HeightLD, Ws2HeightSD, elf3HeightsLD, elf3HeightsSD)

#get range of means
range=sapply(yt, function(i){
  temp=colMeans(i)
  plot(temp)
  c(min(temp), max(temp))
})

#sig diff stdevs?
centredFinal=lapply(yt, function(i){
  i[,143]-mean(i[,143])
})
unadjusted=c(var.test(centredFinal[[1]], centredFinal[[2]], alternative="less")$p.value,
             var.test(centredFinal[[1]], centredFinal[[3]], alternative="less")$p.value,
             var.test(centredFinal[[1]], centredFinal[[4]], alternative="less")$p.value,
             var.test(centredFinal[[2]], centredFinal[[3]], alternative="greater")$p.value,
             var.test(centredFinal[[2]], centredFinal[[4]], alternative="less")$p.value,
             var.test(centredFinal[[3]], centredFinal[[4]], alternative="less")$p.value)

adjusted=p.adjust(unadjusted)



nbasis <- 140
Lfdobj <- 2
lambda <- c(10^(-6), 10^(-4),  10^(-4), 10^(-4))
rng <- c(min(range), max(range))
hgtbasis <- create.bspline.basis(rng, nbasis)
#growfdPar <- fdPar(hgtbasis, Lfdobj, lambda)

##ii) What is g^(-1)(y(t))?
#install.packages("fdapace")
library(fdapace)
invGtofYt=lapply(1:length(yt), function(id){
  i=yt[[id]]
  #temp=colMeans(i)
  
  Ly=lapply(c(1:(dim(i)[1])), function(j){
    i[j,]
  })
  Lt=lapply(c(1:(dim(i)[1])), function(j){
    1:143
  })
  temp=GetMeanCurve(Ly, Lt)$mu
  #print(temp)
  
  rng <- c(min(i), max(i))
  hgtbasis <- create.bspline.basis(rng, nbasis)
  growfdPar <- fdPar(hgtbasis, Lfdobj, lambda[id])
  
  temp2=smooth.basis(temp, c(1:(dim(yt[[1]])[2])), growfdPar)$fd
  plot(temp2)
  points(temp, c(1:(dim(yt[[1]])[2])))
  apply(i, 1, function(j){
    predict.fd(temp2, j)
  })
  
})

##iii) plot g^(-1)(y(t))
#Recall: yt=list(Ws2HeightLD, Ws2HeightSD, elf3HeightsLD, elf3HeightsSD)

sapply(invGtofYt, function(i){
  plot(c(), xlim=c(0, 143), ylim=c(min(i), max(i)))
  apply(i, 2, function(j){
    lines(c(1:143), j)
  })
})


#look at how biological age changes for each unit of chronological age
png(file=paste("plots/02_02_biologicalAgeElf3Arabidopsis.png", sep=""), height=9.2, width=4.5, units="in", res=500)
par(mfcol = c(2,1), oma = c(0.1, 0.1, 0.1, 0.1))

cols_by_cond=c("lightgreen", "forestgreen", "lightblue", "cornflowerblue")
plot(c(), xlim=c(74,(74+143)), ylim=c(74, 300), xlab="Chronological time (h)", ylab="Hypocotyl length scaled by mean growth curve")
sapply(c(1:4), function(id){
  i=invGtofYt[[id]]
  apply(i, 2, function(j){
    lines(c(74:(74+143))[seq(1, 143, 24)]+12, j[seq(1, 143, 24)]+74, col=cols_by_cond[id])
  })
})

legend(80, 300, c("Ws-2 SD", "Ws-2 LD", "elf3-4 SD", "elf3-4 LD"), lwd=1, col=c("forestgreen", "lightgreen", "cornflowerblue", "lightblue"), bty='n')



###colour-code 

#get slopes:
freq=24
slopes=sapply(invGtofYt, function(i){
  apply(i, 2, function(j){
    j[seq(freq+1, 142, freq)]-j[seq(1, 142-freq, freq)]
  })
})

mergedSlopes=as.numeric(unlist(slopes))
hist(mergedSlopes/24, xlab="Rate of maturation (per day)", col='orange', main="")

dev.off()


#################
# Statistical tests to see if the histogram above follows a gamma distribution
#################


sapply(slopes, function(i){
  plot(c(), xlim=c(1, 5), ylim=c(min(i), max(i)))
  apply(i, 2, function(j){
    lines(c(1:5), j)
  })
})

#goodness of fit gamma
library(goft)
gamma_test(mergedSlopes)
gamma_test(as.numeric(slopes[[1]]))
gamma_test(as.numeric(slopes[[2]]))
gamma_test(as.numeric(slopes[[3]]))
gamma_test(as.numeric(slopes[[4]]))

ks.test(as.numeric(slopes[[1]]), as.numeric(slopes[[2]]), "pgamma")

ks.test(as.numeric(slopes[[1]]), as.numeric(slopes[[3]]), "pgamma")

ks.test(as.numeric(slopes[[1]]), as.numeric(slopes[[4]]) , "pgamma")

ks.test(as.numeric(slopes[[2]]), as.numeric(slopes[[3]]) , "pgamma")

ks.test(as.numeric(slopes[[2]]), as.numeric(slopes[[4]]) , "pgamma")

ks.test(as.numeric(slopes[[3]]), as.numeric(slopes[[4]]) , "pgamma")

library(fitdistrplus)
sapply(slopes, function(i){
  gamma_fit_try <- gamma_fit(as.numeric(i))
  print(paste(gamma_fit_try[1,], 1/gamma_fit_try[2,]))
  fit.gamma <- fitdist(as.numeric(i), distr = "gamma", method = "mle")
  print(summary(fit.gamma)$estimate)
})


