##############################
#
# Step 2: Analyse hypocotyl data in Ws-2 and elf3-4

library(fitdistrplus)

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

#########
# Draw the same figures, but draw CV and stdev over time for each condition

heterogeneity_over_time=lapply(c("Ws2", "elf3"), function(geno){
  ids=which(genotypes==geno & condition %in% c("SD", "LD"))
  condition_sub=condition[c(1:3, ids)]
  print(table(condition_sub))
  data_WS2=dataNum[,c(1:3, ids)]
  
  #growth rate
  # plot(data_WS2[,2], data_WS2[,4], type="l")
  
  #height
  height=sapply(c(2:dim(data_WS2)[1]), function(i){
    colSums(data_WS2[1:i, ])
  })
  
  #sd
  h_sd=height[which(condition_sub=="SD"),]
  sd_sd=apply(h_sd, 2, function(val){
    sd(val)
  })
  cv_sd=apply(h_sd, 2, function(val){
    sd(val)/mean(val)
  })
  
  #ld
  h_ld=height[which(condition_sub=="LD"),]
  sd_ld=apply(h_ld, 2, function(val){
    sd(val)
  })
  cv_ld=apply(h_ld, 2, function(val){
    sd(val)/mean(val)
  })
  
  data.frame("St. Dev, SD"=sd_sd,
             "St. Dev, LD"=sd_ld,
             "CV, SD"=cv_sd,
             "CV, LD"=cv_ld)
  })

#st dev figure:
png(filename = "plots/02_03_stDev_cv_overTime.png", height=800, pointsize=20)
par(mfrow = c(2, 1))
font=1
plot(c(), lwd=2, cex.lab=font, cex.axis=font, cex.main=font, cex.sub=font,
       xlim=c(min(data_WS2[2:dim(data_WS2)[1],2]), 
              max(data_WS2[2:dim(data_WS2)[1],2])), 
       ylim=c(0, 1.5), xlab="Time (h)", ylab="St. Dev.")

lines(data_WS2[2:dim(data_WS2)[1],2], 
      heterogeneity_over_time[[1]][,"St..Dev..SD"], col="forestgreen", lwd=2)

lines(data_WS2[2:dim(data_WS2)[1],2], 
      heterogeneity_over_time[[1]][,"St..Dev..LD"], col="lightgreen", lwd=2)

lines(data_WS2[2:dim(data_WS2)[1],2], 
      heterogeneity_over_time[[2]][,"St..Dev..SD"], col="cornflowerblue", lwd=2)

lines(data_WS2[2:dim(data_WS2)[1],2], 
      heterogeneity_over_time[[2]][,"St..Dev..LD"], col="lightblue", lwd=2)


font=1
plot(c(), lwd=2, cex.lab=font, cex.axis=font, cex.main=font, cex.sub=font,
     xlim=c(min(data_WS2[2:dim(data_WS2)[1],2]), 
            max(data_WS2[2:dim(data_WS2)[1],2])), 
     ylim=c(0, 2), xlab="Time (h)", ylab="Coeff. of Var.")

lines(data_WS2[2:dim(data_WS2)[1],2], 
      heterogeneity_over_time[[1]][,"CV..SD"], col="forestgreen", lwd=2)

lines(data_WS2[2:dim(data_WS2)[1],2], 
      heterogeneity_over_time[[1]][,"CV..LD"], col="lightgreen", lwd=2)

lines(data_WS2[2:dim(data_WS2)[1],2], 
      heterogeneity_over_time[[2]][,"CV..SD"], col="cornflowerblue", lwd=2)

lines(data_WS2[2:dim(data_WS2)[1],2], 
      heterogeneity_over_time[[2]][,"CV..LD"], col="lightblue", lwd=2)

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

#in response to reviewer: show that the overlap is still convincing if each light condition is drawn separately
cols_by_cond=c("lightgreen", "forestgreen", "lightblue", "cornflowerblue")
lyt_by_cond=c(1, 2, 1, 2)
plot(c(), xlim=c(74,(74+143)), ylim=c(74, 300), xlab="Chronological time (h)", ylab="Inferred biological age (h)")
sapply(c(1,3), function(id){
  i=invGtofYt[[id]]
  apply(i, 2, function(j){
    lines(c(74:(74+143))[seq(1, 143, 24)]+12, j[seq(1, 143, 24)]+74, col=cols_by_cond[id], lwd=1.5, lty=lyt_by_cond[id])
  })
})

legend(80, 300, c("Ws-2 SD", "Ws-2 LD", "elf3-4 SD", "elf3-4 LD"), lwd=1.5, col=c("forestgreen", "lightgreen", "cornflowerblue", "lightblue"),lty=lyt_by_cond, bty='n')
abline(c(0,1), lty=2)

cols_by_cond=c("lightgreen", "forestgreen", "lightblue", "cornflowerblue")
lyt_by_cond=c(1, 2, 1, 2)
plot(c(), xlim=c(74,(74+143)), ylim=c(74, 300), xlab="Chronological time (h)", ylab="Inferred biological age (h)")
sapply(c(2,4), function(id){
  i=invGtofYt[[id]]
  apply(i, 2, function(j){
    lines(c(74:(74+143))[seq(1, 143, 24)]+12, j[seq(1, 143, 24)]+74, col=cols_by_cond[id], lwd=1.5, lty=lyt_by_cond[id])
  })
})

legend(80, 300, c("Ws-2 SD", "Ws-2 LD", "elf3-4 SD", "elf3-4 LD"), lwd=1.5, col=c("forestgreen", "lightgreen", "cornflowerblue", "lightblue"),lty=lyt_by_cond, bty='n')
abline(c(0,1), lty=2)



#look at how biological age changes for each unit of chronological age
png(file=paste("plots/02_02_biologicalAgeElf3Arabidopsis.png", sep=""), height=9.2, width=4.5, units="in", res=500)
par(mfcol = c(2,1), oma = c(0.1, 0.1, 0.1, 0.1))

cols_by_cond=c("lightgreen", "forestgreen", "lightblue", "cornflowerblue")
lyt_by_cond=c(1, 2, 1, 2)
plot(c(), xlim=c(74,(74+143)), ylim=c(74, 300), xlab="Chronological time (h)", ylab="Inferred biological age (h)")
sapply(c(1:4), function(id){
  i=invGtofYt[[id]]
  apply(i, 2, function(j){
    lines(c(74:(74+143))[seq(1, 143, 24)]+12, j[seq(1, 143, 24)]+74, col=cols_by_cond[id], lwd=1.5, lty=lyt_by_cond[id])
  })
})

legend(80, 300, c("Ws-2 SD", "Ws-2 LD", "elf3-4 SD", "elf3-4 LD"), lwd=1.5, col=c("forestgreen", "lightgreen", "cornflowerblue", "lightblue"),lty=lyt_by_cond, bty='n')
abline(c(0,1), lty=2)


###colour-code 

#get slopes:
freq=24
slopes=sapply(invGtofYt, function(i){
  apply(i, 2, function(j){
    j[seq(freq+1, 142, freq)]-j[seq(1, 142-freq, freq)]
  })
})

mergedSlopes=as.numeric(unlist(slopes))
#hist(mergedSlopes/24, xlab="Rate of maturation (per day)", col='orange', main="")

#draw the curve of the closest fitting gamma distribution
fit_together <- fitdist(mergedSlopes/24, distr = "gamma", method = "mle")
denscomp(fit_together, addlegend=FALSE, xlab="Inferred biological/chronological age (per day)", main="")
legend(1.5, 1.2, "gamma", col="red", lty=1, bty="n")
dev.off()

#as per revision:  have separate gamma dist for each category
png(file=paste("plots/02_04_gamma_split.png", sep=""), height=9.2, width=9.2, units="in", res=500)
par(mfcol = c(2,2), oma = c(0.1, 0.1, 0.1, 0.1))

a=sapply(slopes, function(sl){
mergedSlopes=as.numeric(sl)
print(mergedSlopes)
fit <- fitdist(mergedSlopes/24, distr = "gamma", method = "mle")
denscomp(fit, addlegend=FALSE, xlab="Inferred biological/chronological age (per day)", main="")
fit$loglik
})

#legend(1.5, 1.2, "gamma", col="red", lty=1, bty="n")
dev.off()

#log-ratio test:
LR <- 2 * (sum(a) - fit_together$loglik)
df <- 6  # each gamma has 2 parameters

p_value <- pchisq(LR, df = df, lower.tail = FALSE)
p_value
#small sample sizes make this unrobust

set.seed(123)
B <- 10000     
LR_boot <- numeric(B)

shape0 <- fit_together$estimate["shape"]
rate0  <- fit_together$estimate["rate"]

for (b in seq_len(B)) {
  
  # Simulate under H0
  x_sim <- rgamma(length(mergedSlopes), shape = shape0, rate = rate0)
  
  # Split into groups
  x_sim_list <- list("a"=x_sim[1:35],
                     "b"=x_sim[36:(36+40-1)],
                     "c"=x_sim[(36+40):(36+40+30-1)],
                     "d"=x_sim[(36+40+30):length(mergedSlopes)])
                     #split(x_sim, c(1:4))
  
  # Refit null
  refit <- fitdistr(x_sim, "gamma")$loglik
  
  # Refit alternative
  refit_alt <- 0
  for (g in c(1:4)) {
    print(length(x_sim_list[[g]]))
    refit_alt <- refit_alt + fitdistr(x_sim_list[[g]], "gamma")$loglik
  }
  
  LR_boot[b] <- 2 * (refit_alt - refit)

}

#pvals
df <- 6  # each gamma has 2 parameters
pvals=sapply(LR_boot, function(i){
  pchisq(i, df = df, lower.tail = FALSE)
})

length(which(p_value>pvals))/10000


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


