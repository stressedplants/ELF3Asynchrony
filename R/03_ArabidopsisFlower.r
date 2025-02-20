##################################################
#
# Step 3: establish trait heterogeneity at bolting in Arabidopsis in Ws-2 and Elf3
#
##################################################

#####################################
# establish the relationship between maturation rate and maturity point trait heterogeneity, using photoperiod shifts

#A) relationship between leaf count and bolting time
boltingData=read.csv("data/FloweringTime_Arabidopsis_photoperiodShifts.csv", header=T)
boltingData=boltingData[which(boltingData[,"genotype"]=="Ws-2"),]

rs=c(min(unique(boltingData[,"boltDay"])), max(unique(boltingData[,"boltDay"])))
cs=c(min(unique(boltingData[,"leafcount"])), max(unique(boltingData[,"leafcount"])))
img=matrix(0, nrow=rs[2]-rs[1], ncol=cs[2]-cs[1])
for(i in 1:(dim(boltingData)[1])){
  img[boltingData[i,"boltDay"]-rs[1], boltingData[i, "leafcount"]-cs[1]]=1+img[boltingData[i,"boltDay"]-rs[1], boltingData[i, "leafcount"]-cs[1]]
}
#redo, formatted for ggplot bubble plot
img=do.call(cbind, list("boltDay"=c(0), "leafcount"=c(0), "count"=c(0)))
for(i in 1:(dim(boltingData)[1])){
  if(length(which(img[,"boltDay"]==boltingData[i, "boltDay"] &
                  img[,"leafcount"]==boltingData[i, "leafcount"]))==1){
    print("update")
    img[which(img[,"boltDay"]==boltingData[i, "boltDay"] &
                img[,"leafcount"]==boltingData[i, "leafcount"]),"count"]=1+img[which(img[,"boltDay"]==boltingData[i, "boltDay"] &
                                                                                       img[,"leafcount"]==boltingData[i, "leafcount"]),"count"]
  }else{
    print("add row")
    img=rbind(img, c(boltingData[i,"boltDay"], boltingData[i,"leafcount"], 1))
  }
}
img=img[c(-1),]
img[,"boltDay"]=21+img[,"boltDay"]
library(ggplot2)

ggplot(data.frame(img), aes(x=boltDay, y=leafcount, size = count)) +
  geom_point(colour="orange")+xlab("Days to bolting")+ylab("Leaf count") + theme_minimal() 
ggsave("plots/03_01_bubbleplotLeaf.png", width=5, height=3.5, dpi=500, unit='in')

#supplemental figure
ggplot(data.frame(boltingData), aes(x=daysInLD, y=leafcount))+
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) + ylab("Leaf count") + xlab("Days in LD") + geom_jitter(shape=16, position=position_jitter(0.2)) 
ggsave("plots/03_02_splitByTreatmentLeaf.png", width=5.5, height=4, dpi=500, unit='in')

ggplot(data.frame(boltingData), aes(x=daysInLD, y=boltDay))+
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) + ylab("Day of bolting") + xlab("Days in LD") + geom_jitter(shape=16, position=position_jitter(0.2))
ggsave("plots/03_03_splitByTreatmentDaysToBolt.png", width=5, height=3.5, dpi=500, unit='in')


#B) boxplots of bolting times by genotype
mutantBolt=read.csv("data/FloweringTime_Arabidopsis_elf3Ws2.csv", header=T)
mutantBolt[,"Genotype"]=factor(mutantBolt[,"Genotype"], levels = c("Ws-2", "elf3-4"))
mutantBolt[,"Condition"]=factor(mutantBolt[,"Condition"], levels = c("SD", "LD"))

ggplot(mutantBolt, aes(x=Genotype, y=DaysToFlower, fill=Condition)) +
  geom_boxplot()+ scale_fill_manual(values=c("SD"="snow4", "LD"="snow2")) +theme_minimal()+ scale_x_discrete(labels=expression("Ws-2" = "Ws-2", "elf3-4" = italic("elf3-4"))) + ylab("Days until bolting") 
ggsave("plots/03_04_daysUntilBolting.png", width=3.5, height=3.5 , dpi=500, unit='in')

mutantBolt$merge=paste(mutantBolt$Genotype, mutantBolt$Condition)
TukeyHSD(aov(DaysToFlower ~ merge, mutantBolt))



#C) standard deviation by genotypes bar plot, ordered by mean
table(mutantBolt[,"Genotype"])
mutantBolt$GenoCond=factor(paste(mutantBolt$Genotype, mutantBolt$Condition), levels=c("Ws-2 SD", "Ws-2 LD", "elf3-4 SD", "elf3-4 LD"))
genoconds=unique(mutantBolt$GenoCond)
sds=sapply(genoconds, function(i){sd(mutantBolt[which(mutantBolt$GenoCond==i),"DaysToFlower"])})
temp=data.frame("Genotype"=genoconds, "StDev"=sds)
ggplot(temp, aes(Genotype, StDev))+
  geom_col(fill=c("lightgreen","lightblue", "forestgreen", "cornflowerblue")) + xlab("Genotype/Condition")+ ylab("Standard deviation") + theme_minimal() + scale_x_discrete(labels=expression("Ws-2" = "Ws-2", "elf3-4" = italic("elf3-4")))
ggsave("plots/03_05_daysUntilBoltingStDev.png", width=3.5, height=3.5, dpi=500, unit='in')

#D) bubble plot for mutants:
#formatted for ggplot bubble plot
img=do.call(cbind, list("boltDay"=c(0), "leafcount"=c(0), "count"=c(0)))
for(i in 1:(dim(mutantBolt)[1])){
  if(length(which(img[,"boltDay"]==mutantBolt[i, "DaysToFlower"] &
                  img[,"leafcount"]==mutantBolt[i, "LeafNumber"]))==1){
    print("update")
    img[which(img[,"boltDay"]==mutantBolt[i, "DaysToFlower"] &
                img[,"leafcount"]==mutantBolt[i, "LeafNumber"]),"count"]=1+img[which(img[,"boltDay"]==mutantBolt[i, "DaysToFlower"] &
                                                                                        img[,"leafcount"]==mutantBolt[i, "LeafNumber"]),"count"]
  }else{
    print("add row")
    img=rbind(img, c(mutantBolt[i,"DaysToFlower"], mutantBolt[i,"LeafNumber"], 1))
  }
}
img=img[c(-1),]
library(ggplot2)

ggplot(data.frame(img), aes(x=boltDay, y=leafcount, size = count)) +
  geom_point(colour="orange")+xlab("Days to bolting")+ylab("Leaf count") + theme_minimal() 
ggsave("plots/03_06_bubbleplotLeafMutant.png", width=5, height=3.5, dpi=500, unit='in')

#D) relationship between leaf count and bolting time in mutants
ggplot(mutantBolt, aes(x=Genotype, y=LeafNumber, fill=Condition)) +
  geom_boxplot()+ scale_fill_manual(values=c("SD"="snow4", "LD"="snow2"))+ theme_minimal()+ scale_x_discrete(labels=expression("Ws-2" = "Ws-2", "elf3-4" = italic("elf3-4"))) + ylab("Leaf count at bolting") 

ggsave("plots/03_07_leafCount.png", width=3.5, height=3.5, dpi=500, unit='in')

TukeyHSD(aov(LeafNumber ~ merge, data=mutantBolt))

#E) relationship between standard deviation of leaf count and bolting time

sds=sapply(genoconds, function(i){sd(mutantBolt[which(mutantBolt$GenoCond==i),"LeafNumber"], na.rm=T)})
temp=data.frame("Genotypecondition"=genoconds, "StDev"=sds)
ggplot(temp, aes(Genotypecondition, StDev))+
  geom_col(fill=c("lightgreen","lightblue", "forestgreen", "cornflowerblue")) + xlab("Genotype/Condition") + ylab("Standard deviation") + theme_minimal() + scale_x_discrete(labels=expression("Ws-2" = "Ws-2", "elf3-4" = italic("elf3-4")))                                                                                                          



ggsave("plots/03_08_leafCountStDev.png", width=3.5, height=3.5, dpi=500, unit='in')





