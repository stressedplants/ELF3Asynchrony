################################
#
# Step 5: Look at heterogeneity in traits in barley closer to maturity
#
################################

library(ggpubr)
library(ggplot2)
cols=c(rgb(1, 0.5,0.5,1),rgb(1, 0.7,0.2,1), rgb(0.8, 0.6,1,1))
floweringData=read.csv(file='data/barleyfloweringtimes.csv', header=T)

plot(floweringData$Days.to.flower, floweringData$Shoot.fresh.weight..g.)

plot(floweringData$Days.to.flower, floweringData$Leaf.number, pch=19, col='orange', xlab='Days to flower', ylab='Leaf count')

plot(floweringData$Days.to.flower, floweringData$Tiller.number, pch=19, col='orange', xlab='Days to flower', ylab='Tiller count')

plot(floweringData$Days.to.flower, floweringData$Head.number, pch=19, col='orange', xlab='Days to flower', ylab='Head count')




aG=ggplot(floweringData, aes(x=Days.to.flower, y=Shoot.fresh.weight..g.)) +
  geom_point(colour='orange')+xlab("Days to flowering")+ylab("Shoot fresh weight (g)") + theme_minimal() 

aP=cor.test(floweringData$Days.to.flower, floweringData$Shoot.fresh.weight..g.)
#p-value 3.28e-14, R=0.795

bG=ggplot(floweringData, aes(x=Days.to.flower, y=Leaf.number)) +
  geom_point(colour='orange')+xlab("Days to flowering")+ylab("Leaf count") + theme_minimal() 

bP=cor.test(floweringData$Days.to.flower, floweringData$Leaf.number)
#p-value 2.707e-12, R=0.756

cG=ggplot(floweringData, aes(x=Days.to.flower, y=Head.number)) +
  geom_point(colour='orange')+xlab("Days to flowering")+ylab("Head count") + theme_minimal() 

cP=cor.test(floweringData$Days.to.flower, floweringData$Head.number)
#p-value 0.000502, R=0.436

dG=ggplot(floweringData, aes(x=Days.to.flower, y=Tiller.number)) +
  geom_point(colour='orange')+xlab("Days to flowering")+ylab("Tiller count") + theme_minimal() 

dP=cor.test(floweringData$Days.to.flower, floweringData$Tiller.number)
#P=0.01 R=0.324

pTogether=c(aP$p.value, bP$p.value, cP$p.value, dP$p.value)
p.adjust(pTogether)
#1.314325e-13 8.119701e-12 1.003201e-03 1.163759e-02


#genotypes
floweringData$Genotype=factor(floweringData$Genotype, levels=c('B289', 'B290', 'B284', 'B285', 'Bowman'))

renames=c("Bowman", "eam10.m (Hvlux)", "eam5.x (HvPHYC)", "eam8.k (Hvelf3)", "eam8.w (Hvelf3)", "Antonella")
names(renames)=c("Bowman", "B284", "B285", "B289", "B290", "Antonella")


eG=ggplot(floweringData, aes(x=Genotype, y=Leaf.number)) +
  geom_boxplot()+ geom_jitter(shape=16, color='orange', position=position_jitter(0.2))+theme_minimal()+ theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1)) +
  scale_x_discrete(labels=expression("Bowman"="Bowman", "B284"="eam10.m (Hvlux)", "B285"="eam5.x (HvPHYC)", "B289"="eam8.k (Hvelf3)", "B290"="eam8.w (Hvelf3)")) + ylab("Leaf count") 

fG=ggplot(floweringData, aes(x=Genotype, y=Shoot.fresh.weight..g.)) +
  geom_boxplot()+ geom_jitter(shape=16, color='orange', position=position_jitter(0.2))+theme_minimal()+ theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1)) +
  scale_x_discrete(labels=expression("Bowman"="Bowman", "B284"="eam10.m (Hvlux)", "B285"="eam5.x (HvPHYC)", "B289"="eam8.k (Hvelf3)", "B290"="eam8.w (Hvelf3)")) + ylab("Shoot fresh weight (g)") 

gG=ggplot(floweringData, aes(x=Genotype, y=Head.number)) +
  geom_boxplot()+ geom_jitter(shape=16, color='orange', position=position_jitter(0.2))+theme_minimal()+ theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1)) +
  scale_x_discrete(labels=expression("Bowman"="Bowman", "B284"="eam10.m (Hvlux)", "B285"="eam5.x (HvPHYC)", "B289"="eam8.k (Hvelf3)", "B290"="eam8.w (Hvelf3)")) + ylab("Head count") 


hG=ggplot(floweringData, aes(x=Genotype, y=Tiller.number)) +
  geom_boxplot()+ geom_jitter(shape=16, color='orange', position=position_jitter(0.2))+theme_minimal()+ theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1)) +
  scale_x_discrete(labels=expression("Bowman"="Bowman", "B284"="eam10.m (Hvlux)", "B285"="eam5.x (HvPHYC)", "B289"="eam8.k (Hvelf3)", "B290"="eam8.w (Hvelf3)")) + ylab("Tiller count") 

figure <- ggarrange(aG, fG, bG, eG, cG, gG, dG, hG,
                    labels = c("A", "E", "B", "F", "C", "G", "D", "H"),
                    ncol = 2, nrow = 4)
figure
ggsave("plots/05_02_maturationPtBarley.png", width=10, height=14)

#######################################
#
# Step 5: barley at maturation point

day60=read.csv("data/barley_greenhouse60days.csv", header=T)[1:147,]

sapply(unique(day60$Plant), function(i){
  day60[which(day60$Plant==i), "MaxZS"]
})

renames=c("Bowman", "eam10.m (Hvlux)", "eam5.x (HvPHYC)", "eam8.k (Hvelf3)", "eam8.w (Hvelf3)")
names(renames)=c("b.wt", "b.284", "b.285", "b.289", "b.290")

day60$Plant=factor(renames[day60$Plant], levels=renames)
#
ggplot(day60, aes(x=MaxZS, colour = Plant, fill = Plant)) + 
  geom_histogram(position = 'dodge', binwidth = 0.5, aes(x=MaxZS,y=(0.5)*..density..)) +theme_minimal()+xlab("Zadok's growth stage (Maximum)")+ylab("Proportion")+ylim(0,1)
ggsave("plots/05_01_Zadok.pdf", width=7, height=5)

#day60$Plant=factor(day60$Plant, levels=c("b.wt", "b.284", 'b.289', "b.285", "b.290"))

#add max PH
maxes=apply(day60, 1, function(i){
  max(i[6:19], na.rm=T)
})

day60$maxPh=as.numeric(maxes)

ggplot(day60, aes(x=Plant, y = MaxZS)) + 
  geom_boxplot()

#Plants at each stage (ignoring genotype): relationship between tiller count and maximum tiller height
ggplot(day60, aes(x=T, y=PHAv, colour=Plant))+geom_point()

ggplot(day60[which(day60$MaxZS==6),], aes(x=T, y=PHAv, colour=Plant))+geom_point()

