#######################################
#
# Step 6: Look at the impact of plant size on osmotic stress response 
#
#######################################



library(lme4)
library(lmerTest)
library(ggplot2)
library(ggpubr)


#######################
#
# Make supplemental figures for preliminary experiments
# These were done on young plants, before there were size differences between genotypes

#Load initial stress data, unsynchronised
prelimData=read.csv("data/PrelimStressExperiment.csv", header=T)
#filter ones without CFW measurements
prelimData=prelimData[which(nchar(prelimData[,"CFW"])>0),]
prelimData=prelimData[-1,]
for(i in c(8:27)){
  prelimData[,i]=as.numeric(prelimData[,i])
}


renames=c("Bowman", "eam10.m (Hvlux)", "eam5.x (Hvphyc)", "eam8.k (Hvelf3)", "eam8.w (Hvelf3)")
names(renames)=c("b.wt", "b.284", "b.285", "b.289", "b.290")

prelimData$plant=as.factor(prelimData$plant)
prelimData$plant=factor(prelimData$plant, labels=c("Antonella", "eam10.m (Hvlux)", "eam5.x (Hvphyc)", "eam8.k (Hvelf3)", "eam8.w (Hvelf3)", "Bowman", 'plant'))
prelimData$PEG=as.factor(prelimData$PEG)
prelimData$Temp=as.factor(prelimData$Temp)
prelimData$IFW=as.numeric(prelimData$IFW)
prelimData$CFW=as.numeric(prelimData$CFW)

prelimData=prelimData[which(!is.na(prelimData$IFW) & !is.na(prelimData$CFW)),]

a=ggplot(prelimData[which(prelimData$Temp=="18°C"),], aes(x=IFW, y=CFW, col=PEG))+geom_point()+theme_classic(base_size=15)+xlab("Pre-stress fresh weight (g)")+ylab("Change in fresh weight (g)")+
  geom_smooth(method = "lm", 
              se = F,
              aes(color = factor(PEG))) +ylim(c(-1, 2))


b=ggplot(prelimData[which(prelimData$Temp=="18°C"),], aes(x=plant, y=CFW, col=PEG))+geom_boxplot()+theme_classic(base_size=15)+xlab("Genotype")+ylab("Change in fresh weight (g)") +theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+geom_point(position = position_jitterdodge(), alpha=0.3)

m1=aov(CFW~plant,  data=prelimData[which(prelimData$Temp=="18°C" & prelimData$PEG=="15%"),])
TukeyHSD(m1)
#A AB B B B AB 

m1=aov(CFW~plant,  data=prelimData[which(prelimData$Temp=="18°C" & prelimData$PEG=="0%"),])
TukeyHSD(m1)

c=ggplot(prelimData[which(prelimData$Temp=="18°C"),], aes(x=plant, y=IFW, col="plant"))+geom_boxplot()+theme_classic(base_size=15)+xlab("Genotype")+ylab("Pre-stress fresh weight (g)") +theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+geom_point(position = position_jitterdodge(), alpha=0.3)

m1=aov(IFW~plant,  data=prelimData[which(prelimData$Temp=="18°C"),])
TukeyHSD(m1)
#A-8.w, B-A, 8.w-5.x
#An  A
#10  ABC
#5x  AB
#8k  ABC
#8w  C
#Bo  BC
#A, AB, A, AB, B, AB

figure <- ggarrange(c, b, a,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
figure
ggsave("plots/06_01_unsynchronisedStress.png", width=12, height=5, dpi=600, units= 'in')






#Load synchronised stress data

stressData=read.csv(file='data/StressTestCombined.csv', header=T)
ids=which(stressData$Initalweight>0.45 & stressData$Initalweight<0.55 & stressData$WeightPreStress>0 & stressData$Totalfreshweight>0 & stressData$Totaldryweight>0)
stressData=stressData[ids,]

stressData$Genotype=factor(stressData$Genotype, labels = c("eam8.k (Hvelf3)", "eam8.w (Hvelf3)", "Bowman"))


#weight change
stressData$weightChange=as.numeric(stressData$Totalfreshweight)-as.numeric(stressData$WeightPreStress)
stressData$weightChangePercent=(as.numeric(stressData$Totalfreshweight)-as.numeric(stressData$WeightPreStress))/as.numeric(stressData$WeightPreStress)

scatter=ggplot(data=stressData[which(!is.na(stressData$Green.yellow)),], aes(x=as.numeric(WeightPreStress), y=as.numeric(weightChange), shape=Genotype, colour=Green.yellow))+geom_point()+theme_classic(base_size=15)+xlab("Pre-stress fresh weight (g)")+ylab("Change in fresh weight (g)")+
  geom_smooth(method = "lm", aes(lty=Genotype), colour='grey', 
              se = F) + scale_colour_manual(drop=F, labels=c("green", "yellow-green", "yellow"), values=c("darkgreen", 'limegreen', 'goldenrod'))+labs(colour='Leaf colour')

boxBf=ggplot(data=stressData[which(!is.na(stressData$Green.yellow)),], aes(x=Genotype, y=as.numeric(weightChange)))+geom_boxplot(colour="black")+theme_classic(base_size=15)+xlab("Initial fresh weight (g)")+ylab("Change in fresh weight (g)")+theme_classic(base_size=15)+xlab("Genotype")+ylab("Change in fresh weight (g)") +theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+ geom_jitter(shape=16, color='orange', position=position_jitter(0.2))
#a=ggplot(data=stressData[which(!is.na(stressData$Green.yellow)),], aes(x=Genotype, y=as.numeric(WeightPreStress), col="grey"))+geom_boxplot()+theme_classic(base_size=15)+xlab("Pre-Stress fresh weight (g)")+ylab("Change in fresh weight (g)")+theme_classic(base_size=15)+xlab("Genotype")+ylab("Pre-stress fresh weight (g)") +theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+geom_point(position = position_jitterdodge())
stressData$minusPredicted=stressData$weightChange-predict(lm(weightChange ~ Initalweight, data=stressData))
boxAf=ggplot(data=stressData[which(!is.na(stressData$Green.yellow)),], aes(x=Genotype, y=as.numeric(minusPredicted)))+geom_boxplot(colour="black")+theme_classic(base_size=15)+xlab("Initial fresh weight (g)")+ylab("Change in fresh weight (g)")+theme_classic(base_size=15)+xlab("Genotype")+ylab("Normalised change in fresh weight (g)") +theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+ geom_jitter(shape=16, color='orange', position=position_jitter(0.2))

subTable=stressData[which(!is.na(stressData$Green.yellow)),]
##do TukeyHSD tests
TukeyHSD(aov(weightChange ~ Genotype, data=subTable))
TukeyHSD(aov(minusPredicted ~ Genotype, data=subTable))

figure <- ggarrange(scatter, boxBf, boxAf,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1, widths=c(1, 0.5, 0.5))
figure
ggsave("plots/06_02_synchronisedStress.png", width=12, height=5, dpi=600, units="in")


############
# As per reviewer comments, perform a linear mixed effects model

temp_lmer=lmerTest::lmer('weightChange ~ WeightPreStress + WeightPreStress|Genotype', data=stressData)
temp_lm=lm('weightChange ~ WeightPreStress + Genotype', data=stressData)
summary(temp_lm)
aov=anova(temp_lmer)
show_tests(aov, fractions = TRUE)

stressData$Age_one=factor(stressData$Age_one, levels= c("8", "9", "10", "11", "12"))

#key supplement:
s3=ggplot(data=stressData[which(!is.na(stressData$Green.yellow)),], aes(x=Age_one, y=as.numeric(weightChange), col="grey"))+geom_boxplot()+theme_classic(base_size=15)+xlab("Initial fresh weight (g)")+ylab("Change in fresh weight (g)")+theme_classic(base_size=15)+xlab("Genotype")+ylab("Change in fresh weight (g)") +theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+geom_point(aes(colour=Green.yellow, shape=Genotype), position = position_jitterdodge()) + scale_colour_manual(values=c("darkgreen", 'limegreen', 'grey', 'goldenrod'))
s1=ggplot(data=stressData[which(!is.na(stressData$Green.yellow)),], aes(x=Genotype, y=as.numeric(Initalweight), col="grey"))+geom_boxplot()+theme_classic(base_size=15)+xlab("Initial fresh weight (g)")+ylab("Change in fresh weight (g)")+theme_classic(base_size=15)+xlab("Genotype")+ylab("Initial fresh weight (g)") +theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+geom_point(position = position_jitterdodge())
s2=ggplot(data=stressData[which(!is.na(stressData$Green.yellow)),], aes(x=as.numeric(Initalweight), y=as.numeric(weightChange), shape=Genotype, colour=Green.yellow))+geom_point()+theme_classic(base_size=15)+xlab("Initial fresh weight (g)")+ylab("Change in fresh weight (g)")+
  geom_smooth(method = "lm", aes(colour='grey', lty=Genotype),
              se = F) + scale_colour_manual(values=c("darkgreen", 'limegreen', 'grey', 'goldenrod'))



m1=aov(weightChange~Genotype,  data=stressData[which(!is.na(stressData$Green.yellow)),])
TukeyHSD(m1)
#A A B (0.0275, 0.0177)

m1=aov(as.numeric(Initalweight)~Genotype,  data=stressData[which(!is.na(stressData$Green.yellow)),])
TukeyHSD(m1)
#none significant

m1=aov(as.numeric(WeightPreStress)~Genotype,  data=stressData[which(!is.na(stressData$Green.yellow)),])
TukeyHSD(m1)
#none significant




