#######################################
#
# Step 6: Look at the impact of plant size on osmotic stress response 
#
#######################################



library(lme4)
library(lmerTest)
library(ggplot2)
library(ggpubr)

# files=list.files(path='data/stress/', pattern='*.txt')
# 
# data=sapply(files, function(i){
#   print(i)
#   a=read.table(paste("data/stress/", i, sep=""), header=T, sep='\t')
#   colnames(a)=tolower(colnames(a))
#   colnames(a)[which(colnames(a)=='ageatstress')]='age_two'
#   colnames(a)[which(colnames(a)=='germination.group')]='group'
#   include=which(!is.na(a[,1]))
#   a[include,]
# })
# 
# 
# # try to merge datasets
# colNames=c('genotype', 'replicate', 'germinationdate', 'age_two', 'group', 'set', 'weight_start', 'weight_mid', 'weight_end_total', 'weight_end_root', 'weight_end_shoot')
# 
# 
# 
# 
# data=lapply(data, function(i){
#   i[,colNames]
# })
# 
# data=lapply(c(1:3), function(i){
#   temp=data[[i]]
#   temp$table=rep(i, length(temp[,1]))
#   temp
# })
# 
# dataMerged=do.call(rbind, data)
# 
# #get rid of plants with missing data or that never exceeded weight threshold
# dataMerged_filtered=dataMerged[which(!is.na(dataMerged[,'set']) & dataMerged[,'weight_start']>0.45 & dataMerged[,'weight_mid']>0.45 & dataMerged[,'weight_end_total']>0.45),]
# 
# 
# 
# table(dataMerged_filtered[,"genotype"]) #between 84 and 101 plants of each genotype
# 
# 
# 
# #####Research questions:
# 
# data2018=dataMerged_filtered
# 
# ##############
# # Is growth pre-treatment dependent on chronological age? (and does this vary by genotype)
# 
# ##### visualisation: 
# #all genotypes reach the growth threshold at a similar day
# ggplot(aes(y=age_two, x = as.factor(genotype), fill=as.factor(group)), data=data2018)+geom_boxplot()
# #intriguingly, HvELF3 variant size depends on the day they surpassed the threshold, but this is not the case for Bowman
# ggplot(aes(y=weight_mid, x = as.factor(genotype), fill=as.factor(age_two)), data=data2018)+geom_boxplot()
# #However, this is a bit of an artefact, because the trend is also found in Bowman within each individual replicate
# ggplot(aes(y=weight_mid, x = as.factor(genotype), fill=as.factor(age_two)), data=data2018[which(data2018[,'table']==1),])+geom_boxplot()
# ggplot(aes(y=weight_mid, x = as.factor(genotype), fill=as.factor(age_two)), data=data2018[which(data2018[,'table']==2),])+geom_boxplot()
# ggplot(aes(y=weight_mid, x = as.factor(genotype), fill=as.factor(age_two)), data=data2018[which(dataMerged_filtered[,'table']==3),])+geom_boxplot()
# 
# 
# #####################################
# # Does size at start correlate with start at mid?
# 
# ggplot(aes(x=weight_start, y=weight_mid, colour = as.factor(genotype)), data=dataMerged_filtered)+geom_point()
# 
# ######################
# # Is growth post-treatment dependent on germination time, initial growth rate, later growth rate? (and does this vary by genotype)
# ggplot(aes(x=weight_mid-weight_start, y=weight_end_total, colour = as.factor(genotype)), data=data2018)+geom_point()
# 
# ggplot(aes(x=weight_mid, y=weight_end_total-weight_mid, colour = as.factor(genotype)), data=data2018)+geom_point()
# 
# ###########################################
# #make into a supplement
# #interesting result: the plants that took longer to grow, tend to have less water loss
# ggplot(aes(x=as.factor(genotype), y=weight_end_total-weight_mid, fill = as.factor(age_two)), data=data2018)+geom_boxplot()
# ggplot(aes(x=as.factor(genotype), y=weight_mid, fill = as.factor(age_two)), data=data2018)+geom_boxplot()
# 
# #simpler model
# modelChange=lmer(weight_end_total-weight_mid ~ 1 + age_two + weight_mid + (1 | genotype) + (1| table), data=data2018)
# summary(modelChange)
# ranova(modelChange)
# 
# #if you normalise the weight loss by the initial weight, do you still see a difference by genotype?
# ggplot(data2018, aes(x=as.factor(genotype), y=weight_end_total-weight_mid))+geom_boxplot()
# 
# mod=lmer(weight_end_total-weight_mid ~ weight_mid + (1|table), data=data2018)
# plot(data2018$weight_mid, predict(mod))
# 
# data2018$changeMinusPredict=data2018$weight_end_total-data2018$weight_mid-predict(mod)
# ggplot(data2018, aes(x=as.factor(genotype), y=changeMinusPredict))+geom_boxplot()
# 
# mod=lmer(weight_end_total-weight_mid ~ weight_mid + as.numeric(age_two) + (1|table), data=data2018)
# plot(data2018$weight_mid, predict(mod))
# 
# data2018$changeMinusPredictAge=data2018$weight_end_total-data2018$weight_mid-predict(mod)
# ggplot(data2018, aes(x=as.factor(genotype), y=changeMinusPredictAge))+geom_boxplot()
# 
# 
# ##Do Tukey HSD tests for these two boxplots
# m1=aov(weight_end_total-weight_mid ~ genotype, data=data2018)
# TukeyHSD(m1)
# 
# m1=aov(changeMinusPredict ~ genotype, data=data2018)
# TukeyHSD(m1)
# 
# m1=aov(changeMinusPredictAge ~ genotype, data=data2018)
# TukeyHSD(m1)
# 
# 

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


renames=c("Bowman", "eam10.m (Hvlux)", "eam5.x (HvPHYC)", "eam8.k (Hvelf3)", "eam8.w (Hvelf3)")
names(renames)=c("b.wt", "b.284", "b.285", "b.289", "b.290")

prelimData$plant=as.factor(prelimData$plant)
prelimData$plant=factor(prelimData$plant, labels=c("Antonella", "eam10.m (Hvlux)", "eam5.x (HvPHYC)", "eam8.k (Hvelf3)", "eam8.w (Hvelf3)", "Bowman", 'plant'))
prelimData$PEG=as.factor(prelimData$PEG)
prelimData$Temp=as.factor(prelimData$Temp)
prelimData$IFW=as.numeric(prelimData$IFW)
prelimData$CFW=as.numeric(prelimData$CFW)

prelimData=prelimData[which(!is.na(prelimData$IFW) & !is.na(prelimData$CFW)),]

a=ggplot(prelimData[which(prelimData$Temp=="18°C"),], aes(x=IFW, y=CFW, col=PEG))+geom_point()+theme_minimal()+xlab("Pre-stress fresh weight (g)")+ylab("Change in fresh weight (g)")+
  geom_smooth(method = "lm", 
              se = F,
              aes(color = factor(PEG))) +ylim(c(-1, 2))


b=ggplot(prelimData[which(prelimData$Temp=="18°C"),], aes(x=plant, y=CFW, col=PEG))+geom_boxplot()+theme_minimal()+xlab("Genotype")+ylab("Change in fresh weight (g)") +theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+geom_point(position = position_jitterdodge(), alpha=0.3)

m1=aov(CFW~plant,  data=prelimData[which(prelimData$Temp=="18°C"),])
TukeyHSD(m1)

c=ggplot(prelimData[which(prelimData$Temp=="18°C"),], aes(x=plant, y=IFW, col='grey'))+geom_boxplot()+theme_minimal()+xlab("Genotype")+ylab("Pre-stress fresh weight (g)") +theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+geom_point(position = position_jitterdodge(), alpha=0.3)

m1=aov(IFW~plant,  data=prelimData[which(prelimData$Temp=="18°C"),])
TukeyHSD(m1)
#A-8.w, 8.w-5.x
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

scatter=ggplot(data=stressData[which(!is.na(stressData$Green.yellow)),], aes(x=as.numeric(WeightPreStress), y=as.numeric(weightChange), shape=Genotype, colour=Green.yellow))+geom_point()+theme_minimal()+xlab("Pre-stress fresh weight (g)")+ylab("Change in fresh weight (g)")+
  geom_smooth(method = "lm", aes(lty=Genotype), colour='grey', 
              se = F) + scale_colour_manual(drop=F, labels=c("green", "yellow-green", "yellow"), values=c("darkgreen", 'limegreen', 'goldenrod'))+labs(colour='Leaf colour')

boxBf=ggplot(data=stressData[which(!is.na(stressData$Green.yellow)),], aes(x=Genotype, y=as.numeric(weightChange)))+geom_boxplot(colour="black")+theme_minimal()+xlab("Initial fresh weight (g)")+ylab("Change in fresh weight (g)")+theme_minimal()+xlab("Genotype")+ylab("Change in fresh weight (g)") +theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+ geom_jitter(shape=16, color='orange', position=position_jitter(0.2))
#a=ggplot(data=stressData[which(!is.na(stressData$Green.yellow)),], aes(x=Genotype, y=as.numeric(WeightPreStress), col="grey"))+geom_boxplot()+theme_minimal()+xlab("Pre-Stress fresh weight (g)")+ylab("Change in fresh weight (g)")+theme_minimal()+xlab("Genotype")+ylab("Pre-stress fresh weight (g)") +theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+geom_point(position = position_jitterdodge())
stressData$minusPredicted=stressData$weightChange-predict(lm(weightChange ~ Initalweight, data=stressData))
boxAf=ggplot(data=stressData[which(!is.na(stressData$Green.yellow)),], aes(x=Genotype, y=as.numeric(minusPredicted)))+geom_boxplot(colour="black")+theme_minimal()+xlab("Initial fresh weight (g)")+ylab("Change in fresh weight (g)")+theme_minimal()+xlab("Genotype")+ylab("Normalised change in fresh weight (g)") +theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+ geom_jitter(shape=16, color='orange', position=position_jitter(0.2))

subTable=stressData[which(!is.na(stressData$Green.yellow)),]
##do TukeyHSD tests
TukeyHSD(aov(weightChange ~ Genotype, data=subTable))
TukeyHSD(aov(minusPredicted ~ Genotype, data=subTable))

figure <- ggarrange(scatter, boxBf, boxAf,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1, widths=c(1, 0.5, 0.5))
figure
ggsave("plots/06_02_synchronisedStress.png", width=12, height=5, dpi=600, units="in")




stressData$Age_one=factor(stressData$Age_one, levels= c("8", "9", "10", "11", "12"))

#key supplement:
s3=ggplot(data=stressData[which(!is.na(stressData$Green.yellow)),], aes(x=Age_one, y=as.numeric(weightChange), col="grey"))+geom_boxplot()+theme_minimal()+xlab("Initial fresh weight (g)")+ylab("Change in fresh weight (g)")+theme_minimal()+xlab("Genotype")+ylab("Change in fresh weight (g)") +theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+geom_point(aes(colour=Green.yellow, shape=Genotype), position = position_jitterdodge()) + scale_colour_manual(values=c("darkgreen", 'limegreen', 'grey', 'goldenrod'))
s1=ggplot(data=stressData[which(!is.na(stressData$Green.yellow)),], aes(x=Genotype, y=as.numeric(Initalweight), col="grey"))+geom_boxplot()+theme_minimal()+xlab("Initial fresh weight (g)")+ylab("Change in fresh weight (g)")+theme_minimal()+xlab("Genotype")+ylab("Initial fresh weight (g)") +theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))+geom_point(position = position_jitterdodge())
s2=ggplot(data=stressData[which(!is.na(stressData$Green.yellow)),], aes(x=as.numeric(Initalweight), y=as.numeric(weightChange), shape=Genotype, colour=Green.yellow))+geom_point()+theme_minimal()+xlab("Initial fresh weight (g)")+ylab("Change in fresh weight (g)")+
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



