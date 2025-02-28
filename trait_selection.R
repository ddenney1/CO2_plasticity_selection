---
title: "Trait expression and natural selection in response to a factorial manipulation of temperature and carbon dioxide in the growth chamber"
---

## Preparation and read in data
library(tidyverse)
library(ggplot2)
library(car)
library(lme4)
library(effects)
require(lmerTest)
require(MuMIn)
library(performance)
require(visreg)
library(emmeans)
library(coxme)
library(DHARMa)
library(glmmTMB)
library(ggeffects)
library(ggpubr)

setwd("~/Documents/personnel/Denney/selection")

# Fitness and trait data
fitness <- read.csv("fitness_data.csv", header = T, fileEncoding = "UTF-8-BOM")
sapply(fitness,class)
fitness$Treatment<-as.factor(fitness$Treatment)
fitness$Temperature<-as.factor(fitness$Temperature)
fitness$CO2<-as.factor(fitness$CO2)
fitness$Genotype <-as.factor(fitness$Genotype)
fitness$Population <-as.factor(fitness$Population)
fitness$Tray <-as.factor(fitness$Tray)
fitness$Block<-paste(fitness$Tray, fitness$Treatment, sep="_")
fitness$Block<-as.factor(fitness$Block)
fitness$Cohort<-as.factor(fitness$Cohort)

#Standardize source elevation
fitness$S_elev<- (fitness $Elevation - mean(fitness $Elevation, na.rm = TRUE)) / sd(fitness $Elevation,na.rm = TRUE)

###Calculating phenological traits
#Vernalization ended 4/10/22, ordinal day 100. All dates in the experiment started with day of transplanting (OD_transplant), which was in 2011, so we need to include (365-OD_transplant) to account for the time alive in 2021
fitness$Days_to_Flower<- fitness$Date_PostVern_Flowering-100-(365-fitness$OD_transplant)

#Duration of flowering is the number of days between the beginning and end of flowering
fitness$Flowering_Duration<-fitness$Date_silique-fitness$Date_PostVern_Flowering

#Some plants had multiple stems, which we summed to get total plant height at flowering. If this sum returns a value of 0, we change that to NA because those plants did not have height data recorded
fitness$Height_flowering <- rowSums(fitness[,c("Height1_flowering", "Height2_flowering","Height3_flowering", "Height4_flowering","Height5_flowering", "Height6_flowering","Height7_flowering", "Height8_flowering")], na.rm=TRUE)
fitness["Height_flowering"][fitness["Height_flowering"] == 0] <- NA

##Examining histograms of the traits
hist(fitness$Days_to_Flower)
hist(fitness$Flowering_Duration)
hist(fitness$Height_flowering)
hist(fitness$SLA_rosette_June22)
hist(fitness$LDMC_June22_Rosette_DW.FW_g)
hist(fitness$Root_Shoot)

##Change baseline treatment levels to improve plotting
fitness$Temperature<-factor(fitness$Temperature, levels = c("Low","High"))
fitness$CO2<-factor(fitness$CO2, levels = c("Low","High"))


##### First: Genetic clines and plasticity in functional and phenological traits #####
##############################################
######Root_Shoot ##########
##############################################
RS<-subset(fitness,Root_Shoot!="NA")
hist(fitness$Root_Shoot)
min(fitness$Root_Shoot,na.rm=TRUE)
Rootshoot <- glmmTMB(Root_Shoot ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype) , family=Gamma(link="log"),data = RS)
Anova(Rootshoot,type="III")
emmip(Rootshoot, ~ Temperature, type="response", CIs=TRUE)
treat_means<-emmeans(Rootshoot, specs=c("Temperature","CO2"))
pairs(treat_means)


#Use the DHARMa package to examine the residuals, which are good
simulationOutput <- simulateResiduals(fittedModel= Rootshoot, plot = T, re.form = NULL)
plotResiduals(simulationOutput, RS$S_elev)
plotResiduals(simulationOutput, RS$Temperature)
plotResiduals(simulationOutput, RS$CO2)
performance::check_model(Rootshoot)

##Test the random effects by running models without each of them.
Rootshoot_nogeno <- glmmTMB(Root_Shoot ~ Cohort+Temperature*CO2* S_elev + (1|Block)  , family=Gamma(link="log"),data = fitness)
Rootshoot_noblock <- glmmTMB(Root_Shoot ~ Cohort+Temperature*CO2* S_elev +(1| Genotype) , family=Gamma(link="log"),data = fitness)
anova(Rootshoot,Rootshoot_nogeno)
anova(Rootshoot,Rootshoot_noblock)

##Non-significant
rootshoot_cline<- ggpredict(Rootshoot, terms = c("S_elev[all]"),type="fixed")
plot(rootshoot_cline, show_data=TRUE, show_title =FALSE, show_legend=TRUE, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_y_continuous("Root:Shoot ratio")+ scale_x_continuous("Source elevation (standardized)")

##Box plot
RS_box <-ggplot(fitness, aes(x = CO2, y = Root_Shoot, fill = Temperature,shape= Temperature)) +geom_boxplot(outlier.shape = NA) +xlab("CO2")+ scale_y_continuous("Root to shoot ratio") +geom_point(aes(shape=factor(Temperature)), size = 2,position = position_jitterdodge(0.3))

Fig1a <-RS_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),
                                               axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(),   panel.grid.minor=element_blank(),legend.position = "right")+ scale_x_discrete(labels=c("Low", "High")) +  scale_fill_manual(values = c("#56B4E9","#D55E00"), name = "Temperature", labels = c("Low","High"))+scale_shape_manual(values=c(21,24), name = "Temperature", labels = c("Low","High"))
Fig1a


##############################################
######SLA ##########
##############################################
SLA <- glmmTMB(SLA_rosette_June22 ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype) , family=Gamma(link="log"),data = fitness)
Anova(SLA,type="III")
#Use the DHARMa package to examine the residuals, which are good
simulationOutput <- simulateResiduals(fittedModel= SLA, plot = T, re.form = NULL)
plotResiduals(simulationOutput, fitness$S_elev)
plotResiduals(simulationOutput, fitness$Temperature)
plotResiduals(simulationOutput, fitness$CO2)
performance::check_model(SLA)

##Test the signficance of random effects
SLA_nogeno <- glmmTMB(SLA_rosette_June22 ~ Cohort+Temperature*CO2* S_elev + (1|Block)  , family=Gamma(link="log"),data = fitness)
SLA_noblock <- glmmTMB(SLA_rosette_June22 ~ Cohort+Temperature*CO2* S_elev +(1| Genotype) , family=Gamma(link="log"),data = fitness)
anova(SLA,SLA_nogeno)
anova(SLA,SLA_noblock)


##############################################
######Leaf dry matter content ##########
##############################################
hist(fitness$LDMC_June22_Rosette_DW.FW_g)
min(fitness$LDMC_June22_Rosette_DW.FW_g, na.rm=TRUE)
max(fitness$LDMC_June22_Rosette_DW.FW_g, na.rm=TRUE)

leaf_dry<-subset(fitness,LDMC_June22_Rosette_DW.FW_g!="NA")

##normal distribution
LDMC <- glmmTMB(LDMC_June22_Rosette_DW.FW_g ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype) , data = leaf_dry)
Anova(LDMC,type="III")
treat_means<-emmeans(LDMC, specs=c("Temperature","CO2"))
pairs(treat_means)
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= LDMC, plot = T, re.form = NULL)
testQuantiles(simulationOutput)
plot(simulationOutput, quantreg = T)
plotResiduals(simulationOutput, leaf_dry$S_elev)
plotResiduals(simulationOutput, leaf_dry$Temperature)
plotResiduals(simulationOutput, leaf_dry$CO2)
performance::check_model(LDMC)

## Slope of the relationship between source elevation, CO2, and temperature .
emtrends(LDMC, specs = c("Temperature","CO2"), var = "S_elev")


##Test random effects:
LDMC_nogeno <- glmmTMB(LDMC_June22_Rosette_DW.FW_g ~ Cohort+Temperature*CO2* S_elev + (1|Block) , data = leaf_dry)
LDMC_noblock <- glmmTMB(LDMC_June22_Rosette_DW.FW_g ~ Cohort+Temperature*CO2* S_elev + (1| Genotype) , data = leaf_dry)
anova(LDMC,LDMC_nogeno)
anova(LDMC,LDMC_noblock)

##Alternative statistical distributions

# #Gamma regression
# LDMC_gamma <- glmmTMB(LDMC_June22_Rosette_DW.FW_g ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype) , data = fitness, family=Gamma(link="log"))
# Anova(LDMC_gamma,type="III")
# #Use the DHARMa package to examine the residuals, which are good
# simulationOutput <- simulateResiduals(fittedModel= LDMC_gamma, plot = T, re.form = NULL)
#
# ##beta regression
# LDMC_beta <- glmmTMB(LDMC_June22_Rosette_DW.FW_g ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype) , data = fitness, family=beta_family(link="logit"))
# Anova(LDMC_beta,type="III")
# #Use the DHARMa package to examine the residuals, which are good
# simulationOutput <- simulateResiduals(fittedModel= LDMC_beta, plot = T, re.form = NULL)
#
# model.sel(LDMC,LDMC_beta,LDMC_gamma)

##Box plot
LDMC_box <-ggplot(fitness, aes(x = CO2, y = LDMC_June22_Rosette_DW.FW_g, fill = Temperature,shape= Temperature)) +geom_boxplot(outlier.shape = NA) +xlab("CO2")+ scale_y_continuous("Leaf dry matter content") +geom_point(aes(shape=factor(Temperature)), size = 2,position = position_jitterdodge(0.3))

Fig1b <-LDMC_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),
                                                     axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(),   panel.grid.minor=element_blank(),legend.position = "right")+ scale_x_discrete(labels=c("Low", "High")) +  scale_fill_manual(values = c("#56B4E9","#D55E00"), name = "Temperature", labels = c("Low","High"))+scale_shape_manual(values=c(21,24), name = "Temperature", labels = c("Low","High"))
Fig1b

cols= c("#56B4E9","#D55E00")

  plotLDMC<-visregList(visreg(LDMC,"S_elev", by="Temperature",cond=list("CO2"="Low"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),
                  visreg(LDMC,"S_elev", by="Temperature",cond=list("CO2"="High"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), collapse=TRUE)
  Fig1c<-ggplot(plotLDMC $fit, aes(S_elev, visregFit,group= Temperature, colour= Temperature, fill=factor(Temperature))) +
    geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
    geom_line(aes(group=Temperature)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= leaf_dry, aes(x=S_elev, y=LDMC_June22_Rosette_DW.FW_g, color= Temperature, shape=Temperature), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Source elevation (standardized)")+ scale_y_continuous("Leaf dry matter content",limits=c(0,1))+scale_fill_manual(values = cols, name = "Temperature treatment", labels = c("Low","High"))+scale_colour_manual(values = cols, name = "Temperature treatment", labels = c("Low","High"))+facet_wrap(~CO2)
  Fig1c


##############################################
######Days_to_Flower ##########
##############################################

##Ensure that these data are integers
all.equal(fitness$Days_to_Flower, as.integer(fitness$Days_to_Flower))
min(fitness$Days_to_Flower,na.rm=TRUE)
max(fitness$Days_to_Flower,na.rm=TRUE)

hist(fitness$Days_to_Flower)

flowering<-subset(fitness,Days_to_Flower!="NA")

FT_nb1 <- glmmTMB(Days_to_Flower ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype) , data = flowering,family= nbinom1(link="log"))
Anova(FT_nb1,type="III")
simulationOutput <- simulateResiduals(fittedModel= FT_nb1, plot = T, re.form = NULL)
plotResiduals(simulationOutput, flowering$S_elev)
plotResiduals(simulationOutput, flowering$Temperature)
plotResiduals(simulationOutput, flowering$CO2)

##Test random effects
FT_nb1_nogeno <- glmmTMB(Days_to_Flower ~ Cohort+Temperature*CO2* S_elev + (1|Block) , data = flowering,family= nbinom1(link="log"))
FT_nb1_noBlock <- glmmTMB(Days_to_Flower ~ Cohort+Temperature*CO2* S_elev +(1| Genotype) , data = flowering,family= nbinom1(link="log"))
anova(FT_nb1,FT_nb1_nogeno)
anova(FT_nb1,FT_nb1_noBlock)

plot(predictorEffects(FT_nb1, ~ S_elev), partial.residuals=TRUE, type="response")
emmip(FT_nb1, ~ CO2, type="response", CIs=TRUE)
select_ft<- ggpredict(FT_nb1, terms = c("S_elev[all]","CO2", "Temperature"),type="fixed")
plot(select_ft, show_data=TRUE, show_title =FALSE, show_legend=TRUE, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_y_continuous("Flowering phenology")+ scale_x_continuous("Source elevation (standardized)")


## Here are models with alternative statistical distributions.

# FT_lognormal <- glmmTMB(Days_to_Flower ~ Cohort+Temperature*CO2* S_elev+ (1|Block) +(1| Genotype) , data = flowering,family= lognormal(link="log"))
# simulationOutput <- simulateResiduals(fittedModel= FT_lognormal, plot = T, re.form = NULL)
# plotResiduals(simulationOutput, flowering$S_elev)
# plotResiduals(simulationOutput, flowering$Temperature)
# plotResiduals(simulationOutput, flowering$CO2)
# plotResiduals(simulationOutput, flowering$Cohort)
# Anova(FT_lognormal,type="III")
# plot(predictorEffects(FT_lognormal, ~ S_elev), partial.residuals=TRUE, type="response")
# visreg(FT_lognormal,"S_elev", partial = TRUE, rug = FALSE,scale="response")
#
#
#
# FT_gamma <- glmmTMB(Days_to_Flower ~ Cohort+Temperature*CO2* S_elev+ (1|Block) +(1| Genotype) , data = flowering,family= Gamma(link="log"))
# simulationOutput <- simulateResiduals(fittedModel= FT_gamma, plot = T, re.form = NULL)
# plotResiduals(simulationOutput, flowering$S_elev)
# plotResiduals(simulationOutput, flowering$Temperature)
# plotResiduals(simulationOutput, flowering$CO2)
# Anova(FT_gamma,type="III")
#
#
# FT_normal <- glmmTMB(Days_to_Flower ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype) , data = flowering)
# Anova(FT_normal,type="III")
# simulationOutput <- simulateResiduals(fittedModel= FT_normal, plot = T, re.form = NULL)
# plotResiduals(simulationOutput, flowering$S_elev)
# plotResiduals(simulationOutput, flowering$Temperature)
# plotResiduals(simulationOutput, flowering$CO2)
#
#
# FT_poisson <- glmmTMB(Days_to_Flower ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype) , data = flowering,family= poisson(link="log"))
# Anova(FT_poisson,type="III")
# simulationOutput <- simulateResiduals(fittedModel= FT_poisson, plot = T, re.form = NULL)
# plotResiduals(simulationOutput, flowering$S_elev)
# plotResiduals(simulationOutput, flowering$Temperature)
# plotResiduals(simulationOutput, flowering$CO2)
#
#
# FT <- glmmTMB(Days_to_Flower ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype) , data = flowering,family= nbinom2())
# Anova(FT,type="III")
# #Use the DHARMa package to examine the residuals, which are okay, not great
# simulationOutput <- simulateResiduals(fittedModel= FT, plot = T, re.form = NULL)
# testOutliers(simulationOutput,type = 'bootstrap', nBoot=100)
# plotResiduals(simulationOutput, flowering$S_elev)
# plotResiduals(simulationOutput, flowering$Temperature)
# plotResiduals(simulationOutput, flowering$CO2)
# visreg(FT,"S_elev", partial = FALSE, rug = FALSE,scale="response")
#
# model.sel(FT_lognormal,FT,FT_normal,FT_gamma,FT_poisson,FT_nb1)


##Box plot
FT_box <-ggplot(fitness, aes(x = CO2, y = Days_to_Flower, fill = Temperature,shape= Temperature)) +geom_boxplot(outlier.shape = NA) +xlab("CO2")+ scale_y_continuous("Days until flowering") +geom_point(aes(shape=factor(Temperature)), size = 2,position = position_jitterdodge(0.3))

FT_box_plot <-FT_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),
                                               axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(),   panel.grid.minor=element_blank(),legend.position = "right")+ scale_x_discrete(labels=c("Low", "High")) +  scale_fill_manual(values = c("#56B4E9","#D55E00"), name = "Temperature", labels = c("Low","High"))+scale_shape_manual(values=c(21,24), name = "Temperature", labels = c("Low","High"))
FT_box_plot


##############################################
######Flowering_Duration ##########
##############################################
hist(fitness$Flowering_Duration)
duration <- glmmTMB(Flowering_Duration ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype) , data = fitness, family=Gamma("log"))
Anova(duration,type="III")
#Use the DHARMa package to examine the residuals, which are good
simulationOutput <- simulateResiduals(fittedModel= duration, plot = T, re.form = NULL)
treat_means<-emmeans(duration, specs=c("Temperature","CO2"))
pairs(treat_means)

##Test random effects
duration_nogeno <- glmmTMB(Flowering_Duration ~ Cohort+Temperature*CO2* S_elev + (1|Block) , data = fitness, family=Gamma("log"))
duration_noblock <- glmmTMB(Flowering_Duration ~ Cohort+Temperature*CO2* S_elev + (1| Genotype) , data = fitness, family=Gamma("log"))
anova(duration,duration_nogeno)
anova(duration,duration_noblock)

emmip(duration, ~Temperature , type="response", CIs=TRUE)

##Box plot
FD_box <-ggplot(fitness, aes(x = CO2, y = Flowering_Duration, fill = Temperature,shape= Temperature)) +geom_boxplot(outlier.shape = NA) +xlab("CO2")+ scale_y_continuous("Flowering duration") +geom_point(aes(shape=factor(Temperature)), size = 2,position = position_jitterdodge(0.3))

Fig2a <-FD_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),
                                               axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(),   panel.grid.minor=element_blank(),legend.position = "right")+ scale_x_discrete(labels=c("Low", "High")) +  scale_fill_manual(values = c("#56B4E9","#D55E00"), name = "Temperature", labels = c("Low","High"))+scale_shape_manual(values=c(21,24), name = "Temperature", labels = c("Low","High"))
Fig2a


##############################################
######Height_flowering ##########
##############################################
fitness$height<-fitness$Height_flowering/100
hist(fitness$height)

mod_height <- glmmTMB(height ~ Cohort+Temperature*CO2* S_elev + (1|Block) +(1| Genotype) , data = fitness, family=Gamma(link="log"))
Anova(mod_height,type="III")

#Use the DHARMa package to examine the residuals, which are good
simulationOutput <- simulateResiduals(fittedModel= mod_height, plot = T, re.form = NULL)

##Testing random effects
mod_height_nogeno <- glmmTMB(height ~ Cohort+Temperature*CO2* S_elev + (1|Block)  , data = fitness, family=Gamma(link="log"))
mod_height_noblock <- glmmTMB(height ~ Cohort+Temperature*CO2* S_elev +(1| Genotype) , data = fitness, family=Gamma(link="log"))
anova(mod_height,mod_height_nogeno)
anova(mod_height,mod_height_noblock)

fitness$height<-fitness$Height_flowering/100
hist(fitness$height)

plotheight<-visregList(visreg(mod_height,"S_elev", by="Temperature",cond=list("CO2"="Low"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
                     visreg(mod_height,"S_elev", by="Temperature",cond=list("CO2"="High"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)
Fig2b<-ggplot(plotheight $fit, aes(S_elev, visregFit,group= Temperature, colour= Temperature, fill=factor(Temperature))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Temperature)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= fitness, aes(x=S_elev, y=height, color= Temperature, shape=Temperature), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Source elevation (standardized)")+ scale_y_continuous("Height at flowering (cm)")+scale_fill_manual(values = cols, name = "Temperature treatment", labels = c("Low","High"))+scale_colour_manual(values = cols, name = "Temperature treatment", labels = c("Low","High"))+facet_wrap(~CO2)
Fig2b

plot(predictorEffects(mod_height, ~ S_elev), partial.residuals=TRUE, type="response")

plot(fitness$Mature_silique_number_total~fitness$Mature_length_siliques_total)



##### Adjust p-values using a FDR correction for multiple testing ####
library(stats)
p <- read.csv("pvalue.csv", header = T, fileEncoding = "UTF-8-BOM")
head(p)

##Drop first two columns. Merge back later
pvalue<-p[, -c(1:2)]
pvalue

round(p, 3)
round(p.adjust(pvalue), 3)
p.adj<-round(p.adjust(pvalue, "BH"), 3)

pvalues<-cbind(p,padj)
write.csv(pvalues,"corrected_pvalues.csv")

#### plotting main figures ####
## requires ggpubr to be loaded and code for the models and plots already run

## Figure 1

figure1 <- ggarrange(Fig1a, Fig1b,Fig1c,
                     labels = c("A", "B","C"),
                     ncol = 3, nrow = 1, common.legend = TRUE, legend="right")
figure1

figure2 <- ggarrange(Fig2a, Fig2b,
                     labels = c("A", "B"),
                     ncol = 2, nrow = 1, common.legend = TRUE, legend="right")
figure2



################################################
### Multivariate selection on functional and phenological traits:
################################################
##First, use the vegetative traits to examine the probability of flowering. We cannot includep phenological traits because plants that failed to flower do not have flowering time data.
# Standardize data for selection
fitness$S_SLA<- (fitness $SLA_rosette_June22 - mean(fitness $SLA_rosette_June22, na.rm = TRUE)) / sd(fitness $SLA_rosette_June22,na.rm = TRUE)
fitness$S_RS <- (fitness $Root_Shoot - mean(fitness $Root_Shoot, na.rm = TRUE)) / sd(fitness $Root_Shoot,na.rm = TRUE)
fitness$S_LDMC <- (fitness $LDMC_June22_Rosette_DW.FW_g - mean(fitness $LDMC_June22_Rosette_DW.FW_g, na.rm = TRUE)) / sd(fitness $LDMC_June22_Rosette_DW.FW_g,na.rm = TRUE)

prob_flowered <- glmmTMB(Lifetime_Flowered ~
                                S_SLA*Temperature*CO2+
                                S_RS*Temperature*CO2+
                                S_LDMC*Temperature*CO2+
                                (1|Genotype)+(1|Block),
                              family=binomial(link="logit"),data=fitness)
Anova(prob_flowered,type="III")

##Test random effects by removing them from the model
prob_flowered_nogeno <- glmmTMB(Lifetime_Flowered ~
                           S_SLA*Temperature*CO2+
                           S_RS*Temperature*CO2+
                           S_LDMC*Temperature*CO2+
                           (1|Block),
                         family=binomial(link="logit"),data=fitness)
anova(prob_flowered,prob_flowered_nogeno)

prob_flowered_noblock <- glmmTMB(Lifetime_Flowered ~
                                  S_SLA*Temperature*CO2+
                                  S_RS*Temperature*CO2+
                                  S_LDMC*Temperature*CO2+
                                  (1|Genotype),
                                family=binomial(link="logit"),data=fitness)
anova(prob_flowered,prob_flowered_noblock)


##Combine vegetative and reproductive traits for fecundity selection
flowered<-subset(fitness,Lifetime_Flowered=="1")
#Restandardized trait values with the subset of data
flowered$S_SLA<- (flowered $SLA_rosette_June22 - mean(flowered $SLA_rosette_June22, na.rm = TRUE)) / sd(flowered $SLA_rosette_June22,na.rm = TRUE)
flowered$S_RS <- (flowered $Root_Shoot - mean(flowered $Root_Shoot, na.rm = TRUE)) / sd(flowered $Root_Shoot,na.rm = TRUE)
flowered$S_LDMC <- (flowered $LDMC_June22_Rosette_DW.FW_g - mean(flowered $LDMC_June22_Rosette_DW.FW_g, na.rm = TRUE)) / sd(flowered $LDMC_June22_Rosette_DW.FW_g,na.rm = TRUE)
flowered$S_Days<-(flowered $Days_to_Flower - mean(flowered $Days_to_Flower, na.rm = TRUE)) / sd(flowered $Days_to_Flower,na.rm = TRUE)
flowered$S_height <-(flowered $Height_flowering - mean(flowered $Height_flowering, na.rm = TRUE)) / sd(flowered $Height_flowering,na.rm = TRUE)
flowered$S_duration <-(flowered $Flowering_Duration - mean(flowered $Flowering_Duration, na.rm = TRUE)) / sd(flowered $Flowering_Duration,na.rm = TRUE)


select_veg_repro <- glmmTMB(Mature_silique_number_total ~ Cohort +
                              S_SLA*Temperature*CO2+I(S_SLA^2)*CO2+
                              S_RS*Temperature*CO2+
                              S_LDMC*Temperature*CO2+
                              S_Days*Temperature*CO2+
                              S_duration*Temperature*CO2+I(S_duration^2)+
                              S_height*Temperature*CO2+
                              (1|Genotype)+(1|Block),
                            family=nbinom1(),data=flowered)
Anova(select_veg_repro,type="III")
#Use the DHARMa package to examine the residuals, which are good
simulationOutput <- simulateResiduals(fittedModel= select_veg_repro, plot = T, re.form = NULL)

##To test the random effects
select_veg_repro_nogeno <- glmmTMB(Mature_silique_number_total ~ Cohort +
                              S_SLA*Temperature*CO2+I(S_SLA^2)*CO2+
                              S_RS*Temperature*CO2+
                              S_LDMC*Temperature*CO2+
                              S_Days*Temperature*CO2+
                              S_duration*Temperature*CO2+I(S_duration^2)+
                              S_height*Temperature*CO2+
                             (1|Block),
                            family=nbinom1(),data=flowered)
anova(select_veg_repro,select_veg_repro_nogeno)
select_veg_repro_noblock <- glmmTMB(Mature_silique_number_total ~ Cohort +
                                     S_SLA*Temperature*CO2+I(S_SLA^2)*CO2+
                                     S_RS*Temperature*CO2+
                                     S_LDMC*Temperature*CO2+
                                     S_Days*Temperature*CO2+
                                     S_duration*Temperature*CO2+I(S_duration^2)+
                                     S_height*Temperature*CO2+
                                     (1|Genotype),
                                   family=nbinom1(),data=flowered)
anova(select_veg_repro,select_veg_repro_noblock)


#Plot the selection surfaces
cols= c("#56B4E9","#D55E00")
selectSLA<-visregList(visreg(select_veg_repro,"S_SLA", by="Temperature",cond=list("CO2"="Low"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
                       visreg(select_veg_repro,"S_SLA", by="Temperature",cond=list("CO2"="High"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)
Fig3a<-ggplot(selectSLA $fit, aes(S_SLA, visregFit,group= Temperature, colour= Temperature, fill=factor(Temperature))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Temperature)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= flowered, aes(y=Mature_silique_number_total, x=S_SLA, color= Temperature, shape=Temperature), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Specific leaf area (cm2/g)")+ scale_y_continuous("Fitness (number of mature fruits)")+scale_fill_manual(values = cols, name = "Temperature treatment", labels = c("Low","High"))+scale_colour_manual(values = cols, name = "Temperature treatment", labels = c("Low","High"))+facet_wrap(~CO2)
Fig3a

selectRS<-visregList(visreg(select_veg_repro,"S_RS", by="Temperature",cond=list("CO2"="Low"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
                      visreg(select_veg_repro,"S_RS", by="Temperature",cond=list("CO2"="High"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)
Fig3b<-ggplot(selectRS $fit, aes(S_RS, visregFit,group= Temperature, colour= Temperature, fill=factor(Temperature))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Temperature)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= flowered, aes(y=Mature_silique_number_total, x=S_RS, color= Temperature, shape=Temperature), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Root to shoot ratio (g/g)")+ scale_y_continuous("Fitness (number of mature fruits)")+scale_fill_manual(values = cols, name = "Temperature treatment", labels = c("Low","High"))+scale_colour_manual(values = cols, name = "Temperature treatment", labels = c("Low","High"))+facet_wrap(~CO2)
Fig3b

selectLDMC<-visregList(visreg(select_veg_repro,"S_LDMC", by="Temperature",cond=list("CO2"="Low"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
                       visreg(select_veg_repro,"S_LDMC", by="Temperature",cond=list("CO2"="High"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)
Fig3c<-ggplot(selectLDMC $fit, aes(S_LDMC, visregFit,group= Temperature, colour= Temperature, fill=factor(Temperature))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Temperature)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= flowered, aes(y=Mature_silique_number_total, x=S_LDMC, color= Temperature, shape=Temperature), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Leaf dry matter content (units)")+ scale_y_continuous("Fitness (number of mature fruits)")+scale_fill_manual(values = cols, name = "Temperature treatment", labels = c("Low","High"))+scale_colour_manual(values = cols, name = "Temperature treatment", labels = c("Low","High"))+facet_wrap(~CO2)
Fig3c

selectFT<-visregList(visreg(select_veg_repro,"S_Days", by="Temperature",cond=list("CO2"="Low"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
                     visreg(select_veg_repro,"S_Days", by="Temperature",cond=list("CO2"="High"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)
Fig4a<-ggplot(selectFT $fit, aes(S_Days, visregFit,group= Temperature, colour= Temperature, fill=factor(Temperature))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Temperature)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= flowered, aes(y=Mature_silique_number_total, x=S_Days, color= Temperature, shape=Temperature), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Timing of first flowering (number of days)")+ scale_y_continuous("Fitness (number of mature fruits)")+scale_fill_manual(values = cols, name = "Temperature treatment", labels = c("Low","High"))+scale_colour_manual(values = cols, name = "Temperature treatment", labels = c("Low","High"))+facet_wrap(~CO2)
Fig4a

select_duration<-visregList(visreg(select_veg_repro,"S_duration", by="Temperature",cond=list("CO2"="Low"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
                            visreg(select_veg_repro,"S_duration", by="Temperature",cond=list("CO2"="High"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)
Fig4b<-ggplot(select_duration $fit, aes(S_duration, visregFit,group= Temperature, colour= Temperature, fill=factor(Temperature))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Temperature)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= flowered, aes(y=Mature_silique_number_total, x=S_duration, color= Temperature, shape=Temperature), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Flowering duration (number of days)")+ scale_y_continuous("Fitness (number of mature fruits)")+scale_fill_manual(values = cols, name = "Temperature treatment", labels = c("Low","High"))+scale_colour_manual(values = cols, name = "Temperature treatment", labels = c("Low","High"))+facet_wrap(~CO2)
Fig4b

selectheight<-visregList(visreg(select_veg_repro,"S_height", by="Temperature",cond=list("CO2"="Low"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
                         visreg(select_veg_repro,"S_height", by="Temperature",cond=list("CO2"="High"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)
Fig4c<-ggplot(selectheight $fit, aes(S_height, visregFit,group= Temperature, colour= Temperature, fill=factor(Temperature))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Temperature)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= flowered, aes(y=Mature_silique_number_total, x=S_height, color= Temperature, shape=Temperature), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Total stem height at first flowering (cm)")+ scale_y_continuous("Fitness (number of mature fruits)")+scale_fill_manual(values = cols, name = "Temperature treatment", labels = c("Low","High"))+scale_colour_manual(values = cols, name = "Temperature treatment", labels = c("Low","High"))+facet_wrap(~CO2)
Fig4c

##Assemble the figures
figure3 <- ggarrange(Fig3a, Fig3b,Fig3c,
                     labels = c("A", "B","C"),
                     ncol = 3, nrow = 1, common.legend = TRUE, legend="none")
figure3

figure4 <- ggarrange(Fig4a, Fig4b,Fig4c,
                     labels = c("A", "B","C"),
                     ncol = 3, nrow = 1, common.legend = TRUE, legend="none")
figure4
