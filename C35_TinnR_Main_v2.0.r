######################################################################
## Analysis flow to compare circadian rhythms among Nasonia species ##
######################################################################

## Open entire dataframe and saves as NR (=Nasonia Rhythms)
NR<-read.table("All lines combined_ForR_v2n.txt",header=T)

#### Cleaning the dataset from not interesting columns #####

NRR<-subset(NR,NR$ADRU=='R')

hist(NRR$Tau)
qqnorm(NRR$Tau)

write.table(NRR,'clipboard-128',sep='\t',col.names=NA)

###################################################
## Only main subsets (C35, C52, C71, C84, C86) #### are retained and saved in MainNR
###################################################

NRMain<-subset(NR,NR$Experiment=='C35'|NR$Experiment=='C52'|NR$Experiment=='C71'|NR$Experiment=='C84'|NR$Experiment=='C86')   
NRMain$Experiment=with(NRMain,droplevels(Experiment))   # Drops unused levels (i.e. C61, C76, C83) from the factor 'Experiment'## NRMain = subset of 

## drops AsymC as a control in each experiment
NRMainC<-NRMain[-which(NRMain$Experiment!='C35' & NRMain$Species=='vitripennis'),]  

NRMainCR<-subset(NRMainC,NRMainC$ADRU=='R')

## Prints the dataframe on the clipboard > can be saved ##
write.table(NRMainCR,"clipboard",sep='\t',col.names=NA)

## Rhythmic DD and LL ##
NRMainC_DDr<-subset(NRMainC,NRMainC$Light=='D'& NRMainC$ADRU=='R')        # subset of rhythmic wasps in DD
NRMainC_LLr<-NRMainC[which(NRMainC$Light=='L' & NRMainC$ADRU=='R'),]    # subset of rhythmic wasps in LL

## Boxplots Taus for rhythmic wasps in main experiments ##                                                                                         
boxplot(Tau~Sex*Experiment,data=NRMainC_DDr,plot=TRUE,notch=TRUE,ylim=c(18,33),xlab='Tau in DD',ylab='Tau [h]',main='Tau in DD for rhythmic wasps \n >Species<') ## Boxplot of Tau_DD for every sex and experiment
x11()          # opens another graphic window
boxplot(Tau~Sex*Experiment,data=NRMainC_LLr,plot=TRUE,notch=TRUE,ylim=c(18,33),xlab='Tau in LL',ylab='Tau [h]',main='Tau in LL for rhythmic wasps \n >Species<') ## Boxplot of Tau_LL for every sex and experiment

###############################
#### More specific subsets ####
###############################

# subsets for each experiment and sex in DD #
Dvm<-NRMainC_DDr[NRMainC_DDr$Experiment=='C35' & NRMainC_DDr$Sex=='M','Tau_DD']
Dg52m<-NRMainC_DDr[NRMainC_DDr$Experiment=='C52' & NRMainC_DDr$Sex=='M','Tau_DD']
Dlm<-NRMainC_DDr[NRMainC_DDr$Experiment=='C71' & NRMainC_DDr$Sex=='M','Tau_DD']
Dom<-NRMainC_DDr[NRMainC_DDr$Experiment=='C84' & NRMainC_DDr$Sex=='M','Tau_DD']
Dg86m<-NRMainC_DDr[NRMainC_DDr$Experiment=='C86' & NRMainC_DDr$Sex=='M','Tau_DD']
Dvf<-NRMainC_DDr[NRMainC_DDr$Experiment=='C35' & NRMainC_DDr$Sex=='F','Tau_DD']
Dg52f<-NRMainC_DDr[NRMainC_DDr$Experiment=='C52' & NRMainC_DDr$Sex=='F','Tau_DD']
Dlf<-NRMainC_DDr[NRMainC_DDr$Experiment=='C71' & NRMainC_DDr$Sex=='F','Tau_DD']
Dof<-NRMainC_DDr[NRMainC_DDr$Experiment=='C84' & NRMainC_DDr$Sex=='F','Tau_DD']
Dg86f<-NRMainC_DDr[NRMainC_DDr$Experiment=='C86' & NRMainC_DDr$Sex=='F','Tau_DD']

## subsets for each experiment and sex in LL ##
Lvm<-NRMainC_LLr[NRMainC_LLr$Experiment=='C35' & NRMainC_LLr$Sex=='M','Tau_LL']
Lg52m<-NRMainC_LLr[NRMainC_LLr$Experiment=='C52' & NRMainC_LLr$Sex=='M','Tau_LL']
Llm<-NRMainC_LLr[NRMainC_LLr$Experiment=='C71' & NRMainC_LLr$Sex=='M','Tau_LL']
Lom<-NRMainC_LLr[NRMainC_LLr$Experiment=='C84' & NRMainC_LLr$Sex=='M','Tau_LL']
Lg86m<-NRMainC_LLr[NRMainC_LLr$Experiment=='C86' & NRMainC_LLr$Sex=='M','Tau_LL']
Lvf<-NRMainC_LLr[NRMainC_LLr$Experiment=='C35' & NRMainC_LLr$Sex=='F','Tau_LL']
Lg52f<-NRMainC_LLr[NRMainC_LLr$Experiment=='C52' & NRMainC_LLr$Sex=='F','Tau_LL']
Llf<-NRMainC_LLr[NRMainC_LLr$Experiment=='C71' & NRMainC_LLr$Sex=='F','Tau_LL']
Lof<-NRMainC_LLr[NRMainC_LLr$Experiment=='C84' & NRMainC_LLr$Sex=='F','Tau_LL']
Lg86f<-NRMainC_LLr[NRMainC_LLr$Experiment=='C86' & NRMainC_LLr$Sex=='F','Tau_LL']

## Combining C51 and C86 together ##
Dg5286m<-NRMainC_DDr[NRMainC_DDr$Species=='giraulti'& NRMainC_DDr$Sex=='M','Tau_DD']
Dg5286f<-NRMainC_DDr[NRMainC_DDr$Species=='giraulti'& NRMainC_DDr$Sex=='F','Tau_DD']
Lg5286m<-NRMainC_LLr[NRMainC_LLr$Species=='giraulti'& NRMainC_LLr$Sex=='M','Tau_LL']
Lg5286f<-NRMainC_LLr[NRMainC_LLr$Species=='giraulti'& NRMainC_LLr$Sex=='F','Tau_LL']


##############################################################
## 1. assessing whether possible to mix C52 and C86 (giraulti) : test effect of Experiment##
##############################################################

C5286<-subset(NRR,NRR$Experiment=='C52'|NRR$Experiment=='C86')   ## subset experiments C52 or C86 ##
C5286g<-subset(C5286,C5286$Species=='giraulti')        ## subset giraulti ##
C5286g$Experiment=with(C5286g,droplevels(Experiment))
gir=C5286g[,c('ID','Light','Sex','Experiment','Tau')]
write.table(gir,"clipboard",sep='\t',col.names=NA)
girD=C5286g[C5286g$Light=='D',] ## subset DD ##
girL=C5286g[C5286g$Light=='L',] ## subset LL ##
  
## Let's plot the data ##
par(mfrow=c(2,2))
boxplot(girD$Tau~girD$Experiment+girD$Sex)
boxplot(girL$Tau~girL$Experiment+girL$Sex)
interaction.plot(girD$Experiment,girD$Sex,girD$Tau)
interaction.plot(girL$Experiment,girL$Sex,girL$Tau)
## summary of data ##
tapply(C5286g$Tau,list(C5286g$Experiment,C5286g$Sex,C5286g$Light),mean)
tapply(C5286g$Tau,list(C5286g$Experiment,C5286g$Sex,C5286g$Light),length)
tapply(C5286g$Tau,list(C5286g$Experiment,C5286g$Sex,C5286g$Light),var)

###### TESTS #########
## 1. WITH ANOVA  ##
aovgir<-aov(C5286g$Tau~C5286g$Experiment*C5286g$Sex*C5286g$Light)
summary(aovgir)
aovgir.0<-aov(C5286g$Tau~(C5286g$Experiment+C5286g$Sex+C5286g$Light)^2)  ## model simplification: no 3 terms interaction #
anova(aovgir,aovgir.0)    # > no difference! > drop 3 terms interaction #
aovgir.01<-update(aovgir.0,~.-C5286g$Sex:C5286g$Light)  # model simplification: drop 1st 2-way interaction #
anova(aovgir.0,aovgir.01)   # > no difference! > drop 1st 2-way interaction #
aovgir.02<-update(aovgir.01,~.-C5286g$Experiment:C5286g$Light)    # model simplification: drop 2nd 2-way interaction #
anova(aovgir.01,aovgir.02)  # > no difference! > drop 2nd 2-way interaction #
aovgir.03<-update(aovgir.02,~.-C5286g$Experiment:C5286g$Sex)  # model simpl. drop 3rd 2-way interaction #
anova(aovgir.02,aovgir.03) # > no difference! > drop 3rd 2-way interaction #
aovgir.1<-aov(C5286g$Tau~C5286g$Experiment+C5286g$Sex+C5286g$Light)    ## model simplification: drop all interactions ##
summary(aovgir.1)
anova(aovgir,aovgir.1) # ok
anova3(C5286g$Tau,C5286g$Experiment,C5286g$Sex,C5286g$Light)    ## Tanja test says Experiment is not significant > try leaving out fixed factors #
aovgir.11<-aov(C5286g$Tau~C5286g$Experiment+C5286g$Sex)    ## model simplification: drop Light ##
anova(aovgir.1,aovgir.11)   # > difference! > keep Light
aovgir.12<-aov(C5286g$Tau~C5286g$Experiment+C5286g$Light)    ## model simplification: drop Sex ##
anova(aovgir.1,aovgir.12)   # > difference! > keep Sex
aovgir.13<-aov(C5286g$Tau~C5286g$Light+C5286g$Sex)    ## model simplification: drop Experiment ##
anova(aovgir.1,aovgir.13)   # > no difference! > drop Experiment
# aovgir.13 best model
summary(aovgir.13)
plot(aovgir.13) # diagnostic plots...  more or less ok, but... ##
shapiro.test(resid(aovgir.13)) ## :-(((
bartlett.test(Tau~Experiment*Sex*Light,data=C5286g)   ## well done???

# diagnostic graphs...
plot(aovgir.13,Tau~fitted(.))   # does not work...
qqnorm(aovgir.13,~resid(.)|Light)
qqnorm(aovgir.13,~resid(.)|Sex)
qqnorm(aovgir.13,~resid(.)|ID)

### Split data into L and D + with Tanja's functions?
aovgirD<-aov(girD$Tau~girD$Experiment*girD$Sex)
summary(aovgirD)
shapiro.test(residuals(aovgirD))
anova2(girD$Tau,girD$Experiment,girD$Sex) ## > no difference in experiments
      
aovgirL<-aov(girL$Tau~girL$Experiment*girL$Sex)
summary(aovgirL)
shapiro.test(residuals(aovgirL))
anova2(girL$Tau,girL$Experiment,girL$Sex) ## > no difference in experiments

# with drop function: a seguire... #

# ANOVA with errors... ##
aovgir.2<-aov(C5286g$Tau~C5286g$Experiment*C5286g$Sex*C5286g$Light+Error(C5286g$Light/C5286g$ID))
summary(aovgir.2)
aovgir.3<-aov(C5286g$Tau~C5286g$Experiment*C5286g$Sex*C5286g$Light+Error(C5286g$Experiment/C5286g$Light))    # model simplification #
summary(aovgir.3)
aovgir.4<-aov(C5286g$Tau~C5286g$Experiment*C5286g$Sex*C5286g$Light+Error(C5286g$Sex*C5286g$Light))    # model simplification #
summary(aovgir.4)
anova(aovgir.2,aovgir.3)
## I do not understand why I get different results... Maybe one factor would be better, as so is example in R book (Rats)...##

## 2. GLM ##
library(nlme)
aoverr.1<-lme(Tau~Experiment*Sex*Light,random=~1|ID/Light,data=C5286g,method="ML")
summary(aoverr.1)
aoverr.2<-lme(Tau~Sex*Light,random=~1|ID/Light,data=C5286g,method="ML")      # simplification: drop Experiment #
summary(aoverr.2)
anova(aoverr.1,aoverr.2)
aoverr.3<-lme(Tau~Sex+Light,random=~1|ID/Light,data=C5286g,method="ML")       # simplification: drop interaction #
summary(aoverr.3)
anova(aoverr.3,aoverr.2)
aoverr.4<-lme(Tau~Light,random=~1|ID/Light,data=C5286g,method="ML")           # simplification: drop sex #
aoverr.5<-lme(Tau~Sex,random=~1|ID/Light,data=C5286g,method='ML')             # simplification: drop light #        
anova(aoverr.4,aoverr.3)
anova(aoverr.5,aoverr.3)
aoverr.6<-lme(Tau~Sex+Light,random=~1|ID,data=C5286g,method="ML")
anova(aoverr.3,aoverr.6)
aoverr.7<-update(aoverr.3,random=~1|Light)
anova(aoverr.3,aoverr.7)
## conclusion: aoverr.3 best model > Experiment has no effect BUT ID/Light also has no effect!!! ###

# diagnostic graphs...
plot(aoverr.3)
plot(aoverr.3,Tau~fitted(.))
qqnorm(aoverr.3,~resid(.)|Light)
qqnorm(aoverr.3,~resid(.)|Sex)
qqnorm(aoverr.3,~resid(.)|ID)

# With lmer...
library(lme4)
lmer.1<-lmer(Tau~Sex*Light*Experiment+(1|ID/Light),data=C5286g)  #  Error in function (fr, FL, start, REML, verbose): Number of levels of a grouping factor for the random effects must be less than the number of observations
lmer.1<-lmer(Tau~Sex*Light*Experiment+(1|ID/Light),data=gir)  #  Error in function (fr, FL, start, REML, verbose): Number of levels of a grouping factor for the random effects must be less than the number of observations
lmer.1<-lmer(Tau~Sex*Light*Experiment+(1|Experiment/ID)+(1|Light),data=C5286g)
lmer.1<-lmer(Tau~Sex*Light*Experiment+(1|Experiment/ID),data=C5286g)
summary(lmer.1)
lmer.2<-lmer(Tau~Sex*Light+(1|Experiment/ID),data=C5286g)
summary(lmer.2)
anova(lmer.1,lmer.2)
lmer.3<-lmer(Tau~Sex+Light+(1|Experiment/ID),poisson,data=C5286g)
anova(lmer.2,lmer.3)
summary(lmer.3)
# diagnostic graphs...
plot(lmer.3)   # does not work...
plot(lmer.3,Tau~fitted(.))      # does not work...
qqnorm(lmer.3,~resid(.)|Light)  # does not work...
qqnorm(lmer.3,~resid(.)|Sex)    # does not work...
qqnorm(lmer.3,~resid(.)|ID)     # does not work...

# checking variances...
# dataframe gir: has NA in empty ID levels
# with lme...
lme.1<-lme(Tau~1,random=~1|Light+ID+Sex+Experiment,data=gir,method="ML")    # fail: missing values in object
lme.1<-lme(Tau~1,random=~1|Light+ID+Sex+Experiment,data=C5286g,method="ML")   # Error in getGroups.data.frame(dataMix, groups) : Invalid formula for groups #
lme.1<-lme(Tau~1,random=~1|Light/ID/Sex/Experiment,data=gir,method="ML")    # error: missing values in object
lme.1<-lme(Tau~1,random=~1|Light/ID/Sex/Experiment,data=C5286g,method="ML")    # it works!
summary(lme.1)
# with lmer...
library(lme4)
lme.2<-lmer(Tau~1+(1|Experiment/Sex/ID/Light),data=gir)  # Error in function (fr, FL, start, REML, verbose): Number of levels of a grouping factor for the random effects must be less than the number of observations ##
lme.2<-lmer(Tau~1+(1|Light/ID/Sex/Experiment),data=gir)   #  Error in function (fr, FL, start, REML, verbose): Number of levels of a grouping factor for the random effects must be less than the number of observations ##
lme.2<-lmer(Tau~1+(1|Light/Sex/Experiment),data=gir)    # it works
summary(lme.2)
lme.3<-lmer(Tau~1+(1|Light)+(1|Sex)+(1|Experiment)+(1|ID),data=gir) ##  Warning message: In mer_finalize(ans) : singular convergence (7) ##
lme.3<-lmer(Tau~1+(1|Light)+(1|Sex)+(1|Experiment)+(1|ID),data=C5286g)
summary(lme.3) 
# conclusions:
# lmer better than lme to 'see' variance contribution
# lme and lmer do not like NA values 


#stuff
sds<-c(bla bla bla bla)
vars<-sds^2
100*vars/sum(vars)


tau2<-read.table("X://Data//Desktop//Tau2.txt", header=T) 
cor.test(tau2$Tau[tau2$Light=="L"],tau2$Tau[tau2$Light=="D"])     #correlation test, there's no correlation   :) good

s<-lm(tau2$Tau[tau2$Light=="L"]~tau2$Tau[tau2$Light=="D"]) ##linear regression it does not make any sense
plot(tau2$Tau[tau2$Light=="L"],tau2$Tau[tau2$Light=="D"],ylim=c(20,30))
abline(s)

## my previous trials...

library(nlme)
attach(C5286g)
lmeC5286g.1<-lme(Tau~Experiment*Sex*Light,random=~1|Experiment/Sex/Light,data=C5286g,method="ML")
summary(lmeC5286g.1)
lmeC5286g.2<-lme(Tau~Experiment+Sex+Light,random=~1|Experiment/Sex/Light,data=C5286g,method="ML")
lmeC5286g.3<-lme(Tau~Experiment,random=~1,data=C5286g,method="ML")
anova(lmeC5286g.1,lmeC5286g.2)
summary(lmeC5286g.2)
summary.aov(lmeC5286g.2)
plot(lmeC5286g.2)
plot(lmeC5286g.2,Tau~fitted(.))
qqnorm(lmeC5286g.2,~resid(.)|Species)
qqnorm(lmeC5286g.2,~resid(.)|Sex)
qqnorm(lmeC5286g.2,~resid(.))
                    
#### Conclusion: no difference between C52 and C86 experiments >> data are pooled #####
             


###############################################
### 2. Comparing vitripennis lines  ###
###############################################

## a) Comparing Asym vitripennis lines among experiments ###

asymR<-subset(NR,NR$Line=='AsymC'& NR$ADRU=='R')
asymR$Experiment=with(asymR,droplevels(Experiment))
tapply(asymR$Tau,list(asymR$Sex,asymR$Experiment,asymR$Light),length)     ## to see how many datapoints every group contains ##
boxplot(asymR$Tau~asymR$Experiment*asymR$Light*asymR$Sex,main='Box_AsymR',ylab='Tau',xlab='Asym lines in different experiments')     ## visualize the new subset ##

# let's make subsets, as it is difficult to interpret the whole set...##
### F and M subsets ###
asymRF<-asymR[asymR$Sex=='F',]       # only asymR Females
asymRM<-asymR[asymR$Sex=='M',]       # only asymR Males
# and plot them... #
par(mfrow=c(2,2))
interaction.plot(asymRF$Light,asymRF$Experiment,asymRF$Tau)
interaction.plot(asymRM$Light,asymRM$Experiment,asymRM$Tau)
interaction.plot(asymR$Sex,asymR$Experiment,asymR$Tau)
interaction.plot(asymR$Light,asymR$Experiment,asymR$Tau)

# anova's...#                             
aovasymR<-aov(asymR$Tau~asymR$Experiment*asymR$Light*asymR$Sex)
summary(aovasymR)
shapiro.test(residuals(aovasymR))
anova3(asymR$Tau,asymR$Experiment,asymR$Light,asymR$Sex,nbrep=5000)

## Also according to Tanja, there is no difference among experiments... ###

## model simplification...? ##
          

aovasymRF<-aov(asymRF$Tau~asymRF$Experiment)#Note: no Light as a factor because no data!!!#
summary(aovasymRF)
shapiro.test(residuals(aovasymRF))    
anova1(asymRF$Tau,asymRF$Experiment)
TukeyHSD(aovasymRF)

aovasymRM<-aov(asymRM$Tau~asymRM$Experiment*asymRM$Light)
summary(aovasymRM)
shapiro.test(residuals(aovasymRM))  
anova2(asymRM$Tau,asymRM$Experiment,asymRM$Light)
TukeyHSD(aovasymRM)

## Conclusion: there are no significant differences between Asym data > they cab be pooled ###


### b) Comparing All vitripennis lines ###

vitR<-subset(NR,NR$Species=='vitripennis'&NR$ADRU=='R')
vitR$Experiment=with(vitR,droplevels(Experiment))
vitR$Line=with(vitR,droplevels(Line))
tapply(vitR$Tau,list(vitR$Sex,vitR$Experiment,vitR$Line,vitR$Light),length)     ## to see how many datapoints every group contains ##
boxplot(vitR$Tau~vitR$Line*vitR$Light*vitR$Sex,main='Box_VitripennisR')               ## visualize the new subset ##

### D and L subsets ###
vitRD<-vitR[vitR$Light=='D',]       # only vitR DD
vitRL<-vitR[vitR$Light=='L',]       # only vitR FF

### let's have a look at the data ###
par(mfrow=c(2,2))       ## plot name: Int_vitRDL_Sex&Light##
interaction.plot(vitR$Sex,vitR$Line,vitR$Tau)
interaction.plot(vitR$Light,vitR$Line,vitR$Tau)
interaction.plot(vitRD$Sex,vitRD$Line,vitRD$Tau)
interaction.plot(vitRL$Sex,vitRL$Line,vitRL$Tau)

# ANOVAs #

aovvitR<-aov(vitR$Tau~vitR$Line*vitR$Light*vitR$Sex)
summary(aovvitR)
shapiro.test(residuals(aovvitR))
anova3(vitR$Tau,vitR$Line,vitR$Light,vitR$Sex)

aovvitRD<-aov(vitRD$Tau~vitRD$Line*vitRD$Sex)
summary(aovvitRD)
shapiro.test(residuals(aovvitRD))
anova2(vitRD$Tau,vitRD$Line,vitRD$Sex)
TukeyHSD(aovvitRD)

aovvitRL<-aov(vitRL$Tau~vitRL$Line*vitRL$Sex)
summary(aovvitRL)
shapiro.test(residuals(aovvitRL))
anova2(vitRL$Tau,vitRL$Line,vitRL$Sex)
TukeyHSD(aovvitRL)

#######################################################################
## 3. Comparing main lines (AsymC, IV72, One, i) in main experiments ##
#######################################################################

## Entire set ##
## Note: giraulti in C52 and C86 are pooled and Asym in all experiments are pooled ###

NRMainR<-subset(NRMain, NRMain$ADRU=='R')    # New set for rythmic wasps in main experiments ##
aovall<-aov(NRMainR$Tau~NRMainR$Sex*NRMainR$Species*NRMainR$Light)
summary(aovall)
shapiro.test(residuals(aovall))
anova3(NRMainR$Tau,NRMainR$Sex,NRMainR$Species,NRMainR$Light,nbrep=1000)

## interaction plots ##
par(mfrow=c(2,2))
interaction.plot(NRMainR$Light,NRMainR$Species,NRMainR$Tau)
interaction.plot(NRMainR$Sex,NRMainR$Species,NRMainR$Tau)
interaction.plot(NRMainR$Light,NRMainR$Sex,NRMainR$Tau)

## Subsets L and D ##

# subset rhythmic wasps, DD and LL #
NRMainRD<-subset(NRMainR, NRMainR$Light=='D')
NRMainRL<-subset(NRMainR, NRMainR$Light=='L')

## Boxplots Taus for rhythmic wasps in main experiments ##                                                                                         
boxplot(Tau~Sex*Species,data=NRMainRD,plot=TRUE,notch=TRUE,ylim=c(19,34),xlab='Tau in DD',ylab='Tau [h]',main='Tau in DD for rhythmic wasps \n >Species<') ## Boxplot of Tau_DD for every sex and experiment
x11()          # opens another graphic window
boxplot(Tau~Sex*Species,data=NRMainRL,plot=TRUE,notch=TRUE,ylim=c(19,34),xlab='Tau in LL',ylab='Tau [h]',main='Tau in LL for rhythmic wasps \n >Species<') ## Boxplot of Tau_LL for every sex and experiment

# ANOVA DD #
aovd<-aov(NRMainRD$Tau~NRMainRD$Sex*NRMainRD$Species)
summary(aovd)
shapiro.test(residuals(aovd))
anova2(NRMainRD$Tau,NRMainRD$Sex,NRMainRD$Species,nbrep=5000) ## check with new F distribution
TukeyHSD(aovd)

# ANOVA LL #
aovl<-aov(NRMainRL$Tau~NRMainRL$Sex*NRMainRL$Species)
summary(aovl)
shapiro.test(residuals(aovl))
anova2(NRMainRL$Tau,NRMainRL$Sex,NRMainRL$Species) ## check with new F distribution
TukeyHSD(aovl)

## interaction plots DD ##
par(mfrow=c(2,2))
interaction.plot(NRMainRD$Sex,NRMainRD$Species,NRMainRD$Tau)
interaction.plot(NRMainRD$Species,NRMainRD$Sex,NRMainRD$Tau)

## interaction plots LL ##
interaction.plot(NRMainRL$Sex,NRMainRL$Species,NRMainRL$Tau)
interaction.plot(NRMainRL$Species,NRMainRL$Sex,NRMainRL$Tau)


######################################################
## 4. Comparing within giraulti, longicornis and oneida lines ##
######################################################

### 4.1. GIRAULTI ### 
## giraulti RV2 - the only line used in different experiments (C52 and C86) - are pooled, as no differences were observed

girR<-subset(NR,NR$Species=='giraulti' & NR$ADRU=='R')
girR$Species=with(girR,droplevels(Species))
girR$Line=with(girR,droplevels(Line))
tapply(girR$Tau,list(girR$Sex,girR$Experiment,girR$Line,girR$Light),length)
boxplot(girR$Tau~girR$Line*girR$Sex*girR$Light,ylab='Tau',main='girR \n Line * Sex * Light')

# Light and darkness subsets
girRL<-subset(girR,girR$Light=='L') # subset giraulti rhythmic Light
girRD<-subset(girR,girR$Light=='D') # subset giraulti rhythmic darkness
par(mfrow=c(2,2))
boxplot(girRL$Tau~girRL$Sex*girRL$Line,main='GirRL \n Sex * Line',ylab='Tau')
boxplot(girRD$Tau~girRD$Sex*girRD$Line,main='GirRD \n Sex * Line',ylab='Tau') 
boxplot(girR$Tau~girR$Line*girR$Light,main='GirR \n Line * Light',ylab='Tau')

# ANOVAs #
aovgirRL<-aov(girRL$Tau~girRL$Line*girRL$Sex)
summary(aovgirRL)
shapiro.test(residuals(aovgirRL))
anova2(girRL$Tau,girRL$Line,girRL$Sex,nbrep=1000)
TukeyHSD(aovgirRL)

aovgirRD<-aov(girRD$Tau~girRD$Line*girRD$Sex)
summary(aovgirRD)
shapiro.test(residuals(aovgirRD))
anova2(girRD$Tau,girRD$Line,girRD$Sex,nbrep=5000)
TukeyHSD(aovgirRD)


### 4.2. LONGICORNIS ###

lonR<-subset(NR,NR$Species=='longicornis'& NR$ADRU=='R')
lonR$Line=with(lonR,droplevels(Line))
tapply(lonR$Tau,list(lonR$Sex,lonR$Experiment,lonR$Line,lonR$Light),length)
boxplot(lonR$Tau~lonR$Line*lonR$Sex*lonR$Light,main='lonR \n Line * Sex * Light', ylab='Tau')

lonRD<-subset(lonR,lonR$Light=='D')
lonRL<-subset(lonR,lonR$Light=='L')
par(mfrow=c(2,2))
boxplot(lonRL$Tau~lonRL$Sex*lonRL$Line,main='lonRL \n Sex * Line',ylab='Tau')
boxplot(lonRD$Tau~lonRD$Sex*lonRD$Line,main='lonRD \n Sex * Line',ylab='Tau')
boxplot(lonR$Tau~lonR$Line*lonR$Light,main='lonR \n Line * Light',ylab='Tau')

# ANOVAs #
aovlonR<-aov(lonR$Tau~lonR$Line*lonR$Sex*lonR$Light)
summary(aovlonR)
shapiro.test(residuals(aovlonR))
anova3(lonR$Tau,lonR$Line,lonR$Sex,lonR$Light,nbrep=1000)

aovlonRD<-aov(lonRD$Tau~lonRD$Line*lonRD$Sex)
summary(aovlonRD)
shapiro.test(residuals(aovlonRD))
anova2(lonRD$Tau,lonRD$Line,lonRD$Sex,nbrep=5000)
TukeyHSD(aovlonRD)

aovlonRL<-aov(lonRL$Tau~lonRL$Line*lonRL$Sex)
summary(aovlonRL)
shapiro.test(residuals(aovlonRL))
anova2(lonRL$Tau,lonRL$Line,lonRL$Sex,nbrep=1000)
TukeyHSD(aovlonRL)

### ONEIDA   (C83, C84)

oneR<-subset(NR,NR$Species=='oneida'&NR$ADRU=='R')
tapply(oneR$Tau,list(oneR$Experiment,oneR$Sex,oneR$Light),length)
oneR$Experiment=with(oneR,droplevels(Experiment))
oneR$Line=with(oneR,droplevels(Line))
oneR$Species=with(oneR,droplevels(Species))
oneR$ADRU=with(oneR,droplevels(ADRU))
tapply(oneR$Tau,list(oneR$Experiment,oneR$Sex,oneR$Light),length)
oneRD<-subset(oneR,oneR$Light=='D')
oneRL<-subset(oneR,oneR$Light=='L')

par(mfrow=c(2,2))
boxplot(oneRL$Tau~oneRL$Line*oneRL$Sex,main='oneRL \n Line * Sex',ylab='Tau')
boxplot(oneRD$Tau~oneRD$Line*oneRD$Sex,main='oneRD \n Line * Sex',ylab='Tau')
boxplot(oneR$Tau~oneR$Line*oneR$Light*oneR$Sex,main='oneR \n Line * Light * Sex',ylab='Tau')

# ANOVA 

aovoneR<-aov(oneR$Tau~oneR$Line*oneR$Light*oneR$Sex)
summary(aovoneR)
shapiro.test(residuals(aovoneR))
anova3(oneR$Tau,oneR$Line,oneR$Light,oneR$Sex,nbrep=1000)

aovoneRL<-aov(oneRL$Tau~oneRL$Line*oneRL$Sex)
summary(aovoneRL)
shapiro.test(residuals(aovoneRL))
anova2(oneRL$Tau,oneRL$Line,oneRL$Sex,nbrep=1000)

aovoneRD<-aov(oneRD$Tau~oneRD$Line*oneRD$Sex)
summary(aovoneRD)
shapiro.test(residuals(aovoneRD))
anova2(oneRD$Tau,oneRD$Line,oneRD$Sex,nbrep=1000)
TukeyHSD(aovoneRD)

# both oneida lines together
boxplot(oneR$Tau~oneR$Sex*oneR$Light,main='oneR \n Sex * Light',ylab='Tau',notch=T)
aovoneRn<-aov(oneR$Tau~oneR$Sex*oneR$Light)
summary(aovoneRn)
shapiro.test(residuals(aovoneRn))
anova2(oneR$Tau,oneR$Sex,oneR$Light,nbrep=1000)

#####################################
# 4. Main lines: with pooled lines ##
#####################################

NRmain2<-subset(NR,NR$Line=='AsymC'|NR$Line=='RV2'|NR$Line=='IVFR2'|NR$Species=='oneida')
NRmain2R<-subset(NRmain2,NRmain2$ADRU=='R')
NRmain2RD<-subset(NRmain2R,NRmain2R$Light=='D')
NRmain2RL<-subset(NRmain2R,NRmain2R$Light=='L')
par(mfrow=c(2,2))
boxplot(NRmain2RL$Tau~NRmain2RL$Sex*NRmain2RL$Species,ylim=c(19,34),ylab='Tau',main='Tau in LL for rhythmic wasps \n NRmain2RL',notch=T)
boxplot(NRmain2RD$Tau~NRmain2RD$Sex*NRmain2RD$Species,ylim=c(19,34),ylab='Tau',main='Tau in DD for rhythmic wasps \n NRmain2RD',notch=T)
boxplot(NRmain2R$Tau~NRmain2R$Sex*NRmain2R$Light*NRmain2R$Species,ylab='Tau',main='NRmain2R',notch=T)

# ANOVAs

aovNRmain2RD<-aov(NRmain2RD$Tau~NRmain2RD$Sex*NRmain2RD$Species)
summary(aovNRmain2RD)
shapiro.test(residuals(aovNRmain2RD))
anova2(NRmain2RD$Tau,NRmain2RD$Sex,NRmain2RD$Species,nbrep=1000)
TukeyHSD(aovNRmain2RD)

aovNRmain2RL<-aov(NRmain2RL$Tau~NRmain2RL$Sex*NRmain2RL$Species)
summary(aovNRmain2RL)
shapiro.test(residuals(aovNRmain2RL))
anova2(NRmain2RL$Tau,NRmain2RL$Sex,NRmain2RL$Species,nbrep=1000)
TukeyHSD(aovNRmain2RL)

######################

# analysis

AsymMRD<-subset(NR,NR$Line=='AsymC'&NR$ADRU=='R'&NR$Light=='D'&NR$Sex=='M')
AsymMRL<-subset(NR,NR$Line=='AsymC'&NR$ADRU=='R'&NR$Light=='L'&NR$Sex=='M')
AsymMR<-subset(NR,NR$Line=='AsymC'&NR$ADRU=='R'&NR$Light=='D'&NR$Sex=='F')



#################################
###         CoG          ########
#################################

phi<-subset(NR,NR$Experiment==c('C35','C52') & NR$ADRU!='D' & NR$Light=='L', select=c(Species,Experiment,Line,Sex,Light,ADRU,PhiMax_Deg,CoG))
phiC52F<-subset(phi,phi$Experiment=='C52'&phi$Species=='giraulti'&phi$Sex=='F')
phiC52M<-subset(phi,phi$Experiment=='C52'&phi$Species=='giraulti'&phi$Sex=='M')
phiC35F<-subset(phi,phi$Experiment=='C35'&phi$Species=='vitripennis'&phi$Sex=='F')
phiC35M<-subset(phi,phi$Experiment=='C35'&phi$Species=='vitripennis'&phi$Sex=='M')
phiC71F<-subset(phi,phi$Experiment=='C71'&phi$Species=='longicornis'&phi$Sex=='F')
phiC71M<-subset(phi,phi$Experiment=='C71'&phi$Species=='longicornis'&phi$Sex=='M')
phiC84F<-subset(phi,phi$Experiment=='C84'&phi$Species=='oneida'&phi$Sex=='F')
phiC84M<-subset(phi,phi$Experiment=='C84'&phi$Species=='oneida'&phi$Sex=='M')
set<-c("phiC52F","phiC52M","phiC35F","phiC35M","phiC71F","phiC71M","phiC84F","phiC84M")
set2<-c(phiC52F,phiC52M)
means<-c(mean(phiC52F$PhiMax_Deg),mean(phiC52M$PhiMax_Deg),mean(phiC35F$PhiMax_Deg),mean(phiC35M$PhiMax_Deg),mean(phiC71F$PhiMax_Deg),mean(phiC71M$PhiMax_Deg),mean(phiC84F$PhiMax_Deg),mean(phiC84M$PhiMax_Deg))
stdev<-c(sd(phiC52F$PhiMax_Deg),sd(phiC52M$PhiMax_Deg),sd(phiC35F$PhiMax_Deg),sd(phiC35M$PhiMax_Deg),sd(phiC71F$PhiMax_Deg),sd(phiC71M$PhiMax_Deg),sd(phiC84F$PhiMax_Deg),sd(phiC84M$PhiMax_Deg))
size<-c(length(phiC52F$PhiMax_Deg),length(phiC52M$PhiMax_Deg),length(phiC35F$PhiMax_Deg),length(phiC35M$PhiMax_Deg),length(phiC71F$PhiMax_Deg),length(phiC71M$PhiMax_Deg),length(phiC84F$PhiMax_Deg),length(phiC84M$PhiMax_Deg))
cogtab<-data.frame(set,means,stdev,size)
names(cogtab)<-c("Set","Phi_Max","St_Dev","N")

tapply(phi$PhiMax_Deg,list(phi$Experiment,phi$Sex,phi$Species),length)

write.table(phi,'clipboard-128',sep='\t',col.names=NA)

model<-aov(NRR$Tau~NRR$Light*NRR$Sex*NRR$Species)
model2.aov<-aov(NRR$Tau~(NRR$Light*NRR$Sex*NRR$Species)^2+(NRR$Species/NRR$Experiment),na.action=NULL)  
model2.ln<-lm(NRR$Tau~(NRR$Light*NRR$Sex*NRR$Species)^2+(NRR$Species/NRR$Experiment),na.action=na.fail)


