#This file contains the R code necessary to replicate the analyses performed in this manuscript:
#Pintar, M.R., N.J. Dorn, J.L. Kline, and J.C. Trexler. Hydrology-mediated ecological function of a large wetland threatened by an invasive predator. Science of the Total Environment.

#Load required packages
library(lme4)
library(lmerTest)
library(MuMIn)
library(plyr)


C <- read.csv(file="Pintar-et-al_STOTEN_swamp-eels_data.csv", head=T)

C$TOPPREDS <- C$E.LEPPLA + C$E.MICSAL + C$E.AMICAL #create top predator variable
C$PERIOD2 <- as.factor(C$PERIOD) #Set site as factor
C$PERIOD3 <- mapvalues(C$PERIOD2, from = c("1", "2", "3", "4", "5"), 
                       to = c("4", "5", "1", "2", "3")) #rename original order of sites so 1 = July
C$PERIOD3 <- ordered(C$PERIOD3, levels=c("1", "2", "3", "4", "5"))
C$PERIOD5 <- as.numeric(C$PERIOD3)
C$PERIOD4 <- as.factor(C$PERIOD5)

#create invasion periods (a = 'before', b = 'during', c = 'after')
C$TSLPD2 <- mapvalues(C$YEAR, from = c("1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009",
                                       "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022"), 
                      to = c("a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a",
                             "b", "b", "b", "b", "b", "c", "c", "c", "c", "c", "c", "c", "c"))

#Create length of previous dry season (LDS) variable
C$LASTDRYSEASON <- ifelse(C$DSLDD>365, "0", C$LASTDAYDRY)
C$LASTDRYSEASON <- as.numeric(C$LASTDRYSEASON)

#log-transform hydrologic variables
C$log.LDS <- log(C$LASTDRYSEASON + 1)
C$log.DSD <- log(C$DSLDD + 1) 

#create datasets by region
SRS <- C[!(C$REGION=="TSL" | C$REGION=="WCA" | C$REGION=="PHD"),] #Shark River Slough
TSL <- C[!(C$REGION=="SRS" | C$REGION=="WCA" | C$REGION=="PHD"),] #Taylor Slough
WCA <- C[!(C$REGION=="TSL" | C$REGION=="SRS" | C$REGION=="PHD"),] #Water Conservation Area 3
PHD <- C[!(C$REGION=="TSL" | C$REGION=="SRS" | C$REGION=="WCA"),] #Panhandle


#Create datasets by invasion periods in TSL
TSL1 <- TSL[!(TSL$YEAR>2009),] #TSL prior swamp eel arrival
TSL2 <- TSL[!(TSL$YEAR<2010),] #after swamp eel detection
TSL3 <- TSL[!(TSL$TSLPD2=="b"),] #Drop middle period when swamp eels were spreading in TSL

#Create function to generate standardized coefficients
stdCoef.merMod <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- fixef(object)*sdx/sdy
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))
}



##### Analyses on individual species #####

### Palaemonetes [Palaemon] paludosus ###
#Baseline period
PP1a <- lmer(log(PALPAL + 1) ~ 
               PERIOD4 + log.LDS + DEPTH_AVE_180DAY + log.DSD + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #180-day depth
dredge(PP1a)

PP1b <- lmer(log(PALPAL + 1) ~ 
               PERIOD4 + log.LDS + DEPTH_AVE_30DAY + log.DSD + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #30-day depth
dredge(PP1b)

PP1 <- lmer(log(PALPAL + 1) ~ 
              PERIOD4 + DEPTH_AVE_180DAY + 
              (1|REGSITE/SITEPLOT), data=TSL1, REML=F) #best hydrologic model
anova(PP1)
summary(PP1)
r.squaredGLMM(PP1)
stdCoef.merMod(PP1)

#Generate predicted densities after invasion based on baseline hydrologic model
TSL2$log.PALPAL.p <- predict(PP1, newdata=TSL2, allow.new.levels=TRUE)
TSL2$PALPAL.p <- exp(TSL2$log.PALPAL.p) -1


#Before vs after invasion
PP2a <- lmer(log(PALPAL + 1) ~ 
               PERIOD4 + DEPTH_AVE_180DAY + 
               (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)


PP2b <- lmer(log(PALPAL + 1) ~  TSLPD2 + 
               PERIOD4 + DEPTH_AVE_180DAY + 
               (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

anova(PP2a, PP2b)

#Model of full time period incorporating densities of predators
PP3 <- lmer(log(PALPAL + 1) ~
              PERIOD4 + DEPTH_AVE_180DAY + 
              log(E.MONALB + 1) + log(E.CICURO + 1) + log(TOPPREDS + 1) +
              (1|REGSITE/SITEPLOT), data=TSL, REML=F)
step(PP3, reduce.random=F, alpha.fixed=.2)

PP3 <- lmer(log(PALPAL + 1) ~
              PERIOD4 + DEPTH_AVE_180DAY + 
              log(E.MONALB + 1) + log(E.CICURO + 1) + 
              (1|REGSITE/SITEPLOT), data=TSL, REML=F)
anova(PP3)
summary(PP3)
r.squaredGLMM(PP3)
stdCoef.merMod(PP3)



### Procambarus alleni ###
#Baseline period
PA1a <- lmer(log(PROALL + 1) ~ 
               PERIOD4 + log.LDS + DEPTH_AVE_180DAY + log.DSD + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #180-day depth
dredge(PA1a)

PA1b <- lmer(log(PROALL + 1) ~ 
               PERIOD4 + log.LDS + DEPTH_AVE_30DAY + log.DSD + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #30-day depth
dredge(PA1b)

PA1 <- lmer(log(PROALL + 1) ~ 
              PERIOD4 + log.LDS + log.DSD + DEPTH_AVE_180DAY + PERIOD4:log.LDS + 
              (1|REGSITE/SITEPLOT), data=TSL1, REML=F) #best hydrologic model
anova(PA1)
summary(PA1)
r.squaredGLMM(PA1)
stdCoef.merMod(PA1)

#Generate predicted densities after invasion based on baseline hydrologic model
TSL2$log.PROALL.p <- predict(PA1, newdata=TSL2, allow.new.levels=TRUE)
TSL2$PROALL.p <- exp(TSL2$log.PROALL.p) -1

#Before vs after invasion
PA2a <- lmer(log(PROALL + 1) ~ 
               PERIOD4 + log.LDS + log.DSD + DEPTH_AVE_180DAY + PERIOD4:log.LDS + 
               (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

PA2b <- lmer(log(PROALL + 1) ~ TSLPD2 +
               PERIOD4 + log.LDS + log.DSD + DEPTH_AVE_180DAY + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

anova(PA2a, PA2b)

#Model of full time period incorporating densities of predators
PA1 <- lmer(log(PROALL + 1) ~ 
              PERIOD4 + log.LDS + log.DSD + DEPTH_AVE_180DAY + PERIOD4:log.LDS +
              log(E.MONALB + 1) + log(E.CICURO + 1) + log(TOPPREDS + 1) + 
              (1|REGSITE/SITEPLOT), data=TSL, REML=F)
step(PA1, reduce.random=F, alpha.fixed=.2)

PA1 <- lmer(log(PROALL + 1) ~ 
              PERIOD4 + log.LDS + log.DSD + DEPTH_AVE_180DAY + PERIOD4:log.LDS + 
              log(E.MONALB + 1) +  
              (1|REGSITE/SITEPLOT), data=TSL, REML=F)

anova(PA1)
summary(PA1)
r.squaredGLMM(PA1)
stdCoef.merMod(PA1)



### Procambarus fallax ###
#Baseline period
PF1a <- lmer(log(PROFAL + 1) ~ 
               PERIOD4 + log.LDS + DEPTH_AVE_180DAY + log.DSD + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #180-day depth
dredge(PF1a)

PF1b <- lmer(log(PROFAL + 1) ~ 
               PERIOD4 + log.LDS + DEPTH_AVE_30DAY + log.DSD + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #30-day depth
dredge(PF1b)

PF1 <- lmer(log(PROFAL + 1) ~ 
              PERIOD4 + log.DSD + DEPTH_AVE_30DAY + 
              (1|REGSITE/SITEPLOT), data=TSL1, REML=F) #best hydrologic model
anova(PF1)
summary(PF1)
r.squaredGLMM(PF1)
stdCoef.merMod(PF1)

#Generate predicted densities after invasion based on baseline hydrologic model
TSL2$log.PROFAL.p <- predict(PF1, newdata=TSL2, allow.new.levels=TRUE)
TSL2$PROFAL.p <- exp(TSL2$log.PROFAL.p) -1


#Before vs after invasion
PF2a <- lmer(log(PROFAL + 1) ~
               PERIOD4 + log.DSD + DEPTH_AVE_30DAY + 
               (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

PF2b <- lmer(log(PROFAL + 1) ~ TSLPD2 + 
               PERIOD4 + log.DSD + DEPTH_AVE_30DAY + 
               (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

anova(PF2a, PF2b)

#Model of full time period incorporating densities of predators
PF1 <- lmer(log(PROFAL + 1) ~
              PERIOD4 + log.DSD + DEPTH_AVE_30DAY + 
              log(E.MONALB + 1) + log(E.CICURO + 1) + log(TOPPREDS + 1) + 
              (1|REGSITE/SITEPLOT), data=TSL, REML=F)
step(PF1, reduce.random=F, alpha.fixed=.2)

PF1 <- lmer(log(PROFAL + 1) ~
              PERIOD4 + log.DSD + DEPTH_AVE_30DAY + 
              log(E.MONALB + 1) + log(TOPPREDS + 1) + 
              (1|REGSITE/SITEPLOT), data=TSL, REML=F)
anova(PF1)
summary(PF1)
r.squaredGLMM(PF1)
stdCoef.merMod(PF1)


### Fundulus chrysotus ###
#Baseline period
Funchr1a <- lmer(log(FUNCHR + 1) ~ 
                   PERIOD4 + log.LDS + DEPTH_AVE_180DAY + log.DSD + PERIOD4:log.LDS +
                   (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #180-day depth
dredge(Funchr1a)

Funchr1b <- lmer(log(FUNCHR + 1) ~ 
                   PERIOD4 + log.LDS + DEPTH_AVE_30DAY + log.DSD + PERIOD4:log.LDS +
                   (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #30-day depth
dredge(Funchr1b)


Funchr1 <- lmer(log(FUNCHR + 1) ~ 
                  PERIOD4 + log.LDS + log.DSD + DEPTH_AVE_30DAY + 
                  (1|REGSITE/SITEPLOT), data=TSL1, REML=F) #best hydrologic model
anova(Funchr1)
summary(Funchr1)
r.squaredGLMM(Funchr1)
stdCoef.merMod(Funchr1)

#Generate predicted densities after invasion based on baseline hydrologic model
TSL2$log.FUNCHR.p <- predict(Funchr1, newdata=TSL2, allow.new.levels=TRUE)
TSL2$FUNCHR.p <- exp(TSL2$log.FUNCHR.p) -1


#Before vs after invasion
Funchr2a <- lmer(log(FUNCHR + 0.5) ~ 
                   PERIOD4 + log.LDS + log.DSD + DEPTH_AVE_30DAY +
                   (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

Funchr2b <- lmer(log(FUNCHR + 0.5) ~+ TSLPD2 + 
                   PERIOD4 + log.LDS + log.DSD + DEPTH_AVE_30DAY +
                   (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

anova(Funchr2a, Funchr2b)

#Model of full time period incorporating densities of predators
Funchr1 <- lmer(log(FUNCHR + 1) ~ 
                  PERIOD4 + log.LDS + log.DSD + DEPTH_AVE_30DAY +
                  log(E.MONALB + 1) + log(E.CICURO + 1) + log(TOPPREDS + 1) + 
                  (1|REGSITE/SITEPLOT), data=TSL, REML=F)
step(Funchr1, reduce.random=F, alpha.fixed=.2)

Funchr1 <- lmer(log(FUNCHR + 1) ~ 
                  PERIOD4 + log.LDS + log.DSD + DEPTH_AVE_30DAY +
                  log(E.MONALB + 1) + log(TOPPREDS + 1) + 
                  (1|REGSITE/SITEPLOT), data=TSL, REML=F)
anova(Funchr1)
summary(Funchr1)
r.squaredGLMM(Funchr1)
stdCoef.merMod(Funchr1)



### Fundulus confluentus ###
#Baseline period
Funcon1a <- lmer(log(FUNCON + 1) ~ 
                   PERIOD4 + log.LDS + DEPTH_AVE_180DAY + log.DSD + PERIOD4:log.LDS +
                   (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #180-day depth
dredge(Funcon1a)

Funcon1b <- lmer(log(FUNCON + 1) ~ 
                   PERIOD4 + log.LDS + DEPTH_AVE_30DAY + log.DSD + PERIOD4:log.LDS +
                   (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #30-day depth
dredge(Funcon1b)

Funcon1 <- lmer(log(FUNCON + 1) ~  
                  log(DSLDD + 1) + DEPTH_AVE_30DAY + 
                  (1|REGSITE/SITEPLOT), data=TSL1, REML=F) #best hydrologic model
anova(Funcon1)
summary(Funcon1)
r.squaredGLMM(Funcon1)
stdCoef.merMod(Funcon1)

#Generate predicted densities after invasion based on baseline hydrologic model
TSL2$log.FUNCON.p <- predict(Funcon1, newdata=TSL2, allow.new.levels=TRUE)
TSL2$FUNCON.p <- exp(TSL2$log.FUNCON.p) -1

#Before vs after invasion
Funcon2a <- lmer(log(FUNCON + 1) ~
                   log(DSLDD + 1) + DEPTH_AVE_30DAY + 
                   (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

Funcon2b <- lmer(log(FUNCON + 1) ~  TSLPD2 + 
                   log(DSLDD + 1) + DEPTH_AVE_30DAY + 
                   (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

anova(Funcon2a, Funcon2b)

#Model of full time period incorporating densities of predators
Funcon1 <- lmer(log(FUNCON + 1) ~ 
                  log(DSLDD + 1) + DEPTH_AVE_30DAY + 
                  log(E.MONALB + 1) + log(E.CICURO + 1) + log(TOPPREDS + 1) + 
                  (1|REGSITE/SITEPLOT), data=TSL, REML=F)
step(Funcon1, reduce.random=F, alpha.fixed=.2)

Funcon1 <- lmer(log(FUNCON + 1) ~ 
                  log(DSLDD + 1) + DEPTH_AVE_30DAY + 
                  log(E.MONALB + 1) + log(E.CICURO + 1) + log(TOPPREDS + 1) + 
                  (1|REGSITE/SITEPLOT), data=TSL, REML=F)
anova(Funcon1)
summary(Funcon1)
r.squaredGLMM(Funcon1)
stdCoef.merMod(Funcon1)



### Gambusia holbrooki ###
#Baseline period
GH1a <- lmer(log(GAMHOL + 1) ~ 
               PERIOD4 + log.LDS + DEPTH_AVE_180DAY + log.DSD + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #180-day depth
dredge(GH1a)

GH1b <- lmer(log(GAMHOL + 1) ~ 
               PERIOD4 + log.LDS + DEPTH_AVE_30DAY + log.DSD + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #30-day depth
dredge(GH1b)

GH1 <- lmer(log(GAMHOL + 1) ~ 
              log.LDS + log.DSD + DEPTH_AVE_180DAY +  
              (1|REGSITE/SITEPLOT), data=TSL1, REML=F) #best hydrologic model
anova(GH1)
summary(GH1)
r.squaredGLMM(GH1)
stdCoef.merMod(GH1)

#Generate predicted densities after invasion based on baseline hydrologic model
TSL2$log.GAMHOL.p <- predict(GH1, newdata=TSL2, allow.new.levels=TRUE)
TSL2$GAMHOL.p <- exp(TSL2$log.GAMHOL.p) -1


#Before vs after invasion
GH2a <- lmer(log(GAMHOL + 1) ~ 
               log.LDS + log.DSD + DEPTH_AVE_180DAY +  
               (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

GH2b <- lmer(log(GAMHOL + 1) ~ TSLPD2 + 
               log.LDS + log.DSD + DEPTH_AVE_180DAY +  
               (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

anova(GH2a, GH2b)

#Model of full time period incorporating densities of predators
GH2b <- lmer(log(GAMHOL + 1) ~ 
               log.LDS + log.DSD + DEPTH_AVE_180DAY +  
               log(E.MONALB + 1) + log(E.CICURO + 1) + log(TOPPREDS + 1) +
               (1|REGSITE/SITEPLOT), data=TSL, REML=F, na.action=na.exclude)
step(GH2b, reduce.random=F, alpha.fixed=.2)

GH2b <- lmer(log(GAMHOL + 1) ~ 
               log.LDS + log.DSD + DEPTH_AVE_180DAY +  
               log(E.MONALB + 1) + log(E.CICURO + 1) + log(TOPPREDS + 1) + 
               (1|REGSITE/SITEPLOT), data=TSL, REML=F, na.action=na.exclude)
anova(GH2b)
summary(GH2b)
r.squaredGLMM(GH2b)
stdCoef.merMod(GH2b)


### Heterandria formosa ###
#Baseline period
HF1a <- lmer(log(HETFOR + 1) ~ 
               PERIOD4 + log.LDS + DEPTH_AVE_180DAY + log.DSD + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #180-day depth
dredge(HF1a)

HF1b <- lmer(log(HETFOR + 1) ~ 
               PERIOD4 + log.LDS + DEPTH_AVE_30DAY + log.DSD + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #30-day depth
dredge(HF1b)

HF1 <- lmer(log(HETFOR + 1) ~ 
              PERIOD4 + log.LDS + log.DSD + PERIOD4:log.LDS + 
              (1|REGSITE/SITEPLOT), data=TSL1, REML=F) #best hydrologic model
anova(HF1)
summary(HF1)
r.squaredGLMM(HF1)
stdCoef.merMod(HF1)

#Generate predicted densities after invasion based on baseline hydrologic model
TSL2$log.HETFOR.p <- predict(HF1, newdata=TSL2, allow.new.levels=TRUE)
TSL2$HETFOR.p <- exp(TSL2$log.HETFOR.p) -1


#Before vs after invasion
HF2a <- lmer(log(HETFOR + 1) ~ 
               PERIOD4 + log.LDS + log.DSD + PERIOD4:log.LDS + 
               (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

HF2b <- lmer(log(HETFOR + 1) ~ TSLPD2 + 
               PERIOD4 + log.LDS + log.DSD + PERIOD4:log.LDS + 
               (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

anova(HF2a, HF2b)

#Model of full time period incorporating densities of predators
HF1 <- lmer(log(HETFOR + 1) ~ 
              PERIOD4 + log.LDS + log.DSD + PERIOD4:log.LDS + 
              log(E.MONALB + 1) + log(E.CICURO + 1) + log(TOPPREDS + 1) +
              (1|REGSITE/SITEPLOT), data=TSL, REML=F)
step(HF1, reduce.random=F, alpha.fixed=.2)

HF1 <- lmer(log(HETFOR + 1) ~ 
              PERIOD4 + log.LDS + log.DSD + PERIOD4:log.LDS + 
              log(E.MONALB + 1) + log(TOPPREDS + 1) +
              (1|REGSITE/SITEPLOT), data=TSL, REML=F)
anova(HF1)
summary(HF1)
r.squaredGLMM(HF1)
stdCoef.merMod(HF1)



### Jordanella floridae ###
#Baseline period
JF1a <- lmer(log(JORFLO + 1) ~ 
               PERIOD4 + log.LDS + DEPTH_AVE_180DAY + log.DSD + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #180-day depth
dredge(JF1a)

JF1b <- lmer(log(JORFLO + 1) ~ 
               PERIOD4 + log.LDS + DEPTH_AVE_30DAY + log.DSD + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #30-day depth
dredge(JF1b)

JF1 <- lmer(log(JORFLO + 1) ~ 
              PERIOD4 + log.DSD + 
              (1|REGSITE/SITEPLOT), data=TSL1, REML=F) #best hydrologic model
anova(JF1)
summary(JF1)
r.squaredGLMM(JF1)
stdCoef.merMod(JF1)


#Generate predicted densities after invasion based on baseline hydrologic model
TSL2$log.JORFLO.p <- predict(JF1, newdata=TSL2, allow.new.levels=TRUE)
TSL2$JORFLO.p <- exp(TSL2$log.JORFLO.p) - 1


#Before vs after invasion
JF2a <- lmer(log(JORFLO + 1) ~ 
               PERIOD4 + log.DSD + 
               (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

JF2b <- lmer(log(JORFLO + 1) ~  TSLPD2 + 
               PERIOD4 + log.DSD + 
               (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

anova(JF2a, JF2b)

#Model of full time period incorporating densities of predators
JF1 <- lmer(log(JORFLO + 1) ~ 
              PERIOD4 + log.LDS + 
              log(E.MONALB + 1) + log(E.CICURO + 1) + log(TOPPREDS + 1) + 
              (1|REGSITE/SITEPLOT), data=TSL, REML=F)
step(JF1, reduce.random=F, alpha.fixed=.2)

JF1 <- lmer(log(JORFLO + 1) ~ 
              PERIOD4 + log.LDS + 
              log(E.MONALB + 1) + log(E.CICURO + 1) + log(TOPPREDS + 1) + 
              (1|REGSITE/SITEPLOT), data=TSL, REML=F)
anova(JF1)
summary(JF1)
r.squaredGLMM(JF1)
stdCoef.merMod(JF1)



### Lucania goodei ###
#Baseline period
LG1a <- lmer(log(LUCGOO + 1) ~ 
               PERIOD4 + log.LDS + DEPTH_AVE_180DAY + log.DSD + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #180-day depth
dredge(LG1a)

LG1b <- lmer(log(LUCGOO + 1) ~ 
               PERIOD4 + log.LDS + DEPTH_AVE_30DAY + log.DSD + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL1, REML=F, na.action="na.fail") #30-day depth
dredge(LG1b)

LG1 <- lmer(log(LUCGOO + 1) ~ 
              PERIOD4 + log.LDS + log.DSD + DEPTH_AVE_180DAY + PERIOD4:log.LDS + 
              (1|REGSITE/SITEPLOT), data=TSL1, REML=F) #best hydrologic model
anova(LG1)
summary(LG1)
r.squaredGLMM(LG1)
stdCoef.merMod(LG1)

#Generate predicted densities after invasion based on baseline hydrologic model
TSL2$log.LUCGOO.p <- predict(LG1, newdata=TSL2, allow.new.levels=TRUE)
TSL2$LUCGOO.p <- exp(TSL2$log.LUCGOO.p) -1


#Before vs after invasion
LG2a <- lmer(log(LUCGOO + 1) ~ 
               PERIOD4 + log.LDS + log.DSD + DEPTH_AVE_180DAY + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

LG2b <- lmer(log(LUCGOO + 1) ~  TSLPD2 + 
               PERIOD4 + log.LDS + log.DSD + DEPTH_AVE_180DAY + PERIOD4:log.LDS +
               (1|REGSITE/SITEPLOT), data=TSL3, REML=F, na.action=na.exclude)

anova(LG2a, LG2b)

#Model of full time period incorporating densities of predators
LG2 <- lmer(log(LUCGOO + 0.5) ~ 
              PERIOD4 + log.LDS + log.DSD + DEPTH_AVE_180DAY + PERIOD4:log.LDS + 
              log(E.MONALB + 1) + log(E.CICURO + 1) + log(TOPPREDS + 1) + 
              (1|REGSITE/SITEPLOT), data=TSL, REML=F, na.action=na.exclude)
step(LG2, reduce.random=F, alpha.fixed=.2)

LG2 <- lmer(log(LUCGOO + 0.5) ~ 
              PERIOD4 + log.LDS + log.DSD + DEPTH_AVE_180DAY + PERIOD4:log.LDS + 
              log(E.MONALB + 1) + log(TOPPREDS + 1) + 
              (1|REGSITE/SITEPLOT), data=TSL, REML=F, na.action=na.exclude)
anova(LG2)
summary(LG2)
stdCoef.merMod(LG2)
r.squaredGLMM(LG2)

