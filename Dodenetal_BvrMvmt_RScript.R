##########################################################
# Script for Doden et al. Movement patterns of resident
# and translocated beavers at multiple spationtemporal
# scales in desert rivers
# Analyses
# Last updated January 10, 2022
##########################################################


#########################################################
####Maximum Distance Detected Analysis #### 
#########################################################

library(tidyverse)
library(ggplot2)
library(nlme)
library(MuMIn)

rm(list=ls())
graphics.off()
options(stringsAsFactors = FALSE)

wd <- "C:/Users/Documents"
setwd(wd)

locMDD <-read.csv("./1_MaxDistDet.csv", #name corresponds to name of tab in Dodenetal_BvrMvmtatMultScales Excel file
                header = T, na.strings = c("", "N/A", "NA"))

# check structure
str(locMDD)

locMDD$BvrID <- as.factor(locMDD$BvrID)

locMDD$Grp_Age <- as.factor(locMDD$Grp_Age)

## transform distance detected
locMDD$logMaxDistDet_km <- log(locMDD$MaxDistDet_km)

## restructure variables for stepwise selection
locMDD$RA <- as.numeric(locMDD$Grp_Age == "RA")
locMDD$RS <- as.numeric(locMDD$Grp_Age == "RS")
locMDD$TS <- as.numeric(locMDD$Grp_Age == "TS")
locMDD$TA <- as.numeric(locMDD$Grp_Age == "TA")

str(locMDD)

## Create box-and-whisker plot of raw data
plot_MDD <- ggplot(locMDD) + 
  geom_boxplot(aes(x=Grp_Age, y=MaxDistDet_km, fill=Grp_Age)) +
  theme_classic() +
  scale_fill_manual(values = c("RA"="chocolate4",
                               "RS"="gray50",
                               "TA"="purple4",
                               "TS"="darkolivegreen")) +
  labs(x="State Category",y="Maximum distance detected (km)") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))  

plot_MDD


## Construct model
summary(m_finMDD <- lm(logMaxDistDet_km ~ 0 + RA + RS + TS + TA, data = locMDD))

# check normality
par(mfrow= c(2,2))
plot(m_finMDD)



#########################################################
###Displacement From Release Analysis #### 
#########################################################

library(tidyverse)
library(ggplot2)
library(nlme)
library(MuMIn)

rm(list=ls())
graphics.off()
options(stringsAsFactors = FALSE)

wd <- "C:/Users/Documents"
setwd(wd)

# load data
loc2 <- read.csv("./2_DisplRel.csv",
                header = T, na.strings = c("", "N/A", "NA"))

# check structure
str(loc2)

loc2$BvrID <- as.factor(loc2$BvrID)
loc2$Grp_Age <- as.factor(loc2$Grp_Age)

## transform and add 1 to all step lengths (so none are 0)
loc2$DistRelNo0 <- loc2$DistRel_km + 1
loc2$LogDistRel <- log(loc2$DistRelNo0)

## transform days since release
loc2$LogDaysRel <- log(loc2$DaysSinceRel)

## remove all kits
loc3 <- loc2 %>% filter(!Grp_Age == "RK" & !Grp_Age == "TK") %>% droplevels()
str(loc3$Grp_Age)

## remove all locations > 182 days post-release (6 months)
# b/c very few individuals were detected after this point
loc4 <- loc3 %>% filter(DaysSinceRel < 183)

## remove all individuals w/ only 1 or 2 samples
locDR <- loc4 %>% group_by(BvrID) %>% filter(n() > 2)

## restructure vars for stepwise selection
locDR$RA <- as.numeric(locDR$Grp_Age == "RA")
locDR$RS <- as.numeric(locDR$Grp_Age == "RS")
locDR$TS <- as.numeric(locDR$Grp_Age == "TS")
locDR$TA <- as.numeric(locDR$Grp_Age == "TA")

## Construct model
summary(m_finDR <- lme(LogDistRel ~ 0 + RA + RS + TS + TA + RA:LogDaysRel + 
                        RS:LogDaysRel + TS:LogDaysRel + TA:LogDaysRel, 
                      random = ~ 1 + LogDaysRel | BvrID, data = locDR))


## test whether adding BvrID as a Random Effect improves model fit
require(lmerTest)
require(lme4)
summary(m_fin_DRtest <- lmer(LogDistRel ~ 0 + RA + RS + TS + TA + RA:LogDaysRel + 
  RS:LogDaysRel + TS:LogDaysRel + TA:LogDaysRel + 
 (1 + LogDaysRel | BvrID), data = locDR))

ranova(m_fin_DRtest) # p-value is significant

# look at pseudo R2 values
pseudoR <- r.squaredGLMM(m_finDR)
pseudoR

# check normality
qqnorm(residuals(m_finDR))
hist(residuals(m_finDR)) 

#### Predict LogDistRel x LogDaysRel #####

# figure out values needed to build dataframe
max(locDR$LogDaysRel)
min(locDR$LogDaysRel)
LogDaysRel <- rep(seq(min(locDR$LogDaysRel), max(locDR$LogDaysRel), length.out = 666), 4)

# here is our model:
# #m_finDR <- lme(LogDistRel ~ 0 + RA + RS + TS + TA + RA:LogDaysRel + 
#                          RS:LogDaysRel + TS:LogDaysRel + TA:LogDaysRel, 
#                        random = ~ 1 + LogDaysRel | BvrID, data = locDR))

## build the dataframe
pred_DR <- data.frame(LogDaysRel = LogDaysRel, 
                       RA = rep(c(1, 0, 0, 0), each = 666),
                       RS = rep(c(0, 1, 0, 0), each = 666),
                       TA = rep(c(0, 0, 1, 0), each = 666),
                       TS = rep(c(0, 0, 0, 1), each = 666))


## generate new predicted values
pred_DR$pred <- predict(m_finDR, newdata = pred_DR,
                         level = 0)


## generate confidence intervals around predicted values
# keep in mind this DOES NOT take uncertainty of random
# effects into account!

# create model matrix
ci <- model.matrix(formula(m_finDR)[-2], pred_DR)

# use matrix multiplication on model matrix & var-covar matrix (vcov)
predvar <- diag(ci %*% vcov(m_finDR) %*% t(ci))

# construct CIs w/ SEs
pred_DR$lower <- with(pred_DR, pred - 1.96*sqrt(predvar))
pred_DR$upper <- with(pred_DR, pred + 1.96*sqrt(predvar))


## Plot setup
str(pred_DR)
# convert RA categories back to factors
locDR$RA <- recode(locDR$Grp_Age, `RS` = "RS", `TS` = "TS", 
                  `RA` = "RA", `TA` = "TA") 
locDR$RA <- as.factor(locDR$RA)

pred_DR$RA <- recode(pred_DR$RA, `0` = "RS_TA_TS", `1` = "RA")
pred_DR$RA <- as.factor(pred_DR$RA)
pred_DR$RS <- recode(pred_DR$RS, `1` = "RS", `0` = "RA_RS_TA")
pred_DR$RS <- as.factor(pred_DR$RS)
pred_DR$TA <- recode(pred_DR$TA, `1` = "TA", `0` = "RA_RS_TS")
pred_DR$TA <- as.factor(pred_DR$TA)
pred_DR$TS <- recode(pred_DR$TS, `1` = "TS", `0` = "RA_RS_TA")
pred_DR$TS <- as.factor(pred_DR$TS)

pred_locRA <- pred_DR %>% filter(RA == "RA")
pred_locRS <- pred_DR %>% filter(RS == "RS")
pred_locRS$RA <- pred_locRS$RS
pred_locTA <- pred_DR %>% filter(TA == "TA")
pred_locTA$RA <- pred_locTA$TA
pred_locTS <- pred_DR %>%  filter(TS == "TS")
pred_locTS$RA <- pred_locTS$TS

# merge filtered datasets together
pred_locRARS <- rbind(pred_locRA, pred_locRS)
pred_locRARSTS <- rbind(pred_locRARS, pred_locTS)
pred_loc7 <- rbind(pred_locRARSTS, pred_locTA)

str(pred_loc7$RA)
levels(pred_loc7$RA)

pred_DR <- pred_loc7 %>% filter(!RA == "RS_TA_TS" & !RA == "RA_RS_TA" & 
                                  !RA == "RA_RS_TS") %>% droplevels()

# check to make sure levels are correct
str(pred_DR$RA)
levels(pred_DR$RA)

str(locDR$RA)
levels(locDR$RA)

# back-transform predicted values
pred_DR$ExpPred <- exp(pred_DR$pred) - 1

# back-transform CIs
pred_DR$Expupper <- exp(pred_DR$upper) - 1
pred_DR$Explower <- exp(pred_DR$lower)  - 1

# exponentiate Log(Step Duration)
pred_DR$DaysSinceRel <- exp(pred_DR$LogDaysRel)

# change all negative ExpLower and ExpPred values to very small value
# (can't have negative numbers on log scale for plotting y axis)
pred_DR$Explower[pred_DR$Explower < 0] <- 0.00001
pred_DR$Explower[pred_DR$Expupper < 0] <- 0.00001
pred_DR$ExpPred[pred_DR$ExpPred < 0] <- 0.00001

DR_plot <- ggplot(data = locDR, aes(x = DaysSinceRel, 
                                        y = DistRel_km,
                                        fill = RA)) +
  geom_point(aes(shape = RA, color = RA)) +
  geom_jitter(aes(shape = RA, color = RA), height = 0.75) +
  scale_color_manual(values = c("RA"="chocolate4",
                                "RS"="gray50",
                                "TA"="purple4",
                                "TS"="darkolivegreen")) +
  geom_ribbon(data = pred_DR, # add CIs around line
              aes(y = NULL, 
                  ymin = Explower, 
                  ymax = Expupper,            
                  color = NULL, 
                  fill = RA),
              alpha = .2) + # darkness of fill
  geom_line(data = pred_DR, aes(y = ExpPred,
                                 group = RA, color = RA), 
            size = 0.8) +
  scale_fill_manual(values = c("RA"="chocolate4",
                               "RS"="gray50",
                               "TA"="purple4",
                               "TS"="darkolivegreen")) +
  labs(y="Distance from Release (km)", 
       x= "Time Since Release (days)") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) + 
  scale_x_continuous(breaks = seq(0,160, by = 20)) + 
  scale_y_continuous(trans = 'log10', labels = label_comma())


DR_plot


#########################################################
##Step Length Analysis #### 
#########################################################

library(tidyverse)
library(ggplot2)
library(scales)
library(cowplot)
library(nlme)
library(MuMIn)

rm(list=ls())
graphics.off()
options(stringsAsFactors = FALSE)

wd <- "C:/Users/Documents"
setwd(wd)


# load data
loc <- read.csv("./3_StepLength.csv",
                header = T, na.strings = c("", "N/A", "NA"))

# This dataset has removed step lengths = 0

## remove release rows, step length = NA
# (should already be removed but just in case)
loc2 <- loc %>% filter(!LocType == "release")

# check structure
str(loc2)

loc2$BvrID <- as.factor(loc2$BvrID)
loc2$Grp_Age <- as.factor(loc2$Grp_Age)
loc2$Sex <- as.factor(loc2$Sex)
loc2$Dischg_start <- as.factor(loc2$Dischg_start)
loc2$StepDuration_days <- as.numeric(loc2$StepDuration_days)
# step duration and days since release are in days
# step length is in meters

## transform all step lengths
loc2$LogStepLen <- log(loc2$StepLength)

## transform step duration
loc2$LogStepDur <- log(loc2$StepDuration)

## transform days since rel
loc2$LogDaysRel <- log(loc2$DaysSinceRel)

# add a column coding for whether the detection
# was a PIA detection 
unique(loc2$LocType)
loc2$LocType <- as.factor(loc2$LocType)

loc3 <- loc2 %>% mutate(Is_PIT = if_else(LocType == "pit", 1, 0))

##filter to include only pts w/ step duration < 60.9 days (2 months)
loc4 <- loc3 %>% filter(StepDuration_days <= 60.8)

# keep only steps > 0 length (should already be filtered out)
loc5 <- loc4 %>% filter(StepLength_m > 0)

# remove steps w/ NDVI < 0 (Green River, water)
loc6 <- loc5 %>% filter(NDVI_start >= 0)

## remove all steps > 182 days post-release (6 months)
# b/c very few indiv detected after this point
loc7 <- loc6 %>% filter(DaysSinceRel < 183)

## restructure variables for stepwise selection
loc7$RA <- as.numeric(loc7$Grp_Age == "RA")
loc7$RS <- as.numeric(loc7$Grp_Age == "RS")
loc7$TS <- as.numeric(loc7$Grp_Age == "TS")
loc7$TA <- as.numeric(loc7$Grp_Age == "TA")
loc7$DC_low <- as.numeric(loc7$Dischg_start == "Low")
loc7$DC_high <- as.numeric(loc7$Dischg_start == "High")
loc7$DC_med <- as.numeric(loc7$Dischg_start == "Median")

## remove all indivduals w/ < 2 samples
locSL <- loc7 %>% group_by(BvrID) %>% filter(n() >= 2)

## SCALE continuous variables (NDVI, ln(DaysSinceRel, DaysHeld))
# do NOT need to scale log(StepDur) b/c not interested in the effect size
# of this covariate (b/c always included in all models)
LogDaysRel_mean <- mean(locSL$LogDaysRel)
LogDaysRel_sd <- sd(locSL$LogDaysRel)
locSL$ScaleLogDaysRel <- (locSL$LogDaysRel - LogDaysRel_mean) / LogDaysRel_sd

NDVI_start_mean <- mean(locSL$NDVI_start)
NDVI_start_sd <- sd(locSL$NDVI_start)
locSL$ScaleNDVI_start <- (locSL$NDVI_start - NDVI_start_mean) / NDVI_start_sd

## Construct model 
summary(m_finSL <- lme(LogStepLen ~ 0 + RA + RS + TS + TA +
                        Is_PIT + DC_low + DC_high + ScaleNDVI_start +
                        RA:ScaleLogDaysRel + RS:ScaleLogDaysRel + 
                        TS:ScaleLogDaysRel + TA:ScaleLogDaysRel +
                        RA:LogStepDur + RS:LogStepDur + 
                        TS:LogStepDur + TA:LogStepDur, 
                      random = ~ 1 + LogStepDur | BvrID, data = locSL))


##test whether adding BvrID as Random Effect improves model fit
require(lmerTest)
require(lme4)

summary(m_fin_testSL <- lmer(LogStepLen ~ 0 + RA + RS + TS + TA +
                         Is_PIT + DC_low + DC_high + ScaleNDVI_start +
                         RA:ScaleLogDaysRel + RS:ScaleLogDaysRel + 
                         TS:ScaleLogDaysRel + TA:ScaleLogDaysRel +
                         RA:LogStepDur + RS:LogStepDur + 
                         TS:LogStepDur + TA:LogStepDur + 
                         (1 + LogStepDur | BvrID), data = locSL))

ranova(m_fin_testSL) # p-value is significant

# look at pseudo R2 values
pseudoR <- r.squaredGLMM(m_finSL)
pseudoR 

# check normality
qqnorm(residuals(m_finSL))
hist(residuals(m_finSL))


#### Predict LogStepLength x LogStepDuration ####

# figure out values needed to build dataframe
max(locSL$LogStepDur)
min(locSL$LogStepDur)
LogStepDur <- rep(seq(min(locSL$LogStepDur), max(locSL$LogStepDur), length.out = 666), 2)
mean(locSL$ScaleLogDaysRel)
mean(locSL$ScaleNDVI_start)

# Here is the model:
# m_finSL <- lme(LogStepLen ~ 0 + RA + RS + TS + TA +
#                  Is_PIT + DC_low + DC_high + ScaleNDVI_start +
#                  RA:ScaleLogDaysRel + RS:ScaleLogDaysRel + 
#                  TS:ScaleLogDaysRel + TA:ScaleLogDaysRel +
#                  RA:LogStepDur + RS:LogStepDur + 
#                  TS:LogStepDur + TA:LogStepDur, 
#                random = ~ 1 + LogStepDur | BvrID, data = locSL))

## build the dataframe
pred_SL <- data.frame(LogStepDur = LogStepDur, 
                      RA = rep(c(1, 0, 0, 0), each = 666),
                      RS = rep(c(0, 1, 0, 0), each = 666),
                      TA = rep(c(0, 0, 1, 0), each = 666),
                      TS = rep(c(0, 0, 0, 1), each = 666),
                       Is_PIT = rep(1,2664), # hold constant
                       DC_low = rep(0,2664), # hold constant
                       DC_high = rep(0,2664),
                       ScaleNDVI_start = rep(-2.132716e-16, 2664), # mean value
                       ScaleLogDaysRel = rep(-7.218419e-17, 2664)) # mean value

## generate new predicted values
pred_SL$pred <- predict(m_finSL, newdata = pred_SL,
                         level = 0)

## generate confidence intervals around predicted values

# create model matrix
ci <- model.matrix(formula(m_finSL)[-2], pred_SL) 

# use matrix mult on model matrix & var-covar matrix (vcov)
# pull all values on diagonal - same var as predicted values
predvar <- diag(ci %*% vcov(m_finSL) %*% t(ci))

# construct approx CIs w/ SEs (uses 2 here for t mult, or
# 1.96 as z mult) and add to new dataset
# both assume normal distribution?
pred_SL$lower <- with(pred_SL, pred - 1.96*sqrt(predvar))
pred_SL$upper <- with(pred_SL, pred + 1.96*sqrt(predvar))

## Plot Setup
str(pred_SL)

# convert RA categories back to factors
locSL$RA <- recode(locSL$Grp_Age, `RS` = "RS", `TS` = "TS", 
                   `RA` = "RA", `TA` = "TA") 
locSL$RA <- as.factor(locSL$RA)

pred_SL$RA <- recode(pred_SL$RA, `0` = "RS_TA_TS", `1` = "RA")
pred_SL$RA <- as.factor(pred_SL$RA)
pred_SL$RS <- recode(pred_SL$RS, `1` = "RS", `0` = "RA_RS_TA")
pred_SL$RS <- as.factor(pred_SL$RS)
pred_SL$TA <- recode(pred_SL$TA, `1` = "TA", `0` = "RA_RS_TS")
pred_SL$TA <- as.factor(pred_SL$TA)
pred_SL$TS <- recode(pred_SL$TS, `1` = "TS", `0` = "RA_RS_TA")
pred_SL$TS <- as.factor(pred_SL$TS)

pred_locRA <- pred_SL %>% filter(RA == "RA")
pred_locRS <- pred_SL %>% filter(RS == "RS")
pred_locRS$RA <- pred_locRS$RS
pred_locTA <- pred_SL %>% filter(TA == "TA")
pred_locTA$RA <- pred_locTA$TA
pred_locTS <- pred_SL %>%  filter(TS == "TS")
pred_locTS$RA <- pred_locTS$TS

# merge filtered datasets together
pred_locRARS <- rbind(pred_locRA, pred_locRS)
pred_locRARSTS <- rbind(pred_locRARS, pred_locTS)
pred_loc7 <- rbind(pred_locRARSTS, pred_locTA)

str(pred_loc7$RA)
levels(pred_loc7$RA)

pred_SL <- pred_loc7 %>% filter(!RA == "RS_TA_TS" & !RA == "RA_RS_TA" & 
                                  !RA == "RA_RS_TS") %>% droplevels()

# check to make sure levels are correct
str(pred_SL$RA)
levels(pred_SL$RA)

str(locSL$RA)
levels(locSL$RA)

# back-transform predicted values
pred_SL$ExpPred <- exp(pred_SL$pred)

# back-transform CIs
pred_SL$Expupper <- exp(pred_SL$upper)
pred_SL$Explower <- exp(pred_SL$lower)

# exponentiate Log(Step Duration)
pred_SL$StepDuration_days <- exp(pred_SL$LogStepDur)

## Plot
SL_plot <- ggplot(data = locSL, aes(x = StepDuration_days, 
                                         y = StepLength_m,
                                         fill = RA)) +
  geom_point(aes(shape = RA, color = RA)) +
  geom_jitter(aes(shape = RA, color = RA), height = 0.5, width = 5) +
  scale_color_manual(values = c("RA"="chocolate4",
                                "RS"="gray50",
                                "TA"="purple4",
                                "TS"="darkolivegreen")) +
  geom_ribbon(data = pred_SL, # add CIs around line
              aes(y = NULL, 
                  ymin = Explower, 
                  ymax = Expupper,            
                  color = NULL, 
                  fill = RA),
              alpha = .2) + # darkness of fill
  geom_line(data = pred_SL, aes(y = ExpPred, 
                                 group = RA, color = RA), 
            size = 0.8) +
  scale_fill_manual(values = c("RA"="chocolate4",
                               "RS"="gray50",
                               "TA"="purple4",
                               "TS"="darkolivegreen")) +
  labs(y="Step Length (m)", 
       x= "Step Duration (days)") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) +
  coord_cartesian(xlim = c(0,50))+
  scale_x_continuous(breaks = seq(0,50, by = 5)) + 
  scale_y_continuous(trans = 'log10', labels = label_comma())

SL_plot 


#### Speed at multiple scales ####


#### Hourly Rate ####

# Here is the model:
# m_finSL <- lme(LogStepLen ~ 0 + RA + RS + TS + TA +
#                  Is_PIT + DC_low + DC_high + ScaleNDVI_start +
#                  RA:ScaleLogDaysRel + RS:ScaleLogDaysRel + 
#                  TS:ScaleLogDaysRel + TA:ScaleLogDaysRel +
#                  RA:LogStepDur + RS:LogStepDur + 
#                  TS:LogStepDur + TA:LogStepDur, 
#                random = ~ 1 + LogStepDur | BvrID, data = locSL))

LogDaysRel <- rep(seq(0, max(locSL$LogDaysRel), length.out = 666), 4)

## build the dataframe
pred_SL <- data.frame(LogStepDur = (log(1/24)), # hourly step duration 
                      RA = rep(c(1, 0, 0, 0), each = 666),
                      RS = rep(c(0, 1, 0, 0), each = 666),
                      TA = rep(c(0, 0, 1, 0), each = 666),
                      TS = rep(c(0, 0, 0, 1), each = 666),
                      Is_PIT = rep(1,2664), # hold constant
                      DC_low = rep(0,2664), # hold constant
                      DC_high = rep(0,2664),
                      ScaleNDVI_start = rep(-2.132716e-16, 2664), # mean value
                      ScaleLogDaysRel = (LogDaysRel - LogDaysRel_mean) / LogDaysRel_sd)

## generate new predicted values
pred_SL$pred <- predict(m_finSL, newdata = pred_SL,
                        level = 0)

## generate confidence intervals around predicted values

# create model matrix
ci <- model.matrix(formula(m_finSL)[-2], pred_SL) 

# use matrix mult on model matrix & var-covar matrix (vcov)
# pull all values on diagonal - same var as predicted values
predvar <- diag(ci %*% vcov(m_finSL) %*% t(ci))

# construct approx CIs w/ SEs (uses 2 here for t mult, or
# 1.96 as z mult) and add to new dataset
# both assume normal distribution?
pred_SL$lower <- with(pred_SL, pred - 1.96*sqrt(predvar))
pred_SL$upper <- with(pred_SL, pred + 1.96*sqrt(predvar))

## Plot, held constant at low discharge, PIA detection = Yes, mean NDVI, mean Days Since Release
str(pred_SL)

# convert RA categories back to factors
locSL$RA <- recode(locSL$Grp_Age, `RS` = "RS", `TS` = "TS", 
                   `RA` = "RA", `TA` = "TA") 
locSL$RA <- as.factor(locSL$RA)

pred_SL$RA <- recode(pred_SL$RA, `0` = "RS_TA_TS", `1` = "RA")
pred_SL$RA <- as.factor(pred_SL$RA)
pred_SL$RS <- recode(pred_SL$RS, `1` = "RS", `0` = "RA_RS_TA")
pred_SL$RS <- as.factor(pred_SL$RS)
pred_SL$TA <- recode(pred_SL$TA, `1` = "TA", `0` = "RA_RS_TS")
pred_SL$TA <- as.factor(pred_SL$TA)
pred_SL$TS <- recode(pred_SL$TS, `1` = "TS", `0` = "RA_RS_TA")
pred_SL$TS <- as.factor(pred_SL$TS)

pred_locRA <- pred_SL %>% filter(RA == "RA")
pred_locRS <- pred_SL %>% filter(RS == "RS")
pred_locRS$RA <- pred_locRS$RS
pred_locTA <- pred_SL %>% filter(TA == "TA")
pred_locTA$RA <- pred_locTA$TA
pred_locTS <- pred_SL %>%  filter(TS == "TS")
pred_locTS$RA <- pred_locTS$TS

# merge filtered datasets together
pred_locRARS <- rbind(pred_locRA, pred_locRS)
pred_locRARSTS <- rbind(pred_locRARS, pred_locTS)
pred_loc7 <- rbind(pred_locRARSTS, pred_locTA)

str(pred_loc7$RA)
levels(pred_loc7$RA)

pred_SL <- pred_loc7 %>% filter(!RA == "RS_TA_TS" & !RA == "RA_RS_TA" & 
                                  !RA == "RA_RS_TS") %>% droplevels()

# check to make sure levels are correct
str(pred_SL$RA)
levels(pred_SL$RA)

str(locSL$RA)
levels(locSL$RA)

# back-transform predicted values
pred_SL$ExpPred <- exp(pred_SL$pred)

# back-transform CIs
pred_SL$Expupper <- exp(pred_SL$upper)
pred_SL$Explower <- exp(pred_SL$lower)

# exponentiate Log(Step Duration)
pred_SL$StepDuration_days <- exp(pred_SL$LogStepDur)

# exponentiate Log(DaysRel)
pred_SL$DaysRel <- exp((LogDaysRel * LogDaysRel_sd) - LogDaysRel_mean)


## Plot speed vs. time since release on HOURLY scale
Plot_hr <- ggplot(data = pred_SL, aes(x = DaysRel, 
                                    y = ExpPred,
                                    fill = RA)) +
  scale_color_manual(values = c("RA"="chocolate4",
                                "RS"="gray50",
                                "TA"="purple4",
                                "TS"="darkolivegreen")) +
  geom_line(data = pred_SL, aes(y = ExpPred, 
                                group = RA, color = RA), 
            size = 1.25) +
  scale_fill_manual(values = c("RA"="chocolate4",
                               "RS"="gray50",
                               "TA"="purple4",
                               "TS"="darkolivegreen")) +
  labs(y="Speed (m/hour)", 
       x= "Time Since Release (days)") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))

Plot_hr 


#### Daily Rate ####

# Here is the model:
# m_finSL <- lme(LogStepLen ~ 0 + RA + RS + TS + TA +
#                  Is_PIT + DC_low + DC_high + ScaleNDVI_start +
#                  RA:ScaleLogDaysRel + RS:ScaleLogDaysRel + 
#                  TS:ScaleLogDaysRel + TA:ScaleLogDaysRel +
#                  RA:LogStepDur + RS:LogStepDur + 
#                  TS:LogStepDur + TA:LogStepDur, 
#                random = ~ 1 + LogStepDur | BvrID, data = locSL))

LogDaysRel <- rep(seq(0, max(locSL$LogDaysRel), length.out = 666), 4)

## build the dataframe
pred_SL <- data.frame(LogStepDur = (log(1)), # daily step duration
                      RA = rep(c(1, 0, 0, 0), each = 666),
                      RS = rep(c(0, 1, 0, 0), each = 666),
                      TA = rep(c(0, 0, 1, 0), each = 666),
                      TS = rep(c(0, 0, 0, 1), each = 666),
                      Is_PIT = rep(1,2664), # hold constant
                      DC_low = rep(0,2664), # hold constant
                      DC_high = rep(0,2664),
                      ScaleNDVI_start = rep(-2.132716e-16, 2664), # mean value
                      ScaleLogDaysRel = (LogDaysRel - LogDaysRel_mean) / LogDaysRel_sd)

## generate new predicted values
pred_SL$pred <- predict(m_finSL, newdata = pred_SL,
                        level = 0)

## generate confidence intervals around predicted values

# create model matrix
ci <- model.matrix(formula(m_finSL)[-2], pred_SL) 

# use matrix mult on model matrix & var-covar matrix (vcov)
# pull all values on diagonal - same var as predicted values
predvar <- diag(ci %*% vcov(m_finSL) %*% t(ci))

# construct approx CIs w/ SEs (uses 2 here for t mult, or
# 1.96 as z mult) and add to new dataset
# both assume normal distribution?
pred_SL$lower <- with(pred_SL, pred - 1.96*sqrt(predvar))
pred_SL$upper <- with(pred_SL, pred + 1.96*sqrt(predvar))

## Plot, held constant at low discharge, PIA detection = Yes, mean NDVI, mean Days Since Release
str(pred_SL)

# convert RA categories back to factors
locSL$RA <- recode(locSL$Grp_Age, `RS` = "RS", `TS` = "TS", 
                   `RA` = "RA", `TA` = "TA") 
locSL$RA <- as.factor(locSL$RA)

pred_SL$RA <- recode(pred_SL$RA, `0` = "RS_TA_TS", `1` = "RA")
pred_SL$RA <- as.factor(pred_SL$RA)
pred_SL$RS <- recode(pred_SL$RS, `1` = "RS", `0` = "RA_RS_TA")
pred_SL$RS <- as.factor(pred_SL$RS)
pred_SL$TA <- recode(pred_SL$TA, `1` = "TA", `0` = "RA_RS_TS")
pred_SL$TA <- as.factor(pred_SL$TA)
pred_SL$TS <- recode(pred_SL$TS, `1` = "TS", `0` = "RA_RS_TA")
pred_SL$TS <- as.factor(pred_SL$TS)

pred_locRA <- pred_SL %>% filter(RA == "RA")
pred_locRS <- pred_SL %>% filter(RS == "RS")
pred_locRS$RA <- pred_locRS$RS
pred_locTA <- pred_SL %>% filter(TA == "TA")
pred_locTA$RA <- pred_locTA$TA
pred_locTS <- pred_SL %>%  filter(TS == "TS")
pred_locTS$RA <- pred_locTS$TS

# merge filtered datasets together
pred_locRARS <- rbind(pred_locRA, pred_locRS)
pred_locRARSTS <- rbind(pred_locRARS, pred_locTS)
pred_loc7 <- rbind(pred_locRARSTS, pred_locTA)

str(pred_loc7$RA)
levels(pred_loc7$RA)

pred_SL <- pred_loc7 %>% filter(!RA == "RS_TA_TS" & !RA == "RA_RS_TA" & 
                                  !RA == "RA_RS_TS") %>% droplevels()

# check to make sure levels are correct
str(pred_SL$RA)
levels(pred_SL$RA)

str(locSL$RA)
levels(locSL$RA)

# back-transform predicted values
pred_SL$ExpPred <- exp(pred_SL$pred)

# back-transform CIs
pred_SL$Expupper <- exp(pred_SL$upper)
pred_SL$Explower <- exp(pred_SL$lower)

# exponentiate Log(Step Duration)
pred_SL$StepDuration_days <- exp(pred_SL$LogStepDur)

# exponentiate Log(DaysRel)
pred_SL$DaysRel <- exp((LogDaysRel * LogDaysRel_sd) - LogDaysRel_mean)

## Plot speed vs. time since release on DAILY scale
Plot_day <- ggplot(data = pred_SL, aes(x = DaysRel, 
                                      y = ExpPred,
                                      fill = RA)) +
  scale_color_manual(values = c("RA"="chocolate4",
                                "RS"="gray50",
                                "TA"="purple4",
                                "TS"="darkolivegreen")) +
  geom_line(data = pred_SL, aes(y = ExpPred, 
                                group = RA, color = RA), 
            size = 1.25) +
  scale_fill_manual(values = c("RA"="chocolate4",
                               "RS"="gray50",
                               "TA"="purple4",
                               "TS"="darkolivegreen")) +
  labs(y="Speed (m/day)", 
       x= "Time Since Release (days)") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))

Plot_day 


#### Monthly Rate ####

# Here is the model:
# m_finSL <- lme(LogStepLen ~ 0 + RA + RS + TS + TA +
#                  Is_PIT + DC_low + DC_high + ScaleNDVI_start +
#                  RA:ScaleLogDaysRel + RS:ScaleLogDaysRel + 
#                  TS:ScaleLogDaysRel + TA:ScaleLogDaysRel +
#                  RA:LogStepDur + RS:LogStepDur + 
#                  TS:LogStepDur + TA:LogStepDur, 
#                random = ~ 1 + LogStepDur | BvrID, data = locSL))

LogDaysRel <- rep(seq(0, max(locSL$LogDaysRel), length.out = 666), 4)

## build the dataframe
pred_SL <- data.frame(LogStepDur = (log(30)), #monthly step duration
                      RA = rep(c(1, 0, 0, 0), each = 666),
                      RS = rep(c(0, 1, 0, 0), each = 666),
                      TA = rep(c(0, 0, 1, 0), each = 666),
                      TS = rep(c(0, 0, 0, 1), each = 666),
                      Is_PIT = rep(1,2664), # hold constant
                      DC_low = rep(0,2664), # hold constant
                      DC_high = rep(0,2664),
                      ScaleNDVI_start = rep(-2.132716e-16, 2664), # mean value
                      ScaleLogDaysRel = (LogDaysRel - LogDaysRel_mean) / LogDaysRel_sd)

## generate new predicted values
pred_SL$pred <- predict(m_finSL, newdata = pred_SL,
                        level = 0)

## generate confidence intervals around predicted values

# create model matrix
ci <- model.matrix(formula(m_finSL)[-2], pred_SL) 

# use matrix mult on model matrix & var-covar matrix (vcov)
# pull all values on diagonal - same var as predicted values
predvar <- diag(ci %*% vcov(m_finSL) %*% t(ci))

# construct approx CIs w/ SEs (uses 2 here for t mult, or
# 1.96 as z mult) and add to new dataset
# both assume normal distribution?
pred_SL$lower <- with(pred_SL, pred - 1.96*sqrt(predvar))
pred_SL$upper <- with(pred_SL, pred + 1.96*sqrt(predvar))

## Plot, held constant at low discharge, PIA detection = Yes, mean NDVI, mean Days Since Release
str(pred_SL)

# convert RA categories back to factors
locSL$RA <- recode(locSL$Grp_Age, `RS` = "RS", `TS` = "TS", 
                   `RA` = "RA", `TA` = "TA") 
locSL$RA <- as.factor(locSL$RA)

pred_SL$RA <- recode(pred_SL$RA, `0` = "RS_TA_TS", `1` = "RA")
pred_SL$RA <- as.factor(pred_SL$RA)
pred_SL$RS <- recode(pred_SL$RS, `1` = "RS", `0` = "RA_RS_TA")
pred_SL$RS <- as.factor(pred_SL$RS)
pred_SL$TA <- recode(pred_SL$TA, `1` = "TA", `0` = "RA_RS_TS")
pred_SL$TA <- as.factor(pred_SL$TA)
pred_SL$TS <- recode(pred_SL$TS, `1` = "TS", `0` = "RA_RS_TA")
pred_SL$TS <- as.factor(pred_SL$TS)

pred_locRA <- pred_SL %>% filter(RA == "RA")
pred_locRS <- pred_SL %>% filter(RS == "RS")
pred_locRS$RA <- pred_locRS$RS
pred_locTA <- pred_SL %>% filter(TA == "TA")
pred_locTA$RA <- pred_locTA$TA
pred_locTS <- pred_SL %>%  filter(TS == "TS")
pred_locTS$RA <- pred_locTS$TS

# merge filtered datasets together
pred_locRARS <- rbind(pred_locRA, pred_locRS)
pred_locRARSTS <- rbind(pred_locRARS, pred_locTS)
pred_loc7 <- rbind(pred_locRARSTS, pred_locTA)

str(pred_loc7$RA)
levels(pred_loc7$RA)

pred_SL <- pred_loc7 %>% filter(!RA == "RS_TA_TS" & !RA == "RA_RS_TA" & 
                                  !RA == "RA_RS_TS") %>% droplevels()

# check to make sure levels are correct
str(pred_SL$RA)
levels(pred_SL$RA)

str(locSL$RA)
levels(locSL$RA)

# back-transform predicted values
pred_SL$ExpPred <- exp(pred_SL$pred)

# back-transform CIs
pred_SL$Expupper <- exp(pred_SL$upper)
pred_SL$Explower <- exp(pred_SL$lower)

# exponentiate Log(Step Duration)
pred_SL$StepDuration_days <- exp(pred_SL$LogStepDur)

# exponentiate Log(DaysRel)
pred_SL$DaysRel <- exp((LogDaysRel * LogDaysRel_sd) - LogDaysRel_mean)


## Plot speed vs. time since release on MONTHLY scale
Plot_month <- ggplot(data = pred_SL, aes(x = DaysRel, 
                                      y = ExpPred,
                                      fill = RA)) +
  scale_color_manual(values = c("RA"="chocolate4",
                                "RS"="gray50",
                                "TA"="purple4",
                                "TS"="darkolivegreen")) +
  geom_line(data = pred_SL, aes(y = ExpPred, 
                                group = RA, color = RA), 
            size = 1.25) +
  scale_fill_manual(values = c("RA"="chocolate4",
                               "RS"="gray50",
                               "TA"="purple4",
                               "TS"="darkolivegreen")) +
  labs(y="Speed (m/month)", 
       x= "Time Since Release (days)") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))

Plot_month 



#########################################################
##Fine-scale Movement Rate Analysis #### 
#########################################################

# First, using the RawData_FineScaleMvmt tab (these are the observer's locations) in the
# Dodenetal_BvrMvmtMultScales spreadsheet, load into ArcPro and use Bearing Distance to Line 
# and the Intersect tools to generate estimated beaver location coordinates resembling the EstLocs_FineScaleMvmt tab


#### Calculate the distance between sequential points along a river ####

rm(list=ls())
graphics.off()
options(stringsAsFactors = FALSE)

wd <- "C:/Users/Documents"
setwd(wd)

# load in the estimated beaver locations created in ArcPro (see note above)
intmvmt <- read.csv("./EstLocs_FineScaleMvmt.csv",
                    header = T, na.strings = c("", "N/A", "NA"))

# Add a column counting from 1 for each bvr
intmvmt2 <- intmvmt %>% group_by(Data_ID) %>% mutate(Survey = 1:n())

# export intmvmt2 as an updated CSV, then separate into beavers from study sites
# Price and Moonshine (one file, name "PrMnMerge.csv"), and Cottonwood (second file, name "Cottonwood.csv)


## Price and Moonshine Study Sites

## load data
prmnloc <- read.csv("./PrMnMerge.csv",
                    header = T, na.strings = c("", "N/A", "NA"))

## import river network
prmnriv2 <- line2network(path = "C:/Users/Documents/River_Shapefiles",
                         layer = "Extend_PriceGrnprmn_MergeProj_20210126")

# add vertices for greater precision in calculations:
prmnriv <- addverts(prmnriv2, mindist=0.1) # add a vertex every 0.1 m

# plot to make sure looks ok
plot(prmnriv)

## Convert XY data to snapped river locs
prmnsnap <- xy2segvert(x=prmnloc$POINT_X, y=prmnloc$POINT_Y, 
                       rivers=prmnriv)

head(prmnsnap)
hist(prmnsnap$snapdist,main = "snapping distance (m)")

# plot to make sure looks ok
zoomtoseg(seg=c(1:12), rivers=prmnriv) # zoom to segment of interest
points(prmnloc$POINT_X, prmnloc$POINT_Y, pch=16, col="red") # orig pts
riverpoints(seg=prmnsnap$seg, vert=prmnsnap$vert, rivers = prmnriv,
            pch=15, col = "blue") # pts snapped to river


## compute network distance btw sequential observations of an individual:
dist <- riverdistanceseq(unique=prmnloc$Data_ID, survey=prmnloc$Survey,
                         seg=prmnsnap$seg,vert=prmnsnap$vert,
                         rivers=prmnriv)
# make sure to use Data_ID (not BvrID) for unique or will only keep one sampling session per beaver

dist

# export as a CSV


## Cottonwood Study Site

## load data
cottloc <- read.csv("./Cottonwood.csv",
                    header = T, na.strings = c("", "N/A", "NA"))

## import river network, make sure to use MERGED shapefile so doesn't remove any segments!
cottriv2 <- line2network(path = "C:/Users/Documents/River_Shapefiles",
                         layer = "Cottonwood_MergeProjUpdated_20210326")

# add vertices for greater precision in calculations:
cottriv <- addverts(cottriv2, mindist=0.1) # add a vertex every 0.1 m

# plot to make sure looks ok
plot(cottriv)

## Convert XY data to snapped river locs
cottsnap <- xy2segvert(x=cottloc$POINT_X, y=cottloc$POINT_Y, 
                       rivers=cottriv)

head(cottsnap)
hist(cottsnap$snapdist,main = "snapping distance (m)")

# plot to make sure looks ok
zoomtoseg(seg=c(5), rivers=cottriv) # zoom to segment of interest
points(cottloc$POINT_X, cottloc$POINT_Y, pch=16, col="red") # orig pts
riverpoints(seg=cottsnap$seg, vert=cottsnap$vert, rivers = cottriv,
            pch=15, col = "blue") # pts snapped to river


## compute network distance btw sequential observations of an individual:
dist <- riverdistanceseq(unique=cottloc$Data_ID, survey=cottloc$Survey,
                         seg=cottsnap$seg,vert=cottsnap$vert,
                         rivers=cottriv)
# make sure to use Data_ID (not BvrID) for unique or will only keep one sampling session per beaver

dist

# export as a CSV 

## next, merge all of the distance columns into a single column,
# and lastly merge the two PrMn and Cottonwood CSVs together, and add the dropped variables 
# back to the dataset. This should be the final product needed to model median step length below


#### Modelling Median Step Length ####

#set options
rm(list=ls())
graphics.off()
options(stringsAsFactors = FALSE)

wd <- "C:/Users/Documents"
setwd(wd)

# load data
loc1 <- read.csv("./4_FineScaleMvmt.csv",
                 header = T, na.strings = c("", "N/A", "NA"))

# check structure
str(loc1)

loc1$BvrID <- as.factor(loc1$BvrID)
loc1$Grp_Age <- as.factor(loc1$Grp_Age)
loc1$DayNight <- as.factor(loc1$DayNight)

## transform mean/median dist traveled
loc1$LogMed_Distm <- log(loc1$Median_Distm)

## restructure vars for stepwise selection
loc1$RA <- as.numeric(loc1$Grp_Age == "RA")
loc1$TS <- as.numeric(loc1$Grp_Age == "TS")
loc1$TA <- as.numeric(loc1$Grp_Age == "TA")
loc1$Night <- as.numeric(loc1$DayNight == "Night")


## Create box-and-whisker plot of data
# boxplot of breakdown of mean distance moved in 5 min by bvr state category

# check structure
str(loc1)

loc1$Grp_Age <- as.factor(loc1$Grp_Age)
loc1$DayNight <- as.factor(loc1$DayNight)

#merge Grp Age and DayNight columns
loc1$GrpTOD <- paste(loc1$Grp_Age, loc1$DayNight)

plot_FS <- ggplot(loc1) + 
  geom_boxplot(aes(x=GrpTOD, y=Median_Distm, fill = GrpTOD)) +
  theme_classic() +
  scale_fill_manual(values = c("RA Day"="chocolate4",
                               "RA Night" = "chocolate4",
                               "TA Day"="purple4",
                               "TA Night"="purple4",
                               "TS Day"="darkolivegreen",
                               "TS Night" = "darkolivegreen")) +
  labs(x="State Category",y="Median distance traveled in 5 minutes (m)") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")) 

plot_FS


## Construct model
summary(m_finFS <- lme(LogMed_Distm ~ 0 + RA + TS + TA + Night, 
                      random = ~1|BvrID, 
                      data = loc1)) 

## test to see whether adding BvrID as RE improves model
require(lmerTest)
require(lme4)
summary(m_finFS_test <- lmer(LogMed_Distm ~ 
                             Night +
                             ScaleMeanNDVI + 
                             ScaleLogDaysRel + (1|BvrID), 
                           data = loc1))

ranova(m_finFS_test) # p-value of LRT is NOT significant

# look at pseudo R2 values
pseudoR <- r.squaredGLMM(m_finFS)
pseudoR 

# check normality
qqnorm(residuals(m_finFS))
hist(residuals(m_finFS)) 

