################################INSTALL AND LOAD PACKAGES################################
rm(list=ls())
library(lme4)
library(lmerTest)
library(MuMIn)
library(emmeans)
library(plotly)
library(tidyr)
library(dplyr)

###Import data set###
setwd('~/Dropbox/FBQ Moorea 2020/Distracted Clams/1) Data set/')
DistractedClams <- read.csv("Clam_Multimodal_Datasheet_FINAL.csv")
Exemplartest <- read.csv("exemplar test.csv")

################################checking distributions before analysis################################
###Checking normal distribution of variables using histograms
##Log transform data that is not normally distributed

#Retractions#
hist(DistractedClams$Retractions) #not normal
DistractedClams$log_retractions <- log10(DistractedClams$Retractions+1)
hist(DistractedClams$log_retractions) #normal

hist(Exemplartest$Retractions) #not normal
Exemplartest$log_retractions <- log10(Exemplartest$Retractions+1)
hist(Exemplartest$log_retractions) #normal

#Latency to Close#
hist(DistractedClams$Latency_to_close) #normal
hist(Exemplartest$Latency_to_close) #normal

#Latency to Reemerge#
hist(DistractedClams$Latency_to_reemerge) #not normal
DistractedClams$log_ltr <- log10(DistractedClams$Latency_to_reemerge)
hist(DistractedClams$log_ltr)#normal

hist(Exemplartest$Latency_to_reemerge) #not normal
Exemplartest$log_ltr <- log10(Exemplartest$Latency_to_reemerge)
hist(Exemplartest$log_ltr) #normal

###Testing for correlation between size and depth###
m0 <- lm(Depth~Size, data=DistractedClams)
summary(m0)
#No correlation found between size and depth of clams

###Set following variables as categorical variables###
DistractedClams$Plume_presence <- as.factor(DistractedClams$Plume_presence)
Exemplartest$Sound_exemplars <- as.factor(Exemplartest$Sound_exemplars)

################################Retraction Analyses################################
####Linear Mixed Effects Models####
###Control for fixed effects of Treatment, Trial number, Plume Presence, and Size
###Control for random effect of individual (Subject ID)


###MODEL: Set Control as Reference
DistractedClams$Treatment = factor(DistractedClams$Treatment,levels=c("Control","Flow","Sound", "Flow_and_Sound"))
m1 <- lmer(log_retractions~ Treatment + Trial_number + Plume_presence + Size + (1| Subject_ID), data = DistractedClams)
summary(m1)
#all treatments are significant from control (p < 0.001 for all)

###Histogram, QQ-plot, (f,r) plot ###
##Check residuals of linear mixed model for normality##
r <- residuals(m1)
hist(r)
qqnorm(r)
f <- fitted(m1)
plot(f,r)
#all normal


###MODEL: Set Flow as Reference
DistractedClams$Treatment = factor(DistractedClams$Treatment,levels=c("Flow", "Sound", "Flow_and_Sound", "Control"))
m1.2 <- lmer(log_retractions~ Treatment + Trial_number + Plume_presence + Size + (1| Subject_ID), data = DistractedClams)
summary(m1.2)
#flow is significantly diff from sound (p = 0.024) and control (p < 0.001), more retractions than both

###Histogram, QQ-plot, (f,r) plot###
r <- residuals(m1.2)
hist(r)
qqnorm(r)
f <- fitted(m1.2)
plot(f,r)
#all normal


###MODEL: Set Sound as Reference
DistractedClams$Treatment = factor(DistractedClams$Treatment,levels=c("Sound", "Flow_and_Sound", "Control","Flow")) #Sets sound as reference
m1.3 <- lmer(log_retractions~ Treatment + Trial_number + Plume_presence + Size + (1| Subject_ID), data = DistractedClams)
summary(m1.3)
#sound is significantly different from control, flow, and multimodal flow and sound (p = 0.02985)

###Histogram, QQ-plot, (f,r) plot###
r <- residuals(m1.3)
hist(r)
qqnorm(r)
f <- fitted(m1.3)
plot(f,r)
#all normal


###MODEL: Set Multimodal Flow and Sound as Reference
DistractedClams$Treatment = factor(DistractedClams$Treatment,levels=c("Flow_and_Sound", "Control","Flow","Sound")) #Sets multimodal flow and sound as reference
m1.4 <- lmer(log_retractions~ Treatment + Trial_number + Plume_presence + Size + (1| Subject_ID), data = DistractedClams)
summary(m1.4)
#Multimodal flow and sound is only significant to control and sound

###Histogram, QQ-plot, (f,r) plot###
r <- residuals(m1.4)
hist(r)
qqnorm(r)
f <- fitted(m1.4)
plot(f,r)
#all normal


###Testing if exemplars had a significant effect###
m1.5 <- lm(log_retractions~ Sound_exemplars + Trial_number + Size, data = Exemplartest)
summary(m1.5)
anova(m1.5)
#exemplars are not significant


####ANOVA####
###Compare linear mixed model with linear model without random effect of individual.
###Test if random effect of individual is significant
DistractedClams$Treatment = factor(DistractedClams$Treatment,levels=c("Control","Flow","Sound", "Flow_and_Sound")) #set control as reference
m1.6 <- lm(log_retractions~ Treatment + Trial_number + Plume_presence + Size, data = DistractedClams)
anova(m1, m1.6)
#There is a significant effect of individual (p = 0.007281). Chisquare = 7.2023, Df = 1


####Emmeans###
###Estimated marginal means for each treatment###
retractionemmeans <- emmeans(m1, "Treatment")
retractionemmeans
plot(emmeans(m1, "Treatment"), comparisons = TRUE)
eff_size(retractionemmeans, sigma = sigma(m1), edf = 47) #Cohen's d effect size

###Graphing Emmeans###
##Retractions Emmeans Bar Graph##
se <- c("0.0405", "0.0381","0.0392","0.0387")
xform <- list(title= "Treatment",
              titlefont=list(size=20),
              tickfont=list(size=15),
              categoryorder = "array",
              categoryarray = c ("control", "flow only", "sound only", "flow and sound"))
r <- plot_ly(x = c("control", "flow only", "sound only", "flow and sound"),
             y = c(0.494, 0.774, 0.663, 0.771),
             name = "Retractions",
             type = "bar",
             error_y = list(array = ~ se,color = 'black'),
             marker=list (color='white',
                          line=list(color='black',
                                    width=1.5)))

r2 <- r %>% layout(
  xaxis= xform,
  yaxis=  list(title="Estimated Marginal Means \u00B1 SE",
               titlefont=list(size=20),
               tickvals = list(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,0.7, 0.8),
               tickformat = (".2f"),
               linecolor = toRGB("black"),
               showgrid=FALSE,
               showline=TRUE))

r2


################################Latency to Close Analyses################################
####Linear Mixed Effects Models####
###Control for fixed effects of Treatment, Trial number, Plume Presence, and Size
###Control for random effect of individual (Subject ID)


###MODEL: Set Control as Reference
DistractedClams$Treatment = factor(DistractedClams$Treatment,levels=c("Control","Flow","Sound", "Flow_and_Sound"))
m2 <- lmer(Latency_to_close~ Treatment + Trial_number + Plume_presence + Size + (1| Subject_ID), data = DistractedClams)
summary(m2)
#Multimodal flow and sound is significant, increased latency to close (p = 0.0369)
#Size is significant p = 0.046

###Histogram, QQ-plot, (f,r) plot###
r <- residuals(m2)
hist(r)
qqnorm(r)
f <- fitted(m2)
plot(f,r)
#all normal


###MODEL: Set Flow as Reference
DistractedClams$Treatment = factor(DistractedClams$Treatment,levels=c("Flow", "Sound", "Flow_and_Sound", "Control"))
m2.2 <- lmer(Latency_to_close~ Treatment + Trial_number + Plume_presence + Size + (1| Subject_ID), data = DistractedClams)
summary(m2.2)
#Flow is almost significant to flow and sound (p = 0.0963)

###Histogram, QQ-plot, (f,r) plot###
r <- residuals(m2.2)
hist(r)
qqnorm(r)
f <- fitted(m2.2)
plot(f,r)
#all normal


###MODEL: Set Sound as Reference
DistractedClams$Treatment = factor(DistractedClams$Treatment,levels=c("Sound", "Flow_and_Sound", "Control","Flow"))
m2.3 <- lmer(Latency_to_close~ Treatment + Trial_number + Plume_presence + Size + (1| Subject_ID), data = DistractedClams)
summary(m2.3)
#almost significant to flow and sound (p = 0.0534)

###Histogram, QQ-plot, (f,r) plot###
r <- residuals(m2.3)
hist(r)
qqnorm(r)
f <- fitted(m2.3)
plot(f,r)
#all normal


#MODEL: Set Multimodal Flow and Sound as Reference
DistractedClams$Treatment = factor(DistractedClams$Treatment,levels=c("Flow_and_Sound", "Control","Flow","Sound"))
m2.4 <- lmer(Latency_to_close~ Treatment + Trial_number + Plume_presence + Size + (1| Subject_ID), data = DistractedClams)
summary(m2.4)

###Histogram, QQ-plot, (f,r) plot###
r <- residuals(m2.4)
hist(r)
qqnorm(r)
f <- fitted(m2.4)
plot(f,r)
#all normal


###Testing if exemplars had a significant effect###
DistractedClams$Treatment = factor(DistractedClams$Treatment,levels=c("Control","Flow","Sound", "Flow_and_Sound"))
m2.5 <- lm(Latency_to_close~ Sound_exemplars + Trial_number + Size, data = Exemplartest)
summary(m2.5)
anova(m2.5)
#no significant effect of exemplars


####ANOVA####
###Compare linear mixed model with linear model without random effect of individual.
###Test if random effect of individual is significant
m2.6 <- lm(Latency_to_close~ Treatment + Trial_number + Plume_presence + Size, data = DistractedClams)
anova(m2, m2.6)
#significant (p < 0.001) p = 0.0009221, Chisq = 10.978, df = 1


####Emmeans###
###Estimated marginal means for each treatment###
closeemmeans <- emmeans(m2, "Treatment")
closeemmeans
plot(emmeans(m2, "Treatment"), comparisons = TRUE)
eff_size(closeemmeans, sigma = sigma(m2), edf = 95.8) #Cohen's d effect size

###Graphing Emmeans###
#Latency to Close Emmeans Bar Graph
se2 <- c("0.299", "0.282","0.299","0.291")
xform <- list(title= "Treatment",
              titlefont=list(size=20),
              tickfont=list(size=15),
              categoryorder = "array",
              categoryarray = c ("control", "flow only", "sound only", "flow and sound"))
c <- plot_ly(x = c("control", "flow only", "sound only", "flow and sound"),
             y = c(4.47,4.65,4.52,5.22),
             name = "Latency to Close",
             type = "bar",
             error_y = list(array = ~ se2,color = 'black'),
             marker=list (color='white',
                          line=list(color='black',
                                    width=1.5)))

c2 <- c %>% layout(
  xaxis= xform,
  yaxis=  list(title="Estimated Marginal Means \u00B1 SE",
               titlefont=list(size=20),
               tickvals = list(0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0),
               tickformat = (".1f"),
               linecolor = toRGB("black"),
               showgrid=FALSE,
               showline=TRUE))

c2


################################latency to reemerge analyses################################
####Linear Mixed Effects Models####
###Control for fixed effects of Treatment, Trial number, Plume Presence, and Size
###Control for random effect of individual (Subject ID)


###MODEL: Set Control as Reference
DistractedClams$Treatment = factor(DistractedClams$Treatment,levels=c("Control","Flow","Sound", "Flow_and_Sound"))
m3 <- lmer(log_ltr ~ Treatment + Trial_number + Plume_presence + Size + (1| Subject_ID), data = DistractedClams)
summary(m3)
#There is no significance between treatments and control. Flow is almost significant (p-value = 0.12741)
#Plume presence is significant (p = 0.00155)

###Histogram, QQ-plot, (f,r) plot###
r <- residuals(m3)
hist(r)
qqnorm(r)
f <- fitted(m3)
plot(f,r)
#all normal


#MODEL: Set Flow as Reference
DistractedClams$Treatment = factor(DistractedClams$Treatment,levels=c("Flow", "Sound", "Flow_and_Sound", "Control"))
m3.2 <- lmer(log_ltr~ Treatment + Trial_number + Plume_presence + Size + (1| Subject_ID), data = DistractedClams)
summary(m3.2)
#Flow is almost significant to sound (p = 0.07605) and is signifcant to flow and sound (p = 0.04510)

###Histogram, QQ-plot, (f,r) plot###
r <- residuals(m3.2)
hist(r)
qqnorm(r)
f <- fitted(m3.2)
plot(f,r)
#all normal

#MODEL: Set Sound as Reference
DistractedClams$Treatment = factor(DistractedClams$Treatment,levels=c("Sound", "Flow_and_Sound", "Control","Flow"))
m3.3 <- lmer(log_ltr~ Treatment + Trial_number + Plume_presence + Size + (1| Subject_ID), data = DistractedClams)
summary(m3.3)
#Sound is almost significant to flow.

###Histogram, QQ-plot, (f,r) plot###
r <- residuals(m3.3)
hist(r)
qqnorm(r)
f <- fitted(m3.3)
plot(f,r)
#all normal


#MODEL: Set Multimodal Sound and Flow as Reference
DistractedClams$Treatment = factor(DistractedClams$Treatment,levels=c("Flow_and_Sound", "Control","Flow","Sound"))
m3.4 <- lmer(log_ltr~ Treatment + Trial_number + Plume_presence + Size + (1| Subject_ID), data = DistractedClams)
summary(m3.4)
#Multimodal flow and sound is significant to flow (p = 0.04510)

###Histogram, QQ-plot, (f,r) plot###
r <- residuals(m3.4)
hist(r)
qqnorm(r)
f <- fitted(m3.4)
plot(f,r)
#all normal


###Testing if exemplars had a significant effect###
DistractedClams$Treatment = factor(DistractedClams$Treatment,levels=c("Control","Flow","Sound", "Flow_and_Sound"))
m3.5 <- lm(log_ltr~ Sound_exemplars + Trial_number + Size, data = Exemplartest)
anova(m3.5)
#exemplars are not significant


####ANOVA####
###Compare linear mixed model with linear model without random effect of individual.
###Test if random effect of individual is significant
m3.6 <- lm(log_ltr~ Treatment + Trial_number + Plume_presence + Size, data = DistractedClams)
anova(m3,m3.6)
#significant (p< 0.001) Chisquare = 22.175, df = 1


####Emmeans###
###Estimated marginal means for each treatment###
reemergeemmeans <- emmeans(m3, "Treatment")
reemergeemmeans
plot(emmeans(m3, "Treatment"), comparisons = TRUE)
eff_size(reemergeemmeans, sigma = sigma(m3), edf = 81.1) #Cohen's d effect size

###Graphing Emmeans###
#Latency to Reemerge Emmeans Bar Graph
se3 <- c("0.0704", "0.0627","0.0682","0.0644")
xform <- list(title= "Treatment",
              titlefont=list(size=20),
              tickfont=list(size=15),
              categoryorder = "array",
              categoryarray = c ("control", "flow only", "sound only", "flow and sound"))
re <- plot_ly(x = c("control", "flow only", "sound only", "flow and sound"),
              y = c(1.71,1.60,1.73,1.74),
              name = "Latency to Close",
              type = "bar",
              error_y = list(array = ~ se3,color = 'black'),
              marker=list (color='white',
                           line=list(color='black',
                                     width=1.5)))

re2 <- re %>% layout(
  xaxis= xform,
  yaxis=  list(title="Estimated Marginal Means \u00B1 SE",
               titlefont=list(size=20),
               tickvals = list(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0),
               tickformat = (".1f"),
               linecolor = toRGB("black"),
               showgrid=FALSE,
               showline=TRUE))

re2