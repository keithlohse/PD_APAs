# Title: Larger anticipatory postural adjustments relate to later and larger reactive steps in people with PD
# Script File by: Keith Lohse
# Date: 2018-05-01

## Opening packages ------------------------------------------------------------
library("ggplot2"); library("lme4"); library("car"); library("dplyr"); 
library("lmerTest"); library("influence.ME"); library("boot")

##----------------------- Data Cleaning and QA ---------------------------------
## Setting the Directory -------------------------------------------------------
# Make sure the data files from GitHub are saved locally and set the working
# directory to that location. 
getwd()
setwd("C:/Users/yourname/Documents/GitHub/PD_APAs")
list.files()

# Participant #15 was excluded based on influence statistics. 
# To recreate the final analyses, import data_PETERSON_MASTER_no15
# To load the full dataset without exclusions, import data_PETERSON_MASTER_05012018.csv
DATA<-read.csv("./data_PETERSON_MASTER_no15.csv", header = TRUE, sep=",",  
               na.strings=c("NA","NaN"," ",""))


head(DATA)
DATA$subject<-factor(DATA$subject) # Make sure that SubID is treated as a factor
DATA$FOG.c<-(DATA$FOGstatus*2)-1 # FOGstatus was dummy coded, so I created a 
# coded version of this variable.
# We can also create a categorical label for this variable to be used in graphs
DATA$group=factor(DATA$FOGstatus)
levels(DATA$group) <- list("PD+FOG"="1", "PD-FOG"="0")
summary(DATA$group)



## Removing Outlying APAs ------------------------------------------------------
# Note that these exclusions were based on the full dataset that included Participant 15.
summary(DATA$APA) # mean APA = 5.22
sd(DATA$APA, na.rm=TRUE) # SD APA = 8.36
# Greater than 2SD = 5.22+2*(8.36) = 21.94 ~ 22
# Therefore, we will cut off everything greater than 22
plot(density(DATA$APA, na.rm=TRUE)); abline(col="red", lwd=2, v=22)

DAT2<-subset(DATA, APA <= 22)
# This leaves 1293/1400 rows, but 55 of those were missing values,
# So really its 1293/1345,  really we've excluded ~4% of our data.
summary(DAT2$APA)
#-------------------------------------------------------------------------------

## Removing Outlying COP_pos ---------------------------------------------------
summary(DATA$COP_pos) # mean COP_pos = 3.12
sd(DATA$COP_pos, na.rm=TRUE) # SD APA = 12.60
# UL for COP_pos is thus Mean + 2SD = 3.12+2*(12.60) = 28.32
# LL for COP_pos is thus Mean - 2SD = 3.12-2*(12.60) = -22.08
# Therefore, we will cut off everything greater than 22
plot(density(DATA$COP_pos, na.rm=TRUE)); abline(col="red", lwd=2, v=c(-22,28))

DAT3<-subset(DAT2, COP_pos >= -22) # removes lower than LL
DAT3<-subset(DAT3, COP_pos <= 28) # removes greater than UL
# This leaves 1246/1400 rows, but 55 of those were missing values,
# So really its 1246/1345,  really we've excluded ~7% of our data.
summary(DAT3$COP_pos)
summary(DAT3$APA)
# ------------------------------------------------------------------------------


## More Data organization-----------------------------------------------------------------------

# Next we can break the data into forward and backward steps for analysis
FORWARD<-subset(DAT3, direction == "forward")
BACKWARD<-subset(DAT3, direction =="backward")

# Note that we have two different methods for identifying step number in the 
# data. 
summary(FORWARD$trial) # relative trial numbers for each step type (1-25)
summary(FORWARD$trial_num) # absolute trial numbers (1-50)
# For greatest accuracy in our models, we want to use the absolute trial numbers.

##  --------------------- Analysis of Foward Steps -----------------------------
## Plot Trial and APA for forward steps ----------------------------------------------
g1<-ggplot(FORWARD, aes(x = trial_num, y = APA)) +
  geom_point(aes(fill=as.factor(subject)), pch=21, size=2, stroke=1.25) +
  geom_line(aes(group=subject)) +
  stat_smooth(aes(col=subject), method="lm", lwd=1.5, se=FALSE) +
  facet_wrap(~FOGstatus)
g2<-g1+scale_x_continuous(name = "Trial Number") +
  scale_y_continuous(name = "APA") + labs(title="Forward Steps")
g3 <- g2 + theme_bw() + 
  theme(axis.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        plot.title=element_text(size=16,face="bold", hjust = 0.5)) +
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position="none")

plot(g3)


## Forward step models ------------------------------------------------------------------
#First, we will create a mean centered version of the APA variable
FORWARD$APA.c<-FORWARD$APA-mean(FORWARD$APA)
summary(FORWARD$APA)
FORWARD$subject<-factor(FORWARD$subject)

#  Next, we will build a series of three models
# 1) Trial
# 2) Trial + APA.c
# 3) Trial + APA.c + FOG.c + APC.c*FOG.c

# Fwd Step length ----------------------------------------------------------------------------
# Model #1 Controlling for trial number
f01<-lmer(StepLength~
            # Fixed-effects
            1+
            trial_num+ # adding in interaction
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f01) # during forwad steping, there is an effect of trial # (i.e. time)

# Model #2 Adding the main effect of APA
f02<-lmer(StepLength~
            # Fixed-effects
            1+
            trial_num+APA.c+ 
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f02)

# Model #3 Adding in the Interaction FOG.c
f03<-lmer(StepLength~
            # Fixed-effects
            1+
            trial_num+APA.c*FOG.c+ 
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f03)

#In order to compare between models, we will use the ANOVA function
anova(f01, f02, f03)
# Model #2 is the best fitting model, but we will always present Model #3 as the 
# model of interest. 

# Checking Normality ----------------------------------------------------------
# Given the large sample size (in terms of data points) the shaprio-wilk test
# will likely be statistically significant for even small deviations from normality
# As such, I think it is best is we check normality visually.
# We can do this with the qqnorm plot:
qqnorm(resid(f03))
# ... and a plot of the density of residuals:
plot(density(resid(f03)))
# If these look like a line and a bell-curve, we're probably okay!

     # We can also check heteroscedasticity by plotting the residuals againts the 
# predicted values... ideally this looks like a cloud, without a clear pattern:
plot(x=fitted(f03), y=resid(f03))


# Checking for Influential Participants ---------------------------------------
# Create an object, "x", that stores the different influence statistics,
# this can take a while to run:
x<-influence(f03,"subject")

plot(x) # We can visually see if individual estimates for any of the effects 
# are wildly different from the others

cooks.distance(x) # we can calculte the influence of any given subject on the 
# full model by inspecting the cooks distances (ideally all below 1)


## Bootstrap CIs for model F03 -------------------------------------------------
b_par<-bootMer(x=f03,FUN=fixef, nsim=500, seed=1)
boot.ci(b_par,type="perc", index = 1) #Set index to 1-5 for different effects
# 1= INT; 2 = trial_num; 3 = APA.c; 4 = FOG.c; 5 = APAxFOG Interaction


## Figure 2A -------------------------------------------------------------------
head(FORWARD)
g1<-ggplot(FORWARD, aes(x = APA, y = StepLength)) +
  geom_point(aes(fill=subject), pch=21, size=2, stroke=1.25) +
  scale_fill_grey()+
  stat_smooth(aes(group=subject), col="grey", method="lm", lwd=1, se=FALSE) +
  facet_wrap(~group)
g2<-g1+scale_x_continuous(name = "APA") +
  scale_y_continuous(name = "Step Length (m)") + labs(title="Forward Stepping")
g3 <- g2 + theme_bw() + 
  theme(axis.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        plot.title=element_text(size=16,face="bold", hjust = 0.5)) +
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position="none")

plot(g3)

# Fixed Effect of APA by FOG, calculated from Model 3
group<-c("PD+FOG", "PD-FOG")
Intercepts<-c(0.205715, 0.218065)
Slopes<- c(0.000912, 0.002354)
dd<-data.frame(group,Intercepts,Slopes)  
dd
g4<- g3 + geom_abline(aes(intercept=Intercepts, slope=Slopes, lty=group), lwd=2, data=dd)

plot(g4)
# ------------------------------------------------------------------------------

## Fwd Step latency ------------------------------------------------------------
# Model #1 Controlling for trial number
f01<-lmer(Sla~
            # Fixed-effects
            1+
            trial_num+ # adding in interaction
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f01) # during forwad steping, there is an effect of trial # (i.e. time)

# Model #2 Adding the main effect of APA
f02<-lmer(Sla~
            # Fixed-effects
            1+
            trial_num+APA.c+ 
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f02)

# Model #3 Adding in the Interaction FOG.c
f03<-lmer(Sla~
            # Fixed-effects
            1+
            trial_num+APA.c*FOG.c+ 
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f03)

#In order to compare between models, we will use the ANOVA function
anova(f01, f02, f03)
# Model #2 is the best fitting model, but we will always present Model #3 as the 
# model of interest. 

# Checking Normality ----------------------------------------------------------
# Given the large sample size (in terms of data points) the shaprio-wilk test
# will likely be statistically significant for even small deviations from normality
# As such, I think it is best is we check normality visually.
# We can do this with the qqnorm plot:
qqnorm(resid(f03))
# ... and a plot of the density of residuals:
plot(density(resid(f03)))
# If these look like a line and a bell-curve, we're probably okay!

# We can also check heteroscedasticity by plotting the residuals againts the 
# predicted values... ideally this looks like a cloud, without a clear pattern:
plot(x=fitted(f03), y=resid(f03))


# Checking for Influential Participants ---------------------------------------
# Create an object, "x", that stores the different influence statistics,
# this can take a while to run:
x<-influence(f03,"subject")

plot(x) # We can visually see if individual estimates for any of the effects 
# are wildly different from the others

cooks.distance(x) # we can calculte the influence of any given subject on the 
# full model by inspecting the cooks distances (ideally all below 1)

## Bootstrap CIs for model F03 -------------------------------------------------
b_par<-bootMer(x=f03,FUN=fixef, nsim=500, seed=1)
boot.ci(b_par,type="perc", index = 1) #Set index to 1-5 for different effects
# 1= INT; 2 = trial_num; 3 = APA.c; 4 = FOG.c; 5 = APAxFOG Interaction

# Figure 2B --------------------------------------------------------------------
head(FORWARD)
g1<-ggplot(FORWARD, aes(x = APA, y = Sla)) +
  geom_point(aes(fill=subject), pch=21, size=2, stroke=1.25) +
  scale_fill_grey()+
  stat_smooth(aes(group=subject), col="grey", method="lm", lwd=1, se=FALSE) +
  facet_wrap(~group)
g2<-g1+scale_x_continuous(name = "APA") +
  scale_y_continuous(name = "Step Latency (s)") + labs(title="Forward Stepping")
g3 <- g2 + theme_bw() + 
  theme(axis.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        plot.title=element_text(size=16,face="bold", hjust = 0.5)) +
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position="none")

plot(g3)

# Fixed Effect of APA by FOG, calculated from Model 3
group<-c("PD+FOG", "PD-FOG")
Intercepts<-c(0.259852, 0.255392)
Slopes<- c(0.002198, 0.005606)
dd<-data.frame(group,Intercepts,Slopes)  
dd
g4<- g3 + geom_abline(aes(intercept=Intercepts, slope=Slopes, lty=group), lwd=2, data=dd)

plot(g4)
#-------------------------------------------------------------------------------


# Figure 4A --------------------------------------------------------------------
# Step Latency by Trial
g1<-ggplot(FORWARD, aes(x = trial_num, y = Sla)) +
  geom_point(aes(fill=subject), pch=21, size=2, stroke=1.25) +
  scale_fill_grey()+
  stat_smooth(aes(group=subject), col="grey", method="lm", lwd=1, se=FALSE) +
  facet_wrap(~group)
g2<-g1+scale_x_continuous(name = "Trial Number") +
  scale_y_continuous(name = "Step Latency (s)") + labs(title="Forward Stepping")
g3 <- g2 + theme_bw() + 
  theme(axis.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        plot.title=element_text(size=16,face="bold", hjust = 0.5)) +
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position="none")

plot(g3)

# Fixed Effect of APA by FOG, calculated from Model 3
group<-c("PD+FOG", "PD-FOG")
Intercepts<-c(0.259852, 0.255392)
Slopes<- c(0.0004769, 0.0004769)
dd<-data.frame(group,Intercepts,Slopes)  
dd
g4<- g3 + geom_abline(aes(intercept=Intercepts, slope=Slopes, lty=group), lwd=2, data=dd)

plot(g4)
# ------------------------------------------------------------------------------





# Fwd ML_MOS ----------------------------------------------------------------------------
# Model #1 Controlling for trial number
f01<-lmer(ML_MOS~
            # Fixed-effects
            1+
            trial_num+ # adding in interaction
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f01) # during forwad steping, there is an effect of trial # (i.e. time)

# Model #2 Adding the main effect of APA
f02<-lmer(ML_MOS~
            # Fixed-effects
            1+
            trial_num+APA.c+ 
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f02)

# Model #3 Adding in the Interaction FOG.c
f03<-lmer(ML_MOS~
            # Fixed-effects
            1+
            trial_num+APA.c*FOG.c+ 
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f03)

#In order to compare between models, we will use the ANOVA function
anova(f01, f02, f03)
# Model #2 is the best fitting model, but we will always present Model #3 as the 
# model of interest. 

# Checking Normality ----------------------------------------------------------
# Given the large sample size (in terms of data points) the shaprio-wilk test
# will likely be statistically significant for even small deviations from normality
# As such, I think it is best is we check normality visually.
# We can do this with the qqnorm plot:
qqnorm(resid(f03))
# ... and a plot of the density of residuals:
plot(density(resid(f03)))
# If these look like a line and a bell-curve, we're probably okay!

# We can also check heteroscedasticity by plotting the residuals againts the 
# predicted values... ideally this looks like a cloud, without a clear pattern:
plot(x=fitted(f03), y=resid(f03))


# Checking for Influential Participants ---------------------------------------
# Create an object, "x", that stores the different influence statistics,
# this can take a while to run:
x<-influence(f03,"subject")

plot(x) # We can visually see if individual estimates for any of the effects 
# are wildly different from the others

cooks.distance(x) # we can calculte the influence of any given subject on the 
# full model by inspecting the cooks distances (ideally all below 1)


## Bootstrap CIs for model F03 -------------------------------------------------
b_par<-bootMer(x=f03,FUN=fixef, nsim=500, seed=1)
boot.ci(b_par,type="perc", index = 1) #Set index to 1-5 for different effects
# 1= INT; 2 = trial_num; 3 = APA.c; 4 = FOG.c; 5 = APAxFOG Interaction

# Figure 2C --------------------------------------------------------------------
head(FORWARD)
g1<-ggplot(FORWARD, aes(x = APA, y = ML_MOS)) +
  geom_point(aes(fill=subject), pch=21, size=2, stroke=1.25) +
  scale_fill_grey()+
  stat_smooth(aes(group=subject), col="grey", method="lm", lwd=1, se=FALSE) +
  facet_wrap(~group)
g2<-g1+scale_x_continuous(name = "APA") +
  scale_y_continuous(name = "Margin of Stability - ML (m)") + labs(title="Forward Stepping")
g3 <- g2 + theme_bw() + 
  theme(axis.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        plot.title=element_text(size=16,face="bold", hjust = 0.5)) +
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position="none")

plot(g3)

# Fixed Effect of APA by FOG, calculated from Model 3
group<-c("PD+FOG", "PD-FOG")
Intercepts<-c(0.01431, 0.018009)
Slopes<- c(0.00014, -0.00021)
dd<-data.frame(group,Intercepts,Slopes)  
dd
g4<- g3 + geom_abline(aes(intercept=Intercepts, slope=Slopes, lty=group), lwd=2, data=dd)

plot(g4)
# ------------------------------------------------------------------------------

# Fwd AP_MOS -------------------------------------------------------------------
# Model #1 Controlling for trial number
f01<-lmer(AP_MOS~
            # Fixed-effects
            1+
            trial_num+ # adding in interaction
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f01) # during forwad steping, there is an effect of trial # (i.e. time)

# Model #2 Adding the main effect of APA
f02<-lmer(AP_MOS~
            # Fixed-effects
            1+
            trial_num+APA.c+ 
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f02)

# Model #3 Adding in the Interaction FOG.c
f03<-lmer(AP_MOS~
            # Fixed-effects
            1+
            trial_num+APA.c*FOG.c+ 
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f03)

#In order to compare between models, we will use the ANOVA function
anova(f01, f02, f03)
# Model #2 is the best fitting model, but we will always present Model #3 as the 
# model of interest. 

# Checking Normality ----------------------------------------------------------
# Given the large sample size (in terms of data points) the shaprio-wilk test
# will likely be statistically significant for even small deviations from normality
# As such, I think it is best is we check normality visually.
# We can do this with the qqnorm plot:
qqnorm(resid(f03))
# ... and a plot of the density of residuals:
plot(density(resid(f03)))
# If these look like a line and a bell-curve, we're probably okay!

# We can also check heteroscedasticity by plotting the residuals againts the 
# predicted values... ideally this looks like a cloud, without a clear pattern:
plot(x=fitted(f03), y=resid(f03))


# Checking for Influential Participants ---------------------------------------
# Create an object, "x", that stores the different influence statistics,
# this can take a while to run:
x<-influence(f03,"subject")

plot(x) # We can visually see if individual estimates for any of the effects 
# are wildly different from the others

cooks.distance(x) # we can calculte the influence of any given subject on the 
# full model by inspecting the cooks distances (ideally all below 1)


## Bootstrap CIs for model F03 -------------------------------------------------
b_par<-bootMer(x=f03,FUN=fixef, nsim=500, seed=1)
boot.ci(b_par,type="perc", index = 1) #Set index to 1-5 for different effects
# 1= INT; 2 = trial_num; 3 = APA.c; 4 = FOG.c; 5 = APAxFOG Interaction


# Figure 2D --------------------------------------------------------------------
head(FORWARD)
g1<-ggplot(FORWARD, aes(x = APA, y = AP_MOS)) +
  geom_point(aes(fill=subject), pch=21, size=2, stroke=1.25) +
  scale_fill_grey()+
  stat_smooth(aes(group=subject), col="grey", method="lm", lwd=1, se=FALSE) +
  facet_wrap(~group)
g2<-g1+scale_x_continuous(name = "APA") +
  scale_y_continuous(name = "Margin of Stability - AP (m)") + labs(title="Forward Stepping")
g3 <- g2 + theme_bw() + 
  theme(axis.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        plot.title=element_text(size=16,face="bold", hjust = 0.5)) +
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position="none")

plot(g3)

# Fixed Effect of APA by FOG, calculated from Model 3
group<-c("PD+FOG", "PD-FOG")
Intercepts<-c(0.075872, 0.06788)
Slopes<- c(-0.00149, 0.00117)
dd<-data.frame(group,Intercepts,Slopes)  
dd
g4<- g3 + geom_abline(aes(intercept=Intercepts, slope=Slopes, lty=group), lwd=2, data=dd)

plot(g4)
# ------------------------------------------------------------------------------


# Fwd Step Width ---------------------------------------------------------------
# Model #1 Controlling for trial number
f01<-lmer(stepWidth~
            # Fixed-effects
            1+
            trial_num+ # adding in interaction
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f01) # during forwad steping, there is an effect of trial # (i.e. time)

# Model #2 Adding the main effect of APA
f02<-lmer(stepWidth~
            # Fixed-effects
            1+
            trial_num+APA.c+ 
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f02)

# Model #3 Adding in the Interaction FOG.c
f03<-lmer(stepWidth~
            # Fixed-effects
            1+
            trial_num+APA.c*FOG.c+ 
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f03)

#In order to compare between models, we will use the ANOVA function
anova(f01, f02, f03)
# Model #2 is the best fitting model, but we will always present Model #3 as the 
# model of interest. 

# Checking Normality ----------------------------------------------------------
# Given the large sample size (in terms of data points) the shaprio-wilk test
# will likely be statistically significant for even small deviations from normality
# As such, I think it is best is we check normality visually.
# We can do this with the qqnorm plot:
qqnorm(resid(f03))
# ... and a plot of the density of residuals:
plot(density(resid(f03)))
# If these look like a line and a bell-curve, we're probably okay!

# We can also check heteroscedasticity by plotting the residuals againts the 
# predicted values... ideally this looks like a cloud, without a clear pattern:
plot(x=fitted(f03), y=resid(f03))


# Checking for Influential Participants ---------------------------------------
# Create an object, "x", that stores the different influence statistics,
# this can take a while to run:
x<-influence(f03,"subject")

plot(x) # We can visually see if individual estimates for any of the effects 
# are wildly different from the others

cooks.distance(x) # we can calculte the influence of any given subject on the 
# full model by inspecting the cooks distances (ideally all below 1)


## Bootstrap CIs for model F03 -------------------------------------------------
b_par<-bootMer(x=f03,FUN=fixef, nsim=500, seed=1)
boot.ci(b_par,type="perc", index = 1) #Set index to 1-5 for different effects
# 1= INT; 2 = trial_num; 3 = APA.c; 4 = FOG.c; 5 = APAxFOG Interaction


# Fwd AP_COMatFO ----------------------------------------------------------------------------
# Model #1 Controlling for trial number
f01<-lmer(AP_COMatFO~
            # Fixed-effects
            1+
            trial_num+ # adding in interaction
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f01) # during forwad steping, there is an effect of trial # (i.e. time)

# Model #2 Adding the main effect of APA
f02<-lmer(AP_COMatFO~
            # Fixed-effects
            1+
            trial_num+APA.c+ 
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f02)

# Model #3 Adding in the Interaction FOG.c
f03<-lmer(AP_COMatFO~
            # Fixed-effects
            1+
            trial_num+APA.c*FOG.c+ 
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f03)

#In order to compare between models, we will use the ANOVA function
anova(f01, f02, f03)
# Model #3 is the best fitting model, but we will always present Model #3 as the 
# model of interest. 

# Checking Normality ----------------------------------------------------------
# Given the large sample size (in terms of data points) the shaprio-wilk test
# will likely be statistically significant for even small deviations from normality
# As such, I think it is best is we check normality visually.
# We can do this with the qqnorm plot:
qqnorm(resid(f03))
# ... and a plot of the density of residuals:
plot(density(resid(f03)))
# If these look like a line and a bell-curve, we're probably okay!

# We can also check heteroscedasticity by plotting the residuals againts the 
# predicted values... ideally this looks like a cloud, without a clear pattern:
plot(x=fitted(f03), y=resid(f03))


# Checking for Influential Participants ---------------------------------------
# Create an object, "x", that stores the different influence statistics,
# this can take a while to run:
x<-influence(f03,"subject")

plot(x) # We can visually see if individual estimates for any of the effects 
# are wildly different from the others

cooks.distance(x) # we can calculte the influence of any given subject on the 
# full model by inspecting the cooks distances (ideally all below 1)


## Bootstrap CIs for model F03 -------------------------------------------------
b_par<-bootMer(x=f03,FUN=fixef, nsim=500, seed=1)
boot.ci(b_par,type="perc", index = 1) #Set index to 1-5 for different effects
# 1= INT; 2 = trial_num; 3 = APA.c; 4 = FOG.c; 5 = APAxFOG Interaction


# Figure 2F --------------------------------------------------------------------
head(FORWARD)
g1<-ggplot(FORWARD, aes(x = APA, y = AP_COMatFO)) +
  geom_point(aes(fill=subject), pch=21, size=2, stroke=1.25) +
  scale_fill_grey()+
  stat_smooth(aes(group=subject), col="grey", method="lm", lwd=1, se=FALSE) +
  facet_wrap(~group)
g2<-g1+scale_x_continuous(name = "APA") +
  scale_y_continuous(name = "COM at Foot Off - AP (m)") + labs(title="Forward Stepping")
g3 <- g2 + theme_bw() + 
  theme(axis.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        plot.title=element_text(size=16,face="bold", hjust = 0.5)) +
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position="none")

plot(g3)

# Fixed Effect of APA by FOG
summary(FORWARD$APA)
group<-c("PD+FOG", "PD-FOG")
Intercepts<-c(0.137392, 0.129328)
Slopes<- c(0.0000346, 0.000398)
dd<-data.frame(group,Intercepts,Slopes)  
dd
g4<- g3 + geom_abline(aes(intercept=Intercepts, slope=Slopes, lty=group), lwd=2, data=dd)

plot(g4)
# ------------------------------------------------------------------------------


# Fwd ML_COMatFO ---------------------------------------------------------------
# Model #1 Controlling for trial number
f01<-lmer(ML_COMatFO~
            # Fixed-effects
            1+
            trial_num+ # adding in interaction
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f01) # during forwad steping, there is an effect of trial # (i.e. time)

# Model #2 Adding the main effect of APA
f02<-lmer(ML_COMatFO~
            # Fixed-effects
            1+
            trial_num+APA.c+ 
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f02)

# Model #3 Adding in the Interaction FOG.c
f03<-lmer(ML_COMatFO~
            # Fixed-effects
            1+
            trial_num+APA.c*FOG.c+ 
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f03)

#In order to compare between models, we will use the ANOVA function
anova(f01, f02, f03)
# Model #2 is the best fitting model, but we will always present Model #3 as the 
# model of interest. 

# Checking Normality ----------------------------------------------------------
# Given the large sample size (in terms of data points) the shaprio-wilk test
# will likely be statistically significant for even small deviations from normality
# As such, I think it is best is we check normality visually.
# We can do this with the qqnorm plot:
qqnorm(resid(f03))
# ... and a plot of the density of residuals:
plot(density(resid(f03)))
# If these look like a line and a bell-curve, we're probably okay!

# We can also check heteroscedasticity by plotting the residuals againts the 
# predicted values... ideally this looks like a cloud, without a clear pattern:
plot(x=fitted(f03), y=resid(f03))


 # Checking for Influential Participants ---------------------------------------
# Create an object, "x", that stores the different influence statistics,
# this can take a while to run:
x<-influence(f03,"subject")

plot(x) # We can visually see if individual estimates for any of the effects 
# are wildly different from the others

cooks.distance(x) # we can calculte the influence of any given subject on the 
# full model by inspecting the cooks distances (ideally all below 1)


## Bootstrap CIs for model F03 -------------------------------------------------
b_par<-bootMer(x=f03,FUN=fixef, nsim=500, seed=1)
boot.ci(b_par,type="perc", index = 1) #Set index to 1-5 for different effects
# 1= INT; 2 = trial_num; 3 = APA.c; 4 = FOG.c; 5 = APAxFOG Interaction


# Figure 2E --------------------------------------------------------------------
head(FORWARD)
g1<-ggplot(FORWARD, aes(x = APA, y = ML_COMatFO)) +
  geom_point(aes(fill=subject), pch=21, size=2, stroke=1.25) +
  scale_fill_grey()+
  stat_smooth(aes(group=subject), col="grey", method="lm", lwd=1, se=FALSE) +
  facet_wrap(~group)
g2<-g1+scale_x_continuous(name = "APA") +
  scale_y_continuous(name = "COM at Foot Off - ML (m)") + labs(title="Forward Stepping")
g3 <- g2 + theme_bw() + 
  theme(axis.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        plot.title=element_text(size=16,face="bold", hjust = 0.5)) +
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position="none")

plot(g3)

# Fixed Effect of APA by FOG, calculated from Model 3
group<-c("PD+FOG", "PD-FOG")
Intercepts<-c(0.0000764, 0.001281)
Slopes<- c(-0.00018, -0.00012)
dd<-data.frame(group,Intercepts,Slopes)  
dd
g4<- g3 + geom_abline(aes(intercept=Intercepts, slope=Slopes, lty=group), lwd=2, data=dd)

plot(g4)
# ------------------------------------------------------------------------------

# Fwd COM displacement ---------------------------------------------------------
# Model #1 Controlling for trial number
f01<-lmer(COMdisp~
            # Fixed-effects
            1+
            trial_num+ # adding in interaction
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f01) # during forwad steping, there is an effect of trial # (i.e. time)

# Model #2 Adding the main effect of APA
f02<-lmer(COMdisp~
            # Fixed-effects
            1+
            trial_num+APA.c+ 
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f02)

# Model #3 Adding in the Interaction FOG.c
summary(FORWARD$COMdisp)
f03<-lmer(COMdisp~
            # Fixed-effects
            1+
            trial_num+APA.c*FOG.c+ 
            # Random-effects
            (1+trial_num|subject), data=FORWARD, REML=FALSE)
summary(f03)

#In order to compare between models, we will use the ANOVA function
anova(f01, f02, f03)
# Model #1 is the best fitting model, but we will always present Model #3 as the 
# model of interest. 

# Checking Normality ----------------------------------------------------------
# Given the large sample size (in terms of data points) the shaprio-wilk test
# will likely be statistically significant for even small deviations from normality
# As such, I think it is best is we check normality visually.
# We can do this with the qqnorm plot:
qqnorm(resid(f03))
# ... and a plot of the density of residuals:
plot(density(resid(f03)))
# If these look like a line and a bell-curve, we're probably okay!

# We can also check heteroscedasticity by plotting the residuals againts the 
# predicted values... ideally this looks like a cloud, without a clear pattern:
plot(x=fitted(f03), y=resid(f03))


# Checking for Influential Participants ---------------------------------------
# Create an object, "x", that stores the different influence statistics,
# this can take a while to run:
x<-influence(f03,"subject")

plot(x) # We can visually see if individual estimates for any of the effects 
# are wildly different from the others

cooks.distance(x) # we can calculte the influence of any given subject on the 
# full model by inspecting the cooks distances (ideally all below 1)


## Bootstrap CIs for model F03 -------------------------------------------------
b_par<-bootMer(x=f03,FUN=fixef, nsim=500, seed=1)
boot.ci(b_par,type="perc", index = 1) #Set index to 1-5 for different effects
# 1= INT; 2 = trial_num; 3 = APA.c; 4 = FOG.c; 5 = APAxFOG Interaction





##  --------------------- Analysis of BACKWARD Steps -----------------------------
## Plot Trial and APA for forward steps ----------------------------------------------
g1<-ggplot(BACKWARD, aes(x = trial_num, y = APA)) +
  geom_point(aes(fill=as.factor(subject)), pch=21, size=2, stroke=1.25) +
  geom_line(aes(group=subject)) +
  stat_smooth(aes(col=subject), method="lm", lwd=1.5, se=FALSE) +
  facet_wrap(~FOGstatus)
g2<-g1+scale_x_continuous(name = "Trial Number") +
  scale_y_continuous(name = "APA") + labs(title="Forward Steps")
g3 <- g2 + theme_bw() + 
  theme(axis.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        plot.title=element_text(size=16,face="bold", hjust = 0.5)) +
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position="none")

plot(g3)


## Backward step models ------------------------------------------------------------------
#First, we will create a mean centered version of the APA variable
BACKWARD$APA.c<-BACKWARD$APA-mean(BACKWARD$APA)
summary(BACKWARD$APA)
BACKWARD$subject<-factor(BACKWARD$subject)
summary(BACKWARD$subject)

#  Next, we will build a series of three models
# 1) Trial
# 2) Trial + APA.c
# 3) Trial + APA.c + FOG.c + APC.c*FOG.c

# Bckwrd Step length ----------------------------------------------------------------------------
# Model #1 Controlling for trial number
f01<-lmer(StepLength~
            # Fixed-effects
            1+
            trial_num+ # adding in interaction
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f01) # during forwad steping, there is an effect of trial # (i.e. time)

# Model #2 Adding the main effect of APA
f02<-lmer(StepLength~
            # Fixed-effects
            1+
            trial_num+APA.c+ 
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f02)

# Model #3 Adding in the Interaction FOG.c
summary(BACKWARD$StepLength)
f03<-lmer(StepLength~
            # Fixed-effects
            1+
            trial_num+APA.c*FOG.c+ 
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f03)

#In order to compare between models, we will use the ANOVA function
anova(f01, f02, f03)
# Model #1 is the best fitting model, but we will always present Model #3 as the 
# model of interest. 

# Checking Normality ----------------------------------------------------------
# Given the large sample size (in terms of data points) the shaprio-wilk test
# will likely be statistically significant for even small deviations from normality
# As such, I think it is best is we check normality visually.
# We can do this with the qqnorm plot:
qqnorm(resid(f03))
# ... and a plot of the density of residuals:
plot(density(resid(f03)))
# If these look like a line and a bell-curve, we're probably okay!

# We can also check heteroscedasticity by plotting the residuals againts the 
# predicted values... ideally this looks like a cloud, without a clear pattern:
plot(x=fitted(f03), y=resid(f03))


# Checking for Influential Participants ---------------------------------------
# Create an object, "x", that stores the different influence statistics,
# this can take a while to run:
x<-influence(f03,"subject")

plot(x) # We can visually see if individual estimates for any of the effects 
# are wildly different from the others

cooks.distance(x) # we can calculte the influence of any given subject on the 
# full model by inspecting the cooks distances (ideally all below 1)


## Bootstrap CIs for model F03 -------------------------------------------------
b_par<-bootMer(x=f03,FUN=fixef, nsim=500, seed=1)
boot.ci(b_par,type="perc", index = 1) #Set index to 1-5 for different effects
# 1= INT; 2 = trial_num; 3 = APA.c; 4 = FOG.c; 5 = APAxFOG Interaction


## Bckwrd Step latency ---------------------------------------------------------
# Model #1 Controlling for trial number
f01<-lmer(Sla~
            # Fixed-effects
            1+
            trial_num+ # adding in interaction
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f01) # during forwad steping, there is an effect of trial # (i.e. time)

# Model #2 Adding the main effect of APA
f02<-lmer(Sla~
            # Fixed-effects
            1+
            trial_num+APA.c+ 
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f02)

# Model #3 Adding in the Interaction FOG.c
f03<-lmer(Sla~
            # Fixed-effects
            1+
            trial_num+APA.c*FOG.c+ 
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f03)

#In order to compare between models, we will use the ANOVA function
anova(f01, f02, f03)
# Model #3 is the best fitting model, but we will always present Model #3 as the 
# model of interest. 

# Checking Normality ----------------------------------------------------------
# Given the large sample size (in terms of data points) the shaprio-wilk test
# will likely be statistically significant for even small deviations from normality
# As such, I think it is best is we check normality visually.
# We can do this with the qqnorm plot:
qqnorm(resid(f03))
# ... and a plot of the density of residuals:
plot(density(resid(f03)))
# If these look like a line and a bell-curve, we're probably okay!

# We can also check heteroscedasticity by plotting the residuals againts the 
# predicted values... ideally this looks like a cloud, without a clear pattern:
plot(x=fitted(f03), y=resid(f03))


# Checking for Influential Participants ---------------------------------------
# Create an object, "x", that stores the different influence statistics,
# this can take a while to run:
x<-influence(f03,"subject")

plot(x) # We can visually see if individual estimates for any of the effects 
# are wildly different from the others

cooks.distance(x) # we can calculte the influence of any given subject on the 
# full model by inspecting the cooks distances (ideally all below 1)

## Bootstrap CIs for model F03 -------------------------------------------------
b_par<-bootMer(x=f03,FUN=fixef, nsim=500, seed=1)
boot.ci(b_par,type="perc", index = 1) #Set index to 1-5 for different effects
# 1= INT; 2 = trial_num; 3 = APA.c; 4 = FOG.c; 5 = APAxFOG Interaction

# Figure 3A --------------------------------------------------------------------
head(BACKWARD)
g1<-ggplot(BACKWARD, aes(x = APA, y = Sla)) +
  geom_point(aes(fill=subject), pch=21, size=2, stroke=1.25) +
  scale_fill_grey()+
  stat_smooth(aes(group=subject), col="grey", method="lm", lwd=1, se=FALSE) +
  facet_wrap(~group)
g2<-g1+scale_x_continuous(name = "APA") +
  scale_y_continuous(name = "Step Latency (s)") + labs(title="Backward Stepping")
g3 <- g2 + theme_bw() + 
  theme(axis.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        plot.title=element_text(size=16,face="bold", hjust = 0.5)) +
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position="none")

plot(g3)

# Fixed Effect of APA by FOG
group<-c("PD+FOG", "PD-FOG")
Intercepts<-c(0.281225, 0.278637)
Slopes<- c(0.000821, 0.002875)
dd<-data.frame(group,Intercepts,Slopes)  
dd
g4<- g3 + geom_abline(aes(intercept=Intercepts, slope=Slopes, lty=group), lwd=2, data=dd)

plot(g4)
# ------------------------------------------------------------------------------

# Bckwrd ML_MOS ----------------------------------------------------------------
# Model #1 Controlling for trial number
f01<-lmer(ML_MOS~
            # Fixed-effects
            1+
            trial_num+ # adding in interaction
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f01) # during forwad steping, there is an effect of trial # (i.e. time)

# Model #2 Adding the main effect of APA
f02<-lmer(ML_MOS~
            # Fixed-effects
            1+
            trial_num+APA.c+ 
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f02)

# Model #3 Adding in the Interaction FOG.c
f03<-lmer(ML_MOS~
            # Fixed-effects
            1+
            trial_num+APA.c*FOG.c+ 
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f03)

#In order to compare between models, we will use the ANOVA function
anova(f01, f02, f03)
# Model #1 is the best fitting model, but we will always present Model #3 as the 
# model of interest. 

# Checking Normality ----------------------------------------------------------
# Given the large sample size (in terms of data points) the shaprio-wilk test
# will likely be statistically significant for even small deviations from normality
# As such, I think it is best is we check normality visually.
# We can do this with the qqnorm plot:
qqnorm(resid(f03))
# ... and a plot of the density of residuals:
plot(density(resid(f03)))
# If these look like a line and a bell-curve, we're probably okay!

# We can also check heteroscedasticity by plotting the residuals againts the 
# predicted values... ideally this looks like a cloud, without a clear pattern:
plot(x=fitted(f03), y=resid(f03))


# Checking for Influential Participants ---------------------------------------
# Create an object, "x", that stores the different influence statistics,
# this can take a while to run:
x<-influence(f03,"subject")

plot(x) # We can visually see if individual estimates for any of the effects 
# are wildly different from the others

cooks.distance(x) # we can calculte the influence of any given subject on the 
# full model by inspecting the cooks distances (ideally all below 1)


## Bootstrap CIs for model F03 -------------------------------------------------
b_par<-bootMer(x=f03,FUN=fixef, nsim=500, seed=1)
boot.ci(b_par,type="perc", index = 1) #Set index to 1-5 for different effects
# 1= INT; 2 = trial_num; 3 = APA.c; 4 = FOG.c; 5 = APAxFOG Interaction


# Figure 3B --------------------------------------------------------------------
g1<-ggplot(BACKWARD, aes(x = APA, y = ML_MOS)) +
  geom_point(aes(fill=subject), pch=21, size=2, stroke=1.25) +
  scale_fill_grey()+
  stat_smooth(aes(group=subject), col="grey", method="lm", lwd=1, se=FALSE) +
  facet_wrap(~group)
g2<-g1+scale_x_continuous(name = "APA") +
  scale_y_continuous(name = "Margin of Stability - ML (m)") + labs(title="Backward Stepping")
g3 <- g2 + theme_bw() + 
  theme(axis.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        plot.title=element_text(size=16,face="bold", hjust = 0.5)) +
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position="none")

plot(g3)

# Fixed Effect of APA by FOG, calculated from Model 3
group<-c("PD+FOG", "PD-FOG")
Intercepts<-c(0.03308, 0.01209)
Slopes<- c(0.000115, 0.000398)
dd<-data.frame(group,Intercepts,Slopes)  
dd
g4<- g3 + geom_abline(aes(intercept=Intercepts, slope=Slopes, lty=group), lwd=2, data=dd)

plot(g4)
# ------------------------------------------------------------------------------

# Figure 4B --------------------------------------------------------------------
g1<-ggplot(BACKWARD, aes(x = trial_num, y = ML_MOS)) +
  geom_point(aes(fill=subject), pch=21, size=2, stroke=1.25) +
  scale_fill_grey()+
  stat_smooth(aes(group=subject), col="grey", method="lm", lwd=1, se=FALSE) +
  facet_wrap(~group)
g2<-g1+scale_x_continuous(name = "Trial Number") +
  scale_y_continuous(name = "Margin of Stability - ML (m)") + labs(title="Backward Stepping")
g3 <- g2 + theme_bw() + 
  theme(axis.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        plot.title=element_text(size=16,face="bold", hjust = 0.5)) +
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position="none")

plot(g3)

# Fixed Effect of Trial
group<-c("PD+FOG", "PD-FOG")
Intercepts<-c(0.03308, 0.01209)
Slopes<- c(-0.0001924, -0.0001924)
dd<-data.frame(group,Intercepts,Slopes)  
dd
g4<- g3 + geom_abline(aes(intercept=Intercepts, slope=Slopes, lty=group), lwd=2, data=dd)

plot(g4)
# ------------------------------------------------------------------------------


# Bckwrd AP_MOS ----------------------------------------------------------------
# Model #1 Controlling for trial number
f01<-lmer(AP_MOS~
            # Fixed-effects
            1+
            trial_num+ # adding in interaction
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f01) # during forwad steping, there is an effect of trial # (i.e. time)

# Model #2 Adding the main effect of APA
f02<-lmer(AP_MOS~
            # Fixed-effects
            1+
            trial_num+APA.c+ 
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f02)

# Model #3 Adding in the Interaction FOG.c
f03<-lmer(AP_MOS~
            # Fixed-effects
            1+
            trial_num+APA.c*FOG.c+ 
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f03)

#In order to compare between models, we will use the ANOVA function
anova(f01, f02, f03)
# Model #2 is the best fitting model, but we will always present Model #3 as the 
# model of interest. 

# Checking Normality ----------------------------------------------------------
# Given the large sample size (in terms of data points) the shaprio-wilk test
# will likely be statistically significant for even small deviations from normality
# As such, I think it is best is we check normality visually.
# We can do this with the qqnorm plot:
qqnorm(resid(f03))
# ... and a plot of the density of residuals:
plot(density(resid(f03)))
# If these look like a line and a bell-curve, we're probably okay!

# We can also check heteroscedasticity by plotting the residuals againts the 
# predicted values... ideally this looks like a cloud, without a clear pattern:
plot(x=fitted(f03), y=resid(f03))


# Checking for Influential Participants ---------------------------------------
# Create an object, "x", that stores the different influence statistics,
# this can take a while to run:
x<-influence(f03,"subject")

plot(x) # We can visually see if individual estimates for any of the effects 
# are wildly different from the others

cooks.distance(x) # we can calculte the influence of any given subject on the 
# full model by inspecting the cooks distances (ideally all below 1)


## Bootstrap CIs for model F03 -------------------------------------------------
b_par<-bootMer(x=f03,FUN=fixef, nsim=500, seed=1)
boot.ci(b_par,type="perc", index = 1) #Set index to 1-5 for different effects
# 1= INT; 2 = trial_num; 3 = APA.c; 4 = FOG.c; 5 = APAxFOG Interaction


# Figure 3C --------------------------------------------------------------------
g1<-ggplot(BACKWARD, aes(x = APA, y = AP_MOS)) +
  geom_point(aes(fill=subject), pch=21, size=2, stroke=1.25) +
  scale_fill_grey()+
  stat_smooth(aes(group=subject), col="grey", method="lm", lwd=1, se=FALSE) +
  facet_wrap(~group)
g2<-g1+scale_x_continuous(name = "APA") +
  scale_y_continuous(name = "Margin of Stability - AP (m)") + labs(title="Backward Stepping")
g3 <- g2 + theme_bw() + 
  theme(axis.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        plot.title=element_text(size=16,face="bold", hjust = 0.5)) +
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position="none")

plot(g3)

# Fixed Effect of APA by FOG, calculated from Model 3
group<-c("PD+FOG", "PD-FOG")
Intercepts<-c(0.057266, 0.061485)
Slopes<- c(-0.00157, -0.00269)
dd<-data.frame(group,Intercepts,Slopes)  
dd
g4<- g3 + geom_abline(aes(intercept=Intercepts, slope=Slopes, lty=group), lwd=2, data=dd)

plot(g4)
# ------------------------------------------------------------------------------


# Figure 4C --------------------------------------------------------------------
g1<-ggplot(BACKWARD, aes(x = trial_num, y = AP_MOS)) +
  geom_point(aes(fill=subject), pch=21, size=2, stroke=1.25) +
  scale_fill_grey()+
  stat_smooth(aes(group=subject), col="grey", method="lm", lwd=1, se=FALSE) +
  facet_wrap(~group)
g2<-g1+scale_x_continuous(name = "Trial Number") +
  scale_y_continuous(name = "Margin of Stability - AP (m)") + labs(title="Backward Stepping")
g3 <- g2 + theme_bw() + 
  theme(axis.text=element_text(size=16, colour="black"), 
        axis.title=element_text(size=16,face="bold"),
        plot.title=element_text(size=16,face="bold", hjust = 0.5)) +
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position="none")

plot(g3)

# Fixed Effect of Trial
group<-c("PD+FOG", "PD-FOG")
Intercepts<-c(0.057266, 0.061485)
Slopes<- c(0.000424, 0.000424)
dd<-data.frame(group,Intercepts,Slopes)  
dd
g4<- g3 + geom_abline(aes(intercept=Intercepts, slope=Slopes, lty=group), lwd=2, data=dd)

plot(g4)
# ------------------------------------------------------------------------------

# Bckwrd Step Width ------------------------------------------------------------
# Model #1 Controlling for trial number
f01<-lmer(stepWidth~
            # Fixed-effects
            1+
            trial_num+ # adding in interaction
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f01) # during forwad steping, there is an effect of trial # (i.e. time)

# Model #2 Adding the main effect of APA
f02<-lmer(stepWidth~
            # Fixed-effects
            1+
            trial_num+APA.c+ 
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f02)

# Model #3 Adding in the Interaction FOG.c
summary(BACKWARD$stepWidth)
f03<-lmer(log(stepWidth)~
            # Fixed-effects
            1+
            trial_num+APA.c*FOG.c+ 
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f03)

#In order to compare between models, we will use the ANOVA function
anova(f01, f02, f03)
# Model #1 is the best fitting model, but we will always present Model #3 as the 
# model of interest. 

# Checking Normality ----------------------------------------------------------
# Given the large sample size (in terms of data points) the shaprio-wilk test
# will likely be statistically significant for even small deviations from normality
# As such, I think it is best is we check normality visually.
# We can do this with the qqnorm plot:
qqnorm(resid(f03))
# ... and a plot of the density of residuals:
plot(density(resid(f03)))
# If these look like a line and a bell-curve, we're probably okay!

# We can also check heteroscedasticity by plotting the residuals againts the 
# predicted values... ideally this looks like a cloud, without a clear pattern:
plot(x=fitted(f03), y=resid(f03))


# Checking for Influential Participants ---------------------------------------
# Create an object, "x", that stores the different influence statistics,
# this can take a while to run:
x<-influence(f03,"subject")

plot(x) # We can visually see if individual estimates for any of the effects 
# are wildly different from the others

cooks.distance(x) # we can calculte the influence of any given subject on the 
# full model by inspecting the cooks distances (ideally all below 1)


## Bootstrap CIs for model F03 -------------------------------------------------
b_par<-bootMer(x=f03,FUN=fixef, nsim=500, seed=1)
boot.ci(b_par,type="perc", index = 1) #Set index to 1-5 for different effects
# 1= INT; 2 = trial_num; 3 = APA.c; 4 = FOG.c; 5 = APAxFOG Interaction


# Bckwrd AP_COMatFO ------------------------------------------------------------
# Model #1 Controlling for trial number
f01<-lmer(AP_COMatFO~
            # Fixed-effects
            1+
            trial_num+ # adding in interaction
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f01) # during forwad steping, there is an effect of trial # (i.e. time)

# Model #2 Adding the main effect of APA
f02<-lmer(AP_COMatFO~
            # Fixed-effects
            1+
            trial_num+APA.c+ 
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f02)

# Model #3 Adding in the Interaction FOG.c
f03<-lmer(AP_COMatFO~
            # Fixed-effects
            1+
            trial_num+APA.c*FOG.c+ 
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f03)

#In order to compare between models, we will use the ANOVA function
anova(f01, f02, f03)
# Model #1 is the best fitting model, but we will always present Model #3 as the 
# model of interest. 

# Checking Normality ----------------------------------------------------------
# Given the large sample size (in terms of data points) the shaprio-wilk test
# will likely be statistically significant for even small deviations from normality
# As such, I think it is best is we check normality visually.
# We can do this with the qqnorm plot:
qqnorm(resid(f03))
# ... and a plot of the density of residuals:
plot(density(resid(f03)))
# If these look like a line and a bell-curve, we're probably okay!

# We can also check heteroscedasticity by plotting the residuals againts the 
# predicted values... ideally this looks like a cloud, without a clear pattern:
plot(x=fitted(f03), y=resid(f03))


# Checking for Influential Participants ---------------------------------------
# Create an object, "x", that stores the different influence statistics,
# this can take a while to run:
x<-influence(f03,"subject")

plot(x) # We can visually see if individual estimates for any of the effects 
# are wildly different from the others

cooks.distance(x) # we can calculte the influence of any given subject on the 
# full model by inspecting the cooks distances (ideally all below 1)


## Bootstrap CIs for model F03 -------------------------------------------------
b_par<-bootMer(x=f03,FUN=fixef, nsim=500, seed=1)
boot.ci(b_par,type="perc", index = 1) #Set index to 1-5 for different effects
# 1= INT; 2 = trial_num; 3 = APA.c; 4 = FOG.c; 5 = APAxFOG Interaction



# Bckwrd ML_COMatFO ------------------------------------------------------------
# Model #1 Controlling for trial number
f01<-lmer(ML_COMatFO~
            # Fixed-effects
            1+
            trial_num+ # adding in interaction
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f01) # during forwad steping, there is an effect of trial # (i.e. time)

# Model #2 Adding the main effect of APA
f02<-lmer(ML_COMatFO~
            # Fixed-effects
            1+
            trial_num+APA.c+ 
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f02)

# Model #3 Adding in the Interaction FOG.c
f03<-lmer(ML_COMatFO~
            # Fixed-effects
            1+
            trial_num+APA.c*FOG.c+ 
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f03)

#In order to compare between models, we will use the ANOVA function
anova(f01, f02, f03)
# Model #1 is the best fitting model, but we will always present Model #3 as the 
# model of interest. 

# Checking Normality ----------------------------------------------------------
# Given the large sample size (in terms of data points) the shaprio-wilk test
# will likely be statistically significant for even small deviations from normality
# As such, I think it is best is we check normality visually.
# We can do this with the qqnorm plot:
qqnorm(resid(f03))
# ... and a plot of the density of residuals:
plot(density(resid(f03)))
# If these look like a line and a bell-curve, we're probably okay!

# We can also check heteroscedasticity by plotting the residuals againts the 
# predicted values... ideally this looks like a cloud, without a clear pattern:
plot(x=fitted(f03), y=resid(f03))


# Checking for Influential Participants ---------------------------------------
# Create an object, "x", that stores the different influence statistics,
# this can take a while to run:
x<-influence(f03,"subject")

plot(x) # We can visually see if individual estimates for any of the effects 
# are wildly different from the others

cooks.distance(x) # we can calculte the influence of any given subject on the 
# full model by inspecting the cooks distances (ideally all below 1)


## Bootstrap CIs for model F03 -------------------------------------------------
b_par<-bootMer(x=f03,FUN=fixef, nsim=500, seed=1)
boot.ci(b_par,type="perc", index = 1) #Set index to 1-5 for different effects
# 1= INT; 2 = trial_num; 3 = APA.c; 4 = FOG.c; 5 = APAxFOG Interaction



# Bckwrd COM displacement ----------------------------------------------------------------------------
# Model #1 Controlling for trial number
f01<-lmer(COMdisp~
            # Fixed-effects
            1+
            trial_num+ # adding in interaction
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f01) # during forwad steping, there is an effect of trial # (i.e. time)

# Model #2 Adding the main effect of APA
f02<-lmer(COMdisp~
            # Fixed-effects
            1+
            trial_num+APA.c+ 
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f02)

# Model #3 Adding in the Interaction FOG.c
plot(density(log(BACKWARD$COMdisp), na.rm=TRUE))
f03<-lmer(COMdisp~
            # Fixed-effects
            1+
            trial_num+APA.c*FOG.c+ 
            # Random-effects
            (1+trial_num|subject), data=BACKWARD, REML=FALSE)
summary(f03)

#In order to compare between models, we will use the ANOVA function
anova(f01, f02, f03)
# Model #3 is the best fitting model, but we will always present Model #3 as the 
# model of interest. 

# Checking Normality ----------------------------------------------------------
# Given the large sample size (in terms of data points) the shaprio-wilk test
# will likely be statistically significant for even small deviations from normality
# As such, I think it is best is we check normality visually.
# We can do this with the qqnorm plot:
qqnorm(resid(f03))
# ... and a plot of the density of residuals:
plot(density(resid(f03)))
# If these look like a line and a bell-curve, we're probably okay!

# We can also check heteroscedasticity by plotting the residuals againts the 
# predicted values... ideally this looks like a cloud, without a clear pattern:
plot(x=fitted(f03), y=resid(f03))


# Checking for Influential Participants ---------------------------------------
# Create an object, "x", that stores the different influence statistics,
# this can take a while to run:
x<-influence(f03,"subject")

plot(x) # We can visually see if individual estimates for any of the effects 
# are wildly different from the others

cooks.distance(x) # we can calculte the influence of any given subject on the 
# full model by inspecting the cooks distances (ideally all below 1)

## Bootstrap CIs for model F03 -------------------------------------------------
b_par<-bootMer(x=f03,FUN=fixef, nsim=500, seed=1)
boot.ci(b_par,type="perc", index = 1) #Set index to 1-5 for different effects
# 1= INT; 2 = trial_num; 3 = APA.c; 4 = FOG.c; 5 = APAxFOG Interaction
