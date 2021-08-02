
##### Predation model ####
##########################

dat<-read.csv("Ballard_Predation.csv",header=T)
require(ggplot2)
headobs<-(dat$Observer1)# primary observer (random effect)
date<-(dat$Date) #date as factor 
treat<-as.factor(dat$Treatment)# treatment
pred<-(dat$Seal_Pred)#number of predation events (count data) 
effort<-(dat$dur_hr)# obsrevation time per session
julian<-(dat$Julian_day_from_start)#julian day from start of baseline period, coded as a single random efefcst this will ultimately also be categorical
seals<-(dat$Seal_Freq)
obsloc<-(dat$Observer1_Location)# number of seals present (count data)
obs<-(dat$Observer1)
#####Predation rates##### 

library(emmeans)
require(glmmTMB)
install.packages('TMB', type = 'source') 
install.packages("glmmTMB", type = "source")
library(glmmTMB)
library(ggeffects)

###Model selection process using information criteria (AIC)#####
###Step 1: error distribution and zero inflation term on fully opulated model
m1<-glmmTMB(pred~treat+(1|julian)+(1|obs)+(1|obsloc)+offset(log(effort)),ziformula=~0,family=poisson)
m2<-glmmTMB(pred~treat+(1|julian)+(1|obs)+offset(log(effort)),ziformula=~1,family=nbinom1)# does not converge
m3<-glmmTMB(pred~treat+(1|julian)+(1|obs)+offset(log(effort)),ziformula=~1,family=nbinom2)# does not converge
m4<-glmmTMB(pred~treat+(1|julian)+(1|obs)+offset(log(effort)),ziformula=~0,family=poisson)# lowest AIC
m5<-glmmTMB(pred~treat+(1|julian)+(1|obs)+offset(log(effort)),ziformula=~0,family=nbinom1)# does not converge
m6<-glmmTMB(pred~treat+(1|julian)+(1|obs)+offset(log(effort)),ziformula=~0,family=nbinom2)
summary(m1,m2,m3,m4,m5,m6)

### ramdom effects###
m1<-glmmTMB(pred~treat+(1|julian)+(1|obs)+offset(log(effort)),ziformula=~0,family=poisson)
m2<-glmmTMB(pred~treat+offset(log(effort)),ziformula=~0,family=poisson)
m3<-glmmTMB(pred~treat+(julian+0|treat)+offset(log(effort)),ziformula=~0,family=poisson)#lowest AIC, including observer IDdoes not seem to much improved model fit 
m4<-glmmTMB(pred~treat+(1|julian)+(1|obsloc)+offset(log(effort)),ziformula=~0,family=poisson)
m5<-glmmTMB(pred~treat+(1|obs)+offset(log(effort)),ziformula=~0,family=poisson)
m5<-glmmTMB(pred~treat+(1|obsloc)+offset(log(effort)),ziformula=~0,family=poisson)
AIC(m1,m2,m3,m4,m5,m6)
###Results####
summary(m3)
#final model uses Poisson error distribution and contains treatment as a fixed effect, julin day/date as a random effects (categorical) and observation as an offset (needs to be log(offset) because of log link))
#model with lowest AIC, inclusion of a simple random effect julian/date with an unstructured correlatoion structure is also neccessary address autocorrelation of residuals 
exp(confint(m1))# coefficients and 95% CI

### plot mean predation rates###
contr<-emmeans(
  m1,
  ~ treat,
  type = "response", offset=0)

plot(contr,xlab = "Predation rate (events/hour): model estimate and 95% confidence interval  ",ylab = "Treatment")

## Diagnostics
res<-resid(m1) #somewhat unsure whether response (default) or pearson should be used
fit<-fitted(m1)
acf(res)# some mild violation at lag 16 but mich improved compared to model with no correlation structure. No violations at short lags. 
plot(fit, res)#ok
plot(res~treat)# spread uneven
plot(res~julian)
hist(res,breaks=c(30))# quite skewed
qqnorm(res)
qqline(res)# hmm

### Model diagnostics using Dhaarma
require(DHARMa)
n_sim <- 150
simulationOutput <- simulateResiduals(fittedModel = m3, n = n_sim)
plot(simulationOutput)# looks ok

##### Distance model ####
##########################
dat<-read.csv("Ballard_distance.csv",header=T)
###Baseline period is included here (if it has distance measurements) as it provides useful information. 
# Note that it was only excluded from the fish passage data (GAMM) because I think that we cannot reasonably fit a polynom to a time series which is data deficient (i.e. only data for one tretament level is present).

sess<-dat$session_id 
date<-(dat$date) #date as factor
treat<-as.factor(dat$treatment)# treatment
noseals<-(dat$num_idv)# no of individual seals
obsloc<-(dat$observer_location)#observer location
species<-as.factor(dat$species)#species
bear<-dat$platform_bearing#bearing to seal
dist<-dat$platform_distance#distance to seal, all obs without distance measure excluded, distance zero was set to dist=0.1 (Gamma can only be used for positive non-zero values)
forag<-dat$foraging #foraging behaviour (re-cided from yes or no)
crash<-dat$crash
###Note: observer variable somehow not in the distance data sheet, so not used!###

###Model selection process using information criteria (AIC)#####
####Select optimal random effects combination #####
m1<-glmmTMB(dist~treat*species+(1|date/sess)+(1|obsloc),family=Gamma(link="log"))# lowest AIC
m2<-glmmTMB(dist~treat*species+(1|date/sess)+(1|forag),family=Gamma(link="log"))
m3<-glmmTMB(dist~treat*species+(1|obsloc),family=Gamma(link="log"))
m4<-glmmTMB(dist~treat*species+(1|forag),family=Gamma(link="log"))
m5<-glmmTMB(dist~treat*species+(1|date/sess)+(1|obsloc)+(1|forag),family=Gamma(link="log"))# AIC very similar to m1
m6<-glmmTMB(dist~treat*species+(1|date/sess),family=Gamma(link="log"))
AIC(m1,m2,m3,m4,m5,m6)#model sselection based on AIC

###select optimal fixed effecst combination while keeping previously selected random effect combination (see Zuur et al. 2009###
m1<-glmmTMB(dist~treat*species+(1|date/sess)+(1|obsloc),family=Gamma(link="log"))# lowest AIC
m2<-glmmTMB(dist~treat+species+(1|date/sess)+(1|obsloc),family=Gamma(link="log"))
m3<-glmmTMB(dist~treat+(1|date/sess)+(1|obsloc),family=Gamma(link="log"))
AIC(m1,m2,m3) #model sselection based on AIC
anova(m1,m2,m3)#alternaive model selection methods, just to check, interaction term really important

###Results####
#final model includes session nested within date as a random effect and observers location, the interaction term between treatment efects and species improves the model significantly
summary(m1)
exp(confint(m1))# coefficients and 95% CI

### plot mean distance estimate/predictiosn  with 95% confidence intervals (CI)###
contr<-emmeans(
  m1,
  ~ treat*species,
  type = "response")

plot(contr,xlab = "Mean distance: model prediction and 95% confidence interval  ",ylab = "Treatment")

### contrasts and comparisons between TAST on and off for both species.  
emmeans(m1,poly~treat|species, type="response")# contrasts for predictions and coeffs

## Diagnostics
res<-resid(m1) #somewhat unsure whether response (default) or pearson should be used
fit<-fitted(m1)
acf(res)# some mild violation but overall pretty ok
plot(fit, res)# mild clsuter and somewhat unequal spread across fitted values
plot(res~treat)# slooks good
plot(res~date)#looks good
hist(res,breaks=c(30))#bit skewed towards lower values
qqnorm(res)
qqline(res)

### Model diagnostics using Dhaarma
require(DHARMa)
n_sim <- 200
simulationOutput <- simulateResiduals(fittedModel = m1, n = n_sim)
plot(simulationOutput)# Dharma standardized residuals vs predicted values look  good
### some indication of violation of normality in QQ plot but Dharma tends to be a bit hypersensitve in flagging 'issues',so might ok  


