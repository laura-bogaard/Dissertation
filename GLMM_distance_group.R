dat<-read.csv("All_Trials_reducedhalfhour.csv",header=T)

require(dplyr)
require(mgcv)
require(lme4)
require(ggplot2)

##subsets data
sound<-subset(dat, PB_C=='Playback')
control<-subset(dat, PB_C=='Control')
dist<-dat$Dist
time<-dat$TimeTrial
trial<-(dat$Trial_ID)
TB<-as.factor(dat$Time_bin)#15 min
treat<-as.factor(dat$PB_C)
site<-dat$Site
w<-dat$Wind
ss<-dat$SS

require(ggplot2)
##30min time bins, it seems that the column has been removed from the data sheet
pl2<-(aggregate.data.frame(dist,list(trial,treat,TB,ss,site,trial),FUN=(mean),data=dat))
datq<-pl2
ggplot(datq,aes(y=pl2$x,x=pl2$Group.3,fill=pl2$Group.2))+geom_boxplot(width=0.6)+xlab("Time (15mins)")+ylab("distance (m)")

#Re-assign variables
dist<-pl2$x
treat<-pl2$Group.2
time<-pl2$Group.3
ss<-pl2$Group.4
site<-pl2$Group.5
trial<-pl2$Group.6
group<-pl2$Group.7

require(glmmTMB)
###Model selection process using information criteria (AIC)#####
#### Steps 1: select optimal error distribution and whether zero-inflation is needed

m1<-glmmTMB(dist~treat*time+ar1(time+0|trial)+(1|ss)+(1|site)+(1|trial),ziformula=~0,family=Gamma(link="log"))
m2<-glmmTMB(dist~treat*time+ar1(time+0|trial)+(1|ss)+(1|site)+(1|trial),ziformula=~0,family=nbinom1)
m3<-glmmTMB(dist~treat*time+ar1(time+0|trial)+(1|ss)+(1|site)+(1|trial),ziformula=~0,family=nbinom2)#best
AIC(m1,m2,m3)

#### Step 2:select optimalcomination of fixed effects on the fully populated model
#NA ....because fixed effects are given by theoretical considerations, i.e showing time bin by treatment  


#### Step 3:select optimalcomination of random effects on the fully populated model
###This shouldn't be based purely based fit (AIC) but also on the 
### data type/distribution and model diagnostics (i.e.one might want to choose neg binom over a poisson if address the issues of patterning in the residuals).
m1<-glmmTMB(dist~treat*time+ar1(time+0|trial),ziformula=~0,family=nbinom2)#best
m2<-glmmTMB(dist~treat*time+(1|ss)+ar1(time+0|trial),ziformula=~0,family=nbinom2)
m3<-glmmTMB(dist~treat*time+(1|trial)+ar1(time+0|trial),ziformula=~0,family=nbinom2)
m4<-glmmTMB(dist~treat*time+(1|site)+ar1(time+0|trial),ziformula=~0,family=nbinom2)
m5<-glmmTMB(dist~treat*time+(1|ss)+(1|site)+ar1(time+0|trial),ziformula=~0,family=nbinom2)
m6<-glmmTMB(dist~treat*time+(1|ss)+(1|trial)+ar1(time+0|trial),ziformula=~0,family=nbinom2)
m7<-glmmTMB(dist~treat*time+(1|site)+(1|trial)+ar1(time+0|trial),ziformula=~0,family=nbinom2)
m8<-glmmTMB(dist~treat*time+(1|site)+(1|trial)+(1|ss)+ar1(time+0|trial),ziformula=~0,family=nbinom2)#best
AIC(m1,m2,m3,m4,m5,m6,m7,m8)

#calculate contrasts (pairwise comparisson)
#Use p-values from 'cont'. Use 'EstCi' for coefficients/estimates (which are strangely called 'ratio') and confidence intervals (called 'upper/lower CL')
# As confint allows the argument type="resposne these are already on the scale of the response variable (and don't need to be back-transformed). 
require(emmeans)
em<- emmeans(m1,poly~treat|time)
#em<- emmeans(m1,pairwise~treat|time)
EstCi<-confint(em,type="response")# calculate contrats with Tukey HSD p-value adjustment
plot(em,comparisons = TRUE)

### Model diagnostics (1st choice)###

res<-resid(m1,type="pearson")# somewhat unsure whether response (default) or pearson should be used
fit<-fitted(m1)
acf(res)# fairly good, looks really good with res of type="response"
plot(fit, res)#residuals by fitted, some pattern, prob because of zeros
plot(treat,res)# res by predictor, spread uneven
plot(time,res)#residuals by predictor
plot(trial,res)#residuals by ramdom effect
plot(site,res)#residuals by random effect
hist(res,breaks=c(30))#histogram of resisuals, looks ok(ish)
#qq plot for random effects
#histogram & qqnorm plots of residuals
qqnorm(res)
qqline(res)# ok

### Model diagnostics using Dhaarma
require(DHARMa)
res = simulateResiduals(m1)
plot(res, rank = T)# QQ plot and residuals vs fitted. There are issues with the latter relating to models of type glmmTMB)
n_sim <- 150
simulationOutput <- simulateResiduals(fittedModel = m1, n = n_sim)
plot(simulationOutput)

