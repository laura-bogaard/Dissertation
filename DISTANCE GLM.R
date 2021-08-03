## poisson log GLM with Distance

#load data
distdat <- read.csv("cleandist.csv") 

#load libraries
library(MuMIn)
library(tidyverse)

# select only data we need
distdat <- distdat%>%
  select(tide_dist, dur_hr, treatment, idv_rl, jul_day, start_hr, foraging, obs, tide_ht, 
         fish)
names(distdat)

#subset data

distance <- distdat$tide_dist #tide corrected distance from tast to seal
duration <- distdat$dur_hr #survey duration in hrs
treatment <- as.factor(distdat$treatment) #treatment as factor
rl <- distdat$idv_rl # rl predicted at each distance from emperical TL model 
day <- distdat$jul_day #julian day from deployment
hour <- distdat$start_hr # hour of day when survey began
forage <- as.factor(distdat$foraging) #binary dummy var 1 = foraging, 0 = not foraging
obs <- as.factor(distdat$obs) # session observer
tide <- distdat$tide_ht #height of the tide in m
fish <- distdat$fish  # net fish counts through the locks that day


# EDA
#qplot(data = distdat, x  = tide_dist, y = idv_rl) # decreasing logarhythmic relationship, gonna have to throw this one out dependence duh. 
qplot(data = distdat, x  = tide_dist, y = jul_day) #grouped by fish phases
qplot(data = distdat, x  = tide_dist, y = start_hr) # pretty regular
qplot(data = distdat, x  = tide_dist, y = tide_ht) # pretty random, slightly clumped closer when tide is high
qplot(data = distdat, x  = tide_dist, y = foraging) # more likely to forage within 100m
qplot(data = distdat, x  = tide_dist, y = treatment) # slightly skewed closer when OFF
qplot(data = distdat, x  = tide_dist, y = fish) #kinda random, follows expected detection f()

qplot(data = distdat, x  = tide_dist, binwidth = 25) # looks poisson dist

# mean ≠ var so would be over-dispersed ie quasipoisson 
mean(distance)
var(distance)

#start with full model
m1 <- glm(distance ~ treatment + day + hour + forage + obs+ tide + fish, na.action = "na.fail",
          family = "quasipoisson", data = distdat) 

m1a <- glm(distance ~ treatment + day + hour + forage + obs+ tide + fish, na.action = "na.fail",
          family = "poisson", data = distdat) 

chat <- deviance(m1a) / df.residual(m1a)

# hack dredge function by fitting a poisson model as m1a and then selecting QAIC
dredge(m1a, rank = "QAIC", chat = chat, evaluate = T, fixed = c("treatment"))

# best models is day, fish, forage +hour + tide + treat

# compare with backwards selection 

summary(m1)

plot(m1)
#kick out observer
m2 <- glm(distance ~ treatment + day + hour + forage + tide + fish, na.action = "na.fail",
                family = "quasipoisson", data = distdat) 
summary(m2)

# kick out fish
m3 <- glm(distance ~ treatment + day + hour + forage +tide, na.action = "na.fail",
          family = "quasipoisson", data = distdat) 
summary(m3)

#kickout hour
m4 <- glm(distance ~ treatment + day + forage + tide, na.action = "na.fail",
          family = "quasipoisson", data = distdat) 
summary(m4)
# 
# #kickout day
# m5 <- glm(distance ~ treatment + rl+ forage + fish, na.action = "na.fail",
#           family = "quasipoisson", data = distdat) 
# summary(m5)
# 
# #kickout fish
# m6 <- glm(distance ~ treatment + rl, na.action = "na.fail",
#           family = "quasipoisson", data = distdat) 
# summary(m6)

# final model matched with dredge and has treatment, day, forage, tide, as significant predictors of distance-- but RL is a function of distance... hm 

### Model diagnostics (for 1st choice model)

res<-resid(m4,type="response")# somewhat unsure whether response (default) or pearson should be used
fit<-fitted(m4)
acf(res)# baad
plot(fit, res)#residuals by fitted, some pattern, prob because of zeros
plot(treat,res)# res by predictor, spread uneven
plot(time,res)#residuals by predictor
plot(trial,res)#residuals by ramdom effect
plot(site,res)#residuals by random effect
hist(res,breaks=c(30))#histogram of residuals, looks ok(ish)
#qq plot for random effects
#histogram & qqnorm plots of residuals
qqnorm(res)
qqline(res)# ok

