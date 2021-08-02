## poisson log GLM with Distance

#load data
distdat <- read.csv("cleandist.csv") 

#load libraries
library(MuMIn)
library(tidyverse)

# select only data we need

distdat <- distdat%>%
  select(tide_dist, dur_hr, treatment, rlmid_emp, jul_day, hour, foraging, obs, tide_ht, 
         fish)
names(distdat)


#start with full model
m1 <- glm(tide_dist ~ treatment + s(rlmid_emp) + s(jul_day) + hour 
          + foraging + tide_ht + fish, na.action = "na.fail",
          family = "poisson", offset = log(dur_hr), data = distdat)

chat <- deviance(m1) / df.residual(m1)
# hack dredge function by fitting a poisson model as m1 and then selecting QAIC
dredge(m1, rank = "QAIC", chat = chat, evaluate = T, fixed = c("treatment"))

# best models are tied: 
# 1: julday, obs, rlmemp, treat
# 2: hour, julday, rlmemp treat


# but what if i kickout observer from the beginning, then the other is still the top mod
# interesting if I try a gam with a smooth for rl and day... 


summary(m1)

plot(m1)
#kick out observer
m2 <- glm(tide_dist ~ treatment + rlmid_emp + jul_day + hour 
                + foraging + tide_ht + fish,family = "quasipoisson", offset = log(dur_hr), data = distdat)

summary(m2)
# kick out tide height
m3 <- glm(tide_dist ~ treatment + rlmid_emp + jul_day + hour 
          + foraging + fish, family = "quasipoisson", offset = log(dur_hr), data = distdat)

summary(m3)

#kickout foraging
m4 <- glm(tide_dist ~ treatment + rlmid_emp + jul_day + hour 
          + fish, family = "quasipoisson", offset = log(dur_hr), data = distdat)

summary(m4)

#kickout fish
m5 <- glm(tide_dist ~ treatment + rlmid_emp + jul_day + hour, offset = log(dur_hr), family = "quasipoisson", data = distdat)

summary(m5)

#kickout hour
m6 <- glm(tide_dist ~ treatment + rlmid_emp + jul_day, offset = log(dur_hr), family = "quasipoisson", data = distdat)

summary(m6)



