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
m1 <- glm(tide_dist ~ treatment + rlmid_emp + jul_day + hour 
          + foraging + obs + tide_ht + fish, na.action = "na.fail",
          family = "poisson", offset = log(dur_hr), data = distdat)
chat <- deviance(m1) / df.residual(m1)

dredge(m1, rank = "QAIC", chat = chat, evaluate = T, fixed = c("treatment"))
summary(m1)

plot(m1)
#kick out observer
m2 <- glm(tide_dist ~ treatment + rlmid_emp + jul_day + hour 
                + foraging + tide_ht + fish,family = "poisson", offset = log(dur_hr), data = distdat)

summary(m2)
# kick out tide height
m3 <- glm(tide_dist ~ treatment + rlmid_emp + jul_day + hour 
          + foraging + fish, family = "poisson", offset = log(dur_hr), data = distdat)

summary(m3)

#kickout foraging
m4 <- glm(tide_dist ~ treatment + rlmid_emp + jul_day + hour 
          + fish, family = "poisson", offset = log(dur_hr), data = distdat)

summary(m4)

#kickout fish
m5 <- glm(tide_dist ~ treatment + rlmid_emp + jul_day + hour, offset = log(dur_hr), family = "poisson", data = distdat)

summary(m5)

#kickout hour
m6 <- glm(tide_dist ~ treatment + rlmid_emp + jul_day, offset = log(dur_hr), family = "poisson", data = distdat)

summary(m6)





QAIC(m1, m2, m3, m4, m5, m6, chat = chat)
?QAIC

