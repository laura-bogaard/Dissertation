## This script contains the code for modelling binomial GAMs for TAST data
## Before you run this script, run TAST EDA 2 script first for wrangled variables
# Laura Bogaard


###############################################################################
# PRESENCE ABSENCE
# BINOMIAL GAMs
# set seed for reproducibility
set.seed(12345)

str(final_data1)
#format data
final_data2 <- transform(final_data1, 
                          pa = as.factor(pa),
                          hour = as.integer(hour),
                          treatment = as.factor(treatment),
                          dist2sq = distance,
                          RL = as.numeric(RL_emperical),
                          fish_day = as.integer(fish_day),
                          obs = as.factor(obs),
                          x = as.numeric(x),
                          y = as.numeric(y),
                          day = Julian_day_from_start)

# DESCRIP: model pa with smooth of RL
m1 <- gam(pa ~ treatment + s(RL), select = T, data = final_data1, choose.k = T, 
          family = binomial, method = "REML")
summary(m1)
gam.check(m1)
plot(m1, all.terms = TRUE, trans = plogis, shift = coef(m1)[1], pages = 1)

# INTERP: all terms significant, auto selection of k was k = 10

# DESCRIP: model pa with smooth of RL + all variables
m1.5 <- gam(pa ~ treatment + hour + s(fish_day, k = 3) + s(RL, k = 3), select = T, 
            data = final_data2, choose.k = T, family = binomial, method = "REML")
summary(m1.5)
gam.check(m1.5)
plot.gam(m1, all.terms = TRUE, trans = plogis, shift = coef(m1)[1], pages = 1)
# INTERP: full model doesn't have treatment as sig, but has observer effects, 
# because I don't really care about this I'll take out obs and try again. V2: now 
# treat is v sig and we'll drop day. V3:



# DESCRIP:pa with smooth of xy, 
m2 <- gam(pa ~ treatment + hour + fish_day + s(x,y), 
          data = final_data2, 
          choose.k = T,
          family = binomial, method = "REML")
summary(m2)
gam.check(m2)
# INTERP: all terms highly sig, dropped obs and day covar due to insignificance, 
# note that some obs were sig and others werent-- what does this mean?

# pa w smooth of distance
m3 <- gam(pa ~ treatment + hour + day + obs + fish_day + s(distance), 
          data = final_data2, family = binomial)
summary(m3)
plot(m3)
# model pa as func of smooth int between x and y binomial
m20 <- gam(pa ~ s(x, y) + treatment + s(day, k= 10) + s(RL, k = 10) + obs + s(hour, k= 10) + fish_day, data = final_data2, 
           method = "REML", family = "binomial")

summary(m20)
gam.check(m20)

#plot all terms of the model on one page, with the intercept and uncertainty, on the probability scale
plot(m23, all.terms = TRUE, trans = plogis, shift = coef(m23)[1], seWithMean = T, pages = 1)

# make predictions on the probability scale, with uncertainty 
##** remember to measure uncertainty on link s(rather than response) cale BEFORE transforming to prob scale w plogis
##* explain predictions by terms, type = "terms" , shows the contribution of each term 
##* to each prediction
##* 
predict(m23, type = "link", se.fit = T) #

# model PA with smooth for temporal variables
m22 <- gam(pa ~ s(x, y, k=30) + treatment + day + s(RL, k=2) + obs + hour + fish_day, 
           data = final_data2, 
           method = "REML", family = "binomial")

summary(m22)
gam.check(m22)
plot(m22, all.terms = TRUE, trans = plogis, shift = coef(m20)[1], pages = 1)

AIC(m1, m2, m3, m21, m22)

###############################################################################
# negative binomial models

#  first change pa to numeric for NB models
final_data_nb <- transform(final_data2, 
                         pa = as.numeric(pa))

#DESCRIP:
m4 <- gam(pa ~ treatment + hour + day + obs + fish_day + RL, data = final_data_nb, method = "REML", family = nb)
summary(m4)
gam.check(m4)
#INTERP:

#DESCRIP:
m5 <- gam(pa ~ treatment + hour + day + obs + fish_day + s(x,y), 
          data = final_data_nb, family = nb)
summary(m5)
gam.check(m5)

#DESCRIP:
m6 <- gam(pa ~ + treatment + hour + day + obs + fish_day, 
          data = final_data_nb, family = nb)
summary(m6)
gam.check(m6)
#INTERP:

# model pa as a func of smooth x an dsmooth y,  NB
m21 <- gam(pa ~ s(x, y, k= 20) + treatment, data = final_data_nb, 
           method = "REML", family = "nb")
summary(m21)
gam.check(m21)
plot(m21, all.terms = TRUE, trans = plogis, shift = coef(m21)[1], pages = 1)


AIC(m4, m5, m6, m21)

## COUNTS PER SQUARE (group size)
# try adding spatial smoother, need to code in dummy var for x, and y? losing information here
# m7 <- gam(num_indv ~ treatment + TOD + DOY + observer + fish_day + RL, data = Df1, family = poisson, offset = Area)
# summary(m7)
# 
# m8 <- gam(num_indv ~ treatment + TOD + DOY + observer + fish_day + s(x,y), data = Df1, family = poisson, offset = Area)
# summary(m8)
# 
# m9<- gam(num_idv ~ s(platform_distance, platform_bearing, k = 100),
#               data = valid_bearings,
#               family = poisson, method = "REML")
# 
# m10 <- gam(num_idv ~ s(platform_distance, platform_bearing, k = 100),
#               data = valid_bearings,
#               family = poisson, method = "REML")

# ## PREDATION RATE
# m11 <- glm(pred_rate ~ treatment + TOD + DOY + observer + fish_day + s(x,y), data = final_data1, family = poisson) 

###############################################################################
# FORAGING

#format data

foragemod_dat <- transform(foragemod_dat,
                           session_id = as.character(session_id),
                           foraging = as.factor(foraging),
                           obs = obs.x,
                           fish_day = fish_day.x,
                           hour = as.integer(hour.x),
                           treatment = as.factor(treatment.x),
                           Y = abs(Y),
                           day = Julian_day_from_start,
                           RL = as.numeric(RL_emperical))

m12 <- gam(foraging ~ treatment + hour + obs + fish_day + s(X,Y), 
           data = foragemod_dat, family = binomial, method = "REML")
 summary(m12)

 m13 <- gam(foraging ~ treatment + hour + day + obs + fish_day + s(RL, k = 50), 
            data = foragemod_dat, family = binomial, method = "REML")
summary(m13)


##remember to change foraging to numeric for - binom models
m14 <- gam(as.numeric(foraging) ~ treatment + hour + day + obs + fish_day, 
           data = foragemod_dat, family = binom, method = "REML")
summary(m14)

m15 <- gam(as.numeric(foraging) ~ treatment + hour + day + obs + fish_day + RL, 
           data = foragemod_dat, family = nb, method = "REML")
summary(m15)

AIC(m12, m13, m14, m15)

################################################################## 
## DISTANCE MODEL
data_distmod <- transform(distance, 
                          platform_distance = as.integer(platform_distance),
                          treatment = as.factor(treatment),
                          hour = as.integer(hour),
                          day = Julian_day_from_start.x,
                          rl = RL_emperical,
                          fish_day = as.integer(fish_day),
                          obs = as.factor(obs),
                          foraging = as.factor(foraging))
## Full  distance model GAM
m16 <- gam(platform_distance ~ treatment + foraging, data = data_distmod, 
              family = quasipoisson)
summary(m16)
gam.check(m16)


#FUll distane model glm
m17 <- glm(platform_distance ~ treatment + day, data = data_distmod, 
           family = quasipoisson)
summary(m17)
plot(m17)

m17 <- glm(platform_distance ~ treatment + day, data = data_distmod, 
           family = quasipoisson)


## Full  distance model GAM
m16.5 <- gam(platform_distance ~ treatment + foraging + s(fish_day, k = 25), data = data_distmod, 
           family = nb)
summary(m16.5)
gam.check(m16.5)


#FUll distane model glm
m17.5 <- glm(platform_distance ~ treatment + hour + rl + obs + fish_day, data = data_distmod, 
           family = nb)
summary(m17)
plot(m17)

################################################################################


############################################################
## FISH COUNT
# 
# m18 <- glm(fish_day ~ treatment + DOY + session_id + foraging, family = quasipoisson, data = final_data1)
# m18
# m19 <- glm(fish_day ~ treatment + DOY, family = poisson, data = final_data1)
# m19

#### code for plotting heat maps adapted from pg 147 adv data notes 
############################################################
library(fields) # quilt.plot(.)
par(mfrow=c(2,2))
# Loop across all phases
for (phase in c("A", "B", "C"))
{
  bWant <- df$Phase %in% phase # a boolean specifying phase
  fields::quilt.plot(df$XPos[bWant],
                     df$YPos[bWant],
                     fitted(mdl)[bWant],
                     nrow=25, ncol=60,
                     zlim=range(fitted(mdl)),
                     main=paste0("Phase ", phase))
}
  
### ADAPTED code for SALSA 2D implementation from adv data notes pg 172 
############################################################

##
library(MRSea)

df <- final_data2
# Create the variable "response" needed by SALSA 2D
df$response <- as.numeric(df$pa)
# Create the spatial variables x.pos and y.pos needed by SALSA 2D
df$x.pos <- df$x
df$y.pos <- df$y
# create the successes and failures columns needed for binomial response
df$successes <- as.numeric(df$pa)
df$failures <- ifelse(df$pa == 1, 0, 1)

# Set initial model without spline-based terms 
initialModel <- glm(response ~ as.factor(treatment) + as.factor(obs) + fish_day + day + hour + RL,
                    family=nb,
                    data=df)
# Set SALSA arguments
factorList <- c("obs", "treatment")
varlist <- c("fish_day", "day", "hour", "RL")
salsa1DList <- list(fitnessMeasure="QAIC", minKnots_1d= c(4,4),
                    maxKnots_1d= c(6,6), startKnots_1d= c(3,3),
                    degree=c(2,2), maxIterations=100,
                    gaps=0)
# Run SALSA 1D
salsa1D <- MRSea::runSALSA1D(initialModel, salsa1DList, varlist,
                             factorList, varlist_cyclicSplines=NULL,
                             splineParams=NULL, datain=df,
                             suppress.printout=TRUE, removal=FALSE,
                             panelid=NULL)
## [1] "Phase will be fitted as a factor variable; there are non-zero counts for all levels"
## [1] "Month will be fitted as a factor variable; there are non-zero counts for all levels"
# Store "best" model
bestModel1D <- salsa1D$bestModel






                              