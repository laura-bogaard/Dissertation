## poisson log GLM with Distance

#load data
distdat <- read.csv("cleandist.csv") 

#load libraries
library(MuMIn)
library(tidyverse)

# select only data we need
distdat <- distdat%>%
  select(tide_dist, blockid, session_id,dur_hr, treatment, idv_rl, jul_day, start_hr, foraging, obs, tide_ht, 
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

#########
# EDA
qplot(data = distdat, x  = tide_dist, y = idv_rl) # decreasing logarhythmic relationship,
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
qqline(res)# ok GLM MODELS
##############

#try to run GAMs but weird... cant get a good fit

m1 <- gam(tide_dist ~ treatment + s(jul_day, k = 2) + s(start_hr, k = 3) + foraging + obs+ s(tide_ht, k = 2) + s(fish, k = 3), na.action = "na.fail",
          family = "quasipoisson", data = distdat) 
# drop observer
m2 <- gam(tide_dist ~ treatment + s(jul_day, k = 2) + s(start_hr, k = 3) + foraging + s(tide_ht, k = 2) + s(fish, k = 3), na.action = "na.fail",
          family = "quasipoisson", data = distdat) 
summary(m2)
gam.check(m2)

m3 <- gam(tide_dist ~ treatment + s(jul_day, k = 2) + s(start_hr, k = 2) + foraging + s(tide_ht, k = 2) + s(fish, k = 2), na.action = "na.fail",
                family = "quasipoisson", data = distdat) 
summary(m3)

m4 <- gam(tide_dist ~ treatment + s(jul_day, k = 3) + s(start_hr, k = 3) + foraging + s(tide_ht, k = 3) + s(fish, k = 3), na.action = "na.fail",
                family = "quasipoisson", data = distdat) 
summary(m4)
gam.check(m4)

m5 <- gam(tide_dist ~ treatment + s(jul_day, k = 4) + s(start_hr, k = 4) + foraging + s(tide_ht, k = 4) + s(fish, k = 4), na.action = "na.fail",
          family = "quasipoisson", data = distdat) 
summary(m5)
gam.check(m5)

m6 <- gam(tide_dist ~ treatment + s(jul_day, k = 5) + s(start_hr, k = 5) + foraging + s(tide_ht, k = 5) + s(fish, k = 5), na.action = "na.fail",
          family = "quasipoisson", data = distdat) 

summary(m6)
gam.check(m6)

QAIC(m1, m2, m3, m4, m5, m6)

###############
############# Distance GAM with GEE and SALSA 1D
df <- distdat

df$foraging <- as.factor(df$foraging)

# make splineParams object
splineParams<-makesplineParams(data=df, varlist=c('jul_day'))
str(splineParams)
 
fullmod <- glm(tide_dist ~ treatment + bs(jul_day, knots = splineParams[[2]]$knots)+ start_hr + as.factor(foraging) + obs+ tide_ht + fish, 
         na.action = "na.fail", family = "quasipoisson", data = df) 
# check for collinearity
car::vif(fullmod) ## smoothed Day is SUPER high (18) and Observer (6) so defo collinearity  

# check for correlated residuals
lawstat::runs.test(residuals(fullmod, type = "pearson")) # returns a low p-val indicatinf adn issue with correlated resid
plotRunsProfile(fullmod, varlist = c("jul_day")) #definite correlation


# ACF plot with blocks as survey sesh
distdat$session_id <- factor(distdat$session_id)
runACF(distdat$session_id, fullmod, store = F) #pretty dang ac, use this as our cor block

#Acf plot with blocks as survey area blocks
distdat$blockid <- as.factor(distdat$blockid)
runACF(distdat$session_id, fullmod, store = F) # not super ac

# plot cumulative residuals for model with tide_ht as smooth term
plotCumRes(fullmod, varlist= c("tide_ht"), splineParams)

# try SALSA1D
library(devtools)
devtools::install_github(repo="lindesaysh/MRSea@v1.01")

# Create the variable "response" needed by SALSA
df$response <- df$tide_dist

# assign a numeric id to each session id
#FELIX HELP THIS WONT ALLOW ME TO RUN A GEE
seshlevels <-  levels(as.factor(df$session_id)) 
numericsesh <- data.frame(session_id = levels(as.factor(df$session_id)),
           block_as_session = as.factor(seq(1, length(seshlevels), 1)))

df <- left_join(df, numericsesh, by = "session_id" )
arrange(df, session_id)

df$foldid<-getCVids(data = df, folds = 5, block = "session_id")

# Set initial model without spline-based terms
initialModel <- glm(response ~ treatment, family= quasipoisson, data=df)

# Set SALSA arguments (no RL because ON and OFF)
factorList <- c("obs", "treatment", "foraging")
varList <- c("tide_ht", "fish", "jul_day", "start_hr")
salsa1DList <- list(fitnessMeasure="QAIC", minKnots_1d=rep(1, 4),
                    maxKnots_1d=rep(5, 4), startKnots_1d=rep(1, 4),
                    degree=rep(2, 4), maxIterations=100,gaps=rep(0, 4))

# Load library
library(MRSea) 

set.seed(53195)
# Run SALSA 1D
salsa1D <- MRSea::runSALSA1D(initialModel, salsa1DList, varList,
                           factorList, varlist_cyclicSplines=NULL,
                           splineParams=NULL, datain=df,
                           suppress.printout=TRUE, removal=FALSE,
                           panelid=NULL)

# Pick best model based on fitness QAIC
bestModel1D <- salsa1D$bestModel

# gamMRSea(formula = response ~ treatment + bs(tide_ht, knots = splineParams[[2]]$knots, degree = splineParams[[2]]$degree, Boundary.knots = splineParams[[2]]$bd) + 
#                                           bs(fish, knots = splineParams[[3]]$knots, degree = splineParams[[3]]$degree, Boundary.knots = splineParams[[3]]$bd) + 
#                                           bs(jul_day, knots = splineParams[[4]]$knots, degree = splineParams[[4]]$degree,  Boundary.knots = splineParams[[4]]$bd) + 
#                                           bs(start_hr, knots = splineParams[[5]]$knots,degree = splineParams[[5]]$degree, Boundary.knots = splineParams[[5]]$bd), 
#                                           family = quasipoisson(link = log), data = df, splineParams = splineParams)
# summary
summary(bestModel1D)
stats::anova(bestModel1D, test="F") # the model kept treatment, tide_ht, fish, day, and hour as significant predictors of distance
str(bestModel1D)
salsa1D$modelFits1D

#show patrial plots
MRSea::runPartialPlots(bestModel1D, varlist.in= varList,
                       showKnots=T, type="link", data=df)

# get CV score
cv1 <- getCV_CReSS(datain = df, baseModel = salsa1D$bestModel,
                   salsa1D$splineParams) ## this doesnt work...

########### Diagnostics

# now that we have out best model, we should re-assess runs test for best model
lawstat::runs.test(residuals(bestModel1D, type = "pearson"))
#There is significant positive residual correlation (p << 0:05 and test statistic is
                                                   #negative) so we rerun the model as a GEE.
car::vif(bestModel1D) #better! all 1 ish

runDiagnostics(bestModel1D) # obs vs fit = poor fit, Marginal R-squared and concordance correlation are low
# no violation of mean variance relationship, redis vs fit == patternless, good

## Need to put into a GEE framework
# Re-fit the chosen model as a GEE (based on SALSA knot placement) and
# GEE p-values
install.packages('geepack')

geeModel <- geepack::geeglm(formula(bestModel1D), data = df, family = poisson,
                   id = block_as_session)

summary(geeModel)
class(geeModel)
str(geeModel)

# table of p-values
# specifying varlist and factorlist makes shorter variable names
getPvalues(geeModel, varlist = c("tide_ht", "fish","jul_day", "start_hr"), factorlist = c("treatment"))
dists <- splineParams[[1]]$dist

plotCumRes(geeModel, varlist=c("tide_ht", "fish","jul_day", "start_hr"), splineParams=splineParams, d2k=dists)
plotRunsProfile(geeModel, varlist=c("tide_ht", "fish","jul_day", "start_hr"))



####### Predictions
#subset data into OFF and ON 

distance_ON <- subset(df, df$treatment == "ON")
distance_OFF <- subset(df, df$treatment == "OFF")

# what does the model predect as the expected distance
mean(predict.gamMRSea(object = bestModel1D, newdata = distance_ON)) #99
mean(predict.gamMRSea(object = bestModel1D, newdata = distance_OFF)) #79

#but is it significant?





# # Select initial knot locations using a space-filling design
# knotGrid <- MRSea::getKnotgrid(coordData=cbind(df$x.pos, df$x.pos),
#                                numKnots=300, plot=FALSE)
# 
# # Create distance matrix
# # (i.e. distance between each data points and knots)
# distMatrix <- MRSea::makeDists(datacoords=cbind(df$x.pos, df$y.pos),
#                                knotcoords=knotGrid)
# # Set SALSA parameters # no interaction
# salsa2DList <- list(fitnessMeasure="AIC", knotgrid=knotGrid,
#                     startKnots=6, minKnots=2, maxKnots=20,
#                     gap=0)
# #
# set.seed(53195)
# # Run SALSA 2D
# salsa2D <- MRSea::runSALSA2D(bestModel1D, salsa2DList, d2k=distMatrix$dataDist,
#                              k2k=distMatrix$knotDist, splineParams=NULL,
#                              tol=0, chooserad=F, panels=NULL,
#                              suppress.printout=TRUE)
# 
# # Store "best" model
# bestModel2D <- salsa2D$bestModel
# # Summary of results
# summary(bestModel2D)
# 
# # Make partial plots?? not sure if this is useful for binary response
# #MRSea::runPartialPlots(bestModel2D, varlist.in=" ",
# #showKnots=TRUE, type="link", data=df)
# 
# ## Map model for On and Off
# par(mfrow=c(1,2))
# # Loop across all phases
# for (treat in c("ON", "OFF"))
# {
#   bWant <- df$treat %in% treat # a boolean specifying phase
#   fields::quilt.plot(df$XPos[bWant],
#                      df$YPos[bWant],
#                      fitted(bestModel2D)[bWant],
#                      
#                      nrow=10, ncol=10,
#                      zlim=range(fitted(bestModel2D)),
#                      main=paste0("TAST ", treat))
#   if (treat %in% "ON")
#   {
#     #somethign about turn off legend for one of the 
#   }
#   
# }
# 


