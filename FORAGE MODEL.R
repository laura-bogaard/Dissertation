
#load data
distdat <- read.csv("cleandist.csv") 

#load libraries
library(MuMIn)
library(tidyverse)
library(MRSea)
library(mgcv)
library(splines)

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

############# Distance GAM with GEE and SALSA 1D
df <- distdat

df$foraging <- as.numeric(df$foraging)

# make splineParams object
splineParams<- makesplineParams(data=df, varlist=c("tide_dist", "start_hr", "tide_ht", "fish", 'jul_day'))
str(splineParams)

# full mod w knot at each mean
fullmod <- glm(foraging ~ treatment + obs + 
                 bs(tide_dist, knots = splineParams[[2]]$knots) + 
                 bs(start_hr, knots = splineParams[[3]]$knots)+  
                 bs(tide_ht, knots = splineParams[[4]]$knots)+ 
                 bs(fish, knots = splineParams[[5]]$knots)+  
                 bs(jul_day, knots = splineParams[[6]]$knots),
                 na.action = "na.fail", family = "binomial", data = df) 

# check for collinearity
car::vif(fullmod) #none! 

# check for correlated residuals
lawstat::runs.test(residuals(fullmod, type = "pearson")) 
# returns a low p-val indicatinf adn issue with correlated resid -22, p-value < 2.2e-16
plotRunsProfile(fullmod, varlist = c("tide_dist", "start_hr", "tide_ht", "fish", 'jul_day'))
 #significant positive autocorrelation for all var

# ACF plot with blocks as survey sesh
df$session_id <- factor(distdat$session_id)
runACF(df$session_id, fullmod, store = F) #pretty dang ac, use this as our cor block

# plot cumulative residuals for model with tide_ht as smooth term
plotCumRes(fullmod, varlist = c("tide_dist", "start_hr", "tide_ht", "fish", 'jul_day'), splineParams)

# try SALSA1D
library(devtools)
devtools::install_github(repo="lindesaysh/MRSea@v1.01")

# Create the variable "response" needed by SALSA
df$response <- df$foraging

# assign a numeric id to each session id
# #FELIX HELP THIS WONT ALLOW ME TO RUN A GEE
# seshlevels <-  levels(as.factor(df$session_id)) 
# numericsesh <- data.frame(session_id = levels(as.factor(df$session_id)),
#                           block_as_session = as.factor(seq(1, length(seshlevels), 1)))
# 
# df <- left_join(df, numericsesh, by = "session_id" )
# arrange(df, session_id)
# 
# df$foldid<-getCVids(data = df, folds = 5, block = "session_id")

# Set initial model without spline-based terms
initialModel <- glm(response ~ treatment, family = binomial, data=df)

# Set SALSA arguments (no RL because ON and OFF)
factorList <- c("treatment", "obs")
varList = c("tide_dist", "start_hr", "tide_ht", "fish", 'jul_day')
salsa1DList <- list(fitnessMeasure="AIC", minKnots_1d=rep(1, 5),
                    maxKnots_1d=rep(5, 5), startKnots_1d=rep(1, 5),
                    degree=rep(2, 5), maxIterations=100, gaps=rep(1, 5))

# Load library
library(MRSea) 

set.seed(53195)
# Run SALSA 1D
salsa1D <- MRSea::runSALSA1D(initialModel, salsa1DList, varList,
                             factorList, varlist_cyclicSplines=NULL,
                             splineParams=NULL, datain=df,
                             suppress.printout=TRUE, removal=T,
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
stats::anova(bestModel1D, test="Chisq") # the model kept treatment, tide_ht, fish, day, and hour as significant predictors of distance
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
predslink <- predict.gamMRSea(object = bestModel1D, newdata = df)
preds <- plogis(predslink)

#but is it significant?
library(boot)
# make predictions for foraging probability when tasst on

# make predictions for forage pr when tast OFF

# do this 100 times



# 





boot(df, statistic = )






do.bootstrap.gam(df, df, ddf.obj=NULL, object = bestModel1D, B=20)

do.bootstrap.cress(preds, df, bestModel1D, B=50)



