### Binomial Presence Absence models

## read in data
padat <- read.csv("~/Desktop/TASTBALLARD/TAST WD/Dissertation/pa_data.csv")

## load libraries
library(tidyverse)
library(car)
library(lawstat)
library(devtools)
devtools::install_github(repo="lindesaysh/MRSea@v1.01")
library(MRSea)
library(mgcv)
library(splines)
library(geepack)

## format data
df <- padat 
# rename so its easier to type :)
df$pa <- as.numeric(df$pa)
df$treatment <- as.factor(df$treatment)
df$obs <- as.factor(df$obs)
df$x.pos <- df$x
df$y.pos <- df$y
df$start_hr <- as.factor(df$start_hr)
df$session_grid_block <- as.factor(paste(df$session_id, df$id, sep= ""))

## add numeric ID for each experimental unit (session_id) to cover temporal AC
## in GEE model (potentially)
# seshlevels <-  levels(as.factor(df$session_id)) 
# numericsesh <- data.frame(session_id = levels(as.factor(df$session_id)),
#                           block_as_session = as.factor(seq(1, length(seshlevels), 1)))
# 
# df <- left_join(df, numericsesh, by = "session_id" )
# arrange(df, session_id)

# # TPS Models
# ####### 
# #EDA
# qplot(df$treat, df$pa)
# 
# # start with full linear model
# m0 <- glm(pa ~ treatment + x.pos + y.pos + fish + obs + dist2sq + jul_day + start_hr, data = df, family = binomial)
# m0
# summary(m0)
# # check for collinearity
# vif(m0) # kick out dist to sq
# 
# m1 <- glm(pa ~ treatment + x.pos + y.pos + fish + obs + jul_day + start_hr, data = df, family = binomial)
# summary(m1)
# # check for collinearity
# vif(m1) #  better
# 
# # add spatial interaction to account for spatial autocorrelation
# m2<- gam(pa ~ treatment + s(x.pos, y.pos) + fish + obs + jul_day + start_hr, data = df, family = binomial)
# summary(m2)
# 
# ### try TPS models with and without interactions
# library(mgcv)
# TPS <- mgcv::gam(pa ~ s(fish, k = 2) + s(x.pos, y.pos) +
#                    s(start_hr, k= 2) + s(jul_day, k= 2) + treat + obs,
#                  family=binomial, data=df)
# summary(TPS)
# gam.check(TPS)
# 
# library(mgcv)
# TPS <- mgcv::gam(pa ~ s(fish, k = 2) + s(x.pos, y.pos) +
#                    s(start_hr, k= 2) + s(jul_day, k= 2) + treat + obs,
#                  family=binomial, data=df)
# summary(TPS)
# gam.check(TPS)
# # drop hr cause its linear
# TPS2 <- mgcv::gam(pa ~ s(fish, k = 2) + s(x.pos, y.pos) +
#                    s(start_hr) + s(jul_day, k= 2) + treat + obs,
#                  family=binomial, data=df)
# summary(TPS2)
# gam.check(TPS2)
# 
# # did not 
# 
# 
# 
# TPSInt <- mgcv::gam(pa ~ s(fish, k = 2) + s(x.pos, y.pos, by = treat) +
#                       s(start_hr, k= 2) + s(jul_day, k= 2) + treat +obs,
#                     family= binomial, data=df)
# gam.check(TPSInt)
# 
# summary(TPSInt)
# 
# 
# AIC(TPS, TPS2, TPSInt)
# 
# 
# ## STEP 3: check for correlated residuals in temporal covariates
# runs.test(residuals(m2, type = "pearson", alternative = "two-sided")) 
# 
# plotRunsProfile(m2, varlist = c("x.pos", "y.pos", "jul_day", "start_hr")) 
# #definite residual correlation
# 
# # ACF plot with blocks as survey sesh + gr id
# df$blockid <- as.factor(paste(df$block_as_session, df$id,sep = ""))
# #WARNING THIS RUNS FOREVER
# # runACF(df$blockid, m2, store = F) #very strange output just gives straight red line 
# 
# #Try grid id as a block
# df$id <- as.factor(df$id)
# runACF(df$id, m2, store = F) #interesting.. this is the one
# 
# #try survey sesh as block
# df$session_id <- as.factor(df$session_id) #
# runACF(df$session_id, m2, store = F) 
# # backwards stepwise selection


#########
# Fit GEE CReSS model with SALSA2D 

### STEP 1: check for collinearity between covariates in glm model
fullmod_linear <- glm(pa ~ treatment + jul_day + start_hr + obs + dist2sq + fish + x.pos + y.pos,
               na.action = "na.fail", family = "binomial", data = df)

## check for collinearity
car::vif(fullmod_linear) 
# most between 1 & 2.1 so pretty much good! dist2Sq and x/y pretty collin so throw out dist

## refit without distance 
fullmod_linear <- glm(pa ~ treatment + jul_day + start_hr + obs + fish + x.pos + y.pos + x.pos:y.pos,
               na.action = "na.fail", family = "binomial", data = df)

## recheck for collinearity
car::vif(fullmod_linear) 
# much better, fish and day only slightly collinear but we can keep it. (1-2)

### STEP 2: make splineParams object with knots at the mean of each variable
splineParams<-makesplineParams(data=df, varlist=c("x.pos", "y.pos", "fish","jul_day"))

## fit new full mod with splines binomial dist
fullmod <- glm(pa ~ treat + obs + start_hr + bs(x.pos, knots = splineParams[[2]]$knots) 
                          + bs(y.pos, knots = splineParams[[3]]$knots)
                          + bs(fish,knots = splineParams[[4]]$knots)
                          + bs(jul_day,knots = splineParams[[5]]$knots),
                          na.action = "na.fail", family = "binomial", data = df)

# ## check qbn distribution to see if dispersion para is 0
# fullmod_qbn <- glm(pa ~ treat + obs + bs(x.pos, knots = splineParams[[2]]$knots) 
#                + bs(y.pos, knots = splineParams[[3]]$knots)
#                + bs(fish,knots = splineParams[[4]]$knots)
#                + bs(jul_day,knots = splineParams[[5]]$knots)
#                + bs(start_hr, knots =splineParams[[6]]$knots),
#                na.action = "na.fail", family = "quasibinomial", data = df)
# 
# summary(fullmod_qbn) #nope dispparm = 1.14, not zero
# summary(fullmod)

### STEP 3: check for correlated residuals in descrete covariates using runs tests
varList <- c("y.pos", "y.pos", "fish", "jul_day")
runs.test(residuals(fullmod, type = "pearson", alternative = "two-sided")) 
# neg standardised runs stat and sm p == houston we have correlation
 par(ask = FALSE, mfrow = c(1,2))
runs_plots <- plotRunsProfile(fullmod, varlist = c("jul_day")) 

ps <- par(ask = FALSE, mfrow = c(1,2))
runs_plots <- plotRunsProfile(fullmod, varlist = c("y.pos")) 

ps <- par(ask = FALSE, mfrow = c(1,2))
runs_plots <- plotRunsProfile(fullmod, varlist = c("x.pos")) 

ps <- par(ask = FALSE, mfrow = c(1,2))
runs_plots <- plotRunsProfile(fullmod, varlist = c("fish")) #THis one doesnt work 

# definite correlation

# ACF plot with blocks as survey sesh, day hour
df$session_id <- factor(df$session_id)
df$day_hour_as_block <- as.factor(paste(df$jul_day, df$start_hr,sep = ""))
ps <- par(ask = FALSE, mfrow = c(1,1))
runACF(df$day_hour_as_block, fullmod, store = F) #some positive AC, use this as our cor block

# acf plot w session id as block
ps <- par(ask = FALSE, mfrow = c(1,1))
runACF(df$session_id, fullmod, store = F) #slots of +AC -- unsure which to pick

# try with sessionid + grid id
ps <- par(ask = FALSE, mfrow = c(1,1))
runACF(df$session_grid_block, fullmod, store = F) 
# this takes literally 5 minutes to run. and now we have zero AC--but it feels wrong

### STEP 4: plot cumulative residuals for model to check over or under prediction
par(ask = F, mfrow = c(2,3))
plotCumRes(fullmod, varlist= varList, splineParams) 

# compare with model that just has linear terms
par(ask = F, mfrow = c(2,3))
plotCumRes(fullmod_linear, varlist= c(varList), splineParams) 

#  looks like x and y are linear and all the rest shouldnt be
#   model defo over and under predicting at times, lets run SALSA for more flexibility

### STEP 5: Run SALSA1D
# Create the variable "response" needed by SALSA
df$response <- df$pa

## Set initial model without spline-based terms
initialModel <- glm(response ~ treatment, family=binomial, data=df)

## Set SALSA arguments
factorList <- c("obs", "treatment", "start_hr")
varList <- c("x.pos", "y.pos", "fish", "jul_day")
salsa1DList <- list(fitnessMeasure="cv.gamMRSea", minKnots_1d=rep(1, 4),
                    maxKnots_1d=rep(5, 4), startKnots_1d=rep(1, 4),
                    degree=rep(2, 4), maxIterations=100,
                    gaps=rep(0, 4))




## Run SALSA with removal
set.seed(53195)
salsa1D_RT <- MRSea::runSALSA1D(initialModel, salsa1DList, varList,
                           factorList, varlist_cyclicSplines=NULL,
                           splineParams=NULL, datain=df,
                           suppress.printout=TRUE, removal=T,
                           panelid=NULL)
bestModel1D<- salsa1D_RT$bestModel
summary(bestmod1D)

## Run SALSA without removal
salsa1D_RF <- MRSea::runSALSA1D(initialModel, salsa1DList, varList,
                                factorList, varlist_cyclicSplines=NULL,
                                splineParams=NULL, datain=df,
                                suppress.printout=TRUE, removal=F,
                                panelid=NULL)

bestmod1D_RF<- salsa1D_RF$bestModel
summary(bestmod1D_RF)

### STEP 6: Pick best model based on AIC
AIC(bestmod1D_RT, bestmod1D_RF)

salsa1D_RF$modelFits1D
salsa1D_RT$modelFits1D

# salsa took the initial model with one knot at mean for each var and added 
# 2 knots for fish at 164 and 728; and 3 knots for start hour at 10, 13, 15;
# The result of SALSA is additional flexibility added to the model for the
# relationship between presence absence probability and fish / hr

### STEP 7: Select initial knot locations using a space-filling design

# Next we add a two dimensional smooth of geographic coordinates (s(x.pos,
# y.pos))# To test for a redistribution of animals with a treatment effect we fit an
# interaction term between the smooth of coordinates and impact.
# SALSA is used to determine spatially adaptive knot locations for this
# smooth term.

## create knotgrid
knotgrid <- MRSea::getKnotgrid(coordData=cbind(df$x.pos, df$y.pos),
                               numKnots=300, plot=FALSE)

### STEP 8: Create distance matrix
# i.e. distance between each data points and knots)

## make distance matrices for datatoknots and knottoknots
distMats <- makeDists(cbind(df$x.pos, df$y.pos), na.omit(knotgrid))

## create sequence of radii
r_seq <- getRadiiChoices(numberofradii= 8, distMats$dataDist, basis = "gaussian")



### STEP 9: Set SALSA parameters start knots as 6
salsa2DList6 <- list(fitnessMeasure="AIC", knotgrid=knotgrid,
                    startKnots=6, minKnots=2, maxKnots=20,
                    gap=1, r_seq = r_seq, interactionTerm = "treatment")


# splineParams must be an object in workspace
## update splineParams with the SALSA1D results
splineParams <- salsa1D_RT$splineParams

### STEP 10: RUN SALSA 2D
set.seed(53195)
salsa2D_6 <- MRSea::runSALSA2D(bestModel1D, salsa2DList6, d2k=distMats$dataDist,
                             k2k=distMats$knotDist, splineParams=NULL,
                             tol=0, chooserad=F, panels=NULL,
                             suppress.printout=TRUE)
## Trial start knots as 7
salsa2DList7 <- list(fitnessMeasure="AIC", knotgrid=knotgrid,
                     startKnots=7, minKnots=2, maxKnots=20,
                     gap=1, r_seq = r_seq, interactionTerm = "treatment")
# reset
splineParams <- salsa1D_RT$splineParams
set.seed(53195)
salsa2D_7 <- MRSea::runSALSA2D(bestModel1D, salsa2DList7, d2k=distMats$dataDist,
                                                   k2k=distMats$knotDist, splineParams=NULL,
                                                   tol=0, chooserad=F, panels=NULL,
                                                   suppress.printout=TRUE)

## Trial start knots as 8
salsa2DList8 <- list(fitnessMeasure="AIC", knotgrid=knotgrid,
                    startKnots=8, minKnots=2, maxKnots=20,
                    gap=1, r_seq = r_seq, interactionTerm = "treatment")
#reset
splineParams <- salsa1D_RT$splineParams
set.seed(53195)
salsa2D_8 <- MRSea::runSALSA2D(bestModel1D, salsa2DList8, d2k=distMats$dataDist,
                             k2k=distMats$knotDist, splineParams=NULL,
                             tol=0, chooserad=F, panels=NULL,
                             suppress.printout=TRUE)

## Trial start knots as 9
salsa2DList9 <- list(fitnessMeasure="AIC", knotgrid=knotgrid,
                     startKnots=9, minKnots=2, maxKnots=20,
                     gap=1, r_seq = r_seq, interactionTerm = "treatment")
#reset
splineParams <- salsa1D_RT$splineParams
set.seed(53195)
salsa2D_9 <- MRSea::runSALSA2D(bestModel1D, salsa2DList9, d2k=distMats$dataDist,
                                                   k2k=distMats$knotDist, splineParams=NULL,
                                                   tol=0, chooserad=F, panels=NULL,
                                                   suppress.printout=TRUE)

## Trial start knots as 10
salsa2DList10 <- list(fitnessMeasure="AIC", knotgrid=knotgrid,
                     startKnots=10, minKnots=2, maxKnots=20,
                     gap=1, r_seq = r_seq, interactionTerm = "treatment")
#reset
splineParams <- salsa1D_RT$splineParams
set.seed(53195)
salsa2D_10 <- MRSea::runSALSA2D(bestModel1D, salsa2DList10, d2k=distMats$dataDist,
                                                   k2k=distMats$knotDist, splineParams=NULL,
                                                   tol=0, chooserad=F, panels=NULL,
                                                   suppress.printout=TRUE)


## compare CV Scores
CV1D <- getCV_CReSS(df, salsa1D_RT$bestModel) # 1D model
CV2D6 <- getCV_CReSS(df, salsa2D_6$bestModel) # 6 start knots
CV2D7 <- getCV_CReSS(df, salsa2D_7$bestModel) 
CV2D8 <- getCV_CReSS(df, salsa2D_8$bestModel)
CV2D9 <- getCV_CReSS(df, salsa2D_9$bestModel) 
CV2D10 <- getCV_CReSS(df, salsa2D_10$bestModel)

compareCV <- data.frame(model = as.character(c("CV1D", "CV2D6", "CV2D7", "CV2D8", "CV2D9", "CV2D10")),
                cv_score = c(CV1D, CV2D6, CV2D7, CV2D8, CV2D9, CV2D10)) 
                

## compare AIC scores
AIC(salsa1D_RT$bestModel, salsa2D_6$bestModel, salsa2D_7$bestModel, salsa2D_8$bestModel, 
salsa2D_9$bestModel, salsa2D_10$bestModel)

## Store "best" model
baseModel <- salsa2D_6$bestModel

## Summary of results
summary(baseModel)
stats::anova(baseModel, test="Chisq")

## update spline parameter object
splineParams <- salsa2D_6$splineParams

### STEP 11: check for residual correlation
runs.test(residuals(baseModel, type = "pearson")) #still neg stat, still tiny p

### STEP 12: Significant positive resid correlation therefore refit as GEE
# the data must be ordered by block (which this is)
# and the blockid must be numeric

## specify parameters for local radial:
radiusIndices <- splineParams[[1]]$radiusIndices
dists <- splineParams[[1]]$dist
radii <- splineParams[[1]]$radii
#aR <- splineParams[[1]]$invInd[splineParams[[1]]$knotPos] 
# ^ this comes up as null, aR not within SP, in larger model object (hope ok)
aR <- as.vector(salsa2D_6$aR$aR1) 

## update model in workspace with parameters for spatial smooth (above)
baseModel <- update(baseModel, . ~ .)

## Re-fit the chosen model as a GEE (based on SALSA knot placement) and
geeModel <- geepack::geeglm(formula(baseModel), data = df, family = "binomial",
                   id = session_grid_block)
## GEE p-values
getPvalues(geeModel) 
summary(geeModel)

### STEP 13: Remove impact to compare
geeModel$formula 
# grab model from above ^
noint.model <- response ~  treatment + bs(x.pos, knots = splineParams[[2]]$knots,degree = splineParams[[2]]$degree, Boundary.knots = splineParams[[2]]$bd) +
                           bs(y.pos, knots = splineParams[[3]]$knots, degree = splineParams[[3]]$degree,Boundary.knots = splineParams[[3]]$bd) + 
                           bs(fish, knots = splineParams[[4]]$knots, degree = splineParams[[4]]$degree, Boundary.knots = splineParams[[4]]$bd) + 
                           bs(jul_day, knots = splineParams[[5]]$knots, degree = splineParams[[5]]$degree,Boundary.knots = splineParams[[5]]$bd) + 
                           bs(start_hr,knots = splineParams[[6]]$knots, degree = splineParams[[6]]$degree, Boundary.knots = splineParams[[6]]$bd) + 
                           LRF.g(radiusIndices, dists, radii, aR)
nointModel <- geepack::geeglm(formula(noint.model), data = df, family = "binomial",id = session_grid_block)
## compare p-values
getPvalues(nointModel)

## compare AIC
AIC(geeModel, noint.model)

### Step 14: Partial residual plots doesnt work for noncontinuous variables?
# par(mfrow = c(2, 3))
runPartialPlots(model = geeModel, data = df)
runPartialPlots(geeModel, df, varlist.in = c("fish"))
# runPartialPlots(geeModel, df, varlist = c("start_hr"))
# runPartialPlots(geeModel, df, varlist = c("jul_day"))

### Step 15: DIAGNOSTICS
# create observed vs fitted and fitted vs residual plots
runDiagnostics(geeModel) 
# How to interpret these with binary output

### Step 16: COVRATIO and PRESS Statistic --- This will take 353 minutes or 6 hours... run over night
timeInfluenceCheck(geeModel, df$session_grid_block, dists, splineParams)
# influence plots (covratio and press statistics)
# influence <- runInfluence(geeModel, data$session_grid_block,  dists, splineParams)

### Step 18: residual plots
resids <- fitted(geeModel) - df$pa
dims <- getPlotdimensions(df$x.pos, df$y.pos, 10, 10)
par(mfrow = c(1, 2), mar = c(3, 3, 3, 5))
quilt.plot(df$x.pos[df$treatment == "OFF"], df$y.pos[df$treatment == "OFF"], 
           resids[df$treatment == "OFF"], asp = 1, ncol = dims[2], nrow = dims[1],
           zlim = c(-2.2, 2.2), main = "OFF",xlab = "x", ylab = "y")
quilt.plot(df$x.pos[df$treatment == "ON"], df$y.pos[df$treatment == "ON"], 
           resids[df$treatment == "ON"], asp = 1, ncol = dims[2], nrow = dims[1],
           zlim = c(-2.2, 2.2), main = "ON", xlab = "x", ylab = "y")

## Step 19: Influence plots
timeInfluenceCheck(geeModel, df$block_as_session, dists, splineParams)
# influence plots (covratio and press statistics) this will take 6hrs
# influence <- runInfluence(bestModel2D, df$block_as_session, dists, splineParams)


## Map model for On and Off
par(mfrow=c(1,2))
# Loop across all phases
for (treat in c("ON", "OFF"))
{
  bWant <- df$treat %in% treat # a boolean specifying phase
  fields::quilt.plot(df$XPos[bWant],
                     df$YPos[bWant],
                     fitted(bestModel2D)[bWant],
                     
                     nrow=10, ncol=10,
                     zlim=range(fitted(bestModel2D)),
                     main=paste0("TAST ", treat))
  if (treat %in% "ON")
  {
  #somethign about turn off legend for one of the 
  }

}

#### Predictions 
predict.gamMRSea()
############
