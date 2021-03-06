### Binomial Presence Absence models

## read in data
padat <- read.csv("~/Desktop/TASTBALLARD/TAST WD/Dissertation/pa_data.csv")

## load libraries
library(tidyverse)
library(car)
library(lawstat)
library(MRSea)
library(mgcv)
library(splines)
library(geepack)

## format data
df <- padat %>%
  select(pa, treatment, obs, x, y, session_id, id, jul_day, start_hr, obs, fish)
df$pa <- as.numeric(df$pa)
df$treatment <- as.factor(df$treatment)
df$obs <- as.factor(df$obs)
df$x.pos <- df$x
df$y.pos <- df$y

#########
# Fit GEE CReSS model with SALSA2D 

### STEP 1: check for collinearity between covariates in glm model
fullmod_linear <- glm(pa ~ treatment + jul_day + start_hr + obs + dist2sq + fish + x.pos + y.pos,
               na.action = "na.fail", family = "binomial", data = df)
# check for collinearity
car::vif(fullmod_linear) # most between 1 & 2.1 so pretty much good! dist2Sq and x/y pretty collin so throw out dist

# refit without distance 
fullmod_linear <- glm(pa ~ treatment + jul_day + start_hr + obs + fish + x.pos + y.pos + x.pos:y.pos,
               na.action = "na.fail", family = "binomial", data = df)

# check for collinearity
car::vif(fullmod_linear) # much better, fish and day only slightly collinear but we can keep it. (1-2)

# STEP 2: make splineParams object with knots at the mean of each variable
splineParams<-makesplineParams(data=df, varlist=c("x.pos", "y.pos", "fish","jul_day", "start_hr"))

# fit new full mod with splines binomial dist
fullmod <- glm(pa ~ treatment + obs + bs(x.pos, knots = splineParams[[2]]$knots) 
                          + bs(y.pos, knots = splineParams[[3]]$knots)
                          + bs(fish,knots = splineParams[[4]]$knots)
                          + bs(jul_day,knots = splineParams[[5]]$knots)
                          + bs(start_hr, knots =splineParams[[6]]$knots),
                          na.action = "na.fail", family = "binomial", data = df)

# check qbn distribution to see if dispersion para is 0
fullmod_qbn <- glm(pa ~ treatment + obs + bs(x.pos, knots = splineParams[[2]]$knots) 
               + bs(y.pos, knots = splineParams[[3]]$knots)
               + bs(fish,knots = splineParams[[4]]$knots)
               + bs(jul_day,knots = splineParams[[5]]$knots)
               + bs(start_hr, knots =splineParams[[6]]$knots),
               na.action = "na.fail", family = "quasibinomial", data = df)

summary(fullmod_qbn) #nope dispparm = 1.14, not zero
summary(fullmod)

### STEP 3: check for correlated residuals in numeric covariates
varList <- c("XPos", "YPos", "fish", "jul_day", "start_hr")
runs.test(residuals(fullmod, type = "pearson", alternative = "two-sided")) # neg standardised runs stat and sm p == houston we have correlation
par(ask = FALSE, mfrow = c(2,2))
runs_plots <- plotRunsProfile(fullmod, varlist = c("start_hr")) #definite correlation

# ACF plot with blocks as survey sesh, day hour
df$session_id <- factor(df$session_id)
df$day_hour_as_block <- as.factor(paste(df$jul_day, df$start_hr,sep = ""))
par(ask = FALSE, mfrow = c(1,1))
runACF(df$day_hour_as_block, fullmod, store = F) #some positive AC, use this as our cor block

# acf plot w session id as block
par(ask = FALSE, mfrow = c(1,1))
runACF(df$session_id, fullmod, store = F) #slots of +AC -- unsure which to pick

# try with sessionid+ grid id
df$session_grid_block <- as.factor(paste(df$session_id, df$id, sep= ""))
par(ask = FALSE, mfrow = c(1,1))
runACF(df$session_grid_block, fullmod, store = F) 
# this takes literally 5 minutes to run. and now we have zero AC--but it feels wrong

### STEP 4: plot cumulative residuals for model to check over or under prediction
par(ask = F, mfrow = c(3,3))
plotCumRes(fullmod, varlist= varList, splineParams) 

# compare with model that just has linear terms
par(ask = F, mfrow = c(3,3))
plotCumRes(fullmod_linear, varlist= c(varList), splineParams) 

#looks like x and y are linear and all the rest shouldnt be
# model defo over and under predicting at times, lets run SALSA for more flexibility

# STEP 5: Run SALSA
library(devtools)
devtools::install_github(repo="lindesaysh/MRSea@v1.01")

# Create the variable "response" needed by SALSA
df$response <- df$pa

# create foldid variable # not necessary in new version of salsa
#df$foldid<-getCVids(data = df, folds = 5, block = "day_hour_as_block")

# Set initial model without spline-based terms
initialModel <- glm(response ~ treatment, family=binomial, data=df)

# Set SALSA arguments
factorList <- c("obs", "treatment")
varList <- c("x.pos", "y.pos", "fish", "jul_day", "start_hr")
salsa1DList <- list(fitnessMeasure="AIC", minKnots_1d=rep(1, 5),
                    maxKnots_1d=rep(5, 5), startKnots_1d=rep(1, 5),
                    degree=rep(2, 5), maxIterations=100,
                    gaps=rep(0, 5))

# Load library
library(MRSea)
set.seed(53195)
# Run SALSA with removal
salsa1D_RT <- MRSea::runSALSA1D(initialModel, salsa1DList, varList,
                           factorList, varlist_cyclicSplines=NULL,
                           splineParams=NULL, datain=df,
                           suppress.printout=TRUE, removal=T,
                           panelid=NULL)
bestmod1D_RT<- salsa1D_RT$bestModel
summary(bestmod1D_RT)

# STEP 6: Pick best model based on AIC
AIC(bestmod1D_RT, bestmod1D_RF)

# STEP 7: Select initial knot locations using a space-filling design

#create knotgrid
knotgrid <- MRSea::getKnotgrid(coordData=cbind(df$x.pos, df$y.pos),
                               numKnots=300, plot=FALSE)

# STEP 8: Create distance matrix
# (i.e. distance between each data points and knots)
# make distance matrices for datatoknots and knottoknots

distMats <- makeDists(cbind(df$x.pos, df$y.pos), na.omit(knotgrid))

# create sequence of radii
r_seq <- getRadiiChoices(numberofradii= 8, distMats$dataDist, basis = "gaussian")

###### without Interaction
#######
# STEP 9: Set SALSA parameters start knots as 6
salsa2DList <- list(fitnessMeasure="AIC", knotgrid=knotgrid,
                    startKnots=6, minKnots=2, maxKnots=20,
                    gap=1, r_seq = r_seq, interactionTerm = "treatment")


# splineParams must be an object in workspace
# update splineParams with the SALSA1D results
splineParams <- salsa1D_RT$splineParams

# STEP 10: RUN SALSA 2D
set.seed(53195)
salsa2D <- runSALSA2D(bestModel1D, salsa2DList, d2k=distMats$dataDist,
                             k2k=distMats$knotDist, splineParams=NULL,
                             tol=0, chooserad=F, panels=NULL,
                             suppress.printout=TRUE)
########## start knots as 7
salsa2DList7 <- list(fitnessMeasure="AIC", knotgrid=knotgrid,
                     startKnots=7, minKnots=2, maxKnots=20,
                     gap=1, r_seq = r_seq, interactionTerm = "treatment")
#reset
splineParams <- salsa1D$splineParams

# salsa w knots at 8
set.seed(53195)
t7 <- system.time({ salsa2D_7 <- MRSea::runSALSA2D(bestModel1D, salsa2DList7, d2k=distMats$dataDist,
                                                   k2k=distMats$knotDist, splineParams=NULL,
                                                   tol=0, chooserad=F, panels=NULL,
                                                   suppress.printout=TRUE)
})

########## start knots as 8
salsa2DList8 <- list(fitnessMeasure="AIC", knotgrid=knotgrid,
                    startKnots=8, minKnots=2, maxKnots=20,
                    gap=1, r_seq = r_seq, interactionTerm = "treatment")
#reset
splineParams <- salsa1D$splineParams

# salsa w knots at 8
set.seed(53195)
t8 <- system.time({ salsa2D_8 <- MRSea::runSALSA2D(bestModel1D, salsa2DList8, d2k=distMats$dataDist,
                             k2k=distMats$knotDist, splineParams=NULL,
                             tol=0, chooserad=F, panels=NULL,
                             suppress.printout=TRUE)
})
########## start knots as 9
salsa2DList9 <- list(fitnessMeasure="AIC", knotgrid=knotgrid,
                     startKnots=9, minKnots=2, maxKnots=20,
                     gap=1, r_seq = r_seq, interactionTerm = "treatment")
#reset
splineParams <- salsa1D$splineParams

# salsa w knots at 9
set.seed(53195)
t9 <- system.time({ salsa2D_9 <- MRSea::runSALSA2D(bestModel1D, salsa2DList9, d2k=distMats$dataDist,
                                                   k2k=distMats$knotDist, splineParams=NULL,
                                                   tol=0, chooserad=F, panels=NULL,
                                                   suppress.printout=TRUE)
})

########## start knots as 10
salsa2DList10 <- list(fitnessMeasure="AIC", knotgrid=knotgrid,
                     startKnots=10, minKnots=2, maxKnots=20,
                     gap=1, r_seq = r_seq, interactionTerm = "treatment")
#reset
splineParams <- salsa1D$splineParams

# salsa w knots at 10
set.seed(53195)
t10 <- system.time({ salsa2D_10 <- MRSea::runSALSA2D(bestModel1D, salsa2DList10, d2k=distMats$dataDist,
                                                   k2k=distMats$knotDist, splineParams=NULL,
                                                   tol=0, chooserad=F, panels=NULL,
                                                   suppress.printout=TRUE)
})


# compare CV Scores
getCV_CReSS(df, salsa1D_RT$bestModel)
getCV_CReSS(df, salsa2D$bestModel) #6 start knots
getCV_CReSS(df, salsa2D_7$bestModel) 
getCV_CReSS(df, salsa2D_8$bestModel)
getCV_CReSS(df, salsa2D_9$bestModel) 
getCV_CReSS(df, salsa2D_10$bestModel)

#compare AIC scores
AIC(salsa1D_RT$bestModel, salsa2D$bestModel, salsa2D_7$bestModel, salsa2D_8$bestModel, 
salsa2D_9$bestModel, salsa2D_10$bestModel)


# Store "best" model
baseModel <- salsa2D_9$bestModel

# Summary of results
summary(baseModel)

stats::anova(baseModel, test="Chisq")

# update spline parameter object
splineParams <- salsa2D_9$splineParams

# STEP 11: check for residual correlation
runs.test(residuals(baseModel, type = "pearson")) #still neg stat, still tiny p

## STEP 12: Significant positive resid correlation therefore refit as GEE
########
# N.B. for the GEE formula, the data must be ordered by block (which this is)
# and the blockid must be numeric
# specify parameters for local radial:

radiusIndices <- splineParams[[1]]$radiusIndices
dists <- splineParams[[1]]$dist
radii <- splineParams[[1]]$radii
aR <- splineParams[[1]]$invInd[splineParams[[1]]$knotPos] ## this comes up as null, so had to workaround aR not within SP, in larger model object
aR <- as.vector(salsa2D_9$aR$aR1) 

# update model in workspace with parameters for spatial smooth (above)
baseModel <- update(baseModel, . ~ .)

# Re-fit the chosen model as a GEE (based on SALSA knot placement) and
# GEE p-values #see if model improves with different blocking str
geeModel <- geepack::geeglm(formula(baseModel), data = df, family = "binomial",
                   id = session_grid_block)
getPvalues(geeModel) 
summary(geeModel)

# STEP 13: Remove impact how to remove impact
geeModel$formula

noint.model <- response ~  treatment + bs(x.pos, knots = splineParams[[2]]$knots,degree = splineParams[[2]]$degree, Boundary.knots = splineParams[[2]]$bd) +
                           bs(y.pos, knots = splineParams[[3]]$knots, degree = splineParams[[3]]$degree,Boundary.knots = splineParams[[3]]$bd) + 
                           bs(fish, knots = splineParams[[4]]$knots, degree = splineParams[[4]]$degree, Boundary.knots = splineParams[[4]]$bd) + 
                           bs(jul_day, knots = splineParams[[5]]$knots, degree = splineParams[[5]]$degree,Boundary.knots = splineParams[[5]]$bd) + 
                           bs(start_hr,knots = splineParams[[6]]$knots, degree = splineParams[[6]]$degree, Boundary.knots = splineParams[[6]]$bd) + 
                           LRF.g(radiusIndices, dists, radii, aR)
nointModel <- geepack::geeglm(formula(noint.model), data = df, family = "binomial",id = session_grid_block)
# reshow p-values
getPvalues(nointModel)

AIC(geeModel)

str(noint.model)

## Step 14: Partial residual plots doesnt work for noncontinuous variables?
# par(mfrow = c(2, 3))
# runPartialPlots(geeModel, df)
# runPartialPlots(geeModel, df, varlist.in = c("fish"))
# runPartialPlots(geeModel, df, varlist = c("start_hr"))
# runPartialPlots(geeModel, df, varlist = c("jul_day"))

## Step 15: DIAGNOSTICS
# create observed vs fitted and fitted vs residual plots
runDiagnostics(geeModel, plotting = b, save = T) #How to interpret these with binary output

## Step 16: COVRATIO and PRESS Statistic --- This will take 353 minutes or 6 hours... run over night
timeInfluenceCheck(geeModel, df$session_grid_block, dists, splineParams)
# influence plots (covratio and press statistics)
# influence <- runInfluence(geeModel, data$session_grid_block,  dists, splineParams)

# step 18: residual plots
# residual plot ## ADD LABELS
resids <- fitted(geeModel) - df$pa
dims <- getPlotdimensions(df$x.pos, df$y.pos, 10, 10)
par(mfrow = c(1, 2), mar = c(3, 3, 3, 5))
quilt.plot(df$x.pos[df$treatment == "OFF"], df$y.pos[df$treatment == "OFF"], 
           resids[df$treatment == "OFF"], asp = 1, ncol = dims[2], nrow = dims[1],
           zlim = c(-2.2, 2.2), main = "OFF",xlab = "x", ylab = "y")
quilt.plot(df$x.pos[df$treatment == "ON"], df$y.pos[df$treatment == "ON"], 
           resids[df$treatment == "ON"], asp = 1, ncol = dims[2], nrow = dims[1],
           zlim = c(-2.2, 2.2), main = "ON", xlab = "x", ylab = "y")
#######
############# Fit with interaction between x, y, treat
# Set SALSA parameters 
salsa2DList <- list(fitnessMeasure="AIC", knotgrid=knotgrid,
                    startKnots=6, minKnots=2, maxKnots=20,
                    gap=1, r_seq = r_seq, interactionTerm = "as.factor(treatment)")
# splineParams must be an object in workspace
# update splineParams with the SALSA1D results
splineParams <- salsa1D$splineParams
# STEP 10: RUN SALSA 2D
set.seed(53195)
salsa2Dint <- MRSea::runSALSA2D(bestModel1D, salsa2DList, d2k=distMats$dataDist,
                             k2k=distMats$knotDist, splineParams=NULL,
                             tol=0, chooserad=F, panels=NULL,
                             suppress.printout=TRUE)


## multiple salsa runs with different starting knots ## doesnt work dont know why

# Store "best" model
bestModel2Dint <- salsa2Dint$bestModel

# Summary of results
summary(bestModel2Dint)
# save the model object
baseModelint <- salsa2Dint$bestModel
# update spline parameter object
splineParams <- salsa2Dint$splineParams

# check for residual correlation
runs.test(residuals(bestModel2Dint, type = "pearson")) #still neg stat, still tiny p

## Significant positive resid correlation therefore refit as GEE CANT
########
# N.B. for the GEE formula, the data must be ordered by block (which this is)
# and the blockid must be numeric
# specify parameters for local radial:
radiusIndices <- splineParams[[1]]$radiusIndices
dists <- splineParams[[1]]$dist
radii <- splineParams[[1]]$radii
aR <- as.vector(salsa2Dint$aR$aR1) ## this comes up as null, so had to workaround aR not within SP, in larger model object
# update model in workspace with parameters for spatial smooth (above)
baseModelint <- update(baseModelint, . ~ .)

# Re-fit the chosen model as a GEE (based on SALSA knot placement) and
# GEE p-values
geeModel <- geepack::geeglm(formula(baseModel), data = df, family = "binomial",
                            id = session_grid_block)
getPvalues(geeModel, varlist = varList, factorlist = "treatment")



## Predictions and Interpretations 
splineParams[[1]]$invInd
# divide data into on and off

###############

runDiagnostics(bestModel2D)
# not great but hard to tell with binary

timeInfluenceCheck(bestModel2D, df$block_as_session, dists, splineParams)
# influence plots (covratio and press statistics)
influence <- runInfluence(bestModel2D, df$block_as_session, dists, splineParams)


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








# legend things?
par(oma=c( 0,0,0,4)) # margin of 4 spaces width at right hand side
set.panel( 1,2) # 2X2 matrix of plots

# now draw all your plots using usual image command
for (  k in 1:4){
  data<- matrix( rnorm(150), 10,15)
  image( data, zlim=c(-4,4), col=tim.colors())
  # and just for fun add a contour plot  
  contour( data, add=TRUE)
}

par(oma=c( 0,0,0,1))# reset margin to be much smaller.
image.plot( legend.only=TRUE, zlim=c(-4,4)) 

# image.plot tricked into  plotting in margin of old setting 

set.panel() # reset plotting device

