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

## check for collinearity
vif(fullmod_linear) 
# most between 1 & 2.1 so pretty much good! dist2Sq and x/y pretty collin so throw out dist

## refit without distance 
fullmod_linear <- glm(pa ~ treatment + jul_day + start_hr + obs + fish + x.pos + y.pos + x.pos:y.pos,
               na.action = "na.fail", family = "binomial", data = df)

## recheck for collinearity
vif(fullmod_linear) 
# much better, fish and day only slightly collinear but we can keep it. (1-2)

# STEP 2: make splineParams object with knots at the mean of each variable
splineParams<-makesplineParams(data=df, varlist=c("x.pos", "y.pos", "fish","jul_day", "start_hr"))

# fit new full mod with splines binomial dist
fullmod <- glm(pa ~ treatment + obs + bs(x.pos, knots = splineParams[[2]]$knots) 
               + bs(y.pos, knots = splineParams[[3]]$knots)
               + bs(fish,knots = splineParams[[4]]$knots)
               + bs(jul_day,knots = splineParams[[5]]$knots)
               + bs(start_hr, knots =splineParams[[6]]$knots),
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
varList <- c("y.pos", "y.pos", "fish", "jul_day", "start_hr")
runs.test(residuals(fullmod, type = "pearson", alternative = "two-sided")) 
# neg standardised runs stat and sm p == we have correlation

par(ask = FALSE, mfrow = c(2,2))
runs_plots <- plotRunsProfile(fullmod, varlist = c("start_hr")) 

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
df$session_grid_block <- as.factor(paste(df$session_id, df$id, sep= ""))
ps <- par(ask = FALSE, mfrow = c(1,1))
runACF(df$session_grid_block, fullmod, store = F) 
# this takes literally 5 minutes to run. and now we have zero AC--but it feels wrong

### STEP 4: plot cumulative residuals for model to check over or under prediction
par(ask = F, mfrow = c(3,3))
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
factorList <- c("obs", "treatment")
varList <- c("x.pos", "y.pos", "fish", "jul_day", "start_hr")
#create salsa list 
salsa1DList <- list(fitnessMeasure="AIC", minKnots_1d=rep(1, 5),
                    maxKnots_1d=rep(5, 5), startKnots_1d=rep(1, 5),
                    degree=rep(2, 5), maxIterations=100,
                    gaps=rep(0, 5))

## Run SALSA with removal
set.seed(53195)
salsa1D_RT <- MRSea::runSALSA1D(initialModel, salsa1DList, varList,
                           factorList, varlist_cyclicSplines=NULL,
                           splineParams=NULL, datain=df,
                           suppress.printout=TRUE, removal=T,
                           panelid=NULL)
bestmod1D_RT<- salsa1D_RT$bestModel

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

# model with removal is the best model
bestModel1D <- salsa1D_RT$bestModel

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
#salsa2DList7 <- list(fitnessMeasure="AIC", knotgrid=knotgrid,
                    ## gap=1, r_seq = r_seq, interactionTerm = "treatment")
# reset
#splineParams <- salsa1D_RT$splineParams
# #set.seed(53195)
# #salsa2D_7 <- MRSea::runSALSA2D(bestModel1D, salsa2DList7, d2k=distMats$dataDist,
#                                                    k2k=distMats$knotDist, splineParams=NULL,
#                                                    tol=0, chooserad=F, panels=NULL,
#                                                    suppress.printout=TRUE)

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
# 
# ## Trial start knots as 9
# #salsa2DList9 <- list(fitnessMeasure="AIC", knotgrid=knotgrid,
#                      startKnots=9, minKnots=2, maxKnots=20,
#                      gap=1, r_seq = r_seq, interactionTerm = "treatment")
# #reset
# #splineParams <- salsa1D_RT$splineParams
# #set.seed(53195)
# #salsa2D_9 <- MRSea::runSALSA2D(bestModel1D, salsa2DList9, d2k=distMats$dataDist,
#                                                    k2k=distMats$knotDist, splineParams=NULL,
#                                                    tol=0, chooserad=F, panels=NULL,
#                                                    suppress.printout=TRUE)

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


# ## compare CV Scores
# CV1D <- getCV_CReSS(df, salsa1D_RT$bestModel) # 1D model
# CV2D6 <- getCV_CReSS(df, salsa2D_6$bestModel) # 6 start knots
# CV2D7 <- getCV_CReSS(df, salsa2D_7$bestModel) 
# CV2D8 <- getCV_CReSS(df, salsa2D_8$bestModel)
# CV2D9 <- getCV_CReSS(df, salsa2D_9$bestModel) 
# CV2D10 <- getCV_CReSS(df, salsa2D_10$bestModel)
# 
# compareCV <- data.frame(model = as.character(c("CV1D", "CV2D6", "CV2D7", "CV2D8", "CV2D9", "CV2D10")),
#                 cv_score = c(CV1D, CV2D6, CV2D7, CV2D8, CV2D9, CV2D10)) 
                

## compare AIC scores
AIC(salsa1D_RT$bestModel, salsa2D_6$bestModel, salsa2D_8$bestModel, salsa2D_10$bestModel)
    
#salsa2D_9$bestModel, 
#salsa2D_7$bestModel, 

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
AIC(geeModel, nointModel)

### Step 14: Partial residual plots 
# par(mfrow = c(2, 3))
runPartialPlots(model = geeModel, data = df)
runPartialPlots(geeModel, df, varlist.in = c("fish"))
# runPartialPlots(geeModel, df, varlist = c("start_hr"))
# runPartialPlots(geeModel, df, varlist = c("jul_day"))

### Step 15: DIAGNOSTICS
# create observed vs fitted and fitted vs residual plots
runDiagnostics(geeModel) 
# How to interpret these with binary output

### Step 16: COVRATIO and PRESS Statistic --- This will take 276 minutes or 4.6 hours... run over night
timeInfluenceCheck(geeModel, df$session_grid_block, dists, splineParams)
# influence plots (covratio and press statistics)
influence <- runInfluence(geeModel, df$session_grid_block,  dists, splineParams)

### Step 17: RAW residual plots
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
influence <- runInfluence(geeModel, df$block_as_session, dists, splineParams)

### Step 18: Predictions
# loading the prediction grid data
df$area <- rep(100, length(df$pa)) #add area column
predictData <- df
# generate a prediction data set -- ID differences only works if there are the same numnber of rows for on and off
## so pred data will be a randomselection of 80% OFF data and nrow(0.8*OFF) of on data, then rbound
# #split data to on and off
# dataON <- filter(df, treatment == "ON")
# dataOFF <- filter(df, treatment == "OFF")
# #shuffle row indices
# set.seed(531)
# rowsOFF <- sample(nrow(dataOFF))
# # reorder by random indices
# dataOFF <- dataOFF[rowsOFF, ]
# #choose split point 80% THROUGH THE DATA
# split <- round(nrow(dataOFF) * 0.80)
# # select 80% of randomly shuffled off data
# dataOFF <- dataOFF[1:split, ]
# # now shuffle on data
# set.seed(531)
# rowsON <- sample(nrow(dataON))
# # reorder 
# dataON <- dataON[rowsON, ]
# #select same number of off as on 
# dataON <- dataON[1:split, ]
# # bind theseback together
# predictData <- rbind(dataON, dataOFF)
# # okay dropped my sample size in HALF ask steve about this
# head(predictData)

# create the distance matrix for predictions
dists <- makeDists(cbind(predictData$x.pos, predictData$y.pos),
                   na.omit(knotgrid), knotmat = FALSE)$dataDist
# use baseModel to make predictions to avoid a warning from
# using geeModel (same answers though)
predslink <- predict(baseModel, predictData, type = "link")
# reversing the logit-link to convert predictions back to the response scale
preds <- plogis(predslink)
length(preds)
df$predictions <- preds

### Step 19 visualising predictions
# plotting the predictions for before and after impact
# get the plot dimensions. We know each cell is 10x10m
dims<-getPlotdimensions(x.pos=predictData$x.pos, predictData$y.pos,
                        segmentWidth=10, segmentLength=10)
par(mfrow=c(1,2), mar=c(3,3,3,5))
quilt.plot(predictData$x.pos[predictData$treatment=="OFF"],
           predictData$y.pos[predictData$treatment=="OFF"],
           preds[predictData$treatment=="OFF"], asp=1, nrow=dims[1], ncol=dims[2], 
           zlim=c(0, maxlim = 0.3), main = "TAST OFF", add.legend = F)

quilt.plot(predictData$x.pos[predictData$treatment=="ON"],
           predictData$y.pos[predictData$treatment=="ON"], preds[predictData$treatment=="ON"],
           asp=1,nrow=dims[1], ncol=dims[2], zlim=c(0, maxlim = 0.3), main = "TAST ON")

# add legend label?
str(preds[predictData$treatment=="ON"])

## Step 20: bootstrap CI
# do the bootstrap # do this at 999 when you get the chance
do.bootstrap.cress(df, df, ddf.obj=NULL, baseModel, splineParams,
                   dists, B=20)

# read in bootstrap predictions
load("predictionboot.RData")
str(bootPreds)

bootPreds$id <- as.factor(df$id)


df$id <- as.factor(df$id)

df$group_by(levels(df$id))
# make percentile confidence intervals




cison <- makeBootCIs(bootPreds[predictData$treatment=="ON",])

cisoff <- makeBootCIs(bootPreds[predictData$treatment=="OFF",])

## Step 21: visualising boot CI
par(mfrow=c(2,2), mar=c(3,4,3,3)) 

#low OFFs
quilt.plot(predictData$x.pos[predictData$treatment=="OFF"],
           predictData$y.pos[predictData$treatment=="OFF"],
           cisoff[,1], asp=1, nrow=dims[1], ncol=dims[2], 
           zlim=c(0, maxlim = 0.4), main = " Lower CI  - TAST OFF", add.legend = F)
par(mar=c(3,1,3,6))
# up OFF
quilt.plot(predictData$x.pos[predictData$treatment=="OFF"],
           predictData$y.pos[predictData$treatment=="OFF"], cisoff[,2],
           asp=1,nrow=dims[1], ncol=dims[2], zlim=c(0, maxlim = 0.4), 
           main = "Upper CI - TAST OFF")
par(mar=c(3,4,3,3))
#low ONs
quilt.plot(predictData$x.pos[predictData$treatment=="ON"],
           predictData$y.pos[predictData$treatment=="ON"],
           cison[,1], asp=1, nrow=dims[1], ncol=dims[2], 
           zlim=c(0, maxlim = 0.4), main = " Lower CI  - TAST ON", add.legend = F)
par(mar=c(3,1,3,6))
# up ON
quilt.plot(predictData$x.pos[predictData$treatment=="ON"],
           predictData$y.pos[predictData$treatment=="ON"], cison[,2],
           asp=1,nrow=dims[1], ncol=dims[2], zlim=c(0, maxlim = 0.4), 
           main = "Upper CI - TAST ON")



## Step 22: Identifying Differences 
differences <- getDifferences(beforePreds = bootPreds[predictData$treatment=="OFF", ],
                              afterPreds = bootPreds[predictData$treatment=="ON", ])


bootPreds$id <- as.factor(df$id)

# get average prediction for each grid id  for on 
boot_off <- bootPreds[predictData$treatment=="OFF",]

# group by grid id

# take mean of each column of each group

# store in a new df with single row for each grid id and row for each boot predict


# get average prediction for each id for off 
bootPreds[predictData$treatment=="ON", ]




for(grid in 1:length(df$id)){

}



# Step 23: Visualising differences
# The median for each after - before difference
medianiff <- differences$meandiff

# The marker for each after - before difference:
# positive ('1') and negative ('-1') significant differences
marker <- differences$significanceMarker

#plot positive differences
par(mfrow = c(1, 1), mar = c(3, 6, 3, 3))
quilt.plot(predictData$x.pos[predictData$treatment=="OFF"],
           predictData$y.pos[predictData$treatment=="OFF"],
           mediandiff, asp = 1, nrow = 7, ncol = 9, add.legend = T)
# add + depending on significance of cells. Just
# requires one significance out of all to be allocated
points(predictData$x.pos[predictData$treatment=="OFF"][marker==1],
       predictData$y.pos[predictData$treatment=="OFF"][marker==1], pch="+",
       col="black", cex=1)
# location of observer
points(0,0,pch= "*", col="green", cex=3)
# 
# #plot negative differences
# quilt.plot(predictData$x.pos[predictData$treatment=="OFF"],
#            predictData$y.pos[predictData$treatment=="OFF"],
#            mediandiff, asp = 1, nrow = 7, ncol = 9)
# add - depending on significance of cells. Just
# requires one significance out of all to be allocated
points(predictData$x.pos[predictData$treatment=="OFF"][marker==(-1)],
       predictData$y.pos[predictData$treatment=="OFF"][marker==(-1)], col="white",
       cex=1)

points(predictData$x.pos[predictData$treatment=="OFF"][marker==(0)],
       predictData$y.pos[predictData$treatment=="OFF"][marker==(0)], col="black",
       cex=.5) 
# location of observer
points(0,0,pch= "*", col="green", cex=3)

## interesting -- how can some squares have both significant increases and significant decreases... 
## see if felix can help>







## Map model for On and Off
par(mfrow=c(1,2))
# Loop across all phases
for (treat in c("ON", "OFF"))
{
  bWant <- df$treat %in% treat # a boolean specifying phase
  fields::quilt.plot(df$x.pos[bWant],
                     df$y.pos[bWant],
                     fitted(geeModel)[bWant],
                     
                     nrow=10, ncol=10,
                     zlim=range(fitted(geeModel)),
                     main=paste0("TAST ", treat))
  if (treat %in% "ON")
  #somethign about turn off legend for one of the
  }

}

#### Predictions 
predict.gamMRSea(geeModel)
############
