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
df <- padat
df$pa <- as.numeric(df$pa)
df$treat <- as.factor(df$treatment)
df$obs <- as.factor(df$obs)
df$XPos <- df$x
df$YPos <- df$y
df$x.pos <- df$x
df$y.pos <- df$y

## add numeric ID for each experimental unit (session_id) to cover temporal AC
## in GEE model (potentially)
seshlevels <-  levels(as.factor(df$session_id)) 
numericsesh <- data.frame(session_id = levels(as.factor(df$session_id)),
                          block_as_session = as.factor(seq(1, length(seshlevels), 1)))

df <- left_join(df, numericsesh, by = "session_id" )
arrange(df, session_id)

#######
#EDA
qplot(df$treat, df$pa)

# start with full linear model
m0 <- glm(pa ~ treatment + x.pos + y.pos + fish + obs + dist2sq + jul_day + start_hr, data = df, family = binomial)
m0
summary(m0)
# check for collinearity
vif(m0) # kick out dist to sq

m1 <- glm(pa ~ treatment + x.pos + y.pos + fish + obs + jul_day + start_hr, data = df, family = binomial)
summary(m1)
# check for collinearity
vif(m1) #  better

# add spatial interaction to account for spatial autocorrelation
m2<- gam(pa ~ treatment + s(x.pos, y.pos) + fish + obs + jul_day + start_hr, data = df, family = binomial)
summary(m2)

## STEP 3: check for correlated residuals in temporal covariates
runs.test(residuals(m2, type = "pearson", alternative = "two-sided")) 

plotRunsProfile(m2, varlist = c("x.pos", "y.pos", "jul_day", "start_hr")) 
#definite residual correlation

# ACF plot with blocks as survey sesh + gr id
df$blockid <- as.factor(paste(df$block_as_session, df$id,sep = ""))
#WARNING THIS RUNS FOREVER
runACF(df$blockid, m2, store = F) #very strange output just gives straight red line 

#Try grid id as a block
df$id <- as.factor(df$id)
runACF(df$id, m2, store = F) #interesting.. this is the one

#try survey sesh as block
df$session_id <- as.factor(df$session_id) #
runACF(df$session_id, m2, store = F) 


# backwards stepwise selection


#########
# Fit GEE CReSS model with SALSA2D 

### STEP 1: check for collinearity between covariates in glm model
fullmod <- glm(pa ~ treatment + jul_day + start_hr + obs + dist2sq + fish + x.pos + y.pos,
               na.action = "na.fail", family = "binomial", data = df)
# check for collinearity
car::vif(fullmod) # most between 1 & 2.1 so pretty much good! dist2Sq and x/y pretty collin so throw out dist

# refit without distance
fullmod_linear <- glm(pa ~ treatment + jul_day + start_hr + obs + fish + x.pos + y.pos,
               na.action = "na.fail", family = "binomial", data = df)

# check for collinearity
car::vif(fullmod) # much better, fish and day only slightly collinear but we can keep it. (1-2)

# STEP 2: make splineParams object with knots at the mean of each variable
splineParams<-makesplineParams(data=df, varlist=c("x.pos", "y.pos", "fish", "jul_day", "start_hr"))

# fit new full mod with splines 
fullmod <- glm(pa ~ treat + obs + bs(x.pos, knots = splineParams[[2]]$knots) 
                          + bs(y.pos, knots = splineParams[[3]]$knots)
                          + bs(fish,knots = splineParams[[4]]$knots)
                          + bs(jul_day,knots = splineParams[[5]]$knots)
                          + bs(start_hr, knots =splineParams[[6]]$knots),
                          na.action = "na.fail", family = "binomial", data = df)

# STEP 3: check for correlated residuals in temporal covariates
varList <- c("XPos", "YPos", "fish", "jul_day", "start_hr")
runs.test(residuals(fullmod, type = "pearson", alternative = "two-sided")) # neg standardised runs stat and sm p == houston we have correlation
plotRunsProfile(fullmod, varlist = c("jul_day", "start_hr")) #definite correlation

# ACF plot with blocks as survey sesh
df$session_id <- factor(df$session_id)
runACF(df$session_id, fullmod, store = F) #some positive AC, use this as our cor block

# STEP 4: plot cumulative residuals for model to check over or under prediction
plotCumRes(fullmod, varlist= varList, splineParams) 

# compare with model that just has linear terms
plotCumRes(fullmod_linear, varlist = varList)

#looks like x and y are linear and all the rest shouldnt be
# model defo over and under predicting at times, lets run SALSA for more flexibility

# STEP 5: Run SALSA
library(devtools)
devtools::install_github(repo="lindesaysh/MRSea@v1.01")

# Create the variable "response" needed by SALSA
df$response <- df$pa

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
# Run SALSA
salsa1D <- MRSea::runSALSA1D(initialModel, salsa1DList, varList,
                           factorList, varlist_cyclicSplines=NULL,
                           splineParams=NULL, datain=df,
                           suppress.printout=TRUE, removal=FALSE,
                           panelid=NULL)

# STEP 6: Pick best model based on AIC
bestModel1D <- salsa1D$bestModel

# salsa took the initial model with one knot at mean for each var and added 
# 2 knots for fish at 164 and 728; and 3 knots for start hour at 10, 13, 15;
# The result of SALSA is additional flexibility added to the model for the
# relationship between presence absence probability and fish / hr

# STEP 7: Select initial knot locations using a space-filling design

# Next we add a two dimensional smooth of geographic coordinates (s(x.pos,
# y.pos))# To test for a redistribution of animals with a treatment effect we fit an
# interaction term between the smooth of coordinates and impact.
# SALSA is used to determine spatially adaptive knot locations for this
# smooth term.

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
splineParams <- salsa1D$splineParams

# STEP 10: RUN SALSA 2D
set.seed(53195)
salsa2D <- MRSea::runSALSA2D(bestModel1D, salsa2DList, d2k=distMats$dataDist,
                             k2k=distMats$knotDist, splineParams=NULL,
                             tol=0, chooserad=F, panels=NULL,
                             suppress.printout=TRUE)

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


# Store "best" model
bestModel2D <- salsa2D$bestModel
baseModel <- salsa2D$bestModel

# Summary of results
summary(bestModel2D)
# save the model object
stats::anova(bestModel2D, test="Chisq")

# update spline parameter object
splineParams <- salsa2D$splineParams

# STEP 11: check for residual correlation
runs.test(residuals(bestModel2D, type = "pearson")) #still neg stat, still tiny p

## Significant positive resid correlation therefore refit as GEE
########
# N.B. for the GEE formula, the data must be ordered by block (which this is)
# and the blockid must be numeric
# specify parameters for local radial:

df <- arrange(df, block_as_session)
radiusIndices <- splineParams[[1]]$radiusIndices
dists <- splineParams[[1]]$dist
radii <- splineParams[[1]]$radii
aR <- splineParams[[1]]$invInd[splineParams[[1]]$knotPos] ## this comes up as null, so had to workaround aR not within SP, in larger model object
aR <- as.vector(salsa2D$aR$aR1) 
# update model in workspace with parameters for spatial smooth (above)
baseModel <- update(baseModel, . ~ .)

# Re-fit the chosen model as a GEE (based on SALSA knot placement) and
# GEE p-values
geeModel <- geepack::geeglm(formula(baseModel), data = df, family = "binomial",
                   id = block_as_session)
getPvalues(geeModel) #treatment, lrf, and treat:lrf not significant
summary(geeModel)

# how to remove impact
noint.model<-update(geeModel, .~. - as.factor(treatment):
                      LocalRadialFunction(radiusIndices, dists, radii, aR))
# reshow p-values
getPvalues(noint.model)



#######
############# Fit with interaction between x, y, treat
# STEP 9: Set SALSA parameters 
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

# STEP 11: check for residual correlation
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
                            id = block_as_session)
getPvalues(geeModel, varlist = varList, factorlist = factorList)



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

