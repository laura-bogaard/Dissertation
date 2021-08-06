### Binomial Presence Absence models

# read in data
padat <- read.csv("~/Desktop/TASTBALLARD/TAST WD/Dissertation/pa_data.csv")

# load libraries
library(tidyverse)
library(car)
library(lawstat)
library(MRSea)
library(mgcv)
library(splines)
library(geepack)
#########
df <- padat
df$pa <- as.numeric(df$pa)
df$treat <- as.factor(df$treatment)
df$obs <- as.factor(df$obs)
df$XPos <- df$x
df$YPos <- df$y
df$x.pos <- df$x
df$y.pos <- df$y

# add numeric ID for each experimental unit (sesh id)
seshlevels <-  levels(as.factor(df$session_id)) 
numericsesh <- data.frame(session_id = levels(as.factor(df$session_id)),
                          block_as_session = as.factor(seq(1, length(seshlevels), 1)))

df <- left_join(df, numericsesh, by = "session_id" )
arrange(df, session_id)

str(df)
#######
#EDA
qplot(df$treat, df$pa)



# start with full model 
m0 <- gam(pa ~ treatment + s(x, y) + s(rl) + s(fish) + obs + s(dist) + day + hour, method = "REML", family = binomial)
m0
summary(m0)
gam.check(m0)
# does not converge

# Try smooths for rl
m1 <- gam(pa ~ treatment + s(x, y) + s(rl) + fish + obs + dist + day + hour, method = "REML", family = binomial)
m1
summary(m1)
gam.check(m1)
# not bad, kick out observer  
m1b <- gam(pa ~ treatment + s(x, y) + s(rl) + fish + dist + day + hour, method = "REML", family = binomial)
m1b
summary(m1b)
gam.check(m1b)
# hm worse

# not bad, kick out observer + dist 
m1c <- gam(pa ~ treatment + s(x, y) + s(rl) + fish + day + hour, method = "REML", family = binomial)
m1c
summary(m1c)
gam.check(m1c)
# hm worse

# not bad, kick out observer + dist 
m1d <- gam(pa ~ treatment + s(x, y) + s(rl) + fish + hour, method = "REML", family = binomial)
m1d
summary(m1d)
gam.check(m1d)
# hm worse


AIC(m1b, m1c, m1d)

#########
# Fit GEE CReSS model with SALSA2D 

### STEP 1: check for collinearity between covariates in glm model
fullmod <- glm(pa ~ treatment + jul_day + start_hr + obs + dist2sq + fish + x.pos + y.pos,
               na.action = "na.fail", family = "binomial", data = df)
# check for collinearity
car::vif(fullmod) # most between 1 & 2.1 so pretty much good! dist2Sq and x/y pretty collin so throw out dist
# refit
fullmod <- glm(pa ~ treatment + jul_day + start_hr + obs + fish + x.pos + y.pos,
               na.action = "na.fail", family = "binomial", data = df)
# check for collinearity
car::vif(fullmod) # much better

# STEP 2: make splineParams object with knots at the mean of each var
splineParams<-makesplineParams(data=df, varlist=c("x.pos", "y.pos", "fish", "jul_day", "start_hr"))
str(splineParams)

# fit new full mod with splines 
fullmod <- glm(pa ~ treat + obs + bs(x.pos, knots = splineParams[[2]]$knots) 
                          + bs(y.pos, knots = splineParams[[3]]$knots)
                          + bs(fish,knots = splineParams[[4]]$knots)
                          + bs(jul_day,knots = splineParams[[5]]$knots)
                          + bs(start_hr, knots =splineParams[[6]]$knots),
                          na.action = "na.fail", family = "binomial", data = df)

# STEP 3: check for correlated residuals
varList <- c("XPos", "YPos", "fish", "jul_day", "start_hr")
lawstat::runs.test(residuals(fullmod, type = "pearson")) # neg stat and sm p == houston we have correlation
plotRunsProfile(fullmod, varlist = c("jul_day")) #definite correlation


# ACF plot with blocks as survey sesh
df$session_id <- factor(df$session_id)
runACF(df$session_id, fullmod, store = F) #pretty dang ac, use this as our cor block

# STEP 4: plot cumulative residuals for model to check over or under prediction
plotCumRes(fullmod, varlist= varList, splineParams) 
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
str(df)
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

# STEP 7: Select initial knot locations using a space-filling design
knotgrid <- MRSea::getKnotgrid(coordData=cbind(df$x.pos, df$y.pos),
                               numKnots=300, plot=FALSE)

# STEP 8: Create distance matrix
# (i.e. distance between each data points and knots)

# make distance matrices for datatoknots and knottoknots

distMats <- makeDists(cbind(df$x.pos, df$y.pos), na.omit(knotgrid))
str(distMats)

# create sequence of radii
r_seq <- getRadiiChoices(numberofradii= 8, distMats$dataDist, basis = "gaussian")

# STEP 9: Set SALSA parameters 
salsa2DList <- list(fitnessMeasure="AIC", knotgrid=knotgrid,
                    startKnots=6, minKnots=2, maxKnots=20,
                    gap=1, r_seq = r_seq)
# splineParams must be an object in workspace
# update splineParams with the SALSA1D results
splineParams <- salsa1D$splineParams
# STEP 10: RUN SALSA 2D
set.seed(53195)
salsa2D <- MRSea::runSALSA2D(bestModel1D, salsa2DList, d2k=distMats$dataDist,
                             k2k=distMats$knotDist, splineParams=NULL,
                             tol=0, chooserad=F, panels=NULL,
                             suppress.printout=TRUE)


## multiple salsa runs with different starting knots ## doesnt work dont know why

# Store "best" model
bestModel2D <- salsa2D$bestModel

# Summary of results
summary(bestModel2D)
# save the model object
baseModel <- salsa2D$bestModel
# update spline parameter object
splineParams <- salsa2D$splineParams

# STEP 11: check for residual correlation
runs.test(residuals(bestModel2D, type = "pearson")) #still neg stat, still tiny p

## Significant positive resid correlation therefore refit as GEE CANT
########
# N.B. for the GEE formula, the data must be ordered by block (which this is)
# and the blockid must be numeric
# specify parameters for local radial:
radiusIndices <- splineParams[[1]]$radiusIndices
dists <- splineParams[[1]]$dist
radii <- splineParams[[1]]$radii
aR <- as.vector(salsa2D$aR$aR1) ## this comes up as null and i dont know why cannot continue 2 days wasted yay. 
# update model in workspace with parameters for spatial smooth (above)
baseModel <- update(baseModel, . ~ .)

# Re-fit the chosen model as a GEE (based on SALSA knot placement) and
# GEE p-values
geeModel <- geepack::geeglm(formula(baseModel), data = df, family = "binomial",
                   id = block_as_session)








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

