### Binomial Presence Absence models

# read in data
padat <- read.csv("~/Desktop/TASTBALLARD/TAST WD/Dissertation/pa_data.csv")

# load libraries
library(tidyverse)
#########
# subset/rename data
pa <- as.factor(padat$pa) #binary presence absensce
fish <- padat$fish #fish count total / day
obs <- as.factor(padat$obs) #observer
treatment <- as.factor(padat$treatment) # tast on or off
x <- padat$x
y <- padat$y
dist <- padat$dist2sq #distance to square centre
rl <- padat$rl_emp_sq #emperically modelled rl for center of each square
day <- padat$jul_day #julian day fromstatrt
hour <- padat$start_hr #survey start hour
#########
df <- padat
df$pa <- as.factor(df$pa)
df$treat <- as.factor(df$treatment)
df$obs <- as.factor(df$obs)
df$XPos <- df$x
df$YPos <- df$y
df$x.pos <- df$x
df$y.pos <- df$y

str(df)
#######
#EDA
qplot(treatment, pa)
qplot(rl, pa)


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

# pa ~ treatment + s(x, y) + s(rl) + fish + obs + 
#   dist + day + hour, method = "REML", family = binomial)

# try SALSA
library(devtools)
devtools::install_github(repo="lindesaysh/MRSea@v1.01")

# Create the variable "response" needed by SALSA
df$response <- df$pa

# Set initial model without spline-based terms
initialModel <- glm(response ~ treat, family=binomial, data=df)

# Set SALSA arguments
factorList <- c("obs", "treat")
varList <- c("dist2sq", "XPos", "YPos", "rl_emp_sq", "fish", "jul_day", "start_hr")
salsa1DList <- list(fitnessMeasure="AIC", minKnots_1d=rep(1, 7),
                    maxKnots_1d=rep(5, 7), startKnots_1d=rep(1, 7),
                    degree=rep(2, 7), maxIterations=100,
                    gaps=rep(0, 7))
str(df)
# Load library
library(MRSea)
set.seed(53195)
# Run SALSA
salsa <- MRSea::runSALSA1D(initialModel, salsa1DList, varList,
                           factorList, varlist_cyclicSplines=NULL,
                           splineParams=NULL, datain=df,
                           suppress.printout=TRUE, removal=FALSE,
                           panelid=NULL)

# Pick best model based on fitnessMeasure
bestModel1D <- salsa$bestModel

# Select initial knot locations using a space-filling design
knotGrid <- MRSea::getKnotgrid(coordData=cbind(df$Xpos, df$Ypos),
                               numKnots=300, plot=FALSE)

# Create distance matrix
# (i.e. distance between each data points and knots)
distMatrix <- MRSea::makeDists(datacoords=cbind(df$x.pos, df$y.pos),
                               knotcoords=knotGrid)
# Set SALSA parameters # no interaction
salsa2DList <- list(fitnessMeasure="AIC", knotgrid=knotGrid,
                    startKnots=6, minKnots=2, maxKnots=20,
                    gap=0)
# Run SALSA 2D
salsa2D <- MRSea::runSALSA2D(bestModel1D, salsa2DList, d2k=distMatrix$dataDist,
                             k2k=distMatrix$knotDist, splineParams=NULL,
                             tol=0, chooserad=F, panels=NULL,
                             suppress.printout=TRUE)
