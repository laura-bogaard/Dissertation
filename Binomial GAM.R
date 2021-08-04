### Binomial Presence Absence models

# read in data
padat <- read.csv("~/Desktop/TASTBALLARD/TAST WD/Dissertation/pa_data.csv")

# load libraries
library(tidyverse)

# subset/rename data
pa <- factor(padat$pa) #binary presence absensce
fish <- padat$fish #fish count total / day
obs <- factor(padat$obs) #observer
treatment <- factor(padat$treatment) 

# start with full model 
m0 <- gam(pa ~ treatment + s(x, y)data = clean_dist, method = "REML", family = binomial)
m0
summary(m0)
gam.check(m0)