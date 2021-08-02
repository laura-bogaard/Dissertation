rm(list = ls())
## TAST MSC MULTINOMIAL GAM SCRIPT
## Laura Bogaard

# Question: what is the prob of being in a block?
            # does this change accross treatment?
# load libraries
library(tidyverse)
library(reshape2)
library(readr)
library(useful)
library(nlme)
library(mgcv)
library(plyr)
library(lubridate)

#set wd 
setwd("~/Desktop/TASTBALLARD/TAST WD/Dissertation")

# load data and adjust data types
#### distance refers to the dataset where each row is a surfacing
distance <- read.csv("Ballard_distance+fishcts.csv", header = T)
head(distance)
#### pred refers to the dataset where each row is a surveying session
pred <- read.csv("Ballard_Predation.csv")
head(pred)

### missing a key piece of the puzzle here: TIDE HEIGHT AFFECTS DISTANCE
tide <- read.csv("tides_ballard.csv", header = T)
head(tide)
str(tide)

# convert to date format and creat keyfor each one in tide
tide$Month <- month(ymd(tide$Date)) 
tide$Day <- day(ymd(tide$Date))
tide$Hour <- hour(hms(tide$Time))
tide$Minute <- minute(hms(tide$Time))
tide$MDHM <- paste(tide$Month,tide$Day, tide$Hour, tide$Minute, sep = "_")

# convert to date format and creat keyfor each one in distance
distance$Month <- month(mdy(distance$date)) 
distance$Day <- day(mdy(distance$date))
distance$Hour <- hour(hm(distance$time))
distance$Minute <- minute(hm(distance$time))
distance$MDHM <- paste(distance$Month,distance$Day, distance$Hour, distance$Minute, sep = "_")

# add tide height to distance
distance <- left_join(distance, tide, by = "MDHM")
distance$tide_ht <- distance$Pred

# correct for platform height above MLLW c^2 - b^2 = a^2
distance$tide_dist_sq <- (distance$platform_dist)^2 - ((6.78 + distance$tide_ht)^2) 

# get rid of NAs
na.omit(distance$tide_dist_sq)

#take sqrt and round boooooom!
distance$tide_dist <- round(sqrt(abs(distance$tide_dist_sq)), 0)

###

# rename sesh id so there is a key variable between distance and pred
pred$session_id <- pred$Session_ID

# first create new column for block # and N/S
distance <- unite(distance, col = "blockid", block:block_ns, sep = "", remove = FALSE, na.rm = FALSE)
head(distance)

# now filter for all distances less than truncated 250m
distance <- distance%>% filter(tide_dist <= 250)
head(distance)

#check that we are now limited to just 1-3
which(distance$block == 4) # some from four need to be removed

# remove these rows
distance <- distance%>% filter(block < 4)
head(distance)

#make treatment a factor
distance$treatment <- factor(distance$treatment)

#make foraging binary
distance$foraging <- ifelse(distance$foraging == "Yes", 1, 0) 

# correct date format
distance$date <- as.Date(distance$date,"%m/%d/%y")

#get julian day and survey effort onto distance data
distance<- pred %>%
  select(session_id, dur_hr, Julian_day_from_start) %>%
  left_join(distance, pred, by = "session_id")

# remove pilot surveys
distance <- slice(distance, 6:n())

# fix julian day to account for this, now call it jul_day
distance$jul_day <- distance$Julian_day_from_start - 18

#fix fish_day
distance$fish <- distance$fish_day

#turn block id into numeric factor dummy variables for multinom model
distance$block_fac <- as.numeric(mapvalues(distance$blockid, warn_missing = T, 
                                            from = c("1S", "1N", "2S", "2N", "3S", "3N"), 
                                            to = c("0", "1", "2", "3", "4", "5")))
# get rid of any rows with NA
na.omit(distance$block_fac)

# should i just seals, kick out sea lions? 
which(distance$species == "SL") #127 observations with SL--keep or toss??

max(distance$tide_ht)
min(distance$tide_ht)

# rescue mission for distance data: 
# figure out if any distances do not line up with blocks ⬇︎
################################################################################
# add a row id 
distance$ID <- seq.int(nrow(distance))
#date and time giving grief drop them for now 
distance <- distance %>%
                select( !date) 
 
block1S <- c(0, 50)
# check if any values fall outside distance range for this block
b1s <- filter(distance, blockid == "1S") 
bad1s <- filter(b1s, tide_dist < block1S[1] | tide_dist > block1S[2])
#okay

block1N <- c(33, 78)
# check if any values fall outside distance range for this block
b1n <- filter(distance, blockid == "1N") 
bad1n <- filter(b1n, tide_dist < block1N[1] | tide_dist > block1N[2])
#okay

block2S <- c(37, 112)
# check if any values fall outside distance range for this block
b2s <- filter(distance, blockid == "2S") 
bad2s <- filter(b2s, tide_dist < block2S[1] | tide_dist > block2S[2])
#okay

block2N <- c(51, 129)
# check if any values fall outside distance range for this block
b2n <- filter(distance, blockid == "2N") 
bad2n <- filter(b2n, tide_dist < block2N[1] | tide_dist > block2N[2])
#okay

block3S <- c(108, 250)
# check if any values fall outside distance range for this block
b3s <- filter(distance, blockid == "3S") 
bad3s <- filter(b3s, tide_dist < block3S[1] | tide_dist > block3S[2])
#okay

block3N <- c(112, 255)
# check if any values fall outside distance range for this block
b3n <- filter(distance, blockid == "3N") 
bad3n <- filter(b3n, tide_dist < block3N[1] | tide_dist > block3N[2])
#okay

#bind them all together ** note that rbind h8s d8s :( **
bad_dists <- rbind(bad1s, bad1n, bad2s, bad2n, bad3s, bad3n)
nrow(bad_dists) # cool only drop 159 observations

#now drop them from D3
distance <- anti_join(distance, bad_dists, by = "ID")
str(distance)
# just a note that this took me 3 hours... but i learned things!!

#final clean up to drop some lingering extra varaibles
distance$hour <- distance$Hour.x

# get rid of all observations made from the gate
dist_no_gate <- filter(distance, observer_location != "Gate") 
## use this later with following 3 lines if we decide to drop all gate obs

# or just those from day 42 that looked wrong on the map
lb42 <- filter(distance, jul_day == 42) # this is the day
LB42 <- filter(lb42, observer_location == "Gate") #this is where the bad bearings come from 
byelau <- anti_join(distance, LB42, by = "MDHM") # drop em like its hot
drop <- which(byelau$MDHM == "9_5_10_28" & byelau$tide_dist == 220) #this one is defo an accident
byelau <- filter(byelau, ID != drop) #bye
#another outlier
ff <- which(byelau$x > 150 & byelau$y > 120) 

###### now clean bearings using Katie's ballard measurements
## TODO validate these with protractor

ang_1s <- c(315, 15)
# check if any values fall outside distance range for this block
ang_b1s <- filter(distance, blockid == "1S") 
bad_1s <- filter(ang_b1s, platform_bearing < block1S[1] & platform_bearing > block1S[2])
#okay none 

ang_1n <- c(180, 110)
# check if any values fall outside distance range for this block
ang_b1s <- filter(distance, blockid == "1S") 
bad_1s <- filter(ang_b1s, platform_bearing < block1S[1] & platform_bearing > block1S[2])
#okay


clean_dist2 <- byelau %>% 
  select(session_id, MDHM, Date, jul_day, Time, hour, dur_hr, treatment, species, 
         num_idv, obs, foraging, blockid, block, block_ns, platform_bearing, 
         platform_distance, tide_dist, tide_ht, fish, block_fac, x, y)

##### now rescue mission for bearings if this was the issue all along... i will freaking scream
clean_dist2$platform_bearing <- clean_dist2$platform_bearing + 20

#correct for values over 360
clean_dist2$platform_bearing <- if_else(clean_dist2$platform_bearing >= 360,
                                        true = clean_dist2$platform_bearing -360,
                                          false = clean_dist2$platform_bearing)
#check
which(clean_dist2$platform_bearing >=360)

### now do calculation to convert distance and bearing (polar coord) to x/y (cartesian coord)

# first convert bearings to radians
theta1 <- (clean_dist2$platform_bearing * (pi/180))

# next convert polar to cartesian
xy_coord <- tibble(pol2cart(clean_dist2$platform_distance, theta1, degrees = F))


# plot to check
plot(xy_coord$x, xy_coord$y)

# bind columns to valid bearing df
clean_dist2$x <- xy_coord$x
clean_dist2$y <- -(xy_coord$y)
head(clean_dist2)

################################################################################
 
### now I'll calculate rl for centre of each block 
# RL refers to the dataset with distances and calculated RLs
RL <- data.frame(read_csv("ballard_rl.csv"))

#check                
head(RL)
tail(RL)
RL[19,]

#weird, remove last 6 rows of NAs
RL <- slice(RL, 1:(n()-6))

#rename distance columns
rldist <- as.numeric(RL$platform_distance) #distance from recording location to platform
location <- factor(RL$location) #whether recording was taken at edge of canal or centre
mean_rl <- as.numeric(RL$mean_rl) #mean received level in dB re 1uPa, assuming spherical spreading

# build simple linear model
rlmod <- lm(mean_rl ~ I(20*log10(rldist)), data = RL) 
summary(rlmod)
rlmod$coefficients
summary(rlmod)$adj.r.squared

#add theo and emperical to RL
RL$RL_theoretical <- (181 - (20*log10(rldist)))
RL$RL_emperical <- c(predict(rlmod, newdata = data.frame(rldist = rldist)))

# now do this in gg plot
RL_plot <- ggplot(data = RL) +
  geom_point(mapping = aes(x = rldist, y = mean_rl, color = location)) +
  stat_smooth(mapping = aes(x = rldist, y = mean_rl), method = "lm", 
              formula = y ~ I(20*log(x)), se=T, 
              level=0.95, color = "blue") +
  geom_smooth(aes(x = rldist, y=RL_theoretical), colour="green", method = "loess", se = T) +
  guides()

# now calculate distance to the centre of each block
blockid <- c("1S", "1N", "2S", "2N", "3S", "3N")
dist2midblk <- c(mean(block1S), mean(block1N), mean(block2S), mean(block2N), mean(block3S), mean(block3N))
#model theorhetical and empirical TL for centre distance of each block
rl_blk_mid_emp <- c(predict(rlmod, newdata = data.frame(rldist = dist2midblk)))
rl_blk_mid_the <- as.numeric(182 - (20*log10(dist2midblk)))

block_rl <- data.frame(blockid = blockid, 
                       dist2midblk = as.numeric(dist2midblk),
                       rl_blk_mid_emp = as.numeric(rlmid_emp), 
                       rl_blk_mid_the = as.numeric(rlmid_the))
str(block_rl)

# join the mean blockrl to distance
clean_dist2 <-  left_join(block_rl, clean_dist2, by = "blockid")

# do the same for every tide corrected distance
tide_dist_df <- data.frame(clean_dist2$tide_dist)
clean_dist2$idv_rl <- predict(rlmod, newdata = tide_dist_df, type = "response")
#####ask fe for help

clean_dist <- clean_dist2
View(clean_dist)
write_csv(clean_dist, "~/Desktop/TASTBALLARD/TAST WD/Dissertation/cleandist.csv")

##### RUN TO HERE FIRST THIGNS FIRST! ###################

################################################################################
#  MULTINOMIAL GAM TIME!
################################################################################
# Start by simplifying the dataset for this particular model
str(clean_dist)
clean_dist$treatment <- ifelse(clean_dist$treatment == "ON", 1, 0)


dataON <- filter(clean_dist, treatment == "ON")
dataOFF <- filter(clean_dist, treatment == "OFF")
hist(dataON$block_fac)
hist(dataOFF$block_fac)
# again, we can see the distance shift in the hists, bot is is significant..?


# start with baby model and build up, forward stepwise
m0 <- gam(list(block_fac ~ treatment, ~ treatment, ~ treatment, ~ treatment, ~ treatment),
          data = clean_dist, method = "REML", family = multinom(K=5)) 
m0
summary(m0)
gam.check(m0)

# add smooth for fish
m1<- gam(list(block_fac ~ treatment + s(fish, k = 25), ~ treatment + s(fish, k = 25), ~ treatment+ s(fish, k = 25), ~ treatment+ s(fish, k = 25), ~ treatment+ s(fish, k = 25)),
          data = clean_dist, method = "REML", family = multinom(K=5)) 
m1
summary(m1)
gam.check(m1)

# fish as linear
m2<- gam(list(block_fac ~ treatment + fish, ~ treatment + fish, ~ treatment+ fish, ~ treatment+ fish, ~ treatment+ fish),
         data = clean_dist, method = "REML", family = multinom(K=5)) 
m2
summary(m2)
gam.check(m2)

# add rl ## did not converge
m3<- gam(list(block_fac ~ treatment + rlmid_emp,
              ~ treatment +  rlmid_emp, ~ treatment + rlmid_emp,
              ~ treatment + rlmid_emp, ~ treatment +  rlmid_emp),
         data = clean_dist, method = "REML", family = multinom(K=5))
m3
summary(m3)
gam.check(m3)

# rl as smooth ## did not converge
m4<- gam(list(block_fac ~ treatment + s(rlmid_emp, k = 2),
              ~ treatment + s(rlmid_emp, k = 2), ~ treatment + s(rlmid_emp, k = 2),
              ~ treatment + s(rlmid_emp, k = 2), ~ treatment  + s(rlmid_emp, k = 2)),
         data = clean_dist, method = "REML", family = multinom(K=5))
m4
summary(m4)
gam.check(m4)

# pay so fish and rlmid are shoddy with the smooths, lets try a temporal smooth
m5 <- gam(list(block_fac ~ treatment +  hour + s(jul_day), 
                ~ treatment + hour + s(jul_day),
                ~ treatment + hour + s(jul_day), 
                ~ treatment + hour + s(jul_day), 
                ~ treatment + hour + s(jul_day)), 
           data = clean_dist,  method = "REML", family = multinom(K = 5))

summary(m5) ## hour not significant, k too low
gam.check(m5)
### 
m6 <- gam(list(block_fac ~ treatment + s(jul_day), 
               ~ treatment + s(jul_day),
               ~ treatment + s(jul_day), 
               ~ treatment + s(jul_day), 
               ~ treatment + s(jul_day)), 
          data = clean_dist,  method = "REML", family = multinom(K = 5))

summary(m6)
gam.check(m6)## smooth significant, but k too low
# worked, TODO play around with K

## hour of day really not significant
m7 <- gam(list(block_fac ~ treatment + hour, 
               ~ treatment + hour,
               ~ treatment + hour, 
               ~ treatment + hour, 
               ~ treatment + hour), 
          data = clean_dist,  method = "REML", family = multinom(K = 5))

summary(m7)
gam.check(m7)
### linear model with day, mostly significant except for block 
m8 <- gam(list(block_fac ~ treatment + jul_day, 
               ~ treatment + jul_day,
               ~ treatment + jul_day, 
               ~ treatment + jul_day, 
               ~ treatment + jul_day), 
          data = clean_dist,  method = "REML", family = multinom(K = 5))

summary(m8)
gam.check(m8)

plogis(predict.gam(m8, type = "response", newdata = dataON))


dataON <- filter(clean_dist, treatment == "ON")
?filter
# full model with smooth for fish_day + rl ## model wont converge with smooth for day or fish
m23 <- gam(list(block_fac ~ treatment + s(fish) + rlmid_emp + hour + jul_day, 
                ~ treatment + s(fish) + rlmid_emp+ hour + jul_day,
                ~ treatment + s(fish) + rlmid_emp+ hour + jul_day, 
                ~ treatment + s(fish) + rlmid_emp+ hour + jul_day, 
                ~ treatment + s(fish) + rlmid_emp+ hour + jul_day), 
           data = clean_dist,  method = "REML", family = multinom(K = 5))

summary(m23)
gam.check(m23)

# full model with smooth for rlmid_emp
m24 <- gam(list(block_fac ~ treatment + fish + s(rlmid_emp) + hour + s(jul_day), 
                ~ treatment + fish + s(rlmid_emp)+ hour + s(jul_day),
                ~ treatment + fish + s(rlmid_emp)+ hour + s(jul_day), 
                ~ treatment + fish + s(rlmid_emp)+ hour + s(jul_day), 
                ~ treatment + fish + s(rlmid_emp)+ hour + s(jul_day), 
                data = clean_dist,  method = "REML", family = multinom(K = 5))) 
summary(m24)
gam.check(m24)

#full model w smooth for fish
# full model with smooth for fish + rlmid_emp
m25 <- gam(list(block_fac ~ treatment + s(fish) + rlmid_emp + hour + jul_day, 
                ~ treatment + s(fish) + rlmid_emp + hour + jul_day,
                ~ treatment + s(fish) + rlmid_emp + hour + jul_day, 
                ~ treatment + s(fish) + rlmid_emp + hour + jul_day, 
                ~ treatment + s(fish) + rlmid_emp + hour + jul_day), 
           data = data_blockmod,  method = "REML", link = logit, family = multinom(K = 5))
summary(m25)
gam.check(m25)

AIC(m23, m24, m25) # best model is with smooth rlmid_emp

# start with base model and build up
m26 <- gam(list(block_fac ~ treatment + fish + s(rlmid_emp), 
                ~ treatment + fish + s(rlmid_emp),
                ~ treatment + fish + s(rlmid_emp), 
                ~ treatment + fish + s(rlmid_emp), 
                ~ treatment + fish + s(rlmid_emp)), 
           data = clean_dist,  method = "REML", 
           family = multinom(K = 5))

summary(m26)
# drop fish, try hour
m27 <- gam(list(block_fac ~ treatment + s(rlmid_emp) + hour, 
                ~ treatment + s(rlmid_emp)+ hour,
                ~ treatment + s(rlmid_emp)+ hour, 
                ~ treatment + s(rlmid_emp)+ hour, 
                ~ treatment + s(rlmid_emp)+ hour), 
                data = clean_dist,  method = "REML", family = multinom(K = 5))
summary(m27)


# did not converge