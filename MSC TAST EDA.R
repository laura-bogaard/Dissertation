## TAST MSC CLEANING AND EDA SCRIPT
## Laura Bogaard

# Load data and libraries
library(tidyverse)
library(reshape2)
library(readr)
library(useful)
library(nlme)
library(mgcv)

#set wd 
setwd("~/Desktop/TASTBALLARD/TAST WD")

# load data and adjust data types
#### distance refers to the dataset where each row is a surfacing
distance <- read.csv("Ballard_distance+fishcts.csv", header = T)
distance$time <- parse_time(distance$time, format = "%T")
distance$treatment <- factor(distance$treatment)

# Seperate time into hms so we can use hour of day as a factor variable 
distance <- distance %>%
  separate(time, c("hour", "minute","second"),
         sep = ":", remove = TRUE)
distance$hour <- as.factor(distance$hour)
# sess<-dat$session_id 
# date<-as.factor(dat$date) #date as factor
# HOD <- as.factor(distance$)
# treat<-as.factor(dat$treatment)# treatment
# noseals<-(dat$num_idv)# no of individual seals
# obsloc<-(dat$observer_location)#observer location
# species<-as.factor(dat$species)#species
# bear<-dat$platform_bearing #bearing to seal
# dist<-dat$platform_distance #distance to seal, all obs without distance measure excluded, distance zero was set to dist=0.1 (Gamma can only be used for positive non-zero values)
# forag<-dat$foraging #foraging behaviour (re-cided from yes or no)
# crash<-dat$crash  


#### pred refers to the dataset where each row is a surveying session
pred <- read.csv("Ballard_Predation.csv")
  pred$Session_ID <- factor(pred$Session_ID)

pred$Date <- as.Date(pred$Date, format = "%Y-%m-%d")
pred$Start_Time <- parse_time(pred$Start_Time, format = "%T")
pred$End_Time <- parse_time(pred$End_Time, format = "%T")

# Seperate time into hms so we can use hour of day as a factor variable  hour = Start time 
pred <- pred %>%
    separate(Start_Time, c("hour", "minute","second"), 
             sep = ":", remove = F)
pred$hour <- as.factor(pred$hour)  
  
#now "hour" in this dataset refers to the hour within which the survey started as a factor
             
  
#### fish refers to the dataset where each row is a date which corresponds # of fish
fish <- read_csv("fishpassage_surveyonly.csv", 
                 col_types = cols(date = col_date(format = "%m/%d/%y"),
                             treatment = col_factor(levels = c("on", "off"))))


################################################################################
#### RL data EDA + wrangling 
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
plot(rldist, mean_rl)

# build simple linear model
rlmod <- lm(mean_rl ~ I(20*log(rldist)), data = RL) # TODO is it right here to have distance 2x?
summary(rlmod)
rlmod$coefficients
summary(rlmod)$adj.r.squared

# For each value of distance, I can get the value of RL estimated by the model, 
# and add it to the current plot 
RLpredict <- predict(rlmod) 
ix <- sort(rldist,index.return=T)$ix
lines(rldist[ix], RLpredict[ix], col=2, lwd=2)

# now do this in gg plot
RL_plot <- ggplot(data = RL) +
  geom_point(mapping = aes(x = rldist, y = mean_rl, color = location)) +
  stat_smooth(mapping = aes(x = rldist, y = mean_rl), method = "lm", 
              formula = y ~ I(20*log10(x)), se=T, 
               level=0.95)
RL_plot

# looovely can fix this up a bit more later
# now select every "on" and apply this model to it 

#create empty vector
e <- rep(NA, nrow(distance))

# loop through distance and for each treatment = on, calculate the received level 
# using the TL equation assuming spherical spreading 
for(i in 1: nrow(distance)){
  if(distance$treatment[i] == "ON"){
    e[i] <- RL$mean_rl[1] - 20*log(distance$platform_distance[i])}
    else{e[i] <- 0
    }
}
e

#add to distance df
distance$rl <- e
plot(distance$platform_distance, distance$rl)

#plot to check
distON <- distance %>%
  filter(treatment == "ON")
qplot(data = distON, x = platform_distance, y = rl)

###TODO Okay but how do I do this using the model that is built with Asila's measurements


################################################################################
# DISTANCE DATA clean and EDA ##
mean(distance$platform_distance) #98.45
var(distance$platform_distance) #4354.38
# its going to be a negative-binomial rather than a poisson... but isnt this obvious?

# correct date format
distance$date <- as.Date(distance$date,"%m/%d/%y")
###TODO comeback to this when i sort out august bearings

# select for days with invalid bearings
data_dist <- distance %>%
  select(date, hour, minute, rl, session_id, chin_day, sock_day, coho_day, fish_day, treatment, species, num_idv, obs, behaviour, 
         foraging, block, block_ns, platform_bearing, platform_distance) %>%
  filter(date < "2020-09-03")

#select and filter out dates prior to 3 SEP -- select valid bearings only
valid_bearings <- distance %>%
  select(date, hour, minute, rl,session_id, chin_day, sock_day, coho_day, fish_day, treatment, species, num_idv, obs, behaviour, 
         foraging, block, block_ns, platform_bearing, platform_distance) %>%
                    filter(date > "2020-09-02")

# fix data types
valid_bearings$platform_bearing <- as.double(valid_bearings$platform_bearing)
valid_bearings$obs <- as.factor(valid_bearings$obs)

# correct for values over 360
valid_bearings$platform_bearing <- if_else(valid_bearings$platform_bearing >= 360, 
                                           true = valid_bearings$platform_bearing -360, 
                                           false = valid_bearings$platform_bearing)
# check
which(valid_bearings$platform_bearing >=360)

### now do calculation to convert distance and bearing (polar coord) to x/y (cartesian coord)

# first convert bearings to radians
valid_bearings$radians <- (valid_bearings$platform_bearing * (pi/180))

# next convert polar to cartesian
xy_coord <- tibble(pol2cart(valid_bearings$platform_distance, valid_bearings$radians, degrees = F))

# plot to check
plot(xy_coord$x, xy_coord$y)

# bind columns to valid bearing df
valid_bearings$x <- xy_coord$x
valid_bearings$y <- abs(xy_coord$y)
head(valid_bearings)

# cool 
# map these points--just to see if they lay out "correctly"
rad_plot <- ggplot(data = valid_bearings, aes(x = x, y = y)) +
            geom_text(aes(label = block, color = block_ns))+
            scale_y_continuous(breaks=seq(0,250, 10))+
            scale_x_continuous(breaks=seq(0,250, 10))
rad_plot

### TODO the basic idea is there, some points seem flipped about their axis. WHY??


# now try making two plots with treatment as aes
treatment_plot <- ggplot(data = valid_bearings, aes(x = x, y = -y)) +
  geom_point(aes(color = treatment))

treatment_plot 
#interesting, but not the most informative

# Boxplot for all distances
distance %>% 
  ggplot(mapping = aes(y = platform_distance, x = treatment, fill = treatment, alpha = 0.5)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("blue", "darkgreen"), name = "Treatment", labels = c("On", "Off")) +
  ylab("Distance") +
  xlab("Treatment") +
  theme_classic()+
  theme(legend.position="none")

# # plot using polar coords, interesting but not informative
# library(plotrix)
# plt.lns <- filter(valid_bearings, platform_distance <= 250)
# angles <- valid_bearings$platform_bearing
# polar.plot(plt.lns$platform_distance, polar.pos=angles, labels="", rp.type = "s")

# Overlapped Hist for Distance by Treatment (Final Figure for Report)
distance %>% 
  ggplot(aes(x = platform_distance, fill = treatment)) + 
  geom_histogram(binwidth = 10, color = "white", position = "dodge", alpha = 0.5) + 
  #scale_alpha_discrete(name = "Fish Species", range = c(0.4, 1), labels = c("Coho", "Chinook")) + 
  scale_fill_manual(name = "Treatment", values = c("blue", "red"), labels = c("Off", "On")) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  scale_x_continuous(breaks = seq(0, 450, 50)) +
  theme_classic() +
  ylab("Frequency") + 
  xlab("Distance from TAST (m)") 


################################################################################

### PREDATION DATA EDA ###

species <- pred$Phase != "Baseline"

#Boxplot comparing pred rate and salmon species/phase with baseline
# TODO: get rid of baseline somehow 
pred %>%
  ggplot(aes (y = pred_rate_per_hour, x = Phase =! "Baseline", fill = Phase))+
  geom_boxplot()+
  scale_fill_manual(name = "Salmon Species", values = c("white", "lightgrey", "darkgrey"), labels = c("Baseline", "Chinook", "Coho")) +
  theme_classic() +
  ylab("Predation Rate")+
  xlab("Salmon Species")

################## 
# FISH DATA EDA
#chinook by treatment
fish %>% 
  mutate(time = as.Date(date, format = "%d-%b")) %>% 
  ggplot(mapping = aes(x = time, y = chinook_day, fill = treatment)) + geom_bar(stat = "identity")
#coho by treatment-
fish %>% 
  mutate(time = as.Date(date, format = "%d-%b")) %>% 
  ggplot(mapping = aes(x = time, y = coho_day, fill = treatment)) + geom_bar(stat = "identity")

# Chinook and Coho by Treatment (FINAL for Report)
fish %>% 
  mutate(time = as.Date(date, format = "%d-%b")) %>% 
  melt(id.vars = c("time", "treatment"), measure.vars = c("coho_day", "chinook_day")) %>% 
  ggplot(aes(x = time, y = value, fill = variable, alpha = treatment)) + 
  geom_col() + 
  scale_alpha_discrete(name = "Treatment", range = c(1, 0.4), labels = c("On", "Off")) + 
  scale_fill_discrete(name = "Fish Species", labels = c("Coho", "Chinook")) +
  theme_classic() +
  ylab("Fish count") + 
  xlab("Date") 

###################################### 
# WRANGLE DATA INTO PRESENCE ABSENCE

# get counts of surfacings within 10m by 10m bins for each session
library(dplyr)
grid_counts <- valid_bearings %>%  
  mutate(
    cut_x = cut(x, breaks = seq(0, 250, by = 10), include.lowest = T),
    cut_y = cut(y, breaks = seq(0, 250, by = 10), include.lowest = T),
  ) %>%
  count(cut_x, cut_y) %>% mutate(n_bin = n())

# make this a function 
gridcount <- function(data){
  data %>%  
  mutate(
    cut_x = cut(x, breaks = seq(0, 250, by = 10), include.lowest = T),
    cut_y = cut(y, breaks = seq(0, 250, by = 10), include.lowest = T),
    distance = 
  ) %>%
    count(cut_x, cut_y) %>% 
    mutate(n_bin = n())
}

# plot all data
grid_count_plot <- ggplot(valid_bearings, aes(x, y)) + 
  stat_bin_2d(binwidth = 10) + 
  geom_point(color = "white", cex = 0.25)
grid_count_plot


# plot on 
gridcount_ON <- filter(valid_bearings, treatment == "ON")
gridcount_ON <- ggplot(gridcount_ON, aes(x, y)) + 
  stat_bin_2d(binwidth = 10) + 
  geom_point(color = "white", cex = 0.25) +
  xlim(0, 250) +
  labs(title = "TAST ON")

gridcount_ON

#plot off
gridcount_OFF <- filter(valid_bearings, treatment == "OFF")
gridcount_OFF <- ggplot(gridcount_OFF, aes(x, y)) + 
  stat_bin_2d(binwidth = 10, show.legend = F) + 
  geom_point(color = "white", cex = 0.25)+
  labs(title = "TAST OFF")
gridcount_OFF

# plot on and off side by side
require(gridExtra)
grid.arrange(gridcount_OFF, gridcount_ON, ncol = 2)

### TODO: how can I add counts for multi-individual surfacings?????

# rearrange pred data set so binary presence/absence for each 10x10m square 
library(reshape2)

#list all unique Session IDs to loop through for valid bearings
session_IDs <- unique(valid_bearings$session_id)

# create empty vector to fill with loop output
datalist = list()

# This loop reates a unique ID for each session + x + y combo with an observation
for (ID in session_IDs) {
  # subset the particular survey session
  df_sub <- subset(valid_bearings, session_id == ID)
  # for that session, list x y square for each observation
  a <- gridcount(df_sub)
  # create unique ID for each session + x + y  combo
  # b <- expand.grid(session_IDs, unique(a$cut_x), unique(a$cut_y))
  # save in new column of a
  datalist[[ID]] <- paste(a$cut_x, a$cut_y, ID, sep = "_")
}
datalist
# turn list into a data frame
datalisst <- plyr::ldply(datalist, rbind)
#group by session make each row an observation
uniqueIDs <- datalisst %>%
  pivot_longer(!.id, names_to = "count", values_to = "unique_ID")
##YAYYYYYYY

# expand grid to create list of every possible combo of session ID and square
head(grid_counts)

all_combo <- expand.grid(unique(grid_counts$cut_x), unique(grid_counts$cut_y), session_IDs)
head(all_combo)

#create 4th column with unique ID
all_combo$Var4 <- paste(all_combo$Var1, all_combo$Var2, all_combo$Var3, sep = "_")

# Create PA column
#loop through all combo, if value exists in unique ids, put a 1, if not, put a 0
Df1 <- all_combo
Df2 <- data.frame(uniqueIDs)

Df1$PA <- as.integer(Df1$Var4 %in% Df2$unique_ID)
head(Df1)


# add column for treatment
library(stringr)
Df1$treatment <- ifelse(str_sub(Df1$Var3, 1, 1) == "T", paste("ON"), paste("OFF"))

# add column for square ID
Df1$square_ID <- paste(Df1$Var1, Df1$Var2, sep = "_")

# get rid of NA squares??
#Df1 <- na.omit(Df1)

# fix "Var" columns
colnames(Df1) <- c("cut_x", "cut_y", "session_id", "unique_id", "p_a", "treatment", "square_ID")

# add a group size and pred rate column by matching session ID
pred$session_id <- pred$Session_ID
pred3 <- pred %>%
  select(session_id, Seal_Freq, pred_rate_per_hour, Observer1, Date, Julian_day_from_start, hour, dur_hr)
Df1 <- left_join(Df1, pred3, by = "session_id")

# add fish count 
dist2 <- distance%>%
  select(session_id, fish_day)
Df1 <- left_join(Df1, dist2, by = "session_id")

#combine distance and valid bearings
# vb1 <- left_join(distance, valid_bearings, by = "session_id")
# View(vb1) 
## DF1 is a data frame where P/A is recorded for each square during each survey
head(Df1)
Df1$treatment <- as.factor(Df1$treatment)
head(valid_bearings)

# make x column from cut_x so we are dealing with a single value for each square
xsq <- regmatches(Df1$cut_x ,gregexpr("(?<=[\\[\\(,])[0-9.]+(?=[\\]\\),])", Df1$cut_x, perl = TRUE))
xsq <- sapply(xsq, function(r) r[2])
Df1$xsq <- as.numeric(xsq)

# make y column # using regular expressiongs (YUCK) 
ysq <- regmatches(Df1$cut_y ,gregexpr("(?<=[\\[\\(,])[0-9.]+(?=[\\]\\),])", Df1$cut_y, perl = TRUE))
ysq <- sapply(ysq, function(x) x[2])
Df1$ysq <- as.numeric(ysq)
head(Df1)
View(Df1)

# add rl column for each square 

# okay now for every PA = 1 bring in Distance and Bearing rl x y etc
PA1 <- Df1%>%
  filter(p_a == 1)
head(PA1)
PA1 <- distinct(PA1, unique_id, .keep_all = TRUE)

# calculate distance to each square using pythagorean theorem then 
# calculate received level at each square
rl_c <- NULL
for(i in 1:nrow(PA1)){
  c[i] <- sqrt(PA1$xsq[i]^2 + PA1$ysq[i]^2)
  rl_c[i] <- 180 - 20*log(c[i])
}
rl_c
PA1$rl_c <- rl_c
head(PA1)

# add this new rl information to DF1


