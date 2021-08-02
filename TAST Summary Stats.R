# Cleaning TAST data
library(tidyverse)
tastdata <- read.csv("TAST_R_DATA_07_NOV_LB.csv")
effortdata <- read_csv("TAST EFFORT_NOV_14_LB.csv")
View(TAST_EFFORT_NOV_14_LB)
dist<- tastdata$platform_distance
treatment <- as.factor(tastdata$treatment)
#box plot of distances on vs off
(boxplot(platform_distance ~ treatment, tastdata))


#subset distances as on/off and coerce to numeric vector
dat_on <- subset(tastdata, subset = treatment == "ON", select = platform_distance)
lapply(dat_on, as.numeric)
dat_on <- na.omit(dat_on) #removes every row that contains an NA value
dat_off <- subset(tastdata, subset = treatment == "OFF", select = platform_distance, rm.NA = TRUE)
lapply(dat_off, as.numeric)
dat_off <- na.omit(dat_off)

# find median and mean
median(dat_on$platform_distance) #87
mean(dat_on$platform_distance) # 109.3623
median(dat_off$platform_distance) #64
mean(dat_off$platform_distance) # 83.99595
# boxplot
boxplot(platform_distance ~ treatment, data = tastdata, 
        ylab = "Distance (m)",
        xlab = " ")
##ggplot box plot 
tastdata %>% 
        ggplot(mapping = aes(y = platform_distance, x = treatment, fill = treatment, alpha = 0.5)) + 
        geom_boxplot() +
        scale_fill_manual(values=c("blue", "darkgreen"), name = "Treatment", labels = c("On", "Off")) +
        ylab("Distance") +
        xlab("Treatment") +
        theme_classic()+
        theme(legend.position="none")
#histogram
hist(dat_on$platform_distance, ylim = c(0, 400), 
     main = "Frequency for TAST ON", 
     xlab = "Distance from TAST (m)")
hist(dat_off$platform_distance, ylim = c(0, 400),
     main = "Frequency for TAST OFF",  
     xlab = "Distance from TAST (m)")
#need to standardize distance bins, plot PDF 
plot(density(dat_on$platform_distance))
plot(density(dat_off$platform_distance))

#### Using ggplot to try and get over lapping hist figure for distance
tastdata %>% 
        ggplot(aes(x = platform_distance, fill = treatment)) + 
        geom_histogram(binwidth = 10, color = "white", position = "dodge", alpha = 0.5) + 
        #scale_alpha_discrete(name = "Fish Species", range = c(0.4, 1), labels = c("Coho", "Chinook")) + 
        scale_fill_manual(name = "Treatment", values = c("blue", "red"), labels = c("Off", "On")) +
        scale_y_continuous(breaks = seq(0, 100, 5)) +
        scale_x_continuous(breaks = seq(0, 450, 50)) +
        theme_classic() +
        ylab("Frequency") + 
        xlab("Distance from TAST (m)") 

#Boxplot comparing predrate and salmon species
effortdata %>%
        ggplot(aes (y = seal_pred_rate, x = fish_species, fill = fish_species))+
        geom_boxplot()+
        scale_fill_manual(name = "Salmon Species", values = c("lightgrey", "darkgrey"), labels = c("Chinook", "Coho")) +
        theme_classic() +
        ylab("Predation Rate")+
        xlab("Salmon Species")

# comparing mean predrate for each salmon species
chinook_pred <- as.numeric(subset(effortdata, subset = fish_species == "Chinook", select = seal_pred_rate))
chinook_pred <- as.numeric(na.omit(chinook_pred))
mean(chinook_pred)
