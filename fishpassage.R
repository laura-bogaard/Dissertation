## Visualizing fish passage data 
## LAURA BOGAARD
library(readr)
fishpassageTAST <- read_csv("fishpassageTAST.csv", col_names = T)
View(fishpassageTAST)
library(tidyverse)
library(dplyr)
library(reshape2)

## Select data for on and off

fish <- data_frame(fishpassageTAST)

chin <- data.frame(fish$chinook_day, fish$date, fish$treatment=="on") 

ggplot(fish, aes(x= chin, color=treatment)) +
  geom_histogram(fill="white", alpha=0.5, position="identity")

#Create a histogram 

ggplot(chinook_dat, aes(x=date, y = chinook_day, color = treatment)) +
  scale_x_date(date_labels = "%U")

hist(fish$chinook_day, 
        xlab = "Date", 
        names.arg = fish$date, 
        main=NULL, 
        breaks = fish$date, 
        col = fish$treatment,
        legend = T)


hist(fish$coho_day, 
        add = T,
        xlab = "Date", 
        names.arg = fish$date, 
        main=NULL, 
        breaks = fish$date, 
        col = fish$treatment,
        legend = T)


hist(y$Length, probability = TRUE, add = TRUE, breaks=seq(0,90,5), 
     col = rgb(0, 0, 1, 0.5))


?barplot

geom_histogram(fill = grey, alpha=0.5)

#Chinook by treatment
fish %>% 
  mutate(time = as.Date(date, format = "%d-%b")) %>% 
  ggplot(mapping = aes(x = time, y = chinook_day, fill = treatment)) + geom_bar(stat = "identity")
#coho by treatment
fish %>% 
  mutate(time = as.Date(date, format = "%d-%b")) %>% 
  ggplot(mapping = aes(x = time, y = coho_day, fill = treatment)) + geom_bar(stat = "identity")
# chinook and coho by treatment =======
#version 1
fish %>% 
  mutate(time = as.Date(date, format = "%d-%b")) %>% 
  melt(id.vars = c("time", "treatment"), measure.vars = c("coho_day", "chinook_day")) %>% 
  ggplot(aes(x = time, y = value, fill = treatment, alpha = variable)) + 
  geom_col() + 
  scale_alpha_discrete(name = "Fish Species", range = c(0.4, 1), labels = c("Coho", "Chinook")) + 
  scale_fill_discrete(name = "Treatment", labels = c("On", "Off")) +
  theme_classic() +
  ylab("Fish count") + 
  xlab("Date") 
# version 2
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

#boxplot by species
ggplot(data = fish, mapping = aes(y = 2020, fill = treatment)) + geom_boxplot()
ggplot(data = fish, mapping = aes(y = coho_day, fill = treatment)) + geom_boxplot()
#Combined boxplot coho chinook
fish %>% 
  mutate(total_fish = chinook_day + coho_day) %>% 
  ggplot(mapping = aes(y = total_fish, x = treatment, fill = treatment, alpha = 0.5)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("blue", "darkgreen"), name = "Treatment", labels = c("On", "Off")) +
  ylab("Total Fish Passsage (Coho and Chinook)") +
  xlab("Treatment") +
  theme_classic()+
  theme(legend.position="none")


  
        