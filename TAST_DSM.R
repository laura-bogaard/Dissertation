## This script contains the code for a Density Surface Model
## routine for Ballard Locks TAST data
# Laura Bogaard

# Step 1: load data (this is only the data I have valid bearings for ie post september 2
# will add earlier points with varaince soon)
data <- valid_bearings
  #read.table("Ballard_distance.csv", header = TRUE, sep = ";", stringsAsFactors=FALSE) this is the RAW data
str(data)

# Step 2: package data for DS Function 
DSdata <- data.frame(SurveyID = data$session_id, #this is the id of each survey
                     Xcoord = as.integer(data$x),
                     Ycoord = as.integer(data$y),
                     OBS =  data$obs,
                     Area = , #not sure what this means
                     num_indv = as.factor(data$num_idv),
                     Effort = length(unique(data$session_id)), #this is the number of surveys T and C, right?
                     distance = na.omit(as.numeric(data$platform_distance)))

# look at distance data to assess potential truncation                   
ggplot(DSdata, aes(x=distance))+
  geom_histogram(binwidth=1)+ 
  theme_classic()


# Step 3: run ds model
library("Distance")
# define distance bins
# bin <- c(0, 50, 100, 150, 200, 250, 300, 350, 400)
# specify model 
df_m0 <- ds(DSdata, truncation=250, transect="point",
            formula= ~1, key="hn")
summary(df_m0)
#plot PDF with DF
plot(df_m0, nc = 10, main = "No covariates",
     pch = 1, cex = 0.5, pdf = TRUE)
#goodness of fit
gof_ds(df_m0, chisq = T)

# add groupsize covar
df_m1 <- ds(DSdata, truncation=250, transect="point",
            formula= ~num_indv, key="hn")

#doesnt fit,  try opbserver covar 

df_m2 <- ds(DSdata, truncation=250, transect="point",
            formula= ~OBS, key="hn")
# not as good
# try combo of covar
df_m3 <- ds(DSdata, truncation=250, transect="point",
            formula= ~OBS + num_indv, key="hn")
#doesnt fit



# GIS DATA -- This is a sstruggle 

install.packages("GDAL")
install.packages("OGR")
install.packages("PROJ.4")
install.packages("sp")
library("rgdal")
library("maptools")
library("ggplot2")
library("plyr")

# provide the correct projection for the data
newproj <- "+proj=longlat +datum=WGS84 +no_defs"
  # import shapefile for the survey area
shape <- readShapeSpatial("ballard_study_area_v2-polygon.shp", proj4string = CRS(newproj),
                          repair=TRUE, force_ring=T, verbose=TRUE)

## import shapef ile for the points
obs_pt <- readShapeSpatial("obs_point-point.shp", proj4string = CRS(newproj),
                        repair=TRUE, force_ring=T, verbose=TRUE)
# make the object simpler
survey.area <- data.frame(shape@polygons[[1]]@Polygons[[1]]@coords)
names(survey.area) <- c("x","y")

# area of study area
shape$data$AREA













