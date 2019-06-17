
# load require libraries
library(rJava)
library(maxent)
library(dismo)
library(maptools)#
library(rgdal)
library(rgeos)
library(sp)
library(raster)#raster files
library(vegan)
library(corrplot)
library(spdep) # neighbours
library(classInt)
library(nlme)
library(colorRamps)#for some crispy colors


presence <- readRDS("data/combined_data/rhus_natalensis_combined.rds")
dim(presence)
names(presence)

bioclim <- getData('worldclim', download=T, path='./data', var='bio', res=2.5)
projection(bioclim)

### Same projection but different projection string, fix
data("wrld_simpl")
projection(wrld_simpl) <- projection(bioclim)
proj4string(bioclim)=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 

#plot africa
africa <- wrld_simpl[wrld_simpl$REGION==2,]

#presences only in Samburu
presences <- SpatialPoints(presence[,c("decimalLongitude", "decimalLatitude")], 
                           proj4string=CRS(projection("+proj=longlat +ellps=WGS84 +datum=WGS84")))
samburu <- readRDS("data/samburu_shp.rds")
plot(samburu, col="blue", axes=T)
samburu <- spTransform(samburu, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) 
#africa <- spTransform(africa, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) 
presences_africa <- presences[africa]
presences_samburu <- presences[samburu]

#plot presences_africa
plot(africa)
points(presences_africa, col="blue", cex=0.5, pch=16)

#withholding 25% of the data or sample for testing 
#75% of the sample size
smp_size <- floor(0.75 * nrow(presence))

## set the seed to make your partition reproducible
set.seed(1000)
train_ind <- sample(seq_len(nrow(presence)), size = smp_size)

train <- presence[train_ind, ]
test <- presence[-train_ind, ]

#spatial presences training data set
presences_train <- SpatialPoints(train[,c("decimalLongitude", "decimalLatitude")], 
                                 proj4string=CRS(projection("+proj=longlat +ellps=WGS84 +datum=WGS84")))
africa <- spTransform(africa, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) 
presences_africa_train <- presences_train[africa]

#spatial presences testing data set
presences_test <- SpatialPoints(test[,c("decimalLongitude", "decimalLatitude")], 
                                proj4string=CRS(projection("+proj=longlat +ellps=WGS84 +datum=WGS84")))
presences_africa_test <- presences_test[africa]
plot(presences_africa_train, add=T)
points(presences_africa_test)

#crop the bioclim variables to africa
bioclim_africa <- crop(bioclim, africa)
bioclim_africa <- mask(bioclim_africa, africa)

#crop the bioclim variables to the study area (samburu)
bioclim_samburu <- crop(bioclim, samburu)
bioclim_samburu <- mask(bioclim_samburu, samburu)


#create background points
set.seed(1000)
background <- randomPoints(bioclim_africa, 10000)


#dividing the background into training and testing
# Now Selecting 75% of data as sample from total 'n' rows of the data  
sample <- sample.int(n = nrow(background), size = floor(.75*nrow(background)), replace = F)
background_train <- background[sample, ]
background_test  <-background[-sample, ]


predictors <- bioclim_africa[[c("bio12", "bio16", "bio7", "bio2",  "bio8")]]

mx_best <- dismo::maxent(predictors, presences_africa_train, background_train)
