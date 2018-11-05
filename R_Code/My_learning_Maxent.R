
# installing many packages at a time
install.packages(c("gbm", "rgbif","spdep", "nfc"))
# load the libraries
library(raster)
library(maptools)
library(dismo)
library(rgdal)
library(rJava)
library(gbm)
library(corrplot)
library(hier.part)

file <- paste(system.file(package = "dismo"), "/ex/bradypus.csv", sep="")
bradypus <- read.table(file, header=T, sep=",")
dim(bradypus)
head(bradypus)

#we only need the 2 and 3 rows
bradypus <- bradypus[,2:3]
head(bradypus)

data("wrld_simpl")
x11()
plot(wrld_simpl, col="red", axes=T)

#environmental data or predictors
# Bioclim data provided by the package dismo
files <- list.files(path=paste(system.file(package="dismo"),"/ex", sep=""), pattern="grd", full.names=TRUE )
predictors <- stack(files)
plot(predictors)

plot(predictors, 1)

plot(wrld_simpl, add=T)
points(bradypus, col="blue", cex=0.5, pch=15) 
# create background points
?set.seed
set.seed(8)
backgrd <- randomPoints(predictors, 500)
x11()
plot(wrld_simpl, col="red", axes=T)
points(bradypus, col="blue", pch=16, cex=0.3)
points(backgrd, col="yellow", pch=16, cex=0.5)

# Extracting values from the raster
envalspre <- extract(predictors, bradypus)
envalabse <- extract(predictors, backgrd)
head(envalspre)
head(envalabse)

#create predictor dataset
pb <- c(rep(1, nrow(envalspre)), rep(0, nrow(envalabse)))
sdmdata <- data.frame(cbind(pb, rbind(envalspre, envalabse)))
head(sdmdata)
tail(sdmdata)
names(sdmdata)
sdmdata$biome
sdmdata[,"biome"] <- as.factor(sdmdata[,"biome"])
head(sdmdata)

#correlation among predictors
x11()
pairs(sdmdata[,2:8], cex=0.1, fig=T)
corrplot(cor(sdmdata[,1:8]), type="lower", diag=F)

#Model fitting using model types
#Algorithms require different data structure, i.e. sdmdata, predictors, bradypus.

#logistic regression should be used with presence-absence data
m1 <- glm(pb ~ bio1 + bio5 + bio12, data = sdmdata, family = "binomial")
summary(m1)

m2 <- step(glm(pb ~ ., data = sdmdata))
summary(m2)

#climatic envelope model only uses only the presence data
bc <- bioclim(envalspre[,c("bio1", "bio12", "bio5")])
summary(bc)
 #maxent uses only presence data
me <- maxent(predictors, bradypus)

#Boosted regression tree withb presence-absence data
brt <- gbm(pb ~ bio1+ bio12 + bio5, data= sdmdata, distribution= "bernoulli")







