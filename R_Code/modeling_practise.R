library(raster)
library(dismo)
library(maptools)
library(FWDselect)
library(rgdal)
library(rJava)

presence <- readRDS("../data/Myrsine_africana_combined.rds")
head(presence)
tail(presence)
dim(presence)

#Bioclim data at coarse resolution can be dowloaded within R
bioclim <- getData("worldclim", download = F, path = "../", var="bio", res=5)
names(bioclim)
projection(bioclim)

### ALWAYS plot the data to see if everything is ok!
x11()
plot(bioclim)
plot(bioclim, 1)
plot(bioclim, 12)

x11()
data("wrld_simpl")
plot(wrld_simpl, col="blue", axes=T)


### Same projection but different projection string, fix
projection(wrld_simpl) <- projection(bioclim)
#help("projection")
help("SpatialPoints")
# create spatial points dataframe 
### Create spatial points dataframe
presences <- SpatialPoints(presence[,c("decimalLongitude", "decimalLatitude")], proj4string=CRS(projection(bioclim)))
x11()
plot(wrld_simpl)
points(decimalLatitude ~ decimalLongitude, data=presences, pch=16, cex=0.5, col="red")

#How about in Kenya 
#kenya <- wrld_simpl[wrld_simpl$NAME=="Kenya",]
#points(presences, kenya)
#OR 

#presenceval <- extract(bioclim, presence)
#data("wrld_simpl")
#x11()
#plot(wrld_simpl, col="blue", axes=T)
#


### Use only europe; decide for your project, what to to!
#europe <- wrld_simpl[wrld_simpl$REGION==150,]
#europe@data 
#europe <- europe[europe$NAME!="Russia",]
#plot(europe)
names(presence)
View(presence)
#use Africa
africa <- wrld_simpl[wrld_simpl$REGION== 2,]
africa@data
plot(africa, col="blue", axes=T)
points(decimalLatitude ~ decimalLongitude, data=presence, col="red", pch=16, cex=0.3)

saveRDS(africa, "Data/Africa_shp.rds")

#### Alternative: Use function extent() to define a region and crop with a rectangle.

# exclude presence data outside Africa (u can use spatial presence data in this case)
presences_africa <- presences[is.na(over(africa, presences)),]
plot(africa, col="blue", axes=T)
points(presences_africa, pch=16, cex=0.3, col="red")

#crop the bioclim variables to the study area (Africa)
bioclim_africa <- crop(bioclim, africa)
bioclim_africa <- mask(bioclim_africa, africa)


saveRDS(bioclim_africa, "Data/bioclim_africa.rds")

# create background points
background <- randomPoints(bioclim_africa, 2000)


#envval <- extract(bio, presence)

#plot(bio,1)
#plot(wrld_simpl, add=T)
#points(presence$decimalLongitude, presence$decimalLatitude, pch=20, cex=0.5, col="red")

#set.seed(0)
#backg <- randomPoints(bio, 500)
#absencevals <- extract(bio, backg)
#pb <- c(rep(1, nrow(envval)), rep(0, nrow(absencevals)))
#sdmdata <- data.frame(cbind(pb, rbind(envval, absencevals))) 
#dim(sdmdata)

me1 <- maxent(bioclim_africa, presences_africa@coords, background)
me <- maxent(bioclim_africa[[c("bio10", "bio18")]], p=presences_africa@coords, a=background) 


install.packages("slam", lib="C:/Users/Dikko/Documents/R/R-3.3.1/library", repos =)


### Arguments of maxent: https://groups.google.com/forum/#!topic/maxent/yRBlvZ1_9rQ

#model evaluation
par(mfrow(c(1,2)))
plot(me1)
plot(me)

e1 <- evaluate(presences_africa, background, me1, bioclim_africa)
e <- evaluate(presences_africa, background, e, bioclim_africa)

e1
e

#assume you have done the variable selection then you evaluate the final model only


thr <- threshold(e)

### ROC-curve and density plot
par(mfrow=c(1,2)) 
plot(e, "ROC")
density(e)

### Response curves ####
response(me)

### 2D response curves, nicer but also more code ####
np <- 30
newdata <- expand.grid(bio10=seq(145, 200, len=np), bio18=seq(0, 240, len=np))
newdata$pred <- predict(me, newdata)

### Use threshold to show distribution
newdata$pred[newdata$pred<thr$sensitivity] <- NA

### Create classes of site suitability
cInt <- classIntervals((newdata$pred))

xdiff <-diff(unique(newdata$bio10))[1]
ydiff <-diff(unique(newdata$bio18))[1]

mypalette <- colorRampPalette(c("lightgreen", "darkgreen"))
newdata$colors <- findColours(cInt, mypalette(length(cInt$brks)))

par(mfrow=c(1,1), mar=c(5,5,1,1))
symbols(x=newdata$bio10, y=newdata$bio18, rectangles=matrix(rep(c(xdiff, ydiff), nrow(newdata)), ncol=2, byrow=T), bg=newdata$colors, fg="white", inches=F, xlab="Temperature of warmest quarter (°dC)", ylab="Precipitation of warmest quarter (mm)")
contour(x=unique(newdata$bio10), y=unique(newdata$bio18), z=matrix(newdata$pred, nrow=np), add=T, levels=unique(round(cInt$brks,1)), labcex = 1.3)
mtext("Cicendia filiformis", side=3, line=-1.3, font=3)
mtext(paste0("AUC = " , round(e@auc, 2), " "), side=1, line=-2.3, adj=1)
mtext(paste0("Pearson r = " , round(e@cor, 2), " "), side=1, line=-1.3, adj=1)

### Q: Is this response curve reasonable? Why?

### Plot distribution map ####
pred <- predict(me, bioclim_europe)
plot(pred)
distr <- pred
distr[distr < thr$sensitivity] <- NA
# cInt <- classIntervals((newdata$pred))

plot(distr, col=mypalette(10), breaks=cInt$brks, legend=F)
points(presences_europe, pch=16, cex=0.1, col="black")
plot(europe, add=T)
mtext("Cicendia filiformis", side=3, line=-1.3, font=3)
mtext(paste0("AUC = " , round(e@auc, 2), " "), side=1, line=-2.3, adj=1)
mtext(paste0("Pearson r = " , round(e@cor, 2), " "), side=1, line=-1.3, adj=1)
### Add legend!

### Climate change projection ####
cc <- getData('CMIP5', var="bio", res=5, rcp=85, model='HD', year=70, download=TRUE, path="data")
cc <- crop(cc, bioclim_europe)
cc <- mask(cc, bioclim_europe)
names(cc) <- names(bioclim_europe)

pred_cc <- predict(me, cc)

distr_cc <- pred_cc
distr_cc[distr_cc[] < thr$sensitivity] <- NA

### Create figure for distribution for current and climate change projection ####
#x11(width=12)
pdf("figures/maps.pdf", width=12, pointsize = 16)
par(mfrow=c(1,2))
plot(distr, col=mypalette(length(cInt$brks)), breaks=cInt$brks, legend=F, xlab="Longitude", ylab="Latitude")
plot(europe, add=T)
mtext("Current climate ", 3,-2.2)
mtext("Cicendia filiformis ", 3,-1.2, font=3)
mtext("WGS84 ", 1,-1.2, adj=1)
plot(distr_cc, col=mypalette(length(cInt$brks)), breaks=cInt$brks, legend=F, xlab="Longitude", ylab="Latitude")
plot(europe, add=T)
mtext("Climate scenario RCP85 HD ", 3,-2.2)
mtext("Cicendia filiformis ", 3,-1.2, font=3)
mtext("WGS84 ", 1,-1.2, adj=1)
dev.off()

### Q: Is this expected change of distribution reasonably from an ecological
### perspective?

### Q: Schould we add a legend? If yes, what's the problem and how does a
### legend make sense?

### Q: Prediction uncertainty: Derive from validation, e.g. k-fold 
### cross-validation. Where is the prediction uncertainty higher/lower?

### Q: Calculate future spatial distribution under the assumption of full 
### dispersal and no dispersal. Which dispersal scenario will be more 
### likely for your species?

### Q: How looks the spatial autocorrelaiton in your data? Is is an issue 
### to be dealt with?












#pairs plot of values of environmental data
pairs(sdmdata[,2:5], cex=0.5, fig=T)


# Model fitting
model1 <- glm(pb ~ bio1 + bio12 + bio5, data = sdmdata)
plot(model1)# learn and check the intepretation of these curves

summary(model1)#learn and check the intepretation of the summary

model2 <- glm(pb ~ ., data=sdmdata)
plot(model2)
summary(model2)#learn and check the intepretation of the summary

#Models that are implemented in dismo do not use a formula (and most models only take presence points). 
#Bioclim is an example. It only uses presence data, so we use ’presvals’ instead of ’sdmdata’

bc <- bioclim(envval[,c("bio1", "bio12", "bio5")])
class(bc)
bc
pairs(bc)

#making predictions
bio1 = c(40,150,200)
bio5 = c(60,115,290)
bio12 = c(600,1600,1700)
pd <- data.frame(cbind(bio1,bio5,bio12))
pd
predict(model1, pd)

predict(bc, pd)

response(bc)


names(bio)
p <- predict(bio, model1)
plot(p)

#making prediction of my data (Myr afr) example using bioclim

myafr_pre <- predict(bc,envval)
plot(myafr_pre)
response(bc)

names(bio)
p1 <- predict(bio, bc)
plot(p1)


#NB TRy making prediction of my data (Myr afr) example using glm


