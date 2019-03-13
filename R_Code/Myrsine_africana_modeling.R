library(raster)
library(dismo)
library(maptools)
library(FWDselect)
library(rgdal)
library(rJava)
library(classInt)
library(corrplot)
library(ENMeval)
library(MaxentVariableSelection)# for variable selection
library(virtualspecies) #for collinearity
library(ISLR)# for cross validation
library(boot)#for cross validation
#library(ggplot2)
#library(ggpubr)
library(rgeos)
library(sp)
presence <- readRDS("data/Myrsine_africana_combined.rds")
head(presence)
tail(presence)
dim(presence)
names(presence)

#Bioclim data at coarse resolution can be dowloaded within R
bioclim <- getData("worldclim", download = F, path = "./data", var="bio", res=2.5)
names(bioclim)
projection(bioclim)

### ALWAYS plot the data to see if everything is ok!
x11()
plot(bioclim, 1)


### Same projection but different projection string, fix
data("wrld_simpl")
projection(wrld_simpl) <- projection(bioclim)

proj4string(bioclim)=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 

#spatial presences
presences <- SpatialPoints(presence[,c("decimalLongitude", "decimalLatitude")], proj4string=CRS(projection("+proj=longlat +ellps=WGS84 +datum=WGS84")))
samburu <- readRDS("data/samburu_shp.rds")
plot(samburu, col="blue", axes=T)
#str(samburu)
samburu <- spTransform(samburu, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) 
presences_samburu <- presences[samburu]
#presences_samburu <- presences[is.na(over(samburu, presences)),]

plot(samburu)
points(presences_samburu, pch=16, cex=0.3, col="red")

#withholding 20% of the data or sample for testing 
fold <- kfold(presences_samburu, k=5)

occtest <- presences_samburu[fold==1,]
occtrain <- presences_samburu[fold !=1,]


#crop the bioclim variables to the study area (samburu)
bioclim_samburu <- crop(bioclim, samburu)
bioclim_samburu <- mask(bioclim_samburu, samburu)

plot(bioclim_samburu[[1]])
plot(bioclim_samburu[[12]])
points(presences_samburu, pch=16, cex=0.3, col="blue")

saveRDS(bioclim_samburu, "Data/bioclim_samburu.rds")

me1 <- dismo::maxent(bioclim_samburu[[c("bio2", "bio4", "bio14", "bio12", "bio6", "bio18")]], occtrain, path=paste0(getwd(), "/figures/myrsine_africana/maxent_outputs"), args=c("-J", "-P"), userfeatures="LQ") # Boden als Faktor

plot(me1)

#response curves
response(me1)

# predict to entire dataset
r <- predict(me1, bioclim_samburu) 
plot(r)

# with some options:
 r1 <- predict(me1, bioclim_samburu, args=c("outputformat=raw"), progress='text', 
      filename='maxent_prediction.grd', overwrite=T)

plot(r1)
points(presences_samburu)


#testing
 #create background points
set.seed(1)
background <- randomPoints(bioclim_samburu, 1000)


#simplest way to evaluate
#e1 <- evaluate(me1, p=occtest, a=background, x=bioclim_samburu)
#e1
#alternative 1
#extracting values from raster
presval_test <- (extract(bioclim_samburu, occtest))
absval_test <- (extract(bioclim_samburu, background))

#e2 <- evaluate(me1, presval_test, absval_test)

# alternative 2 
# predict to testing points 
testp <- predict(me1, presval_test) 
head(testp)
testa <- predict(me1, absval_test) 

e3 <- evaluate(p=testp, a=testa)
e3
threshold(e3)

plot(e3, 'ROC')

#MODELING WITH SOIL

# Bodendaten des African Soil Atlas einlesen
soil_shp <- readOGR("data/Africa_soil_WRB" , layer = "afticasoilmap")

# Gebiet auf Untersuchungsgebiet begrenzen
soil_shp <- crop(soil_shp, samburu)
plot(soil_shp)
#soil_shp <- crop(soil_shp, UG)
#soil_shp <- crop(soil_shp, africa)
soil_shp@data
soil_leg <- soil_shp@data

### remove dublicaded to create a legend

soil <- rasterize(soil_shp, bioclim_samburu, field = "GRIDCODE")

soil <- crop(soil, samburu)
soil <- mask(soil, samburu)

plot(soil)
writeRaster(soil, "soil_wrb.tif", overwrite = TRUE)
#write.table(soil_leg, "wrb_leg.txt")
table(soil[])

# Bodendaten werden in einem separaten Skript auf Hauptbodentypen reduziert
soil <- stack("data/soil_mst.tif")
plot(soil)


#rasterstack aus bioclim_UG und soil bauen mit Befehl stack
bioclim_soil_samburu <- stack(bioclim_samburu, soil_shp)
writeRaster(bioclim_soil_samburu, "bioclim_soil_samburu.tif", overwrite = TRUE)

bioclim_soil_samburu <- stack("bioclim_soil_samburu.tif")
bioclim_soil_samburu@layers


names(bioclim_soil_samburu)
# SDM mit Software Maxent laufen lassen
#me <- maxent(bioclim_soil_UG, presences@coords, background)
#plot(me) # Variable contribution anzeigen lassen, relevante Variablen aussuchen (Vorsicht bei gegenseitiger Beeinflussung), maximal 4 bis 5 Variablen

#quick check if the predictors are missing for half of occurence points
x <- extract(bioclim_soil_samburu, presence)
sum(is.na(x))

#maxent has an inbuilt program that conducts cross validation
# Best model possible with the highest AUC and highest correlation 
me2 <- dismo::maxent(bioclim_soil_samburu[[c("bioclim_soil_samburu.1", "bioclim_soil_samburu.8", "bioclim_soil_samburu.14", "bioclim_soil_samburu.12", "bioclim_soil_samburu.6", "bioclim_soil_samburu.18")]], p = occtrain, args=c("-J", "-P"), path=paste0(getwd(), "/figures/myrsine_africana/maxent_outputs.2"), userfeatures="LQ") 

### Add b

### The soil layer is included as layer 19, I assume. Needs to be checked! Tell maxent that this layer is a factor, not a numeric variable.
me2_soil <- dismo::maxent(bioclim_soil_samburu[[c("bioclim_soil_samburu.1", "bioclim_soil_samburu.8", "bioclim_soil_samburu.14", "bioclim_soil_samburu.12", "bioclim_soil_samburu.6", "bioclim_soil_samburu.18", "bioclim_soil_samburu.19")]], p = occtrain, factors = "bioclim_soil_samburu.19", args=c("-J", "-P"), path=paste0(getwd(), "/figures/myrsine_africana/maxent_outputs.2"), userfeatures="LQ") # Boden als Faktor




plot(me2)

pdf("figures/Model_without_soil_Variable_contribution.pdf", width=7, height=5)

par(mfrow = c(1, 1))
plot(me2)
dev.off()


plot(me2)

pdf("figures/Model_with_soil_Variable_contribution.pdf", width=7, height=5)

par(mfrow = c(1, 1))
plot(me2_soil)
dev.off()

#predict 
r3 <- predict(me2, bioclim_soil_samburu)
plot(r3)

r3_soil <- predict(me2_soil, bioclim_soil_samburu)
plot(r3_soil)

### threshold and threshold dependent model quality measures
### evaluate, threshold

# predict with some options:
r4 <- predict(me2, bioclim_soil_samburu, args=c("outputformat=raw"), progress='text', 
             filename='maxent_prediction_withsoil.grd', overwrite=T)
plot(r4)

#withholding 20% of the data or sample for testing 
fold <- kfold(presences_samburu, k=5)

occtest <- presences_samburu[fold==1,]
occtrain <- presences_samburu[fold !=1,]

#crop the bioclim variables to the study area (samburu)
bioclim_soil_samburu <- crop(bioclim_soil_samburu, samburu)
bioclim_soil_samburu <- mask(bioclim_soil_samburu, samburu)


#testing
#Evaluation
#extracting values from raster
presval_test_me2 <- (extract(bioclim_soil_samburu, occtest))
absval_test_me2 <- (extract(bioclim_soil_samburu, background))

 
# predict to testing points 
testp_me2 <- predict(me2, presval_test_me2) 
head(testp_me2)
testa_me2 <- predict(me2, absval_test_me2) 

e4 <- evaluate(p=testp_me2, a=testa_me2)
e4
thr <- threshold(e4)
thr5 <- threshold(e4, stat="no_omission")# 0% omission rate
thr6 <- threshold(e4, stat="spec_sens")# highest TSS
thr7 <- threshold(e4, stat="sensitivity", sensitivity=0.9)
thr8 <- threshold(e4, stat="sensitivity", sensitivity=0.95)
plot(thr8)
plot(e4, 'ROC')

#e4 <- evaluate(me2,  presval_test, absval_test)
### ROC-curve and density plot
x11()
par(mfrow = c(1, 2)) 
pdf("figures/ROC_density_samburu_mysrine_africana_with_soil_bioclim.pdf", width=7, height=5)
plot(e4, "ROC")
density(e4)
dev.off()

### Response curves ####
pdf("figures/Response_curve_mysrine_africana_with_soil_bioclim_.pdf", width=7, height=5)
response(me2)
dev.off()

# Graphical representation of the distribution map
### 2D response curve
#Values are pred column
np = 10

bio6 <- seq(min(bioclim_soil_samburu[["bioclim_soil_samburu.6"]][], na.rm=T), max(bioclim_soil_samburu[["bioclim_soil_samburu.6"]][], na.rm=T), len=np)
bio1 <- seq(min(bioclim_soil_samburu[["bioclim_soil_samburu.1"]][], na.rm=T), max(bioclim_soil_samburu[["bioclim_soil_samburu.1"]][], na.rm=T), len=np)
bio8 <- seq(min(bioclim_soil_samburu[["bioclim_soil_samburu.8"]][], na.rm=T), max(bioclim_soil_samburu[["bioclim_soil_samburu.8"]][], na.rm=T), len=np)
bio18 <- seq(min(bioclim_soil_samburu[["bioclim_soil_samburu.18"]][], na.rm=T), max(bioclim_soil_samburu[["bioclim_soil_samburu.18"]][], na.rm=T), len=np)
bio14 <- seq(min(bioclim_soil_samburu[["bioclim_soil_samburu.14"]][], na.rm=T), max(bioclim_soil_samburu[["bioclim_soil_samburu.14"]][], na.rm=T), len=np)
bio12 <- seq(min(bioclim_soil_samburu[["bioclim_soil_samburu.12"]][], na.rm=T), max(bioclim_soil_samburu[["bioclim_soil_samburu.12"]][], na.rm=T), len=np)


hist(bioclim_soil_samburu[["bioclim_soil_samburu.6"]][])
quantile(bioclim_soil_samburu[["bioclim_soil_samburu.6"]][], na.rm=T)

bio6 <- seq(quantile(bioclim_soil_samburu[["bioclim_soil_samburu.6"]][], probs=0.05, na.rm=T), quantile(bioclim_soil_samburu[["bioclim_soil_samburu.6"]][], probs=0.95, na.rm=T), len=np)
bio1 <- seq(quantile(bioclim_soil_samburu[["bioclim_soil_samburu.1"]][], probs=0.05, na.rm=T), quantile(bioclim_soil_samburu[["bioclim_soil_samburu.1"]][], probs=0.95, na.rm=T), len=np)
bio8 <- seq(quantile(bioclim_soil_samburu[["bioclim_soil_samburu.8"]][], probs=0.05, na.rm=T), quantile(bioclim_soil_samburu[["bioclim_soil_samburu.8"]][], probs=0.95, na.rm=T), len=np)
bio18 <- seq(quantile(bioclim_soil_samburu[["bioclim_soil_samburu.18"]][], probs=0.05, na.rm=T), quantile(bioclim_soil_samburu[["bioclim_soil_samburu.18"]][], probs=0.95, na.rm=T), len=np)
bio14 <- seq(quantile(bioclim_soil_samburu[["bioclim_soil_samburu.14"]][], probs=0.05, na.rm=T), quantile(bioclim_soil_samburu[["bioclim_soil_samburu.14"]][], probs=0.95, na.rm=T), len=np)
bio12 <- seq(quantile(bioclim_soil_samburu[["bioclim_soil_samburu.12"]][], probs=0.05, na.rm=T), quantile(bioclim_soil_samburu[["bioclim_soil_samburu.12"]][], probs=0.95, na.rm=T), len=np)


#most common type of soil
soilpresences <- extract(bioclim_soil_samburu[["bioclim_soil_samburu.12"]], presence)
soilpresences <- na.omit(soilpresences)
#write.csv(soilpresences, "soilpresences.csv")
unique(soilpresences)
table(soilpresences)
which.max(table(soilpresences))

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

v <- soilpresences
result <- getmode(v)
print(result)

# specify the most common soil type for soil (mode value)
newdata <- expand.grid(bioclim_soil_samburu.6 = bio6, bioclim_soil_samburu.1 = bio1, bioclim_soil_samburu.8 = bio8, bioclim_soil_samburu.18 = bio18, bioclim_soil_samburu.14 = bio14, bioclim_soil_samburu.12 = bio12) #h?ufigster Hauptbodentyp mit Nummer 5
newdata$pred <- predict(me2, newdata)
plot(newdata$pred)

# Presentation of the response curve

### Use threshold to show distribution
newdata$pred[newdata$pred<thr$sensitivity] <- NA

#help(classIntervals)
#help("diff")
#help("unique")
### Create classes of site suitability
cInt <- classIntervals((newdata$pred))

xdiff <- diff(unique(newdata$bioclim_soil_samburu.1))[1]
ydiff <- diff(unique(newdata$bioclim_soil_samburu.12))[12]

mypalette <- colorRampPalette(c("lightgreen", "darkgreen"))
newdata$colors <- findColours(cInt, mypalette(length(cInt$brks)))

pdf("figures/Model_with_soil.pdf", width=7, height=5)
par(mfrow = c(1,1), mar = c(5,5,1,1))
symbols(x = newdata$bioclim_soil_samburu.1, y = newdata$bioclim_soil_samburu.12, rectangles = matrix(rep(c(xdiff, ydiff), nrow(newdata)), ncol = 2, byrow = T), 
        bg = newdata$colors, 
        fg = "white", inches = F, 
        xlab = "Annual Mean Temperature (?C)", 
        ylab = "Precipitation of Driest Quarter (mm)", 
        ylim = c(0, 1200), axes = F)
#contour(x = unique(newdata$bioclim_soil_samburu.1), y = unique(newdata$bioclim_soil_samburu.12), z = matrix(newdata$pred, nrow = np), add = T, levels = unique(round(cInt$brks,1)), labcex = 1.3)
axis(1, at = c(-100, 0, 100, 200, 300, 350), labels = c(-100/10, 0, 100/10, 200/10, 300/10, 350/10))
axis(2)
mtext("Myrsine africana", side = 3, line = -2.3, font = 3)
mtext(paste0("AUC = " , round(e3@auc, 2), " "), side = 3, line = -5, adj = 1)
mtext(paste0("Pearson r? = " , round(e3@cor, 2), " "), side = 3, line = -6.3, adj = 1)
dev.off()

### Plot distribution map ####

x11(width = 12)
par(mfrow = c(1, 1))
pred <- predict(me2, bioclim_soil_samburu)

#plot(pred, xlab = "Longitude", ylab = "Latitude", xlim = c(-30,60), cex.axis = 1.2, cex.lab = 1.2, 
     #legend.args = list(text='Probability of occurence', side = 4, font = 1, line = 2.8, cex = 1.2), 
     #legend.width = 1, legend.shrink = 1)

distr <- pred
plot(distr)
distr[distr < thr$sensitivity] <- NA
cInt <- classIntervals((newdata$pred))


# Distribution map 
# plot(distr, col = mypalette(10), breaks = cInt$brks, legend = F, xlab = "Longitude", ylab = "Latitude")
x11(width = 12)
plot(distr, xlab = "Longitude", ylab = "Latitude", xlim = c(-30,50), cex.axis = 1.2, cex.lab = 1.2, 
     legend.args = list(text='Probability of occurence', side = 4, font = 1, line = 2.8, cex = 1.2), 
     legend.width = 1, legend.shrink = 1)
plot(samburu, add = T)
points(presences, pch = 17, cex = 0.6, col = "red")

mtext("Myrsine africana", side = 1, line = -1.3, font = 3)
mtext(paste0("AUC = " , round(e3@auc, 2), " "), side = 1, line = -2.3, adj = 1)
mtext(paste0("Pearson r? = " , round(e3@cor, 2), " "), side = 1, line = -1.3, adj = 1)
legend("topleft", legend = c("presence data"), pch = 17, col = c("red"), border = "white", cex = 0.8, bty = "n", bg = "white")


### Climate change projection ####
cc <- getData('CMIP5', var="bio", res=5, rcp=85, model='HD', year=70, download=F, path="data")

cc <- crop(cc, samburu)
cc <- mask(cc, samburu)
extent(cc)
extent(bioclim_soil_samburu)
plot(cc[[1]])
plot(bioclim_soil_samburu[[1]]) 

#make the extent of cc and that of bioclim_soil_samburu similar
cc2<- resample(cc, bioclim_soil_samburu, method="bilinear")
cc2 <- crop(cc2, bioclim_soil_samburu)
cc2 <- mask(cc2, bioclim_soil_samburu)
names(cc2) <- names(bioclim_soil_samburu)
plot(cc2[[1]])
plot(bioclim_soil_samburu[[1]]) 

writeRaster(cc2, "cc._samburu_tif", overwrite = TRUE)

cc <- stack("cc._samburu_tif")

pred_cc <- predict(me2, cc)

distr_cc <- pred_cc
distr_cc[distr_cc[] < thr$sensitivity] <- NA

plot(pred_cc)
hist(distr_cc)
hist(distr_cc)

# Graphical representation of the different climate scenarios

### Create figure for distribution for current and climate change projection ####
x11(width = 12)

#pdf("figures/maps_myrsine africana.pdf", width = 12, pointsize = 16)

plot(distr_cc, col = mypalette(length(cInt$brks)), breaks = cInt$brks, legend = F, xlab = "Longitude", ylab = "Latitude")
plot(samburu, add = T)
mtext("RCP85 HD ", 3, -1.2, adj = 0)
mtext("Myrsine africana", 1, -1.2, font = 3)
mtext("WGS84 ", 1, -1.2, adj = 1)

dev.off()

### Create figure for distribution for current and climate change projection ####
#x11(width=12)
pdf("figures/maps_myrsine_afri.pdf", width=12, pointsize = 16)
par(mfrow=c(1,2))
plot(distr, col=mypalette(length(cInt$brks)), breaks=cInt$brks, legend=F, xlab="Longitude", ylab="Latitude")
plot(samburu, add=T)
mtext("Current climate ", 3,-2.2)
mtext("Myrsine_africana", 3,-1.2, font=3)
mtext("WGS84 ", 1,-1.2, adj=1)
plot(distr_cc, col=mypalette(length(cInt$brks)), breaks=cInt$brks, legend=F, xlab="Longitude", ylab="Latitude")
plot(samburu, add=T)
mtext("Climate scenario RCP85 HD ", 3,-2.2)
mtext("Myrsine_africana ", 3,-1.2, font=3)
mtext("WGS84 ", 1,-1.2, adj=1)
dev.off()

