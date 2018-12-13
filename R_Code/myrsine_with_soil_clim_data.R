
library(raster)
library(maptools)
library(gbm)
library(dismo)
library(rgdal)
library(corrplot)
library(rJava)
library(hier.part)
library(colorRamps)
library(classInt)


#read presence data
presence <- readRDS("data/combined_data/myrsine_africana_combined.rds")
head(presence)
names(presence)
dim(presence)
#read bioclim data
bioclim <- getData("worldclim", download = F, path = "./data", var="bio", res=5)
projection(bioclim)

data("wrld_simpl")
plot(wrld_simpl, col="blue")

#same projection but different different projection string, fix
projection(wrld_simpl) <- projection(bioclim)

#plot bioclim

plot(bioclim)
plot(bioclim, 1)
plot(bioclim, 12)

#creating spatial points dataframe

presence <- SpatialPoints(presence[,c("decimalLatitude", "decimalLongitude")], proj4string = CRS(projection(bioclim)))

#use only data from africa
x11()
africa <- wrld_simpl[wrld_simpl$REGION==2,]
plot(africa, col="blue")

points(decimalLatitude ~ decimalLongitude, data= presence, col="red", cex=0.3, pch=16)

saveRDS(africa, "data/africa_shp.rds")
# use data from Kenya only
kenya <- wrld_simpl[wrld_simpl$NAME=="Kenya",]
plot(kenya, col="red")
points(decimalLatitude ~ decimalLongitude, data=presence, pch=16, cex=0.3, col="blue")

#presence africa
presence_africa <- presence[africa]
#presences_africa <- presence[is.na(over(africa, presence)),]


#crop predictor raster for the study region only
bioclim_africa <- crop(bioclim, africa)
bioclim_africa <- mask(bioclim_africa, africa)

saveRDS(bioclim_africa, "data/bioclim_africa.rds")


#background
set.seed(8)
background <- randomPoints(bioclim_africa, 500)

#extracting values from raster
presvals <- extract(bioclim_africa, presence_africa)
absvals <- extract(bioclim_africa, background)

#creating predictor data set
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))

names(sdmdata)
ncol(sdmdata)
dim(sdmdata)
#correlation among the predictor variables
x11()
pairs(sdmdata[,2:10], cex=0.1, fig=T)
#corrplot(cor(sdmdata[,1:19]), type = "lower", diag=F)

corr <- cor(na.omit(sdmdata[,2:7]), method = "spearman")
corrplot(corr, col = NULL, diag = T, type = "lower")
#modeling and variable selection 


me1 <- maxent(bioclim_africa, presence_africa@coords, background)
me <- maxent(bioclim_africa[[c("bio10", "bio18")]], p=presence_africa@coords, a=background)

me
me1





#model evaluation 
par(mfrow=c(1,2))
plot(me1)
plot(me)

e1 <- evaluate(presence_africa, background, me1, bioclim_africa)
e <- evaluate(presence_africa, background, me, bioclim_africa)

e1
e

#assume you have done the variable selection then you evaluate the final model only

thr <- threshold(e)


#ROC-curve and density plots

par(mfrow=c(1,2))
plot(e, "ROC")
density(e)

#response curve
response(me)

#2D response curves are nicer but of course require more codes
np <- 30
newdata <- expand.grid(bio10=seq(145, 200, len=np), bio18=seq(0, 240, len=np))
newdata$predict <- predict(me, newdata)

#use threshold to show distribution
newdata$predict[newdata$predict<thr$sensitivity] <- NA

#create classes of site suitability
cInt <- classIntervals((newdata$predict))

xdiff <- diff(unique(newdata$bio10))[1]
ydiff <- diff(unique(newdata$bio18))[1]

mypalatte <- colorRampPalette(c("lightgreen", "darkgreen"))
newdata$colors <- findColours(cInt, mypalatte(length(cInt$brks)))
par(mfrow=c(1,1), mar=c(5,5,1,1)) 

symbols(x=newdata$bio10, y=newdata$bio18, rectangles=matrix(rep(c(xdiff, ydiff), nrow(newdata)), ncol=2, byrow=T), bg=newdata$colors, fg="white", inches=F, xlab="Temperature of warmest quarter (°dC)", ylab="Precipitation of warmest quarter (mm)")
contour(x=unique(newdata$bio10), y=unique(newdata$bio18), z=matrix(newdata$pred, nrow=np), add=T, levels=unique(round(cInt$brks,1)), labcex = 1.3)
mtext("Myrsine africana", side=3, line=-1.3, font=3)
mtext(paste0("AUC = " , round(e@auc, 2), " "), side=1, line=-2.3, adj=1)
mtext(paste0("Pearson r = " , round(e@cor, 2), " "), side=1, line=-1.3, adj=1)

#plot distribution of the map
pred <- predict(me, bioclim_africa)
plot(pred)
distr <- pred
distr[distr < thr$sensitivity] <- NA

#plot(distr, col=mypalatte(10), breaks = cInt$brks, legend=F)


#climate change projection 
cc <- getData("CMIP5", var="bio", res=5, rcp=85, model="HD", year=70, download = T, path = "./data")
cc <- crop(cc, bioclim_africa)
cc <- mask(cc, bioclim_africa)
names(cc) <- names(bioclim_africa)

predict_cc <- predict(me, cc)
distr_cc <- predict_cc
distr_cc[distr_cc[] < thr$sensitivity] <- NA

#create a figure for the distribution of current and climate change projection
dev.off()
pdf("figures/myrsine_africana_cur_cc_map.pdf", width = 12, pointsize = 16)
par(mfrow=c(1,2))
plot(distr, col=mypalatte(length(cInt$brks)), breaks=cInt$brks, legend=F, xlab="Longitude", ylab="Latitude")
plot(africa, add=T)
mtext("Current climate", 3, -0.1)
mtext("Myrsine africana", 3, -1.1, font = 3)
mtext("wgs84", 1, -1-2, adj = 1)
plot(distr_cc, col=mypalatte(length(cInt$brks)), breaks=cInt$brks, legend=F, xlab="Longitude", ylab="Latitude")
plot(africa, add=T)
mtext("Climate Scenaria RCP85 HD", 3, -0.1)
mtext("Myrsine africana", 3, -1.1, font = 3)
mtext("wgs84", 1, -1.2, adj=1)
dev.off()

#ANALYSIS WITH SOIL

# Bodendaten des African Soil Atlas einlesen
soil_shp <- readOGR("data/Africa_soil_WRB" , layer = "afticasoilmap")

# Gebiet auf Untersuchungsgebiet begrenzen
#soil_shp <- crop(soil_shp, UG)
soil_shp <- crop(soil_shp, africa)
plot(soil_shp)
soil_shp@data
soil_leg <- soil_shp@data

soil <- rasterize(soil_shp, bioclim_africa, field = "GRIDCODE")

#soil <- crop(soil, UG)
#soil <- mask(soil, UG)

#plot(soil)
#writeRaster(soil, "soil_wrb.tif", overwrite = TRUE)
#write.table(soil_leg, "wrb_leg.txt")



# Bodendaten werden in einem separaten Skript auf Hauptbodentypen reduziert
soil <- stack("data/soil_mst.tif")
plot(soil)


#rasterstack aus bioclim_UG und soil bauen mit Befehl stack
bioclim_soil_africa <- stack(bioclim_africa, soil_shp)
writeRaster(bioclim_soil_africa, "bioclim_soil_africa.tif", overwrite = TRUE)

bioclim_soil_africa <- stack("bioclim_soil_africa.tif")
bioclim_soil_africa@layers


names(bioclim_soil_africa)
# SDM mit Software Maxent laufen lassen
#me <- maxent(bioclim_soil_UG, presences@coords, background)
#plot(me) # Variable contribution anzeigen lassen, relevante Variablen aussuchen (Vorsicht bei gegenseitiger Beeinflussung), maximal 4 bis 5 Variablen


# Auswahl der Variablenkombination
me1 <- maxent(bioclim_soil_africa[[c("bioclim_soil_africa.10", "bioclim_soil_africa.18", "bioclim_soil_africa.20")]], p=presence_africa@coords, a=background) # Boden als Faktor
me2 <- maxent(bioclim_soil_africa[[c("bioclim_soil_africa.1", "bioclim_soil_africa.12", "bioclim_soil_africa.20")]], factors = "bioclim_soil_africa.20", p=presences@coords, a=background)
me3 <- maxent(bioclim_soil_africa[[c("bioclim_soil_africa.1", "bioclim_soil_africa.17", "bioclim_soil_africa.20")]], factors = "bioclim_soil_africa.20", p=presences@coords, a=background)
me4 <- maxent(bioclim_soil_africa[[c("bioclim_soil_africa.11", "bioclim_soil_africa.19", "bioclim_soil_africa.20")]], factors = "bioclim_soil_africa.20", p=presences@coords, a=background)
me5 <- maxent(bioclim_soil_UG[[c("bioclim_soil_africa.1", "bioclim_soil_africa.18", "bioclim_soil_africa.20")]], factors = "bioclim_soil_africa.20", p=presences@coords, a=background)





# Modellevaluierung
pdf("Ergebnisse/Modell_mit_Boden/Variable_contribution.pdf", width=7, height=5)
par(mfrow = c(1, 1))
plot(me1)
dev.off()

e1 <- evaluate(presence_africa, background, me1, bioclim_soil_africa)
e2 <- evaluate(presences, background, me2, bioclim_soil_UG)
e3 <- evaluate(presences, background, me3, bioclim_soil_UG)
e4 <- evaluate(presences, background, me4, bioclim_soil_UG)
e5 <- evaluate(presences, background, me5, bioclim_soil_UG)

# G?tema?e anschauen, welches Modell ist besser?, im Buch Franklin beschrieben
# das Modell mit dem h?chsten AUC und der h?chsten Korrelation ist am besten
e1 # AUC: 0.8823984, cor: 0.1735644 
e2 # AUC: 0.8462127, cor: 0.1627853  
e3 # AUC: 0.8858223, cor: 0.1765484
e4 # AUC: 0.8726011, cor: 0.1865184
e5 # AUC: 0.8728477, cor: 0.1699013 

# Modell me3 hat die besten Werte mit bio1 und bio17
# bio1: Annual Mean Temperature (in 10er Grad)
# bio17: Precipitation of Driest Quarter


### Let's assume we have done a valid variable selection and evaluate the final model only
thr <- threshold(e1)
thr


### ROC-curve and density plot
x11()
par(mfrow = c(1, 2)) 
#pdf("Ergebnisse/Modell_mit_Boden/ROC_density_mit Boden.pdf", width=7, height=5)
plot(e1, "ROC")
density(e1)
dev.off()

### Response curves ####
pdf("Ergebnisse/Modell_mit_Boden/Response curve_mit Boden.pdf", width=7, height=5)
x11()
response(me1)
dev.off()



# Graphische Darstellung der Verbreitungskarte

### 2D response curve
#Werte sind pred Spalte
np = 20

#table(is.na(bioclim_soil_UG[["bioclim_soil_UG.1"]][]))
#min(bioclim_soil_UG[["bioclim_soil_UG.1"]][], na.rm=T)
#max(bioclim_soil_UG[["bioclim_soil_UG.1"]][], na.rm=T)

bio1 <- seq(min(bioclim_soil_UG[["bioclim_soil_UG.1"]][], na.rm=T), max(bioclim_soil_UG[["bioclim_soil_UG.1"]][], na.rm=T), len=np)
bio17 <- seq(min(bioclim_soil_UG[["bioclim_soil_UG.17"]][], na.rm=T), max(bioclim_soil_UG[["bioclim_soil_UG.17"]][], na.rm=T), len=np)

hist(bioclim_soil_UG[["bioclim_soil_UG.1"]][])
quantile(bioclim_soil_UG[["bioclim_soil_UG.1"]][])

bio1 <- seq(quantile(bioclim_soil_UG[["bioclim_soil_UG.1"]][], probs=0.05, na.rm=T), quantile(bioclim_soil_UG[["bioclim_soil_UG.1"]][], probs=0.95, na.rm=T), len=np)
bio17 <- seq(quantile(bioclim_soil_UG[["bioclim_soil_UG.17"]][], probs=0.05, na.rm=T), quantile(bioclim_soil_UG[["bioclim_soil_UG.17"]][], probs=0.95, na.rm=T), len=np)


#h?ufigsten Bodentyp herausfinden
soilpresences <- extract(bioclim_soil_UG[["bioclim_soil_UG.20"]], presences)
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
#Modalwert: 5



# bei Boden h?ufigsten Bodentyp festlegen (Modalwert)
newdata <- expand.grid(bioclim_soil_UG.1 = bio1, bioclim_soil_UG.17 = bio17, bioclim_soil_UG.20 = 5) #h?ufigster Hauptbodentyp mit Nummer 5
newdata$pred <- predict(me3, newdata)
plot(newdata$pred)


# Darstellung der Antwortkurve 

### Use threshold to show distribution
newdata$pred[newdata$pred<thr$sensitivity] <- NA


### Create classes of site suitability
cInt <- classIntervals((newdata$pred))

xdiff <- diff(unique(newdata$bioclim_soil_UG.1))[1]
ydiff <- diff(unique(newdata$bioclim_soil_UG.17))[1]

mypalette <- colorRampPalette(c("lightgreen", "darkgreen"))
newdata$colors <- findColours(cInt, mypalette(length(cInt$brks)))

pdf("Ergebnisse/Modell_mit_Boden/Umweltraum.pdf", width=7, height=5)
par(mfrow = c(1,1), mar = c(5,5,1,1))
symbols(x = newdata$bioclim_soil_UG.1, y = newdata$bioclim_soil_UG.17, rectangles = matrix(rep(c(xdiff, ydiff), nrow(newdata)), ncol = 2, byrow = T), 
        bg = newdata$colors, 
        fg = "white", inches = F, 
        xlab = "Annual Mean Temperature (?C)", 
        ylab = "Precipitation of Driest Quarter (mm)", 
        ylim = c(0, 1200), axes = F)
contour(x = unique(newdata$bioclim_soil_UG.1), y = unique(newdata$bioclim_soil_UG.17), z = matrix(newdata$pred, nrow = np), add = T, levels = unique(round(cInt$brks,1)), labcex = 1.3)
axis(1, at = c(-100, 0, 100, 200, 300, 350), labels = c(-100/10, 0, 100/10, 200/10, 300/10, 350/10))
axis(2)
mtext("Carissa edulis", side = 3, line = -2.3, font = 3)
mtext(paste0("AUC = " , round(e3@auc, 2), " "), side = 3, line = -5, adj = 1)
mtext(paste0("Pearson r? = " , round(e3@cor, 2), " "), side = 3, line = -6.3, adj = 1)
dev.off()


### Q: Is this response curve reasonable? Why?


### Plot distribution map ####
par(mfrow = c(1, 1))
pred <- predict(me3, bioclim_soil_UG)

x11(width = 12)

plot(pred, xlab = "Longitude", ylab = "Latitude", xlim = c(-30,60), cex.axis = 1.2, cex.lab = 1.2, 
     legend.args = list(text='Probability of occurence', side = 4, font = 1, line = 2.8, cex = 1.2), 
     legend.width = 1, legend.shrink = 1)

distr <- pred
distr[distr < thr$sensitivity] <- NA
cInt <- classIntervals((newdata$pred))


# Distribution map mit Vorkommen
# plot(distr, col = mypalette(10), breaks = cInt$brks, legend = F, xlab = "Longitude", ylab = "Latitude")
x11(width = 12)
plot(distr, xlab = "Longitude", ylab = "Latitude", xlim = c(-30,50), cex.axis = 1.2, cex.lab = 1.2, 
     legend.args = list(text='Probability of occurence', side = 4, font = 1, line = 2.8, cex = 1.2), 
     legend.width = 1, legend.shrink = 1)
plot(UG, add = T)
points(presences, pch = 17, cex = 0.6, col = "red")

mtext("Carissa edulis", side = 1, line = -1.3, font = 3)
mtext(paste0("AUC = " , round(e3@auc, 2), " "), side = 1, line = -2.3, adj = 1)
mtext(paste0("Pearson r? = " , round(e3@cor, 2), " "), side = 1, line = -1.3, adj = 1)
legend("topleft", legend = c("presence data"), pch = 17, col = c("red"), border = "white", cex = 0.8, bty = "n", bg = "white")






# Modellierung der Reaktion der Artverbreitung auf verschiedene klimatische Ver?nderungen

# ?ber die Variablen rcp, model und year kann ein bestimmtes Klimaszenario erstellt werden

# 'model' should be one of "AC", "BC", "CC", "CE", "CN", "GF", "GD", "GS", "HD", "HG", "HE", "IN", "IP", "MI", "MR", "MC", "MP", "MG", or "NO"
# Globale Klimamodelle - GCM (General Circulation Model)

# 'rcp' should be one of 26, 45, 60, or 85
# Die RCPs repr?sentieren verschiedene Entwicklungspfade der Konzentrationen von Treibhausgasen, Aerosolen und zugeh?rige Emissionen
# Projektionen m?glicher Klima?nderungen im 21. Jahrhundert und dar?ber hinaus 
# RCP2.6 (RF relativ niedrig), RCP4.5 (RF mittel), RCP6.0 (RF hoch) und RCP8.5 (RF sehr hoch)
# Ver?nderung des gesamten anthropogenen Strahlungsantriebs (,radiative forcing' RF in W/m?)

# 'year' should be 50 or 70
# Time periods: 2050 (average for 2041-2060) and 2070 (average for 2061-2080)


# Klimaszenario f?r das Jahr 2070, RCP85 HD
#cc <- getData('CMIP5', var = "bio", res = 5, rcp = 85, model = 'HD', year = 70, download = TRUE)
#list.files()
cc <- stack(x = c("hd85bi701.tif", "hd85bi702.tif", "hd85bi703.tif", "hd85bi704.tif", "hd85bi705.tif", "hd85bi706.tif", 
                  "hd85bi707.tif", "hd85bi708.tif", "hd85bi709.tif", "hd85bi7010.tif", "hd85bi7011.tif", "hd85bi7012.tif",
                  "hd85bi7013.tif", "hd85bi7014.tif", "hd85bi7015.tif", "hd85bi7016.tif", "hd85bi7017.tif", "hd85bi7018.tif",
                  "hd85bi7019.tif"))
cc <- crop(cc, bioclim_UG)
cc <- mask(cc, bioclim_UG)
names(cc) <- names(bioclim_soil_UG)
writeRaster(cc, "cc.tif", overwrite = TRUE)

cc <- stack("cc.tif")

pred_cc <- predict(me3, cc)

distr_cc <- pred_cc
distr_cc[distr_cc[] < thr$sensitivity] <- NA



# Klimaszenario f?r das Jahr 2070, RCP26 HD
list.files()
ccc <- stack(x = c("hd26bi701.tif", "hd26bi702.tif", "hd26bi703.tif", "hd26bi704.tif", "hd26bi705.tif", "hd26bi706.tif", 
                   "hd26bi707.tif", "hd26bi708.tif", "hd26bi709.tif", "hd26bi7010.tif", "hd26bi7011.tif", "hd26bi7012.tif",
                   "hd26bi7013.tif", "hd26bi7014.tif", "hd26bi7015.tif", "hd26bi7016.tif", "hd26bi7017.tif", "hd26bi7018.tif",
                   "hd26bi7019.tif"))
ccc <- crop(ccc, bioclim_UG)
ccc <- mask(ccc, bioclim_UG)

names(ccc) <- names(bioclim_soil_UG)
writeRaster(ccc, "ccc.tif", overwrite = TRUE)
ccc <- stack("ccc.tif")

pred_ccc <- predict(me3, ccc)

distr_ccc <- pred_ccc
distr_ccc[distr_ccc[] < thr$sensitivity] <- NA

hist(distr_ccc)
hist(distr_cc)





# Graphische Darstellung der verschiedenen Klimaszenarien

### Create figure for distribution for current and climate change projection ####
x11(width = 12)

#pdf("figures/maps_Saxifraga_bryoides.pdf", width = 12, pointsize = 16)
par(mfrow = c(1, 2))

#plot(distr, col = mypalette(length(cInt$brks)), breaks = cInt$brks, legend = F, xlab = "Longitude", ylab = "Latitude")
#plot(europe, add = T)
#mtext("Current climate ", 3, -1.2, adj = 0)
#mtext("Saxifraga bryoides", 1, -1.2, font = 3)
#mtext("WGS84 ", 1, -1.2, adj = 1)
#?mtext

plot(distr_ccc, col = mypalette(length(cInt$brks)), breaks = cInt$brks, legend = F, xlab = "Longitude", ylab = "Latitude")
plot(UG, add = T)
mtext("RCP26 HD ", 3, -1.2, adj = 0)
mtext("Carissa edulis", 1, -1.2, font = 3)
mtext("WGS84 ", 1, -1.2, adj = 1)

plot(distr_cc, col = mypalette(length(cInt$brks)), breaks = cInt$brks, legend = F, xlab = "Longitude", ylab = "Latitude")
plot(UG, add = T)
mtext("RCP85 HD ", 3, -1.2, adj = 0)
mtext("Carissa edulis", 1, -1.2, font = 3)
mtext("WGS84 ", 1, -1.2, adj = 1)

dev.off()


