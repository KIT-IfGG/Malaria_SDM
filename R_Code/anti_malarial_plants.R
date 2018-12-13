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
presence1 <- readRDS("data/combined_data/carissa_edulis_combined.rds")
presence2 <- readRDS("data/combined_data/lippia javanica_combined.rds")
presence3 <- readRDS("data/combined_data/lycium_europaeum_combined.rds")
presence4 <- readRDS("data/combined_data/rhus_natalensis_combined.rds")
presence5 <- readRDS("data/combined_data/solanum_incanum_combined.rds")
presence6 <- readRDS("data/combined_data/myrsine_africana_combined.rds")
presence7 <- readRDS("data/combined_data/acacia xanthophloea_combined.rds")
presence8 <- readRDS("data/combined_data/acacia_mellifera_combined.rds")
presence9 <- readRDS("data/combined_data/acacia_tortilis_combined.rds")
presence10 <- readRDS("data/combined_data/ajuga_remota_combined.rds")
presence11 <- readRDS("data/combined_data/balanites_aegyptiaca_combined.rds")
presence12 <- readRDS("data/combined_data/boscia_coriacea_combined.rds")
presence13 <- readRDS("data/combined_data/croton_dichogamus_combined.rds")
presence14 <- readRDS("data/combined_data/croton_megalocarpus_combined.rds")
presence15 <- readRDS("data/combined_data/euclea_divinorum_combined.rds")
presence16 <- readRDS("data/combined_data/euphorbia_heterochroma_combined.rds")
presence17 <- readRDS("data/combined_data/gardenia_jovis-tonantis_combined.rds")
presence18 <- readRDS("data/combined_data/kedrostis_pseudogijef_combined.rds")
presence19 <- readRDS("data/combined_data/rhamnus_stado.rds")
presence20 <- readRDS("data/combined_data/salvadora_persica_combined.rds")
presence21 <- readRDS("data/combined_data/senna_didymobotrya_combined.rds")
presence22 <- readRDS("data/combined_data/teclea_simplicifolia.rds")
presence23 <- readRDS("data/combined_data/zanthoxylum_usambarense_combined.rds")


presence <- rbind(presence1, presence2, presence3, presence4, presence5, presence6, presence7, presence8, presence9, presence10, presence11, presence12, presence13, presence14, presence15, presence16, presence17, presence19, presence20, presence21, presence22, presence23)

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

#correlation among the predictor variables
x11()
pairs(sdmdata[,2:10], cex=0.1, fig=T)
#corrplot(cor(sdmdata[,1:19]), type = "lower", diag=F)

corr <- cor(na.omit(sdmdata[,2:14], method = "spearman"))
corrplot(corr, diag = T, type = "lower")

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


help(mtext)
#create a figure for the distribution of current and climate change projection
dev.off()
pdf("figures/anti_malarial_plants_cur_cc_map.pdf", width = 12, pointsize = 16)
par(mfrow=c(1,2))
plot(distr, col=mypalatte(length(cInt$brks)), breaks=cInt$brks, legend=F, xlab="Longitude", ylab="Latitude")
plot(africa, add=T)
mtext("Current climate", 3, -0.1)
mtext("Anti-malarial Plants", 3, -1.1, font = 3)
mtext("wgs84", 1, -1-2, adj = 1)
plot(distr_cc, col=mypalatte(length(cInt$brks)), breaks=cInt$brks, legend=F, xlab="Longitude", ylab="Latitude")
plot(africa, add=T)
mtext("Climate Scenaria RCP85 HD", 3, -0.1)
mtext("Anti-malarial Plants", 3, -1.1, font = 3)
mtext("wgs84", 1, -1.2, adj=1)
dev.off()
