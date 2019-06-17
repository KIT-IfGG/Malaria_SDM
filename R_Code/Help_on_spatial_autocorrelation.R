
# set java number of cores to use
options(java.parameters = "-Xmx30g")#  

# load require libraries
library(rJava)
library(maxent)
library(dismo)
library(maptools)#
library(rgdal)
library(rgeos)
library(sp)
library(ENMeval)#evaluating and selecting best candidate model
library(knitr)#needed for kable function 
library(raster)#raster files
library(grDevices)#coloured map
library(ape)  # MoranÂ´s I
library(vegan)
library(corrplot)
library(ncf) # correlogramm (there are alternatives!)
library(spdep) # neighbours
library(classInt)
library(nlme)
library(spdep)# neighbors
library(colorRamps)#for some crispy colors

# workstation
options(java.parameters = "-Xmx30g")  
#my computer
options(java.parameters = "-Xmx30g")  

presence <- readRDS("data/combined_data/rhus_natalensis_combined.rds")
dim(presence)
names(presence)

#Bioclim data at coarse resolution can be dowloaded within R
bioclim <- getData("worldclim", download = F, path = "./data", var="bio", res=2.5)
projection(bioclim)

### Same projection but different projection string, fix
data("wrld_simpl")
projection(wrld_simpl) <- projection(bioclim)
proj4string(bioclim)=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
plot(wrld_simpl)
points(presence, col="blue", cex=0.5)

#presences only in Samburu and Africa
presences <- SpatialPoints(presence[,c("decimalLongitude", "decimalLatitude")], 
                           proj4string=CRS(projection("+proj=longlat +ellps=WGS84 +datum=WGS84")))
samburu <- readRDS("data/samburu_shp.rds")
plot(samburu, col="blue", axes=T)
samburu <- spTransform(samburu, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) 
presences_africa <- presences[africa]
presences_samburu <- presences[samburu]

#plot africa
africa <- wrld_simpl[wrld_simpl$REGION==2,]
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

#############################################
####MODELING WITH SOI

# Read data of the African Soil Atlas
soil_shp <- readOGR("data/Africa_soil_WRB" , layer = "afticasoilmap")

# limit area to study area
soil_shp <- crop(soil_shp, africa)
plot(soil_shp)
soil_shp@data
soil_leg <- soil_shp@data

### remove dublicaded to create a legend
dups <- duplicated(soil_leg["GRIDCODE", "SU_WRB1_PH"])

#Only keep record not duplicated
soil_leg <- soil_leg[!dups,]
summary(soil_leg)

#crop and mask to africa
soil_africa <- rasterize(soil_shp, bioclim_africa, field = "GRIDCODE")
soil_africa  <- crop(soil_africa, africa)
soil_africa   <- mask(soil_africa, africa)
plot(soil_africa)
writeRaster(soil_africa, "soil_africa_wrb.tif", overwrite = TRUE)
table(soil_africa[])

#soil data is reduced to main soil types 
soil <- stack("soil_africa_wrb.tif")
plot(soil_africa)

#rasterstack of bioclim_samburu and soil built with command or function stack
bioclim_soil_africa <- stack(bioclim_africa, soil)
writeRaster(bioclim_soil_africa, "bioclim_soil_africa.tif", overwrite = TRUE)

bioclim_soil_africa <- stack("bioclim_soil_africa.tif")
bioclim_soil_africa@layers
names(bioclim_soil_africa)
plot(bioclim_soil_africa[[20]])

#crop and mask to study area
soil_samb <- rasterize(soil_shp, bioclim_samburu, field = "GRIDCODE")
soil_samb  <- crop(soil_samb , samburu)
soil_samb  <- mask(soil_samb , samburu)
plot(soil)
writeRaster(soil_samb, "soil_samb_wrb.tif", overwrite = TRUE)
table(soil_samb[])

#soil data is reduced to main soil types 
soil <- stack("soil_samb_wrb.tif")
plot(soil_samb)

#rasterstack of bioclim_samburu and soil built with command or function stack
bioclim_soil_samburu <- stack(bioclim_samburu, soil)
writeRaster(bioclim_soil_samburu, "bioclim_soil_samburu.tif", overwrite = TRUE)

bioclim_soil_samburu <- stack("bioclim_soil_samburu.tif")
bioclim_soil_samburu@layers
names(bioclim_soil_samburu)
plot(bioclim_soil_samburu[[20]])

#create background points
set.seed(1000)
background <- randomPoints(bioclim_soil_africa, 20000)

#extracting values from raster
presvals <- extract(bioclim_soil_africa, presences_africa)
absvals <- extract(bioclim_soil_africa, background)

#creating predictor data set
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))
names(sdmdata)

#remove NA and check the correlation
sdmdata <- na.omit(sdmdata)
corr <- cor(na.omit(sdmdata[,2:21]), method = "pearson")
corrplot(corr, col = NULL, diag = T, type = "lower") 
pairs(sdmdata[,2:21], cex=0.1, fig=T)

#checking for collinearity of the variables
#retain one variable with higher ecological significance within the pairs
pairs(bioclim_soil_africa, hist=TRUE, cor=TRUE, use="pairwise.complete.obs", maxpixels=100000, method = "pearson")

#dividing the background into training and testing
# Now Selecting 75% of data as sample from total 'n' rows of the data  
sample <- sample.int(n = nrow(background), size = floor(.75*nrow(background)), replace = F)
background_train <- background[sample, ]
background_test  <-background[-sample, ]

predictors <- bioclim_soil_africa[[c("bioclim_soil_africa.12", "bioclim_soil_africa.16",  
                                     "bioclim_soil_africa.7", "bioclim_soil_africa.2", "bioclim_soil_africa.8",  "bioclim_soil_africa.20")]]
#use ENMeval function to assess candidate models
#this may take 10 to 15 minutes 
eval.results2 <- ENMevaluate(occ=presences_africa_train, env=predictors, categoricals = "bioclim_soil_africa.20", algorithm = "maxent.jar", background_train, n.bg = 20000, RMvalues=seq(0.5, 4, 0.5), fc = c("LQ", "L", "H", "LQHP", "LQH", "LQHPT"), method="block", clamp = T, "fadebyclamping=TRUE", numCores = 4, parallel = T)

#assessing the results
results2<- eval.results2@results
kable(results2)#the table of evaluation metrics. 
head(eval.results2@results)

#selecting settings for the more parsimonious model# with lowest AICc
aicmods2 <- which(eval.results2@results$AICc == min(na.omit(eval.results2@results$AICc)))
eval.results2@results[aicmods2,]#the most parsimonious method used the FCs =  LQ and RM=1.5

#creating logistic prediction based on first model
eval.results2@predictions[[46]]
pred <- predict(eval.results2@models[[46]], predictors)
plot(pred)

#Use the setting from the most parsimonious model to build a new Maxent model.
aicmods2 <- which(eval.results2@results$AICc == min(na.omit(eval.results2@results$AICc)))[1]
aicmods2 <- eval.results2@results[aicmods2,]
FC_best2 <- as.character(aicmods2$features[1]) # Get FCs from the AIC model
rm_best2 <- aicmods2$rm # Get RM from the AIC model

# code from dismo to get Maxent version
dismo.vs <- packageVersion('dismo')
v <- maxentJARversion()
alg <- paste("Maxent", v, "via dismo", dismo.vs)


maxent.args[[1]]

# code from dismo to get Maxent version
dismo.vs <- packageVersion('dismo')
v <- maxentJARversion()
alg <- paste("Maxent", v, "via dismo", dismo.vs)

predictors_samburu <- bioclim_soil_samburu[[c("bioclim_soil_samburu.12", "bioclim_soil_samburu.16",  
                                              "bioclim_soil_samburu.7", "bioclim_soil_samburu.2", "bioclim_soil_samburu.8",  "bioclim_soil_samburu.20")]]

#############################################
### Modeling
# Build the model_with_soil
maxent.args <- make.args(RMvalues = rm_best2, fc = FC_best2) # Define arguments
mx_best <- dismo::maxent(predictors, presences_africa_train, background_train, args=maxent.args[[1]], 
                         factors = "bioclim_soil_africa.20",  path = paste0(getwd(), "/Outputs/modelling_outputs/model_results/with_soil/rhus_natalensis"), removeDuplicates=F, overight=T)

#prediction to the study area
names(predictors_samburu) <- names(predictors)
distr_samb <- predict(mx_best, predictors_samburu, overwrite=TRUE, progress = 'text')
distr <- r_best

#save raster file of best model
writeRaster(distr_samb, filename="Outputs/modelling_outputs/raster_final/rhus_natalensis", format="GTiff", overwrite=TRUE)

e_model <- evaluate(p=presences_africa_test, a=background_test, model =mx_best,  x=bioclim_soil_africa)
e_model

kappa <- max(e_model@kappa); kappa2
TPR <- e_model@TPR[which.max( e_model@kappa )]; TPR2 # sensitivity
TNR <- e_model@TNR[which.max( e_model@kappa )]; TNR2#specificity
mat <- e_model@confusion[which.max( e_model@kappa ),];mat
tss <- num/den; tss# TSS = Sensitivity + specificity -1 
AUC <- e_model@auc;AUC

######## Make a list of all of the results obtained #############
accuracies <- data.frame(model=c(kappa, TPR, TNR, tss, AUC)); rownames(accuracies)<-c("Kappa", "TPR", "TNR", "TSS", "AUC"); accuracies

############CHECKING FOR SPATIAL AUTOCORRELATION
#I used a function already developed to check for the residuals and later calculate the MOrans I
#I also checked this uing correlogram in another script
presences_africa_train_sa <- as.data.frame(presences_africa_train)
library(spatstat)
presences_africa_train_sa <- as.data.frame(presences_africa_train)
# - ppm.model : Fitted IPP model
ipp <- function(species.data, predictors){
  # Function to convert raster to im
  raster.to.im <- function(raster){
    raster.mat <- as.matrix(raster)
    
    raster.mat <- apply(raster.mat, 2, rev)
    raster.im <- im(raster.mat,
                    xcol = seq(extent(raster)@xmin,
                               extent(raster)@xmax,
                               length=dim(raster)[2]),
                    yrow = seq(extent(raster)@ymin,
                               extent(raster)@ymax,
                               length=dim(raster)[1]))
    return(raster.im)
  }
  
  # Convert predictors from raster to im, place in covariates list
  covariates <- NULL
  for(i in 1:length(names(predictors))){
    covariates[[i]] <- raster.to.im(subset(predictors, i))
  }
  names(covariates) <- names(predictors)
  # Create trend formula
  trend <- "~"
  for(i in 1:length(names(covariates))){
    trend <- paste(trend, paste(names(covariates)[i], "+", sep=""),
                   sep="")
  }
  trend <- substr(trend, 1, nchar(trend)-1)
  trend <- as.formula(trend)
  # Rasterize species data
  species.grid <- rasterize(species.data, predictors)
  species.grid <- rasterToPoints(species.grid)[,-3]
  # Define ppp object
  W <- levelset(covariates[[1]], -100, ">")
  species.ppp <- ppp(species.grid[,1], species.grid[,2], window=W)
  print(bw.diggle(species.ppp))
  # Define quadrature scheme
  G <- gridcentres(W, nx=dim(predictors)[2], ny=dim(predictors)[1])
  dummy.ppp <- ppp(G$x, G$y, window=W)
  Q <- quadscheme(data=species.ppp, dummy=dummy.ppp)
  # Fit point process model
  ppm.model <- ppm(Q, trend, covariates=covariates)
  return(ppm.model)
}

# run functionb
ipp2 <- ipp(presences_africa_train_sa, predictors)
# - field : Smoothed residual field, raster
residuals <- function(maxent.model, species.data, predictors, type,
                      bandwidth){
  temp.mod <- ipp(species.data, predictors)
  temp.res <- diagnose.ppm(temp.mod,
                           which="smooth",
                           type=type,
                           sigma=bandwidth, plot.it=F)
  field <- temp.res$smooth$Z$v
  field <- apply(field, 2, rev)
  field <- raster(field)
  crs(field) <- crs(predictors)
  extent(field) <- extent(predictors)
  return(field)
}

resids2 <- residuals(mx_best, presences_africa_train_sa, predictors, type = "inverse", bandwidth=4)
plot(resids2)
points(presences_africa_train_sa)

res2_crop <- crop(resids2, samburu)
plot(res2_crop)
plot(samburu, add=T)

# - Cumulative residual plots
cum.res <- function(maxent.model, species.data, predictors){
  maxent.ipp <- ipp(species.data, predictors)
  diagnose <- diagnose.ppm(maxent.ipp)
  x <- diagnose$xcumul$empirical$covariate
  vx <- diagnose$xcumul$empirical$value
  y <- diagnose$ycumul$empirical$covariate
  vy <- diagnose$ycumul$empirical$value
  par(mfrow=c(1,2))
  plot(x, vx, type="l", xlab="Longitude",
       ylab="Cumulative Residual")
  plot(y, vy, type="l", xlab="Latitude",
       ylab="Cumulative Residual")
}

cum.res1<- cum.res(mx_best, presences_africa_train_sa, predictors)

# Calculation of Morans I
xy <- presences_africa_train
plot(xy)
X <- xy@coords[,1]
Y <- xy@coords[,2]

res <- extract(resids2, xy)

# delete NA values
res_na <- which(is.na(res)); res_na
res2 <- res[-res_na]
X <- X[-res_na]
Y <- Y[-res_na]

# distance matrix
xy.dists <- as.matrix(dist(cbind(X, Y)))
xy.dists.inv <- 1/xy.dists # invers
diag(xy.dists.inv) <- 0
xy.dists.inv[is.infinite(xy.dists.inv)] <- 0
xy.dists.inv[1:5, 1:5] # check

#calculate the morans I 
moran <- Moran.I(res2, xy.dists.inv); moran


#######REMOVING SPATIAL AUTOCORRELATION using using spatial Aigenvector Maps

matr_pred <- as.matrix(predictors)
matr_pred[] <- matr_pred

# Calculate response variable based on env
coeffs <- c(3, 0.9, 0.2)
z.value_lin <- coeffs[1] + coeffs[2] * matr_pred + coeffs[3] * matr_pred^2 
  my.mat[] <- z.value_lin
r.response <- raster(my.mat)

### Linear regression ####
dat_lin <- data.frame(matr_pred, z.value_lin)
m <- lm(z.value_lin ~ matr_pred + I(matr_pred^2), dat_lin)
round(coef(m), 1)

### -> Parameters used to create the data can be found by lm.

### Spatial effect ####
# all paiwise euclidean distances between the cells
xy.dist <- dist(presences_africa_train_sa)
# PCNM axes of the dist. matrix (from 'vegan' package)
dat <- data.frame(z.value = z.value, env.value = matr_pred)

pcnm.axes <- pcnm(xy.dist)$vectors
# using 8th PCNM axis as my atificial z variable
z.space <- pcnm.axes[,8]*100 
#z.space <- scale(z.space, center=TRUE, scale=FALSE)

# plotting the artificial spatial effect
my.mat[] <- z.space
r.space <- raster(my.mat)
plot(r.space, axes=F, col=matlab.like(20))

### Create data with environmental and spatial effects ####
z.value <- z.space + z.value_lin
my.mat[] <- z.value
r.z <- raster(my.mat)

# Identify neighbors
xy.neigh <- dnearneigh(as.matrix(presences_africa_train_sa), d1=0, d2=10, longlat = FALSE)
#Check if any point is not tethered
summary(xy.neigh)
xy.neigh_w <- nb2listw(xy.neigh)
options(error=recover)

res <- ME(z.value ~ matr_pred, dat, family = gaussian, listw = xy.neigh_w)

