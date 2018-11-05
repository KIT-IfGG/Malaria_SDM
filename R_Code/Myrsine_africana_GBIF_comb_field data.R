library(dismo)
library(methods)
library(raster)
library(sp)
library(rgdal)
library(rJava)
library(maptools)

#read in species data
species <- data.frame(read.csv("E:/medical_plant_samburu/data/field_data/species.csv", header = T, sep=",", stringsAsFactors = FALSE))

# only select Myrsine africana from the data set
Species_maf <- subset(species, species == "Myrsine africana")


# leere Variable sites_car einf?hren
sites_maf<-c(NULL)
# Die Schleife beginnt in Spalte 2
i=2
while (i <= ncol(Species_maf)){
  if(Species_maf[i] > 0){
    sites_maf<-cbind(sites_maf, colnames(Species_maf[i])) #enth?lt Spalten?berschrift(Site_id)
  }
  i=i+1 #eine Spalte nach rechts r?cken
}


sites <- data.frame(read.csv("E:/medical_plant_samburu/data/field_data/sites.csv", header = T, sep=",", dec=".", stringsAsFactors = F))

# ohne "X" angegeben sind und sonst  %in% kein Ergebnis liefert)
sites_maf <- gsub("X","", sites_maf) # number of sites with Myrsine African 
result <- subset(sites, sites$site_id  %in%  sites_maf)
names(result)


# x,y Koordinaten aus result auslesen und in data frame presences_DikkoJeffGafna zusammenfassen
# und ?berschriften identisch zu den ?berschriften von presences_GBIF festlegen
presences_DikkoJeffGafna <- data.frame(cbind(result$y, result$x))
colnames(presences_DikkoJeffGafna) <- c("decimalLatitude","decimalLongitude")


# Pr?senzdaten anzeigen (nur f?r Testbetrieb,kann sp?ter enfallen)
head(presences_DikkoJeffGafna)
# Number of points in the field presence records
nrow(presences_DikkoJeffGafna)


#dowloading species data from gbif

mam <- gbif('Myrsine', 'africana*', geo = T, removeZeros = T, download = T)

dim(mam)
names(mam)

#Duplicate records differentiating by sub species

dups <- duplicated(mam["lon", "lat"])

#Only keep record not duplicated
mafin <- mam[!dups,]
View(mafin)

#select lat and lon
names(mafin)
mafin <- mafin[,c(82,88)]

mafin <- na.omit(mafin)#not sure of this 

head(mafin)

#ensure the names of both data frames are the same
names(mafin)==names(presences_DikkoJeffGafna)
names(mafin)
names(presences_DikkoJeffGafna)
names(mafin)[1]=c("decimalLatitude")
names(mafin)[2]=c("decimalLongitude")


#total presence data of field data & gbif data
total_presence <- rbind(mafin,presences_DikkoJeffGafna)

presences <- rbind(presences_DikkoJeffGafna, mafin)
nrow(presences)

#plotting the presence data in map to check
library(maptools)
data("wrld_simpl")
plot(wrld_simpl, xlim=c(-80,70), ylim=c(-60,10), axes=T, col="red")
points(presences$decimalLongitude, presences$decimalLatitude, pch=20, cex=0.20, col="blue")

#Checking the number of points
nrow(presences)
head(presences)


#save the total presence data
saveRDS(presences, "E:/medical_plant_samburu/data/field_data/Myrsine_africana_combined.rds")

