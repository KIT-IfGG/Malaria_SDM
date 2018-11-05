
### Read the field data ####
sites <- read.table("data/field_data/sites_add.csv", sep=",", header=T, dec=".")
species <- read.table("data/field_data/species.csv", sep=",", header=T, dec=".", row.names = 1)

species <- na.omit(species)

### Coords for species presences ####
colnames(species)
sites$site_id

### Get rid of the "X" in colnames(species)
site_id_species <- substr(colnames(species), start = 2, stop = 1000000L)

sites <- sites[,c("site_id", "x", "y")]

medical_plans_samburu <- data.frame(species=NULL, site_id = NULL, pa = NULL)
for(i in 1:nrow(species)){
  print(rownames(species)[i])
  medical_plans_samburu <- rbind.data.frame(medical_plans_samburu, data.frame(species = rownames(species)[i], site_id = site_id_species[species[i,] > 0], pa = 1))
  
  medical_plans_samburu <- rbind.data.frame(medical_plans_samburu, data.frame(species = rownames(species)[i], site_id = site_id_species[species[i,] == 0], pa = 0))
  
}

m <- match(medical_plans_samburu$site_id, sites$site_id)

medical_plans_samburu$x <- sites$x[m]
medical_plans_samburu$y <- sites$y[m]
medical_plans_samburu <- medical_plans_samburu[complete.cases(medical_plans_samburu),]

### Check by plotitng ####
single_species <- medical_plans_samburu[medical_plans_samburu$species == "Myrsine africana",]
plot(single_species[,c("x", "y")], pch=16, col=ifelse(single_species$pa, "green", "grey"))

### Write pa dataset
write.csv(medical_plans_samburu, file = "data/medical_plants_samburu_ma.csv")


