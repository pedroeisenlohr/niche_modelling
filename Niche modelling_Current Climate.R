####################################################################################
####################################################################################
#######          NICHE MODELLING WITH BIOMOD2 USING       #######
#######    40 ENVIRONMENTAL VARIABLES (10-km RESOLUTION)  ####### 
#######               SUMMARIZED IN PCA AXES              #######
####################################################################################
####################################################################################


# Contact: Pedro V. Eisenlohr (pedrov.eisenlohr@gmail.com)


###################### Acknowledgments ##########################
### Dr. Guarino Colli's team of Universidade de Brasília.
### Dr. Diogo Souza Bezerra Rocha (Instituto de Pesquisas Jardim Botânico/RJ).
### Drª Marinez Ferreira de Siqueira (Instituto de Pesquisas Jardim Botânico/RJ).
### My students of Ecology Lab (http://pedroeisenlohr.webnode.com.br/laboratorio-de-ecologia/).




###########################################################################################################################################################
####### Environmental data source 
## Temperature and precipitation, solar radiation, water vapor pressure and wind speed: WorldClim 2.0 (http://worldclim.org/version2).
## PET and Aridity Index: Global Aridity and PET Database (http://www.cgiar-csi.org/data/global-aridity-and-pet-database)
## AET and Soil Water Stress: Global High-Resolution Soil-Water Balance (http://www.cgiar-csi.org/data/global-high-resolution-soil-water-balance#download).
## Relative Humidity: Climond (https://www.climond.org/RawClimateData.aspx).
###########################################################################################################################################################



###########################
## SET WORKING DIRECTORY ##
###########################

# Each user should adjust this!
setwd(choose.dir())
getwd()

############################################
## INSTALL AND LOAD THE REQUIRED PACKAGES ##
############################################

install.packages("biomod2", dep=T)
install.packages("colorRamps", dep=T)
install.packages("dismo", dep=T)
install.packages("dplyr", dep=T)
install.packages("maps", dep=T)
install.packages("maptools", dep=T)
install.packages("plotKML", dep=T)
install.packages("raster", dep=T)
install.packages("rgdal", dep=T)
install.packages("RStoolbox", dep=T)
install.packages("foreach", dep=T)
install.packages("doParallel", dep=T)

library(biomod2)
library(colorRamps)
library(dismo)
library(dplyr)
library(maps)
library(maptools)
library(plotKML)
library(raster)
library(rgdal)
library(RStoolbox)
library(foreach)
library(doParallel)


### Loading species occurrence data:
spp<-read.table(file.choose(),row.names=1,header=T,sep=",")
dim(spp)
edit(spp)

###################################################################
### Loading WorldClim 2.0 layers ###
###################################################################

bioclim <- list.files("./Environmental layers/Temperature and Precipitation", full.names=TRUE)
bio <-stack(bioclim)
plot(bio[[1]])


## Crop WorldClim layers
#  *********************
neotrop <- readOGR("./Shapefiles/ShapeNeo/neotropic.shp")
bio.wc2 <- mask(crop(bio,neotrop),neotrop)
bio.wc2
res(bio.wc2)
plot(bio.wc2[[1]])
names(bio.wc2)




##################################################################################
##################################################################################

######### The steps below are required only at the first time ####################

##################################################################################
##################################################################################



####################################################################
### Compiling other rasters to stack ###
####################################################################

#Solar Radiation:
solar.radiation <- list.files("./Environmental layers/Solar Radiation", pattern=".tif", full.names=TRUE)
solar.radiation <- stack(solar.radiation)
solar.radiation.mean <- mean(solar.radiation)
solar.radiation.max <- max(solar.radiation)
solar.radiation.min <- min(solar.radiation)
solar.radiation.mean <- mask(crop(solar.radiation.mean, neotrop),neotrop)
writeRaster(solar.radiation.mean,"SolarRadiationMean")
solar.radiation.max <- mask(crop(solar.radiation.max, neotrop),neotrop)
writeRaster(solar.radiation.max,"SolarRadiationMax")
solar.radiation.min <- mask(crop(solar.radiation.min, neotrop),neotrop)
writeRaster(solar.radiation.min,"SolarRadiationMin")
res(solar.radiation.mean)
plot(solar.radiation.mean)
res(solar.radiation.max)
plot(solar.radiation.max)
res(solar.radiation.min)
plot(solar.radiation.min)


#Water Vapor Pressure:
water.vapor.pressure <- list.files("./Environmental layers/Water Vapor Pressure", pattern=".tif", full.names=TRUE)
water.vapor.pressure <-stack(water.vapor.pressure)
water.vapor.pressure.mean <-mean(water.vapor.pressure)
water.vapor.pressure.max <-max(water.vapor.pressure)
water.vapor.pressure.min <-min(water.vapor.pressure)
water.vapor.pressure.mean <- mask(crop(water.vapor.pressure.mean, neotrop),neotrop)
writeRaster(water.vapor.pressure.mean,"WaterVaporPressureMean")
water.vapor.pressure.max <- mask(crop(water.vapor.pressure.max, neotrop),neotrop)
writeRaster(water.vapor.pressure.max,"WaterVaporPressureMax")
water.vapor.pressure.min <- mask(crop(water.vapor.pressure.min, neotrop),neotrop)
writeRaster(water.vapor.pressure.min,"WaterVaporPressureMin")
res(water.vapor.pressure.mean)
plot(water.vapor.pressure.mean)
res(water.vapor.pressure.max)
plot(water.vapor.pressure.max)
res(water.vapor.pressure.min)
plot(water.vapor.pressure.min)


#Wind Speed:
wind.speed <- list.files("./Environmental layers/Wind Speed", pattern=".tif", full.names=TRUE)
wind.speed <- stack(wind.speed)
wind.speed.mean <-mean(wind.speed)
wind.speed.max <-max(wind.speed)
wind.speed.min <-min(wind.speed)
wind.speed.mean <-mask(crop(wind.speed.mean, neotrop),neotrop)
writeRaster(wind.speed.mean, "WindSpeedMean")
wind.speed.max <-mask(crop(wind.speed.max, neotrop),neotrop)
writeRaster(wind.speed.max, "WindSpeedMax")
wind.speed.min <-mask(crop(wind.speed.min, neotrop),neotrop)
writeRaster(wind.speed.min, "WindSpeedMin")
res(wind.speed.mean)
plot(wind.speed.mean)
res(wind.speed.max)
plot(wind.speed.max)
res(wind.speed.min)
plot(wind.speed.min)


#Potential Evapotranspiration:
PET.1km <- raster("./Environmental layers/Potential Evapotranspiration/Global PET - Annual/PET_he_annual/pet_he_yr/w001001.adf")
PET.1km <- mask(crop(PET.1km,neotrop),neotrop)
PET.10km <- resample(PET.1km,bio.wc2)
writeRaster(PET.10km, "PET10km")
res(PET.10km)
plot(PET.10km)


#Aridity Index:
Aridity.1km <- raster("./Environmental layers/Global Aridity and PET database/Global Aridity - Annual/AI_annual/ai_yr/w001001.adf")
Aridity.1km <- mask(crop(Aridity.1km,neotrop),neotrop)
Aridity.10km <- resample(Aridity.1km,bio.wc2)
writeRaster(Aridity.10km, "Aridity10km")
res(Aridity.10km)
plot(Aridity.10km)


#Actual Evapotranspiration:
AET.1km <- raster("./Environmental layers/Global Soil Water Balance and AET/Mean Annual AET/AET_YR/aet_yr/w001001.adf")
AET.1km <- mask(crop(AET.1km,neotrop),neotrop)
AET.10km <- resample(AET.1km,bio.wc2)
writeRaster(AET.10km, "AET10km")
res(AET.10km)
plot(AET.10km)


#Soil Water Stress:
SWS.jan <-raster("./Environmental layers/Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_1/w001001.adf")
SWS.feb <-raster("./Environmental layers/Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_2/w001001.adf")
SWS.mar <-raster("./Environmental layers/Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_3/w001001.adf")
SWS.apr <-raster("./Environmental layers/Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_4/w001001.adf")
SWS.may <-raster("./Environmental layers/Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_5/w001001.adf")
SWS.jun <-raster("./Environmental layers/Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_6/w001001.adf")
SWS.jul <-raster("./Environmental layers/Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_7/w001001.adf")
SWS.aug <-raster("./Environmental layers/Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_8/w001001.adf")
SWS.sep <-raster("./Environmental layers/Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_9/w001001.adf")
SWS.oct <-raster("./Environmental layers/Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_10/w001001.adf")
SWS.nov <-raster("./Environmental layers/Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_11/w001001.adf")
SWS.dec <-raster("./Environmental layers/Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_12/w001001.adf")
SWS.stack <-stack(SWS.jan,SWS.feb,SWS.mar,SWS.apr,SWS.may,SWS.jun,SWS.jul,
				SWS.aug,SWS.sep,SWS.oct,SWS.nov,SWS.dec)

SWS.mean.1km <-mean(SWS.stack)
SWS.mean.1km <-mask(crop(SWS.mean.1km,neotrop),neotrop)
SWS.mean.10km <-resample(SWS.mean.1km, bio.wc2)
writeRaster(SWS.mean.10km,"SWSmean10km")
res(SWS.mean.10km)
plot(SWS.mean.10km)

SWS.max.1km <-max(SWS.stack)
SWS.max.1km <-mask(crop(SWS.max.1km,neotrop),neotrop)
SWS.max.10km <-resample(SWS.max.1km, bio.wc2)
writeRaster(SWS.max.10km,"SWSmax10km")
res(SWS.max.10km)
plot(SWS.max.10km)

SWS.min.1km <-min(SWS.stack)
SWS.min.1km <-mask(crop(SWS.min.1km,neotrop),neotrop)
SWS.min.10km <-resample(SWS.min.1km, bio.wc2)
writeRaster(SWS.min.10km,"SWSmin10km")
res(SWS.min.10km)
plot(SWS.min.10km)


#Relative Humidity at 3pm:
Humidity.3pm.jan <-raster("./Environmental layers/Relative Humidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm01/w001001.adf")
Humidity.3pm.feb <-raster("./Environmental layers/Relative Humidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm02/w001001.adf")
Humidity.3pm.mar <-raster("./Environmental layers/Relative Humidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm03/w001001.adf")
Humidity.3pm.apr <-raster("./Environmental layers/Relative Humidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm04/w001001.adf")
Humidity.3pm.may <-raster("./Environmental layers/Relative Humidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm05/w001001.adf")
Humidity.3pm.jun <-raster("./Environmental layers/Relative Humidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm06/w001001.adf")
Humidity.3pm.jul <-raster("./Environmental layers/Relative Humidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm07/w001001.adf")
Humidity.3pm.aug <-raster("./Environmental layers/Relative Humidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm08/w001001.adf")
Humidity.3pm.sep <-raster("./Environmental layers/Relative Humidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm09/w001001.adf")
Humidity.3pm.oct <-raster("./Environmental layers/Relative Humidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm10/w001001.adf")
Humidity.3pm.nov <-raster("./Environmental layers/Relative Humidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm11/w001001.adf")
Humidity.3pm.dec <-raster("./Environmental layers/Relative Humidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm12/w001001.adf")
Humidity.3pm.stack <-stack(Humidity.3pm.jan, Humidity.3pm.feb, Humidity.3pm.mar, Humidity.3pm.apr, Humidity.3pm.may, Humidity.3pm.jun, Humidity.3pm.jul, 
				Humidity.3pm.aug, Humidity.3pm.sep, Humidity.3pm.oct, Humidity.3pm.nov, Humidity.3pm.dec)

Humidity.3pm.mean.20km <-mean(Humidity.3pm.stack)
Humidity.3pm.mean.20km <-mask(crop(Humidity.3pm.mean.20km,neotrop),neotrop)
Humidity.3pm.mean.10km <-resample(Humidity.3pm.mean.20km, bio.wc2)
writeRaster(Humidity.3pm.mean.10km,"Humidity3pmMean10km")
res(Humidity.3pm.mean.10km)
plot(Humidity.3pm.mean.10km)

Humidity.3pm.max.20km <-max(Humidity.3pm.stack)
Humidity.3pm.max.20km <-mask(crop(Humidity.3pm.max.20km,neotrop),neotrop)
Humidity.3pm.max.10km <-resample(Humidity.3pm.max.20km, bio.wc2)
writeRaster(Humidity.3pm.max.10km,"Humidity3pmMax10km")
res(Humidity.3pm.max.10km)
plot(Humidity.3pm.max.10km)

Humidity.3pm.min.20km <-min(Humidity.3pm.stack)
Humidity.3pm.min.20km <-mask(crop(Humidity.3pm.min.20km,neotrop),neotrop)
Humidity.3pm.min.10km <-resample(Humidity.3pm.min.20km, bio.wc2)
writeRaster(Humidity.3pm.min.10km,"Humidity3pmMin10km")
res(Humidity.3pm.min.10km)
plot(Humidity.3pm.min.10km)


#Relative Humidity at 9am:
Humidity.9am.jan <-raster("./Environmental layers/Relative Humidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham01/w001001.adf")
Humidity.9am.feb <-raster("./Environmental layers/Relative Humidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham02/w001001.adf")
Humidity.9am.mar <-raster("./Environmental layers/Relative Humidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham03/w001001.adf")
Humidity.9am.apr <-raster("./Environmental layers/Relative Humidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham04/w001001.adf")
Humidity.9am.may <-raster("./Environmental layers/Relative Humidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham05/w001001.adf")
Humidity.9am.jun <-raster("./Environmental layers/Relative Humidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham06/w001001.adf")
Humidity.9am.jul <-raster("./Environmental layers/Relative Humidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham07/w001001.adf")
Humidity.9am.aug <-raster("./Environmental layers/Relative Humidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham08/w001001.adf")
Humidity.9am.sep <-raster("./Environmental layers/Relative Humidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham09/w001001.adf")
Humidity.9am.oct <-raster("./Environmental layers/Relative Humidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham10/w001001.adf")
Humidity.9am.nov <-raster("./Environmental layers/Relative Humidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham11/w001001.adf")
Humidity.9am.dec <-raster("./Environmental layers/Relative Humidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham12/w001001.adf")
Humidity.9am.stack <-stack(Humidity.9am.jan, Humidity.9am.feb, Humidity.9am.mar, Humidity.9am.apr, Humidity.9am.may, Humidity.9am.jun, Humidity.9am.jul, 
				Humidity.9am.aug, Humidity.9am.sep, Humidity.9am.oct, Humidity.9am.nov, Humidity.9am.dec)

Humidity.9am.mean.20km <-mean(Humidity.9am.stack)
Humidity.9am.mean.20km <-mask(crop(Humidity.9am.mean.20km,neotrop),neotrop)
Humidity.9am.mean.10km <-resample(Humidity.9am.mean.20km, bio.wc2)
writeRaster(Humidity.9am.mean.10km,"Humidity9amMean10km")
res(Humidity.9am.mean.10km)
plot(Humidity.9am.mean.10km)

Humidity.9am.max.20km <-max(Humidity.9am.stack)
Humidity.9am.max.20km <-mask(crop(Humidity.9am.max.20km,neotrop),neotrop)
Humidity.9am.max.10km <-resample(Humidity.9am.max.20km, bio.wc2)
writeRaster(Humidity.9am.max.10km,"Humidity9amMax10km")
res(Humidity.9am.max.10km)
plot(Humidity.9am.max.10km)

Humidity.9am.min.20km <-min(Humidity.9am.stack)
Humidity.9am.min.20km <-mask(crop(Humidity.9am.min.20km,neotrop),neotrop)
Humidity.9am.min.10km <-resample(Humidity.9am.min.20km, bio.wc2)
writeRaster(Humidity.9am.min.10km,"Humidity9amMin10km")
res(Humidity.9am.min.10km)
plot(Humidity.9am.min.10km)


##################################################################################
##################################################################################

######### The steps above are required only at the first time ####################

##################################################################################
##################################################################################




######################################################################
#################### Automatizing the process above ##################
######################################################################
### Only in case if you have already processed the steps above:
solar.radiation.mean <-raster("./Environmental layers/Solar Radiation/SolarRadiationMean.grd")
solar.radiation.max <-raster("./Environmental layers/Solar Radiation/SolarRadiationMax.grd")
solar.radiation.min <-raster("./Environmental layers/Solar Radiation/SolarRadiationMin.grd")
water.vapor.pressure.mean<-raster("./Environmental layers/Water Vapor Pressure/WaterVaporPressureMean.grd")
water.vapor.pressure.max <-raster("./Environmental layers/Water Vapor Pressure/WaterVaporPressureMax.grd")
water.vapor.pressure.min <-raster("./Environmental layers/Water Vapor Pressure/WaterVaporPressureMin.grd")
wind.speed.mean <-raster("./Environmental layers/Wind Speed/WindSpeedMean.grd")
wind.speed.max <-raster("./Environmental layers/Wind Speed/WindSpeedMax.grd")
wind.speed.min <-raster("./Environmental layers/Wind Speed/WindSpeedMin.grd")
PET.10km <-raster("./Environmental layers/Potential Evapotranspiration/Global PET - Annual/PET10km.grd")
Aridity.10km <-raster("./Environmental layers/Global Aridity/Global Aridity - Annual/Aridity10km")
AET.10km <-raster("./Environmental layers/Actual Evapotranspiration/Mean Annual AET/AET10km.grd")
SWS.mean.10km <-raster("./Environmental layers/Soil Water Stress/Monthly Soil Water Stress/SWSmean10km.grd")
SWS.max.10km <-raster("./Environmental layers/Soil Water Stress/Monthly Soil Water Stress/SWSmax10km.grd")
SWS.min.10km <-raster("./Environmental layers/Soil Water Stress/Monthly Soil Water Stress/SWSmin10km.grd")
Humidity.3pm.mean.10km <-raster("./Environmental layers/Relative Humidity 3pm/Humidity3pmMean10km.grd")
Humidity.3pm.min.10km <-raster("./Environmental layers/Relative Humidity 3pm/Humidity3pmMin10km.grd")
Humidity.3pm.max.10km <-raster("./Environmental layers/Relative Humidity 3pm/Humidity3pmMax10km.grd")
Humidity.9am.mean.10km <-raster("./Environmental layers/Relative Humidity 9am/Humidity9amMean10km.grd")
Humidity.9am.max.10km <-raster("./Environmental layers/Relative Humidity 9am/Humidity9amMax10km.grd")
Humidity.9am.min.10km <-raster("./Environmental layers/Relative Humidity 9am/Humidity9amMin10km.grd")



#########################################################################
############### Stacking all environmental layers #######################
#########################################################################
bio.crop <- stack(bio.wc2, 
	solar.radiation.mean, solar.radiation.max, solar.radiation.min, 
	water.vapor.pressure.mean, water.vapor.pressure.max, water.vapor.pressure.min, 
	wind.speed.mean, wind.speed.max, wind.speed.min, 
	PET.10km, Aridity.10km, AET.10km, 
	SWS.mean.10km, SWS.min.10km, SWS.max.10km,
	Humidity.3pm.mean.10km, Humidity.3pm.min.10km, Humidity.3pm.max.10km, Humidity.9am.mean.10km, Humidity.9am.max.10km, Humidity.9am.min.10km)
bio.crop
res(bio.crop)


#################################################################
##################### PCA #######################################
#################################################################

env.selected1 <- rasterPCA(bio.crop, nComp=8, spca = TRUE)
# Here I selected the first 8 components because they account for about 95% 
# of the total variance considering the 40 predictors of this routine for the 
# entire Neotropical Region.

summary(env.selected1$model)
plot(env.selected1$model)
env.selected <-stack(env.selected1$map)
env.selected
plot(env.selected)
names(env.selected)


# Selecting spatially unique records #
mask <-bio.crop[[1]]
cell <- cellFromXY(mask, spp[,1:2]) # get the cell number for each point
cell
dup <- duplicated(cbind(spp[,1:2],cell))
spp <- spp[!dup, ]# select the records that are not duplicated
dim(spp)


### Creating other required objects for BIOMOD_Formating Data:
# Convert dataset to SpatialPointsDataFrame (only presences)
myRespXY <- spp[,c("long","lat")]
# Creating occurrence data object
occurrence.resp <-  rep(1, length(myRespXY$long))



#################################################
## BUILDING THE MODELS ##
#################################################

setwd("C:/Models")

### for example, number of species occurrence records = 93
# Prepare data

sppBiomodData.PA.equal <- BIOMOD_FormatingData(
	resp.var = occurrence.resp,
	expl.var = env.selected,
	resp.xy = myRespXY,
	resp.name = "Occurrence",
	PA.nb.rep = 10,
	PA.nb.absences = 93,
	PA.strategy = "sre",
	PA.sre.quant = 0.025)
sppBiomodData.PA.equal

sppBiomodData.PA.10000 <- BIOMOD_FormatingData(
	resp.var = occurrence.resp,
	expl.var = env.selected,
	resp.xy = myRespXY,
	resp.name = "Occurrence",
	PA.nb.rep = 10,
	PA.nb.absences = 10000,
	PA.strategy = "sre",
	PA.sre.quant = 0.025)
sppBiomodData.PA.10000


# MaxEnt .jar
  jar <- paste0(system.file(package = "dismo"), "/java/maxent.jar")
  if (file.exists(jar) != T) {
    url = "http://biodiversityinformatics.amnh.org/open_source/maxent/maxent.php?op=download"
    download.file(url, dest = "maxent.zip", mode = "wb")
    unzip("maxent.zip", files = "maxent.jar", exdir = system.file("java", package = "dismo"))
    unlink("maxent.zip")
    warning("Maxent foi colocado no diretório")
  } 
system.file("java", package = "dismo")

# Path to MAXENT
myBiomodOption <- BIOMOD_ModelingOptions(MAXENT.Phillips = list(path_to_maxent.jar="C:/Users/R/win-library/3.3/dismo/java"))


# In this case, 70% of data will be used to train and 30% will be used to test the model:
### Note that keeping VarImport = 100 implies in a very time-consuming process! Unless obtaining the relative importance of
### predictors is a priority for you, I suggest you changing '100' to 'FALSE'.
sppModelOut.PA.equal <- BIOMOD_Modeling(sppBiomodData.PA.equal, 
	models = c("GBM", "CTA", "RF"), 
	NbRunEval = 10,
	DataSplit = 70, 
	Prevalence = NULL, 
	VarImport = 100,
	models.eval.meth = c("KAPPA","TSS","ROC","SR","POD","ACCURACY","BIAS"),
	SaveObj = TRUE,
	rescal.all.models = FALSE,
	do.full.models = FALSE,
	modeling.id = "spp")
sppModelOut.PA.equal


# Parallel processing
cl <- makeCluster(detectCores()) # number of cores in computer
registerDoParallel(cl)
getDoParWorkers()

### Note that keeping VarImport = 100 implies in a very time-consuming process! Unless obtaining the relative importance of
### predictors is a priority for you, I suggest you changing '100' to 'FALSE'.
sppModelOut.PA.10000 <- BIOMOD_Modeling( 
	sppBiomodData.PA.10000, 
	models = c("GLM","GAM","ANN","SRE","FDA","MARS","MAXENT.Phillips"), 
	models.options = myBiomodOption, 
	NbRunEval = 10,
	DataSplit = 70, 
	Prevalence = NULL, 
	VarImport = 100,
	models.eval.meth = c("KAPPA","TSS","ROC","SR","POD","ACCURACY","BIAS"),
	SaveObj = TRUE,
	rescal.all.models = FALSE,
	do.full.models = FALSE,
	modeling.id = "spp")
sppModelOut.PA.10000


###################################
####### EVALUATING MODELS #########
###################################

# Get evaluations
sppModelEval.PA.equal <- get_evaluations(sppModelOut.PA.equal)
sppModelEval.PA.equal
sppVarsEval.PA.equal <- get_variables_importance(sppModelOut.PA.equal)
sppVarsEval.PA.equal
write.table(sppVarsEval.PA.equal, "VariablesImportance1.csv")

sppModelEval.PA.10000 <- get_evaluations(sppModelOut.PA.10000)
sppModelEval.PA.10000
sppVarsEval.PA.10000 <- get_variables_importance(sppModelOut.PA.10000)
sppVarsEval.PA.10000
write.table(sppVarsEval.PA.10000, "VariablesImportance2.csv")

# Get summaries (mean and std.dev.) of model evaluation - 1
sdm.models1 <- c("GBM","CTA","RF") #3 models
sdm.models1
eval.methods1 <- c("KAPPA","TSS","ROC","SR","FAR","ACCURACY","BIAS") #7 evaluation methods
eval.methods1

means.i <- numeric(0)
means.j <- numeric(7)
for (i in 1:3){
	for (j in 1:7){
	means.j[j] <- mean(sppModelEval.PA.equal[paste(eval.methods1[j]),"Testing.data",paste(sdm.models1[i]),,])
	}
	means.i <- c(means.i, means.j)
}

summary.eval.equal <- data.frame(rep(sdm.models1,each=7), rep(eval.methods1,3), means.i)
names(summary.eval.equal) <- c("Model", "Method", "Mean")
summary.eval.equal
write.table(summary.eval.equal,"Models1_Evaluation_Mean.csv")

sd.i <- numeric(0)
sd.j <- numeric(7)
for (i in 1:3){
	for (j in 1:7){
	sd.j[j] <- sd(sppModelEval.PA.equal[paste(eval.methods1[j]),"Testing.data",paste(sdm.models1[i]),,])
	}
	sd.i <- c(sd.i, sd.j)
}

summary.eval.equal <- data.frame(rep(sdm.models1,each=7), rep(eval.methods1,3), sd.i)
names(summary.eval.equal) <- c("Model", "Method", "SD")
summary.eval.equal
write.table(summary.eval.equal,"Models1_Evaluation_SD.csv")


# Get summaries (mean and std.dev.) of model evaluation - 2
sdm.models2 <- c("GLM","GAM","ANN","SRE","MARS","MAXENT.Phillips","FDA") #7 models
sdm.models2
eval.methods2 <- c("KAPPA","TSS","ROC","SR","FAR","ACCURACY","BIAS") #7 evaluation methods
eval.methods2

means.i <- numeric(0)
means.j <- numeric(7)
for (i in 1:7){
	for (j in 1:7){
	means.j[j] <- mean(sppModelEval.PA.10000[paste(eval.methods2[j]),"Testing.data",paste(sdm.models2[i]),,], na.rm=T)
	}
	means.i <- c(means.i, means.j)
}

summary.eval.10000 <- data.frame(rep(sdm.models2,each=7), rep(eval.methods2,7), means.i)
names(summary.eval.10000) <- c("Model", "Method", "Mean")
summary.eval.10000
write.table(summary.eval.10000,"Models2_Evaluation_Mean.csv")

sd.i <- numeric(0)
sd.j <- numeric(7)
for (i in 1:7){
	for (j in 1:7){
	sd.j[j] <- sd(sppModelEval.PA.10000[paste(eval.methods2[j]),"Testing.data",paste(sdm.models2[i]),,])
	}
	sd.i <- c(sd.i, sd.j)
}

summary.eval.10000 <- data.frame(rep(sdm.models2,each=7), rep(eval.methods2,7), sd.i)
names(summary.eval.10000) <- c("Model", "Method", "SD")
summary.eval.10000
write.table(summary.eval.10000,"Models2_Evaluation_SD.csv")



### Which algorithms should be retained for the ensemble model?



#################################
## PRODUCING MODEL PROJECTIONS ##
#################################

spp.projections_1 <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.equal,
	new.env = env.selected,
	proj.name = "Cur1",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")

spp.projections_2 <- BIOMOD_Projection(
	modeling.output = sppModelOut.PA.10000,
	new.env = env.selected,
	proj.name = "Cur2",
	selected.models = "all",
	binary.meth = "ROC",
	output.format = ".grd")


# Stack projections
### You should define where is the file 'proj_Cur1_Occurrence.grd'
projections_cont1 <-stack("./Occurrence/proj_Cur1/proj_Cur1_Occurrence.grd") #onde estão os modelos produzidos por GBM, CTA e RF
names(projections_cont1)

### You should define where is the file 'proj_Cur1_Occurrence_ROCbin.grd'
projections_bin1 <-stack("./Occurrence/proj_Cur1/proj_Cur1_Occurrence_ROCbin.grd") #onde estão os modelos produzidos por GBM, CTA e RF
names(projections_bin1)

names(projections_bin1)=names(projections_cont1)
names(projections_bin1)
plot(projections_bin1)


### You should define the path to the file 'proj_Cur2_Occurrence.grd'
projections_cont2 <-stack("./Occurrence/proj_Cur2/proj_Cur2_Occurrence.grd") #onde estão os modelos produzidos pelos demais algoritmos
names(projections_cont2)

### You should define the path to the file 'proj_Cur2_Occurrence_ROCbin.grd'
projections_bin2 <-stack("./Occurrence/proj_Cur2/proj_Cur2_Occurrence_ROCbin.grd") #onde estão os modelos produzidos pelos demais algoritmos
names(projections_bin2)

names(projections_bin2)=names(projections_cont2)
names(projections_bin2)
plot(projections_bin2)


### Apply the steps below only for retained algorithms
### Ensemble mean model for each algorithm:
projections.RF.all <- subset(projections_bin1, grep("RF", names(projections_bin1)))
projections.RF.mean <- mean(projections.RF.all)
windows(w=6, h=6)
plot(projections.RF.mean, col = matlab.like(100), main = "RF", las = 1)
writeRaster(projections.RF.mean,"RF")

projections.GBM.all <-subset(projections_bin1, grep("GBM", names(projections_bin1)))
projections.GBM.mean <- mean(projections.GBM.all)
windows(w=6, h=6)
plot(projections.GBM.mean, col = matlab.like(100), main = "GBM", las = 1)

projections.CTA.all <-subset(projections_bin1,grep("CTA", names(projections_bin1)))
projections.CTA.mean <- mean(projections.CTA.all)
windows(w=6, h=6)
plot(projections.CTA.mean, col = matlab.like(100), main = "CTA", las = 1)

projections.GLM.all <-subset(projections_bin2,grep("GLM", names(projections_bin2)))
projections.GLM.mean <- mean(projections.GLM.all)
windows(w=6, h=6)
plot(projections.GLM.mean, col = matlab.like(100), main = "GLM", las = 1)

projections.GAM.all <-subset(projections_bin2,grep("GAM", names(projections_bin2)))
projections.GAM.mean <- mean(projections.GAM.all)
windows(w=6, h=6)
plot(projections.GAM.mean, col = matlab.like(100), main = "GAM", las = 1)

projections.ANN.all <- subset(projections_bin2,grep("ANN", names(projections_bin2)))
projections.ANN.mean <- mean(projections.ANN.all)
windows(w=6, h=6)
plot(projections.ANN.mean, col = matlab.like(100), main = "ANN", las = 1)

projections.SRE.all <- subset(projections_bin2,grep("SRE", names(projections_bin2)))
projections.SRE.mean <- mean(projections.SRE.all)
windows(w=6, h=6)
plot(projections.SRE.mean, col = matlab.like(100), main = "SRE", las = 1)

projections.MARS.all <- subset(projections_bin2,grep("MARS", names(projections_bin2)))
projections.MARS.mean <- mean(projections.MARS.all)
windows(w=6, h=6)
plot(projections.MARS.mean, col = matlab.like(100), main = "MARS", las = 1)

projections.FDA.all <- subset(projections_bin2,grep("FDA", names(projections_bin2)))
projections.FDA.mean <- mean(projections.FDA.all)
windows(w=6, h=6)
plot(projections.FDA.mean, col = matlab.like(100), main = "FDA", las = 1)

projections.MAXENT.all <- subset(projections_bin2,grep("MAXENT.Phillips", names(projections_bin2)))
projections.MAXENT.mean <- mean(projections.MAXENT.all)
windows(w=6, h=6)
plot(projections.MAXENT.mean, col = matlab.like(100), main = "MAXENT", las = 1)


############################################################
#### IF YOU PREFER TO OBTAIN ONLY CONTINUOUS PROJECTIONS ###
############################################################

### Apply the steps below only for retained algorithms
### Ensemble mean model for each algorithm:
projections.RF.all <- subset(projections_cont1, grep("RF", names(projections_cont1)))
projections.RF.mean <- mean(projections.RF.all)
plot(projections.RF.mean, col = matlab.like(100), main = "RF", las = 1)

projections.GBM.all <-subset(projections_cont1, grep("GBM", names(projections_cont1)))
projections.GBM.mean <- mean(projections.GBM.all)
plot(projections.GBM.mean, col = matlab.like(100), main = "GBM", las = 1)

projections.CTA.all <-subset(projections_cont1,grep("CTA", names(projections_cont1)))
projections.CTA.mean <- mean(projections.CTA.all)
plot(projections.CTA.mean, col = matlab.like(100), main = "CTA", las = 1)

projections.GLM.all <-subset(projections_cont2,grep("GLM", names(projections_cont2)))
projections.GLM.mean <- mean(projections.GLM.all)
plot(projections.GLM.mean, col = matlab.like(100), main = "GLM", las = 1)

projections.GAM.all <-subset(projections_cont2,grep("GAM", names(projections_cont2)))
projections.GAM.mean <- mean(projections.GAM.all)
plot(projections.GAM.mean, col = matlab.like(100), main = "GAM", las = 1)

projections.ANN.all <- subset(projections_cont2,grep("ANN", names(projections_cont2)))
projections.ANN.mean <- mean(projections.ANN.all)
plot(projections.ANN.mean, col = matlab.like(100), main = "ANN", las = 1)

projections.SRE.all <- subset(projections_cont2,grep("SRE", names(projections_cont2)))
projections.SRE.mean <- mean(projections.SRE.all)
plot(projections.SRE.mean, col = matlab.like(100), main = "SRE", las = 1)

projections.MARS.all <- subset(projections_cont2,grep("MARS", names(projections_cont2)))
projections.MARS.mean <- mean(projections.MARS.all)
plot(projections.MARS.mean, col = matlab.like(100), main = "MARS", las = 1)

projections.FDA.all <- subset(projections_cont2,grep("FDA", names(projections_cont2)))
projections.FDA.mean <- mean(projections.FDA.all)
plot(projections.FDA.mean, col = matlab.like(100), main = "FDA", las = 1)

projections.MAXENT.all <- subset(projections_cont2,grep("MAXENT.Phillips", names(projections_cont2)))
projections.MAXENT.mean <- mean(projections.MAXENT.all)
plot(projections.MAXENT.mean, col = matlab.like(100), main = "MAXENT", las = 1)



##############################################################
################## ENSEMBLE MODEL ############################
##############################################################

projections.all.mean <- mean(projections.RF.mean + projections.GBM.mean +
	projections.CTA.mean + projections.GLM.mean + projections.GAM.mean +
	projections.ANN.mean + projections.SRE.mean + projections.MARS.mean + 
	projections.FDA.mean + projections.MAXENT.mean)
windows(w=6, h=6)
plot(projections.all.mean, col = matlab.like(100), main = "Ensemble - Current Climate", las = 1)


# If you decide to include species occurrence records on the map:
df.pts <- data.frame(Occurrence = rep(1, length(myRespXY$long)), myRespXY)
dim(df.pts)
spp.spdf <- SpatialPointsDataFrame(coords = df.pts[, 2:3], data = data.frame(df.pts[1]),
	proj4string = crs(env.selected))
spp.spdf
windows(w=6, h=6)
plot(projections.all.mean, col = matlab.like(100), main = "Ensemble - Current Climate", las = 1)
plot(spp.spdf, pch = 21, cex = 0.95, col = "gray20", bg = "green", add = T)
plot(spp.spdf.e, pch = 21, cex = 0.95, col = "gray20", bg = "red", add = T)


# You can also add shapefiles/rasters of protected areas, biogeographic regions, vegetation coverage etc.


## Convert final consensus map to binary
## *****************

projections.mean.binary <- BinaryTransformation(projections.all.mean, 0.5) # Or replace 0.5 by another value.
class(projections.mean.binary)
str(values(projections.mean.binary))
summary(values(projections.mean.binary))
windows(w=6, h=6)
plot(projections.mean.binary, col = matlab.like(100), main = "Ensemble - Current Climate (Binary)", las = 1)
