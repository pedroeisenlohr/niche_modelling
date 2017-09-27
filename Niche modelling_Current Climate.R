####################################################################################
####################################################################################
#######          NICHE MODELLING WITH BIOMOD2 USING       #######
#######    40 ENVIRONMENTAL VARIABLES (10-km RESOLUTION)  ####### 
#######                 SUMMARIZED IN PCA AXES            #######
####################################################################################
####################################################################################



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
## Relative Umidity: Climond (https://www.climond.org/RawClimateData.aspx).
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
install.packages("sdm", dep=T)
install.packages("SDMTools", dep=T)
install.packages("sqldf", dep=T)
install.packages("testthat", dep=T)
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
library(sdm)
library(SDMTools)
library(sqldf)
library(testthat)
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

bioclim <- list.files("./wc5 2.0", full.names=TRUE)
bio <-stack(bioclim)
plot(bio[[1]])

#Extracting values from species occurrence records:
presvals <- extract(bio, spp)
dim(presvals)
edit(presvals)

## Crop WorldClim layers
#  *********************
neotrop <- readOGR("./ShapeNeo/neotrpic_mex_contorno.shp")
bio.wc2 <- mask(crop(bio,neotrop),neotrop)
bio.wc2
res(bio.wc2)
plot(bio.wc2[[1]])
names(bio.wc2)


####################################################################
### Compiling other rasters to stack ###
####################################################################

#Solar Radiation:
solar.radiation <- list.files("./Solar Radiation", pattern=".tif", full.names=TRUE)
solar.radiation <- stack(solar.radiation)
solar.radiation.mean <- mean(solar.radiation)
solar.radiation.max <- max(solar.radiation)
solar.radiation.min <- min(solar.radiation)
solar.radiation.mean <- mask(crop(solar.radiation.mean, neotrop),neotrop)
solar.radiation.max <- mask(crop(solar.radiation.max, neotrop),neotrop)
solar.radiation.min <- mask(crop(solar.radiation.min, neotrop),neotrop)
res(solar.radiation.mean)
plot(solar.radiation.mean)
res(solar.radiation.max)
plot(solar.radiation.max)
res(solar.radiation.min)
plot(solar.radiation.min)


#Water Vapor Pressure:
water.vapor.pressure <- list.files("./Water Vapor Pressure", pattern=".tif", full.names=TRUE)
water.vapor.pressure <-stack(water.vapor.pressure)
water.vapor.pressure.mean <-mean(water.vapor.pressure)
water.vapor.pressure.max <-max(water.vapor.pressure)
water.vapor.pressure.min <-min(water.vapor.pressure)
water.vapor.pressure.mean <- mask(crop(water.vapor.pressure.mean, neotrop),neotrop)
water.vapor.pressure.max <- mask(crop(water.vapor.pressure.max, neotrop),neotrop)
water.vapor.pressure.min <- mask(crop(water.vapor.pressure.min, neotrop),neotrop)
res(water.vapor.pressure.mean)
plot(water.vapor.pressure.mean)
res(water.vapor.pressure.max)
plot(water.vapor.pressure.max)
res(water.vapor.pressure.min)
plot(water.vapor.pressure.min)


#Wind Speed:
wind.speed <- list.files("./Wind Speed", pattern=".tif", full.names=TRUE)
wind.speed <- stack(wind.speed)
wind.speed.mean <-mean(wind.speed)
wind.speed.max <-max(wind.speed)
wind.speed.min <-min(wind.speed)
wind.speed.mean <-mask(crop(wind.speed.mean, neotrop),neotrop)
wind.speed.max <-mask(crop(wind.speed.max, neotrop),neotrop)
wind.speed.min <-mask(crop(wind.speed.min, neotrop),neotrop)
res(wind.speed.mean)
plot(wind.speed.mean)
res(wind.speed.max)
plot(wind.speed.max)
res(wind.speed.min)
plot(wind.speed.min)


#Potential Evapotranspiration:
PET.1km <- raster("./Global Aridity and PET database/Global PET - Annual/PET_he_annual/pet_he_yr/w001001.adf")
PET.1km <- mask(crop(PET.1km,neotrop),neotrop)
PET.10km <- resample(PET.1km,bio.wc2)
res(PET.10km)
plot(PET.10km)


#Aridity Index:
Aridity.1km <- raster("./Global Aridity and PET database/Global Aridity - Annual/AI_annual/ai_yr/w001001.adf")
Aridity.1km <- mask(crop(Aridity.1km,neotrop),neotrop)
Aridity.10km <- resample(Aridity.1km,bio.wc2)
res(Aridity.10km)
plot(Aridity.10km)


#Actual Evapotranspiration:
AET.1km <- raster("./Global Soil Water Balance and AET/Mean Annual AET/AET_YR/aet_yr/w001001.adf")
AET.1km <- mask(crop(AET.1km,neotrop),neotrop)
AET.10km <- resample(AET.1km,bio.wc2)
res(AET.10km)
plot(AET.10km)


#Soil Water Stress:
SWS.jan <-raster("./Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_1/w001001.adf")
SWS.feb <-raster("./Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_2/w001001.adf")
SWS.mar <-raster("./Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_3/w001001.adf")
SWS.apr <-raster("./Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_4/w001001.adf")
SWS.may <-raster("./Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_5/w001001.adf")
SWS.jun <-raster("./Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_6/w001001.adf")
SWS.jul <-raster("./Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_7/w001001.adf")
SWS.aug <-raster("./Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_8/w001001.adf")
SWS.sep <-raster("./Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_9/w001001.adf")
SWS.oct <-raster("./Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_10/w001001.adf")
SWS.nov <-raster("./Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_11/w001001.adf")
SWS.dec <-raster("./Global Soil Water Balance and AET/Monthly Soil Water Stress/swc_fr/swc_fr_12/w001001.adf")
SWS.stack <-stack(SWS.jan,SWS.feb,SWS.mar,SWS.apr,SWS.may,SWS.jun,SWS.jul,
				SWS.aug,SWS.sep,SWS.oct,SWS.nov,SWS.dec)

SWS.mean.1km <-mean(SWS.stack)
SWS.mean.1km <-mask(crop(SWS.mean.1km,neotrop),neotrop)
SWS.mean.10km <-resample(SWS.mean.1km, bio.wc2)
res(SWS.mean.10km)
plot(SWS.mean.10km)

SWS.min.1km <-min(SWS.stack)
SWS.min.1km <-mask(crop(SWS.min.1km,neotrop),neotrop)
SWS.min.10km <-resample(SWS.min.1km, bio.wc2)
res(SWS.min.10km)
plot(SWS.min.10km)

SWS.max.1km <-max(SWS.stack)
SWS.max.1km <-mask(crop(SWS.max.1km,neotrop),neotrop)
SWS.max.10km <-resample(SWS.max.1km, bio.wc2)
res(SWS.max.10km)
plot(SWS.max.10km)


#Relative Umidity at 3pm:
Umidity.3pm.jan <-raster("./Relative Umidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm01/w001001.adf")
Umidity.3pm.feb <-raster("./Relative Umidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm02/w001001.adf")
Umidity.3pm.mar <-raster("./Relative Umidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm03/w001001.adf")
Umidity.3pm.apr <-raster("./Relative Umidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm04/w001001.adf")
Umidity.3pm.may <-raster("./Relative Umidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm05/w001001.adf")
Umidity.3pm.jun <-raster("./Relative Umidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm06/w001001.adf")
Umidity.3pm.jul <-raster("./Relative Umidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm07/w001001.adf")
Umidity.3pm.aug <-raster("./Relative Umidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm08/w001001.adf")
Umidity.3pm.sep <-raster("./Relative Umidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm09/w001001.adf")
Umidity.3pm.oct <-raster("./Relative Umidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm10/w001001.adf")
Umidity.3pm.nov <-raster("./Relative Umidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm11/w001001.adf")
Umidity.3pm.dec <-raster("./Relative Umidity at 3 pm/CM10_1975H_Raw_ESRI_RHpm_V1.2/CM10_1975H_Raw_ESRI_RHpm_V1.2/rhpm12/w001001.adf")
Umidity.3pm.stack <-stack(Umidity.3pm.jan, Umidity.3pm.feb, Umidity.3pm.mar, Umidity.3pm.apr, Umidity.3pm.may, Umidity.3pm.jun, Umidity.3pm.jul, 
				Umidity.3pm.aug, Umidity.3pm.sep, Umidity.3pm.oct, Umidity.3pm.nov, Umidity.3pm.dec)

Umidity.3pm.mean.20km <-mean(Umidity.3pm.stack)
Umidity.3pm.mean.20km <-mask(crop(Umidity.3pm.mean.20km,neotrop),neotrop)
Umidity.3pm.mean.10km <-resample(Umidity.3pm.mean.20km, bio.wc2)
res(Umidity.3pm.mean.10km)
plot(Umidity.3pm.mean.10km)

Umidity.3pm.max.20km <-max(Umidity.3pm.stack)
Umidity.3pm.max.20km <-mask(crop(Umidity.3pm.max.20km,neotrop),neotrop)
Umidity.3pm.max.10km <-resample(Umidity.3pm.max.20km, bio.wc2)
res(Umidity.3pm.max.10km)
plot(Umidity.3pm.max.10km)

Umidity.3pm.min.20km <-min(Umidity.3pm.stack)
Umidity.3pm.min.20km <-mask(crop(Umidity.3pm.min.20km,neotrop),neotrop)
Umidity.3pm.min.10km <-resample(Umidity.3pm.min.20km, bio.wc2)
res(Umidity.3pm.min.10km)
plot(Umidity.3pm.min.10km)


#Relative Umidity at 9am:
Umidity.9am.jan <-raster("./Relative Umidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham01/w001001.adf")
Umidity.9am.feb <-raster("./Relative Umidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham02/w001001.adf")
Umidity.9am.mar <-raster("./Relative Umidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham03/w001001.adf")
Umidity.9am.apr <-raster("./Relative Umidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham04/w001001.adf")
Umidity.9am.may <-raster("./Relative Umidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham05/w001001.adf")
Umidity.9am.jun <-raster("./Relative Umidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham06/w001001.adf")
Umidity.9am.jul <-raster("./Relative Umidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham07/w001001.adf")
Umidity.9am.aug <-raster("./Relative Umidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham08/w001001.adf")
Umidity.9am.sep <-raster("./Relative Umidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham09/w001001.adf")
Umidity.9am.oct <-raster("./Relative Umidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham10/w001001.adf")
Umidity.9am.nov <-raster("./Relative Umidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham11/w001001.adf")
Umidity.9am.dec <-raster("./Relative Umidity at 9 am/CM10_1975H_Raw_ESRI_RHam_V1.2/CM10_1975H_Raw_ESRI_RHam_V1.2/rham12/w001001.adf")

Umidity.9am.mean.20km <-stack(Umidity.9am.jan, Umidity.9am.feb, Umidity.9am.mar, Umidity.9am.apr, Umidity.9am.may, Umidity.9am.jun, Umidity.9am.jul, 
				Umidity.9am.aug, Umidity.9am.sep, Umidity.9am.oct, Umidity.9am.nov, Umidity.9am.dec)
Umidity.9am.mean.20km <-mask(crop(Umidity.9am.mean.20km,neotrop),neotrop)
Umidity.9am.mean.10km <-resample(Umidity.9am.mean.20km, bio.wc2)
res(Umidity.9am.mean.10km)
plot(Umidity.9am.mean.10km)

Umidity.9am.max.20km <-max(Umidity.9am.jan, Umidity.9am.feb, Umidity.9am.mar, Umidity.9am.apr, Umidity.9am.may, Umidity.9am.jun, Umidity.9am.jul, 
				Umidity.9am.aug, Umidity.9am.sep, Umidity.9am.oct, Umidity.9am.nov, Umidity.9am.dec)
Umidity.9am.max.20km <-mask(crop(Umidity.9am.max.20km,neotrop),neotrop)
Umidity.9am.max.10km <-resample(Umidity.9am.max.20km, bio.wc2)
res(Umidity.9am.max.10km)
plot(Umidity.9am.max.10km)

Umidity.9am.min.20km <-min(Umidity.9am.jan, Umidity.9am.feb, Umidity.9am.mar, Umidity.9am.apr, Umidity.9am.may, Umidity.9am.jun, Umidity.9am.jul, 
				Umidity.9am.aug, Umidity.9am.sep, Umidity.9am.oct, Umidity.9am.nov, Umidity.9am.dec)
Umidity.9am.min.20km <-mask(crop(Umidity.9am.min.20km,neotrop),neotrop)
Umidity.9am.min.10km <-resample(Umidity.9am.min.20km, bio.wc2)
res(Umidity.9am.min.10km)
plot(Umidity.9am.min.10km)


#########################################################################
############### Stacking all environmental layers #######################
#########################################################################
bio.crop <- stack(bio.wc2, 
	solar.radiation.mean, solar.radiation.max, solar.radiation.min, 
	water.vapor.pressure.mean, water.vapor.pressure.max, water.vapor.pressure.min, 
	wind.speed.mean, wind.speed.max, wind.speed.min, 
	PET.10km, Aridity.10km, AET.10km, 
	SWS.mean.10km, SWS.min.10km, SWS.max.10km,
	Umidity.3pm.mean.10km, Umidity.3pm.min.10km, Umidity.3pm.max.10km, Umidity.9am.mean.10km, Umidity.9am.max.10km, Umidity.9am.min.10km)
bio.crop


#################################################################
##################### PCA #######################################
#################################################################

env.selected1 <- rasterPCA(bio.crop, nComp=7, spca = TRUE)
# Here I selected the first 7 components because they account for about 95% 
# of the total variance of the 40 predictors for Neotropical Region.

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
## BUILDING SPECIES DISTRIBUTION MODELS - SDMs ##
#################################################

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
sppModelOut.PA.equal <- BIOMOD_Modeling(sppBiomodData.PA.equal, 
	models = c("GBM", "CTA", "RF"), 
	NbRunEval = 10,
	DataSplit = 70, 
	Prevalence = NULL, 
	VarImport = 1000,
	models.eval.meth = c("KAPPA","TSS","ROC","SR","FAR","ACCURACY","BIAS"),
	SaveObj = TRUE,
	rescal.all.models = FALSE,
	do.full.models = FALSE,
	modeling.id = "spp")
sppModelOut.PA.equal


# Parallel processing
cl <- makeCluster(detectCores()) # number of cores in computer
registerDoParallel(cl)
getDoParWorkers()
                      
sppModelOut.PA.10000 <- BIOMOD_Modeling( 
	sppBiomodData.PA.10000, 
	models = c("GLM","GAM","ANN","SRE","FDA","MARS","MAXENT.Phillips"), 
	models.options = myBiomodOption, 
	NbRunEval = 10,
	DataSplit = 70, 
	Prevalence = NULL, 
	VarImport = 1000,
	models.eval.meth = c("KAPPA","TSS","ROC","SR","FAR","ACCURACY","BIAS"),
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


### You should define where is the file 'proj_Cur2_Occurrence.grd'
projections_cont2 <-stack("./Occurrence/proj_Cur2/proj_Cur2_Occurrence.grd") #onde estão os modelos produzidos pelos demais algoritmos
names(projections_cont2)

### You should define where is the file 'proj_Cur2_Occurrence_ROCbin.grd'
projections_cont2 <-stack("./Occurrence/proj_Cur2/proj_Cur2_Occurrence_ROCbin.grd") #onde estão os modelos produzidos pelos demais algoritmos
names(projections_bin2)

names(projections_cont2)=names(projections_bin2)
plot(projections_cont2)


### Apply the steps below only for retained algorithms
### Ensemble mean model for each algorithm:
projections.RF.all <- subset(projections_bin1, grep("RF", names(projections_bin1)))
projections.RF.mean <- mean(projections.RF.all)
windows(w=6, h=6)
plot(projections.RF.mean, col = matlab.like(100), main = "RF", las = 1)

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


# You can also add shapefiles of protected areas, biogeographic regions etc.


## Convert final consensus map to binary
## *****************

projections.mean.binary <- BinaryTransformation(projections.all.mean, 0.5) # Or replace 0.5 by another value.
class(projections.mean.binary)
str(values(projections.mean.binary))
summary(values(projections.mean.binary))
windows(w=6, h=6)
plot(projections.mean.binary, col = matlab.like(100), main = "Ensemble - Current Climate (Binary)", las = 1)


