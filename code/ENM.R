##########################################
#### Ecological Niche Modeling############
##########################################

# 1- Libraries required ####
library("biomod2")
library("dismo")
library("dplyr")
library("fuzzySim")
library("magrittr")
library("raster")
library("rgdal")
library("rgeos")
library("readr")
library("sp")
library("tools")


# 2- Ecorregions that will be used to define M ####
ecorregions<-raster::shapefile("World'sEcorregions.shp")
crs.wgs84 <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
ecorregions <- spTransform(ecorregions,crs.wgs84)

# 3- Climatic variables #### 
#Current climatic conditions
bioclima<-stack(list.files(path="/Present", pattern = "*.tif$", full.names=TRUE)) 
#Future climatic variables 
#--------------------------------- GreenControl
covarDataFolder_GreenControl_T1<-stack(list.files(path = "/T1", pattern="*.tif$", full.names=TRUE)) 
covarDataFolder_GreenControl_T2<-stack(list.files(path = "/T2", pattern="*.tif$", full.names=TRUE))   
covarDataFolder_GreenControl_T3<-stack(list.files(path = "/T3", pattern="*.tif$", full.names=TRUE)) 
#---------------------------------
#Repeat with different climatic scenarios and time horizons

# 4-Species occurrences #### 
args <- list.files("species/", pattern = "*.csv$",full.names = TRUE)

inputDataFile <- args[1]
outputFolder <- inputDataFile %>%
  basename %>%
  file_path_sans_ext
outputFolder

#In case the folder name needs to be changed
outputFolder<- gsub(" ", ".", outputFolder1) 
outputFolder
#outputFolder2<- gsub("_", " ", outputFolder1)
#outputFolder2

if (!dir.exists(outputFolder)) {
  dir.create(outputFolder, recursive = TRUE)
}

# 5-Extract envorimental varibales with species occurrences ####
occsData<-readr::read_csv(inputDataFile)
if (dim(occsData)[1] > 500) {occsData <- occsData[sample(1:dim(occsData)[1], 500, replace=FALSE),] }
names(occsData)[names(occsData)== outputFolder] <- outputFolder

names(occsData)[3]<-outputFolder 
sp::coordinates(occsData) <- c("x", "y")
sp::proj4string(occsData) <- crs.wgs84

covarData <- raster::extract(bioclima, occsData)
covarData <- cbind(occsData, covarData)

completeDataCases <- covarData@data %>%
  dplyr::select(.dots=names(bioclima)) %>%
  complete.cases 
covarData <- covarData[completeDataCases, ]


# 6- Selecting variables by eliminating those with high correlation ####
speciesCol <- raster::match(outputFolder, names(occsData))
varCols <- ncol(occsData) + 1

correlacion <- corSelect(
  data = covarData@data,
  sp.cols = speciesCol,
  var.cols = varCols:ncol(covarData),
  cor.thresh = 0.8,
  use = "pairwise.complete.obs"
)
select_var <- correlacion$selected.vars
write(select_var, file = file.path(outputFolder, "selected_variables.csv"))

enviromentalVariables <- bioclima[[select_var]]
selectedVariables <- enviromentalVariables[[select_var]]

#--------------------------------- GreenControl
env_GreenControl_T1t<-covarDataFolder_GreenControl_T1[[select_var]]
env_GreenControl_T2t<-covarDataFolder_GreenControl_T2[[select_var]]
env_GreenControl_T3t<-covarDataFolder_GreenControl_T3[[select_var]]
#---------------------------------
#Repeat with different climatic scenarios and time horizons

# 7- M selection ####
ecoregionsOfInterest <- sp::over(occsData, ecorregions) %>%
  filter(!is.na(ECO_ID))

idsEcoRegions <- unique(ecoregionsOfInterest$ECO_ID)
polygonsOfInterest <- ecorregions[ecorregions$ECO_ID %in% idsEcoRegions, ]
pts_b <- gBuffer(occsData, width=4)
pts_b <- spTransform(pts_b,crs.wgs84)
pts_b <- as(pts_b, 'SpatialPolygonsDataFrame')

# 7.1 Calibration area
polygonsOfInterest<-gIntersection(pts_b, polygonsOfInterest, drop_lower_td = T)

# 7.2 Transference area
polyTransferencia<-gBuffer(polygonsOfInterest, width=2)
polyTransferencia <- as(polyTransferencia, 'SpatialPolygonsDataFrame')

# 7.3 Mask current climatic variables to calibration area
selectedVariablesCrop <- raster::crop(selectedVariables, polygonsOfInterest)
myExpl <- raster::mask(selectedVariablesCrop,polygonsOfInterest) 
myExpl<-stack(myExpl)

# 8- Mask future climatic variables to transference area ####
#GrenControl-----------------------
#GreenControl_T1
env_GreenControl_T1a <- raster::crop(env_GreenControl_T1t, polyTransferencia)
env_GreenControl_T1 <- raster::mask(env_GreenControl_T1a ,  polyTransferencia)
env_GreenControl_T1<-stack(env_GreenControl_T1)
rm(env_GreenControl_T1t, env_GreenControl_T1a)
#GreenControl_T3
env_GreenControl_T3a <- raster::crop(env_GreenControl_T3t, polyTransferencia)
env_GreenControl_T3 <- raster::mask(env_GreenControl_T3a ,  polyTransferencia)
env_GreenControl_T3<-stack(env_GreenControl_T3)
rm(env_GreenControl_T3t, env_GreenControl_T3a)
#GreenControl_T2
env_GreenControl_T2a <- raster::crop(env_GreenControl_T2t, polyTransferencia)
env_GreenControl_T2 <- raster::mask(env_GreenControl_T2a ,  polyTransferencia)
env_GreenControl_T2<-stack(env_GreenControl_T2)
rm(env_GreenControl_T2t, env_GreenControl_T2a)
#---------------------------------
#Repeat with different climatic scenarios and time horizons

# 9- Data formatting ####
occurrences<-data.frame(occsData)
names(occurrences)[names(occurrences)=="outputFolder"] <- outputFolder
occurrences<-dplyr::select(occurrences, c("x", "y",outputFolder))
names(occurrences)

# 10- Modeling ####
DataSpecies<-occurrences
myRespName <- outputFolder
myResp <- as.numeric(DataSpecies[,myRespName])
myRespCoord = DataSpecies[c("x", "y")]
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespCoord,
                                     resp.name = myRespName,
                                     PA.nb.rep = 1,
                                     PA.nb.absences = 10000,
                                     PA.strategy = 'random')

# 11- Parametrization and evaluation ####
myBiomodOption <-BIOMOD_ModelingOptions()
myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('GLM', 'GBM', 'GAM', 'CTA', 'ANN','RF','MAXENT.Phillips'), #SRE, FDA, MARS
  models.options = myBiomodOption,
  NbRunEval=10,
  DataSplit=70,
  models.eval.meth = c('KAPPA','TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = FALSE,
  do.full.models = FALSE,
  modeling.id = paste(outputFolder))
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
write.csv(myBiomodModelEval, file = file.path(outputFolder, "myBiomodModelEval.csv"),
          row.names = FALSE)

# 11- Projection ####
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl,
  proj.name = outputFolder,
  selected.models = 'all',
  binary.meth = "TSS",
  compress = 'xz',
  build.clamping.mask = TRUE,
  output.format = '.grd')

# 12- Ensemble ####
myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('ROC'),
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  eval.metric.quality.threshold = c(0.7),
  prob.mean.weight = T, 
  VarImport = 1)

myVarImportEM<-data.frame(get_variables_importance(myBiomodEM))
myVarImportEM<-myVarImportEM[5]
write.csv(myVarImportEM, file = file.path(outputFolder, "myVarImportEM.csv"),
          row.names = T)


# 13- Ensemble validation ####
myBiomodEMEval<-get_evaluations(myBiomodEM)
write.csv(myBiomodEMEval, file = file.path(outputFolder, "myBiomodEMEval.csv"),
          row.names = FALSE)

# 14- Projection under current climatic conditions ####
myBiomodEM_proj <-BIOMOD_EnsembleForecasting(EM.output  = myBiomodEM,
                                             projection.output = myBiomodProj,
                                             selected.models = 'all',
                                             proj.name = outputFolder,
                                             binary.meth = "TSS")

outpath<-file.path(outputFolder,paste0("proj_",outputFolder))
currentPred <- stack(file.path(outputFolder,paste("proj_",outputFolder,"/proj_",outputFolder,"_",outputFolder,"_ensemble_TSSbin.grd",sep="")))
writeRaster(currentPred,
            file.path(outpath,"TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

currentPred_c <- stack(file.path(outputFolder,paste("proj_",outputFolder,"/proj_",outputFolder,"_",outputFolder,"_ensemble.grd",sep="")))
writeRaster(currentPred_c,
            file.path(outpath,"c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

ClampMask <- stack(file.path(outputFolder,paste("proj_",outputFolder,"/proj_",outputFolder,"_ClampingMask.grd",sep="")))
writeRaster(ClampMask,file.path(outpath,paste0(outputFolder,"_ClampingMask.tif")),overwrite= TRUE)

# 15- Projections under current climatic conditions ####
#GreenControl_T1#
myBiomodEM_GreenControl_T1<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_GreenControl_T1,
                                                        selected.models = 'all',
                                                        proj.name = "GreenControl_T1",
                                                        binary.meth = "TSS")
myBiomodEM_GreenControl_T1
FutureProj <- stack(file.path(outputFolder, paste("proj_GreenControl_T1/proj_GreenControl_T1_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_GreenControl_T1/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)
FutureProj_c <- stack(file.path(outputFolder, paste("proj_GreenControl_T1/proj_GreenControl_T1_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_GreenControl_T1/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c, myBiomodEM_GreenControl_T1)
#---------------------------------
#Repeat with different climatic scenarios and time horizons

