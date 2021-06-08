
library("sp")
library("raster")
library("dplyr")
library("rgdal")
library("fuzzySim")
library("rgeos")
library("magrittr")
library("tools")
library("readr")
library("dismo")
library("biomod2")

start_time<-Sys.time()

setwd("C:/Users/CAROLINA URETA/Desktop/prueba")

####DataFormating ####
shapePath <- 'C:/Users/CAROLINA URETA/Desktop/prueba'
shapeLayer <- "wwf_terr_ecos"
regionalizacion <- rgdal::readOGR(shapePath, shapeLayer)

#
bioclima<-stack(list.files(path="C:/Users/CAROLINA URETA/Desktop/prueba/Presente", pattern = "*.tif$", full.names=TRUE)) #presente

#--------------------------------- BCC_SSP245
covarDataFolder_GreenControl_T1<-stack(list.files(path = "C:/Users/CAROLINA URETA/Desktop/prueba/T1", pattern="*.tif$", full.names=TRUE))
covarDataFolder_GreenControl_T2<-stack(list.files(path = "/LUSTRE/users/see/acuervo/GreenControl/T2", pattern="*.tif$", full.names=TRUE))
covarDataFolder_GreenControl_T3<-stack(list.files(path = "/LUSTRE/users/see/acuervo/GreenControl/T2", pattern="*.tif$", full.names=TRUE))

#--------------------------------- Can_SSP245
covarDataFolder_Green0.5_T1<-stack(list.files(path = "/LUSTRE/users/see/acuervo/Green0.5/T1", pattern="*.tif$", full.names=TRUE))
covarDataFolder_Green0.5_T2<-stack(list.files(path = "/LUSTRE/users/see/acuervo/Green0.5/T2", pattern="*.tif$", full.names=TRUE))
covarDataFolder_Green0.5_T3<-stack(list.files(path = "/LUSTRE/users/see/acuervo/Green0.5/T3", pattern="*.tif$", full.names=TRUE))

#--------------------------------- BCC_SSP585
covarDataFolder_GreenControl_T1<-stack(list.files(path = "/LUSTRE/users/see/acuervo/GreenControl/T1", pattern="*.tif$", full.names=TRUE))
covarDataFolder_GreenControl_T2<-stack(list.files(path = "/LUSTRE/users/see/acuervo/GreenControl/T2", pattern="*.tif$", full.names=TRUE))
covarDataFolder_GreenControl_T3<-stack(list.files(path = "/LUSTRE/users/see/acuervo/GreenControl/T2", pattern="*.tif$", full.names=TRUE))

#--------------------------------- Can_SSP585
covarDataFolder_Green0.5_T1<-stack(list.files(path = "/LUSTRE/users/see/acuervo/Green0.5/T1", pattern="*.tif$", full.names=TRUE))
covarDataFolder_Green0.5_T2<-stack(list.files(path = "/LUSTRE/users/see/acuervo/Green0.5/T2", pattern="*.tif$", full.names=TRUE))
covarDataFolder_Green0.5_T3<-stack(list.files(path = "/LUSTRE/users/see/acuervo/Green0.5/T3", pattern="*.tif$", full.names=TRUE))

#
args <- list.files("C:/Users/CAROLINA URETA/Desktop/prueba", pattern = "*.csv$",full.names = TRUE)


inputDataFile <- args[2]
outputFolder1 <- inputDataFile %>%
  basename %>%
  file_path_sans_ext
outputFolder1
outputFolder<- gsub(" ", ".", outputFolder1)
outputFolder

outputFolder2<- gsub("_", " ", outputFolder1)
outputFolder2

if (!dir.exists(outputFolder)) {
  dir.create(outputFolder, recursive = TRUE)
}

# Extract envorimental varibales with species occurrences
crs.wgs84 <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
occsData <- readr::read_csv(inputDataFile)

if (dim(occsData)[1] > 250) {occsData <- occsData[sample(1:dim(occsData)[1], 250, replace=FALSE),] }
names(occsData)[names(occsData)== outputFolder2] <- outputFolder

sp::coordinates(occsData) <- c("x", "y")
sp::proj4string(occsData) <- crs.wgs84

covarData <- raster::extract(bioclima, occsData)
covarData <- cbind(occsData, covarData)

completeDataCases <- covarData@data %>%
  dplyr::select(.dots=names(bioclima)) %>%
  complete.cases
covarData <- covarData[completeDataCases, ]

####Variables selection####
speciesCol <- match(outputFolder, names(occsData))
varCols <- ncol(occsData) + 1

correlacion <- corSelect(
  data = covarData@data,
  sp.cols = speciesCol,
  var.cols = varCols:ncol(covarData),
  cor.thresh = 0.8,
  use = "pairwise.complete.obs"
)

select_var <- correlacion$selected.vars
write(select_var, file = file.path(outputFolder, "selected_variables.txt"))

#Raster covariables selected for model calibration
enviromentalVariables <- bioclima[[select_var]]

# Raster covariables selected for model calibration
selectedVariables <- enviromentalVariables[[select_var]]

# Raster covariables selected for model transfer
#seleccion variables CC

#--------------------------------- GreenControl
env_GreenControl_T1t<-covarDataFolder_GreenControl_T1[[select_var]]
env_GreenControl_T2t<-covarDataFolder_GreenControl_T2[[select_var]]
env_GreenControl_T3t<-covarDataFolder_GreenControl_T3[[select_var]]
#--------------------------------- Green0.5
env_Green0.5_T1t<-covarDataFolder_Green0.5_T1[[select_var]]
env_Green0.5_T2t<-covarDataFolder_Green0.5_T2[[select_var]]
env_Green0.5_T3t<-covarDataFolder_Green0.5_T3[[select_var]]
#--------------------------------- Green1
env_Green1_T1t<-covarDataFolder_Green1_T1[[select_var]]
env_Green1_T2t<-covarDataFolder_Green1_T2[[select_var]]
env_Green1_T3t<-covarDataFolder_Green1_T3[[select_var]]
#--------------------------------- Green1.5
env_Green1.5_T1t<-covarDataFolder_Green1.5_T1[[select_var]]
env_Green1.5_T2t<-covarDataFolder_Green1.5_T2[[select_var]]
env_Green1.5_T3t<-covarDataFolder_Green1.5_T3[[select_var]]
#--------------------------------- Green1.5
env_Green3_T1t<-covarDataFolder_Green3_T1[[select_var]]
env_Green3_T2t<-covarDataFolder_Green3_T2[[select_var]]
env_Green3_T3t<-covarDataFolder_Green3_T3[[select_var]]

# M####
# Intersects the occurrence data with polygons
ecoregionsOfInterest <- sp::over(occsData, regionalizacion) %>%
  filter(!is.na(ECO_ID))

idsEcoRegions <- unique(ecoregionsOfInterest$ECO_ID)
polygonsOfInterest <- regionalizacion[regionalizacion$ECO_ID %in% idsEcoRegions, ]
pts_b <- gBuffer(occsData, width=4)
pts_b <- as(pts_b, 'SpatialPolygonsDataFrame')

#Poligono area de calibracion
polygonsOfInterest<-gIntersection(pts_b, polygonsOfInterest, drop_lower_td = T)

#Poligono area  transferencia
polyTransferencia<-gBuffer(polygonsOfInterest, width=2)
polyTransferencia <- as(polyTransferencia, 'SpatialPolygonsDataFrame')

# Variables ambientales
# Cortar bioclimas con mascara M
selectedVariablesCrop <- raster::crop(selectedVariables, polygonsOfInterest)
myExpl <- raster::mask(selectedVariablesCrop,polygonsOfInterest) #Species variables delimited by M
myExpl<-stack(myExpl)

# Mask future raster with ecoregions of interest
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

##Green0.5---------------------------
#Green0.5_T1
env_Green0.5_T1a <- raster::crop(env_Green0.5_T1t, polyTransferencia)
env_Green0.5_T1 <- raster::mask(env_Green0.5_T1a ,  polyTransferencia)
env_Green0.5_T1<-stack(env_Green0.5_T1)
rm(env_Green0.5_T1t, env_Green0.5_T1a)
#Green0.5_T2
env_Green0.5_T2a <- raster::crop(env_Green0.5_T2t, polyTransferencia)
env_Green0.5_T2 <- raster::mask(env_Green0.5_T2a ,  polyTransferencia)
env_Green0.5_T2<-stack(env_Green0.5_T2)
rm(env_Green0.5_T2t, env_Green0.5_T2a)
#Green0.5_T3
env_Green0.5_T3a <- raster::crop(env_Green0.5_T3t, polyTransferencia)
env_Green0.5_T3 <- raster::mask(env_Green0.5_T3a ,  polyTransferencia)
env_Green0.5_T3<-stack(env_Green0.5_T3)
rm(env_Green0.5_T3t, env_Green0.5_T3a)

##Green1-------------------------------
#Green1_T1
env_Green1_T1a <- raster::crop(env_Green1_T1t, polyTransferencia)
env_Green1_T1 <- raster::mask(env_Green1_T1a ,  polyTransferencia)
env_Green1_T1<-stack(env_Green1_T1)
rm(env_Green1_T1t, env_Green1_T1a)
#Green1_T2
env_Green1_T2a <- raster::crop(env_Green1_T2t, polyTransferencia)
env_Green1_T2 <- raster::mask(env_Green1_T2a ,  polyTransferencia)
env_Green1_T2<-stack(env_Green1_T2)
rm(env_Green1_T2t, env_Green1_T2a)
#Green1_T3
env_Green1_T3a <- raster::crop(env_Green1_T3t, polyTransferencia)
env_Green1_T3 <- raster::mask(env_Green1_T3a ,  polyTransferencia)
env_Green1_T3<-stack(env_Green1_T3)
rm(env_Green1_T3t, env_Green1_T3a)

##Green1.5-------------------------------
#Green1.5_T1
env_Green1.5_T1a <- raster::crop(env_Green1.5_T1t, polyTransferencia)
env_Green1.5_T1 <- raster::mask(env_Green1.5_T1a ,  polyTransferencia)
env_Green1.5_T1<-stack(env_Green1.5_T1)
rm(env_Green1.5_T1t, env_Green1.5_T1a)
#Green1.5_T2
env_Green1.5_T2a <- raster::crop(env_Green1.5_T2t, polyTransferencia)
env_Green1.5_T2 <- raster::mask(env_Green1.5_T2a ,  polyTransferencia)
env_Green1.5_T2<-stack(env_Green1.5_T2)
rm(env_Green1.5_T2t, env_Green1.5_T2a)
#Green1.5_T3
env_Green1.5_T3a <- raster::crop(env_Green1.5_T3t, polyTransferencia)
env_Green1.5_T3 <- raster::mask(env_Green1.5_T3a ,  polyTransferencia)
env_Green1.5_T3<-stack(env_Green1.5_T3)
rm(env_Green1.5_T3t, env_Green1.5_T3a)

##Green3-------------------------------
#Green3_T1
env_Green3_T1a <- raster::crop(env_Green3_T1t, polyTransferencia)
env_Green3_T1 <- raster::mask(env_Green3_T1a ,  polyTransferencia)
env_Green3_T1<-stack(env_Green3_T1)
rm(env_Green3_T1t, env_Green3_T1a)
#Green3_T2
env_Green3_T2a <- raster::crop(env_Green3_T2t, polyTransferencia)
env_Green3_T2 <- raster::mask(env_Green3_T2a ,  polyTransferencia)
env_Green3_T2<-stack(env_Green3_T2)
rm(env_Green3_T2t, env_Green3_T2a)
#Green3_T3
env_Green3_T3a <- raster::crop(env_Green3_T3t, polyTransferencia)
env_Green3_T3 <- raster::mask(env_Green3_T3a ,  polyTransferencia)
env_Green3_T3<-stack(env_Green3_T3)
rm(env_Green3_T3t, env_Green3_T3a)

#DataFormating####
presencias<-data.frame(occsData)
names(presencias)[names(presencias)=="outputFolder"] <- outputFolder
presencias<-dplyr::select(presencias, c("x", "y",outputFolder))
names(presencias)

## BIOMOD####
DataSpecies<-presencias

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

myBiomodData
#plot(myBiomodData)

### Parametrizacion
myBiomodOption <-BIOMOD_ModelingOptions()

### Modelling
myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('GLM', 'GBM', 'GAM', 'CTA', 'ANN','RF', "MAXENT.Phillips"),
  models.options = myBiomodOption,
  NbRunEval=10,
  DataSplit=70,
  models.eval.meth = c('KAPPA','TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = FALSE,
  do.full.models = FALSE,
  modeling.id = paste(outputFolder))

#myBiomodModelOut
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
write.csv(myBiomodModelEval, file = file.path(outputFolder, "myBiomodModelEval.csv"),
          row.names = FALSE)

### Hacer predicciones sobre el raster
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl,
  proj.name = outputFolder,
  selected.models = 'all',
  binary.meth = "TSS",
  compress = 'xz',
  build.clamping.mask = TRUE,
  output.format = '.grd')

#plot(myBiomodProj)

### ensemble_modeling#####
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

# print summary
#myBiomodEM

# get evaluation scores
myBiomodEMEval<-get_evaluations(myBiomodEM)
write.csv(myBiomodEMEval, file = file.path(outputFolder, "myBiomodEMEval.csv"),
          row.names = FALSE)

####Current####
# Creating ensembles projections
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


## TRANSFERENCIA AL FUTURO ####

###GreenControl_T1####
### Hacer predicciones sobre el raster futuros
myBiomodEM_GreenControl_T1<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_GreenControl_T1,
                                                        selected.models = 'all',
                                                        proj.name = "GreenControl_T1",
                                                        binary.meth = "TSS")

myBiomodEM_GreenControl_T1
#plot(myBiomodEM_GreenControl_T1)

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

###GreenControl_T2####
### Hacer predicciones sobre el raster futuros
myBiomodEM_GreenControl_T2<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_GreenControl_T2,
                                                        selected.models = 'all',
                                                        proj.name = "GreenControl_T2",
                                                        binary.meth = "TSS"
)

#myBiomodEM_GreenControl_T2
#plot(myBiomodEM_GreenControl_T2)

FutureProj <- stack(file.path(outputFolder, paste("proj_GreenControl_T2/proj_GreenControl_T2_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_GreenControl_T2/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_GreenControl_T2/proj_GreenControl_T2_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_GreenControl_T2/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c,myBiomodEM_GreenControl_T2)

###GreenControl_T3####
### Hacer predicciones sobre el raster futuros
myBiomodEM_GreenControl_T3<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_GreenControl_T3,
                                                        selected.models = 'all',
                                                        proj.name = "GreenControl_T3",
                                                        binary.meth = "TSS"
)

#myBiomodEM_GreenControl_T3
#plot(myBiomodEM_GreenControl_T3)

#myBiomodEFPred
FutureProj <- stack(file.path(outputFolder, paste("proj_GreenControl_T3/proj_GreenControl_T3_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_GreenControl_T3/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_GreenControl_T3/proj_GreenControl_T3_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_GreenControl_T3/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj,FutureProj_c, myBiomodEM_GreenControl_T3)

###Green0.5_T1####
### Hacer predicciones sobre el raster futuros
myBiomodEM_Green0.5_T1<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_Green0.5_T1,
                                                    selected.models = 'all',
                                                    proj.name = "Green0.5_T1",
                                                    binary.meth = "TSS")

#myBiomodEM_Green0.5_T1
#plot(myBiomodEM_Green0.5_T1)

FutureProj <- stack(file.path(outputFolder, paste("proj_Green0.5_T1/proj_Green0.5_T1_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_Green0.5_T1/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_Green0.5_T1/proj_Green0.5_T1_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_Green0.5_T1/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c, myBiomodEM_Green0.5_T1)

###Green0.5_T2####
### Hacer predicciones sobre el raster futuros
myBiomodEM_Green0.5_T2<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_Green0.5_T2,
                                                    selected.models = 'all',
                                                    proj.name = "Green0.5_T2",
                                                    binary.meth = "TSS"
)

#myBiomodEM_Green0.5_T2
#plot(myBiomodEM_Green0.5_T2)

FutureProj <- stack(file.path(outputFolder, paste("proj_Green0.5_T2/proj_Green0.5_T2_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_Green0.5_T2/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_Green0.5_T2/proj_Green0.5_T2_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_Green0.5_T2/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)


rm(FutureProj,FutureProj_c, myBiomodEM_Green0.5_T2)

###Green0.5_T3####
### Hacer predicciones sobre el raster futuros
myBiomodEM_Green0.5_T3<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_Green0.5_T3,
                                                    selected.models = 'all',
                                                    proj.name = "Green0.5_T3",
                                                    binary.meth = "TSS"
)

#myBiomodEM_Green0.5_T3
#plot(myBiomodEM_Green0.5_T3)

#myBiomodEFPred
FutureProj <- stack(file.path(outputFolder, paste("proj_Green0.5_T3/proj_Green0.5_T3_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_Green0.5_T3/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_Green0.5_T3/proj_Green0.5_T3_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_Green0.5_T3/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)


rm(FutureProj,FutureProj_c, myBiomodEM_Green0.5_T3)

###Green1_T1####
### Hacer predicciones sobre el raster futuros
myBiomodEM_Green1_T1<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                  new.env = env_Green1_T1,
                                                  selected.models = 'all',
                                                  proj.name = "Green1_T1",
                                                  binary.meth = "TSS")

#myBiomodEM_Green1_T1
#plot(myBiomodEM_Green1_T1)

FutureProj <- stack(file.path(outputFolder, paste("proj_Green1_T1/proj_Green1_T1_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_Green1_T1/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)


FutureProj_c <- stack(file.path(outputFolder, paste("proj_Green1_T1/proj_Green1_T1_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_Green1_T1/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c, myBiomodEM_Green1_T1)

###Green1_T2####
### Hacer predicciones sobre el raster futuros
myBiomodEM_Green1_T2<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                  new.env = env_Green1_T2,
                                                  selected.models = 'all',
                                                  proj.name = "Green1_T2",
                                                  binary.meth = "TSS"
)

#myBiomodEM_Green1_T2
#plot(myBiomodEM_Green1_T2)

FutureProj <- stack(file.path(outputFolder, paste("proj_Green1_T2/proj_Green1_T2_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_Green1_T2/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_Green1_T2/proj_Green1_T2_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_Green1_T2/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c, myBiomodEM_Green1_T2)

###Green1_T3####
### Hacer predicciones sobre el raster futuros
myBiomodEM_Green1_T3<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                  new.env = env_Green1_T3,
                                                  selected.models = 'all',
                                                  proj.name = "Green1_T3",
                                                  binary.meth = "TSS"
)

#myBiomodEM_Green1_T3
#plot(myBiomodEM_Green1_T3)

#myBiomodEFPred
FutureProj <- stack(file.path(outputFolder, paste("proj_Green1_T3/proj_Green1_T3_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_Green1_T3/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_Green1_T3/proj_Green1_T3_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_Green1_T3/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)


rm(FutureProj,FutureProj_c,  myBiomodEM_Green1_T3)

###Green1.5_T1####
### Hacer predicciones sobre el raster futuros
myBiomodEM_Green1.5_T1<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_Green1.5_T1,
                                                    selected.models = 'all',
                                                    proj.name = "Green1.5_T1",
                                                    binary.meth = "TSS")

myBiomodEM_Green1.5_T1
#plot(myBiomodEM_Green1.5_T1)

FutureProj <- stack(file.path(outputFolder, paste("proj_Green1.5_T1/proj_Green1.5_T1_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_Green1.5_T1/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_Green1.5_T1/proj_Green1.5_T1_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_Green1.5_T1/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj,FutureProj_c, myBiomodEM_Green1.5_T1)

###Green1.5_T2####
### Hacer predicciones sobre el raster futuros
myBiomodEM_Green1.5_T2<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_Green1.5_T2,
                                                    selected.models = 'all',
                                                    proj.name = "Green1.5_T2",
                                                    binary.meth = "TSS"
)

#myBiomodEM_Green1.5_T2
#plot(myBiomodEM_Green1.5_T2)

FutureProj <- stack(file.path(outputFolder, paste("proj_Green1.5_T2/proj_Green1.5_T2_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_Green1.5_T2/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_Green1.5_T2/proj_Green1.5_T2_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_Green1.5_T2/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj,FutureProj_c, myBiomodEM_Green1.5_T2)

###Green1.5_T3####
### Hacer predicciones sobre el raster futuros
myBiomodEM_Green1.5_T3<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_Green1.5_T3,
                                                    selected.models = 'all',
                                                    proj.name = "Green1.5_T3",
                                                    binary.meth = "TSS"
)

#myBiomodEM_Green1.5_T3
#plot(myBiomodEM_Green1.5_T3)

#myBiomodEFPred
FutureProj <- stack(file.path(outputFolder, paste("proj_Green1.5_T3/proj_Green1.5_T3_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_Green1.5_T3/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_Green1.5_T3/proj_Green1.5_T3_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_Green1.5_T3/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)
rm(FutureProj,FutureProj_c, myBiomodEM_Green1.5_T3)

###Green3_T1####
### Hacer predicciones sobre el raster futuros
myBiomodEM_Green3_T1<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                  new.env = env_Green3_T1,
                                                  selected.models = 'all',
                                                  proj.name = "Green3_T1",
                                                  binary.meth = "TSS")

myBiomodEM_Green3_T1
#plot(myBiomodEM_Green3_T1)

FutureProj <- stack(file.path(outputFolder, paste("proj_Green3_T1/proj_Green3_T1_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_Green3_T1/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_Green3_T1/proj_Green3_T1_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_Green3_T1/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj,FutureProj_c, myBiomodEM_Green3_T1)

###Green3_T2####
### Hacer predicciones sobre el raster futuros
myBiomodEM_Green3_T2<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                  new.env = env_Green3_T2,
                                                  selected.models = 'all',
                                                  proj.name = "Green3_T2",
                                                  binary.meth = "TSS"
)

#myBiomodEM_Green3_T2
#plot(myBiomodEM_Green3_T2)

FutureProj <- stack(file.path(outputFolder, paste("proj_Green3_T2/proj_Green3_T2_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_Green3_T2/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_Green3_T2/proj_Green3_T2_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_Green3_T2/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj,FutureProj_c, myBiomodEM_Green3_T2)

###Green3_T3####
### Hacer predicciones sobre el raster futuros
myBiomodEM_Green3_T3<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                  new.env = env_Green3_T3,
                                                  selected.models = 'all',
                                                  proj.name = "Green3_T3",
                                                  binary.meth = "TSS"
)

#myBiomodEM_Green3_T3
#plot(myBiomodEM_Green3_T3)
#myBiomodEFPred
FutureProj <- stack(file.path(outputFolder, paste("proj_Green3_T3/proj_Green3_T3_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_Green3_T3/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_Green3_T3/proj_Green3_T3_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_Green3_T3/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj,FutureProj_c, myBiomodEM_Green3_T3)


  
unlink(file.path(outputFolder, paste("models/")), recursive = T)
unlink(list.files(outputFolder, pattern = "*.grd$", full.names=TRUE, recursive = T))
unlink(list.files(outputFolder, pattern = "*.gri$", full.names=TRUE, recursive = T))
#unlink(list.files(outputFolder, pattern = "*.out$", full.names=TRUE, recursive = T))


end_time<-Sys.time()

x<-data.frame(end_time-start_time)
write.csv(x, file = file.path(outputFolder, "time.csv"),
          row.names = FALSE)
