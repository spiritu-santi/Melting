##########################################
#### Covert polygons to raster points ####
##########################################

# 1- Libraries required ####
  library(maptools)
  library(raster)
  library(rgeos)
  library(sp)
  library(tidyverse)
  
# 2- Dissolve to identify all polygons with same ID ####
  GroupofInterest<-readOGR("IUCN/GroupofInterest.shp")
  IDs <- GroupofInterest$id_no
  GroupofInterest_2<-rgeos::gBuffer(GroupofInterest, byid=TRUE, width=0)
  GroupofInterest_3 <- spTransform(GroupofInterest_2,crs.wgs84)
    GroupofInterest_new<-unionSpatialPolygons(GroupofInterest_3, IDs)
  writeOGR(GroupofInterest_new, ".", "GroupofInterest_dissolve", driver="ESRI Shapefile")
  
  # We obtained our endemic species #### 
  # We obtained the endemic species with the "select by location" tool with
  # the spatial selection method for targe layer feature was "are completely within 
  # the source layer feature".
  
# 3- Once we have our endemic species polygon we covert them into coordinates ####
  unique <- unique(endemic@data$binomial,endemic@data$binomial[!is.na(endemic@data$binomial)])
  endemic_limpia <- subset(endemic, !is.na(endemic@data$binomial))
  unique <- unique(endemic_limpia@data$binomial)
  
  bind<-NULL
  total<-NULL
  
  for (i in 1:length(endemic)){
    print(i)
    poly <- endemic_limpia[which(data$binomial == unique[i]), ] 
    r <- raster("MegadiverseCountries.tif") 
    extent(r) <- extent(poly)  
    r <- raster(extent(r), res=0.083333333) # Depending on the distance between points 
    r <- rasterize(poly, r, field=1)  
    xy<-rasterToPoints(r)  
    outname<-poly$binomial  
    bind<-as_tibble(xy) %>%  
      mutate(binomial=as.character(outname)[1])  
    total<-rbind(total, bind) 
    write.csv(xy, paste0("RasterPoints/", as.character(outname)[1], ".csv"),row.names = FALSE)# Create a folder named "RasterPoints"
  }
  
# 4- Adding a third column with the species name ####
  setwd("RasterPoints") # Create a raster points directory
  species<- list.files(pattern="*.csv")
  nombres<-gsub(".csv", "", species) 
  
  speciesNew<- lapply(species, read_csv) 
  for (i in 1:length(species)){          # Rename the third column as the species name
    names(speciesNew)[i]<-names[i]
    names(speciesNew[[i]])[3]<-names[i]
  }
  
# 5- Saving the new files ####
  for(i in 1:length(species)){
    write.csv(speciesNew[[i]], paste0("records/", species[i]), row.names = FALSE)
  }
  
# 6- Filtering species with up to 25 records ####
  speciesNew25<- speciesNew25[lapply(speciesNew, nrow)>=25]
  
  for(i in 1:length(speciesNew25)){
    write.csv(speciesNew25[[i]], 
              paste0("RasterPoints25/", 
                     names(speciesNew25)[i], ".csv"), 
              row.names = FALSE)
  }
  
  