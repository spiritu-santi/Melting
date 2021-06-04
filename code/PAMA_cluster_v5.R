if(!dir.exists("PAM")) dir.create("PAM")
lista <- list.files(pattern=".tif",recursive = T) %>% .[grep("EMwmeanByROC_",.)] %>% .[grep("TSSbin_",.)]
head(lista)
lista2 <- list.files(pattern=".tif",recursive = T) %>% .[grep("EMcvByROC_",.)] %>% .[grep("/c_",.)]
head(lista2)
m <- strsplit(lista,"/")
group <- unlist(lapply(m,"[",1))
country <- unlist(lapply(m,"[",2))
species <- unlist(lapply(m,"[",3))
modelos <- lapply(m,"[",4) %>% mapply(function(x,y) sub(paste("proj_",y,sep=""),"Present",x),.,species) %>% 
  sub("proj_","",.) %>%  mapply(function(x,y) paste(y,x,sep="_"),.,species) %>% unname(.)
m <- strsplit(lista2,"/")
group2 <- unlist(lapply(m,"[",1))
country2 <- unlist(lapply(m,"[",2))
species2 <- unlist(lapply(m,"[",3))
modelos2 <- lapply(m,"[",4) %>% mapply(function(x,y) sub(paste("proj_",y,sep=""),"Present",x),.,species2) %>% 
  sub("proj_","",.) %>%  mapply(function(x,y) paste(y,x,sep="_"),.,species2) %>% unname(.)
unique(country==country2);unique(group==group2);unique(species==species2);unique(modelos==modelos2) ##### just checking the order because I'm obsesive
g <- raster(nrows=180*12,ncols=360*12,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1,crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% as(., 'SpatialPixels')
plan(multisession)

for (i in 1:length(unique(country))){
  iden <- which(country==unique(country)[i])
  lista_cty <- lista[iden]
  lista_cty2 <- lista2[iden]
  cat("Processing - ",unique(country)[i],"\n")
  master_tib <- 1:length(lista_cty) %>% future_map_dfr(function(lis){
    r0 <- lista_cty[lis] %>% raster(.) 
    r1 <- r0 %>% rasterToPoints(.,fun=function(x){x==1}) %>% as.data.frame(.) %>% dplyr::select(.,x,y)
    rv <- lista_cty2[lis] %>% raster(.) %>% rasterToPoints(.) %>% as.data.frame(.) %>% rename(CV=names(.)[3])
    rv_data <- dplyr::select(rv,CV) %>% as_vector(.)
    rv %<>% dplyr::select(.,x,y) %>% SpatialPoints(.,proj4string= CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% 
      sp::over(.,g) %>% enframe(.,name="name") %>% mutate(.,name = modelos[iden][lis]) %>% rename(., CellID=value) %>% 
      mutate (.,CV=rv_data)
    if(dim(r1)[1]>1) { r1 %<>% SpatialPoints(.,proj4string= CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% 
        sp::over(.,g) %>% enframe(.,name="name") %>% mutate(.,name = modelos[iden][lis]) %>% rename(., CellID=value) %>% 
        mutate(.,CV=rv$CV[match(.$CellID,rv$CellID)]); return(r1)}
    if(dim(r1)[1] == 0) { r1 <- tibble(name=modelos[iden][lis],CellID=NA,CV=NA); return(r1)}
    if(dim(r1)[1] == 1) { r1 %<>% rbind(.,.) %>%  SpatialPoints(.,proj4string= CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% 
        sp::over(.,g) %>% enframe(.,name="name") %>% mutate(.,name = modelos[iden][lis]) %>% .[1,] %>% rename(., CellID=value) %>% 
        mutate(.,CV=rv$CV[match(.$CellID,rv$CellID)]); return(r1)}
  },.progress = TRUE)
  save(master_tib,file=paste("PAM/MasterTib",unique(country)[i],".Rdata",sep=""))
  save(lista_cty,file=paste("PAM/MasterList",unique(country)[i],".Rdata",sep=""))
  mastermod <- cbind(country[iden],species[iden],modelos[iden])
  colnames(mastermod) <- c("Group.Country","Species","Model")
  save(mastermod,file=paste("PAM/MasterModels",unique(country)[i],".Rdata",sep=""))
  cat(unique(country)[i],":",dim(master_tib)[1],"entries","\n")
  cat(unique(country)[i],":",dim(mastermod)[1],"models","\n")
  cat(unique(country)[i],":",length(unique(mastermod[,2])),"species","\n")
  rm(list=c("master_tib"))
}