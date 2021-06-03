library(rgdal)
library(raster)
library(sp)
library(tidyverse)
library(dplyr)
library(furrr)
library(purrr)
library(cowplot)
library(stringr)
library(GGally)
library(viridis)
library(crayon)
library(reshape2)
library(data.table)
library(rasterVis)
library(patchwork)
library("rnaturalearth")
library("rnaturalearthdata")
library("BBmisc")
library(sf)


#### Process the new IUCN download ####
#### This includes filtering and writing individual shapes
setwd("~/Desktop/3.Termohalina/redlist_species_data_4b96d3b0-d089-40e4-a6dd-1b0b13f37e21/SHPs")
data <- readOGR("../data_0.shp") 
unique <- unique(data@data$BINOMIAL)
paises <- readOGR("../../PaisesMegadiversos.shp")
counter <- length(unique)
for (i in 1:length(unique)) {
  counter <- counter - 1
  cat("processing species", counter, "\n")
  tmp <- data[data$BINOMIAL == unique[i], ] 
  overlay <- over(tmp,paises) %>% select(.,PAÍS) %>% as_vector(.)
  if(is.na(overlay)) next
  if(length(overlay)!=1) next
  cat("Found species!!","\n")
  writeOGR(tmp, dsn=getwd(), unique[i], driver="ESRI Shapefile",
           overwrite_layer=TRUE)
}

#### Estimate number of new shapes from IUCN by country ####
lista_shps <- list.files(pattern = "*.shp$",path="../redlist_species_data_4b96d3b0-d089-40e4-a6dd-1b0b13f37e21/SHPs/")
paises <- readOGR("../PaisesMegadiversos.shp")
plan(multisession)
por_pais <- lista_shps %>% future_map(function(i) { 
  sp_tmp <- readOGR(paste("../redlist_species_data_4b96d3b0-d089-40e4-a6dd-1b0b13f37e21/SHPs/",i,sep=""))
  overlay <- over(sp_tmp,paises) %>% dplyr::select(.,PAÍS) %>% as_vector(.)
  return(overlay)
},.progress = TRUE)
plan(sequential)

#### These data are for the same IUCN search
taxonomy <- readxl::read_xlsx("../redlist_species_data_5bda9f17-b941-4750-b55a-23e663041813/taxonomy.xlsx")
lista_shps <- lista_shps %>% substr(.,1,nchar(.)-4)
lista_shps <- bind_cols(name=lista_shps,country=unlist(por_pais))
lista_shps <- bind_cols(lista_shps,group=taxonomy$className[match(lista_shps$name,taxonomy$scientificName)])
lista_shps$group <- as_factor(lista_shps$group)
#### Check levels first!!!
levels(lista_shps$group)
levels(lista_shps$group) <- c("Magnoliophyta","Magnoliophyta","Monilophyta","Lycophyta")
lista_shps <- table(lista_shps$country,lista_shps$group) %>% as.data.frame.matrix(.) %>% mutate(country=rownames(.)) %>% pivot_longer(.,1:3)
lista_shps <- bind_cols(lista_shps,code=apply(lista_shps[,1:2],1,paste,collapse="_"))
lista_shps %>% arrange(name) %>% view()

#### Get taxonomy for IUCN species (shape files processed by Angela) and incorporate into general list ####
setwd("~/Desktop/3.Termohalina/REGISTROS_TODAS_v1/UICN_Plantas/")
lista_IUCN <- list.files(pattern="*.csv$",recursive=T)
x.x <- strsplit(lista_IUCN,split="/")
pais <- sapply(x.x,"[",1)
especie <- sapply(x.x,"[",2)
especie <- especie %>% substr(.,1,nchar(.)-4)
lista_taxonomy <- list.files(path="../../IUCN/",pattern="taxonomy.csv",recursive = T)
lista_taxonomy <- lapply(lista_taxonomy, function(x) fread(paste("../../IUCN/",x,sep=""))) %>% bind_rows(.) %>% as_tibble(.)
lista_taxonomy$scientificName <- sub(" ","_",lista_taxonomy$scientificName)
grupo <- lista_taxonomy$className[match(especie,lista_taxonomy$scientificName)]
length(grupo) == length(pais);length(grupo) == length(especie)
grupo <- grupo %>% as_factor(.)
levels(grupo) <- c("Magnoliophyta","Magnoliophyta","Lycophyta","Monilophyta")
tabla_IUCN <- bind_cols(file=lista_IUCN,Group=grupo,Country=pais,Species=especie)
tabla_IUCN$Species <- sub("_"," ",tabla_IUCN$Species)
dups <- c()
for (i in 1:nrow(tabla_IUCN)){ 
  cat(i,"\r")
  test <- read.table(tabla_IUCN$file[i],header=T,sep=",")
  ruta <- paste("~/Desktop/3.Termohalina/REGISTROS_TODAS_v1/",tabla_IUCN$Group[i],
                "/",tabla_IUCN$Country[i],
                "/",tabla_IUCN$Species[i],".csv",sep="")
  
  parent <- paste("~/Desktop/3.Termohalina/REGISTROS_TODAS_v1/",tabla_IUCN$Group[i],
                  "/",tabla_IUCN$Country[i],sep="")
  if(!dir.exists(parent)) dir.create(parent)
  if(file.exists(ruta)) {dups <- c(dups,i); next}
  write.table(test,file=ruta,quote=F,row.names = F,col.names = T,sep=",")
}


#### List occurrence files across all groups ####
#### Need to process before continuing: this might end up being deprecated with the full data.
setwd("~/Desktop/3.Termohalina/REGISTROS_TODAS_v1/")
lista <- list.files(pattern=".csv$",recursive = T)
m <- strsplit(lista,"/")
group <- unlist(lapply(m,"[",1))
country <- unlist(lapply(m,"[",2))
species <- unlist(lapply(m,"[",3)) %>% substr(.,1,nchar(.)-4)
data <- cbind(group,country,species) %>% as_tibble(.)
data$species[which(duplicated(data$species))] ### Checking if there are no duplicates
write.table(data,file="output/ListSpecies_v1.txt",sep=",",row.names = F,col.names = T,quote = F)
data %>% dplyr::select(.,species) %>% unique(.) %>% dim(.)
data %>% dplyr::select(.,group) %>% unique(.) %>% dim(.)
data <- data %>% group_by(.,country,group) %>% summarise(spp = n()) %>%  ungroup(.) %>% 
  mutate(group=factor(group)) %>% arrange(.,group)
data <- bind_cols(data,code=apply(data[,1:2],1,paste,collapse="_"))
richness <- read.table("output/SpeciesbyGroup.txt",header=T,sep=",") %>% 
  pivot_longer(.,-Country) %>% group_by(.,Country,name) %>% 
  mutate(group=factor(name)) %>% 
  arrange(.,name,Country)
richness <- bind_cols(richness,code=apply(richness[,1:2],1,paste,collapse="_"))
data <- bind_cols(data[match(richness$code,data$code),],richness[,c(3,5)])
names(data) <- c("country","group","spp","code_1","total","code_2")
data$code_1==data$code_2
#### Continue from above....
write.table(data,file = "output/TOTAL_spp_v1.txt",quote = F,row.names = F,col.names = F,sep=",")

##### Format data to create circular barplot by group ####
data2 <- data
x <- unlist(lapply(strsplit(data$code_2,"_"),"[[",1))
y <- unlist(lapply(strsplit(data$code_2,"_"),"[[",2))
data$group <- as.character(data$group)
data[which(is.na(data$group)),2] <- y[which(is.na(data$group))]
data[which(is.na(data$country)),1] <- x[which(is.na(data$country))]
which(is.na(data$country))
data.sum <- data %>% group_by(group) %>% summarise(Modeladas=sum(spp,na.rm=T)) %>% 
  mutate(Linaje=c("animal","animal","planta","planta","animal","planta","planta","animal")) %>% 
  mutate(Total=c(8181,11147,1340,291000,6495,11916,1058,11136)) %>% 
  mutate(Proporcion=Modeladas/Total * 100)

data.sum %>% group_by(Linaje) %>% summarise(Total_Mod=sum(Modeladas,na.rm=T),Total_Tot=sum(Total,na.rm=T)) %>% 
  mutate(Total_Mod/Total_Tot * 100)


filling = as.factor(data$group)
data$group <- factor(data$group,levels=c("Mammalia","Squamata","Aves","Amphibia","Lycophyta","Monilophyta","Pinophyta","Magnoliophyta"))
empty_bar <- 3
data <- data.frame( matrix(NA,empty_bar*nlevels(filling), ncol(data)) ) %>% 
  'colnames<-' (colnames(data)) %>% mutate(.,group=rep(levels(filling), each=empty_bar)) %>% 
  rbind(data, .) %>% arrange(group)  %>% mutate(.,id=seq(1, nrow(.))) %>% 
  mutate(spp=log(spp)+.1)
data[which(!is.finite(data$spp)&!is.na(data$spp)),"spp"] <- NA
#mutate(spp=(spp/total)*100)
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
base_data <- data %>% 
  group_by(group) %>%  ##### here!!
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]
data
grouping = data$country
y.limits <- floor(range(data$spp,na.rm = T)[2] + 3)
ys <- seq(2,12,2)
p <- ggplot(data, aes(x=as.factor(id), y=spp, fill=grouping)) +  ### here!!
  geom_bar(aes(x=as.factor(id), y=spp, fill=grouping), ### here!!
           stat="identity", alpha=0.5) +
  scale_fill_viridis_d(option="D") +
  geom_segment(data=grid_data, aes(x = end, y = ys[length(ys)-3], xend = start, yend = ys[length(ys)-3]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = ys[length(ys)-4], xend = start, yend = ys[length(ys)-4]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = ys[length(ys)-5], xend = start, yend = ys[length(ys)-5]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  annotate("text", x = rep(max(data$id),length(ys)-2), y = ys[1:4], 
           label = as.character(ys[1:4]), color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  geom_bar(aes(x=as.factor(id), y=spp, fill=grouping), stat="identity", alpha=0.5) + ### here!!
  ylim(-10,12) + theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")) + coord_polar() + 
  geom_text(data=label_data, aes(x=id, y = spp + 0.2, label=country, hjust=hjust), 
            color="black", fontface="bold",alpha=0.6, size=2, 
            angle= label_data$angle, inherit.aes = FALSE ) +
  geom_segment(data=base_data, aes(x = start, y = -0.5, xend = end, yend = -.5), 
               colour = "black", alpha=0.8, size=0.6, inherit.aes = FALSE )
ggdraw() + 
  draw_image("../Grupos/Siluetas/MammaliaAustralia_Tarsipes_rostratus.png",
             scale=.06,x=.06,y=0.135) + 
  draw_image("../Grupos/Siluetas/SquamataVenezuela_Gonatodes_rozei.png",
             scale=.07,x=.135,y=0.072) +
  draw_image("../Grupos/Siluetas/AvesMexico_Phaetornis_mexicanus.png",
             scale=.06,x=.125,y=-0.05) +
  draw_image("../Grupos/Siluetas/AmphibiaColombia_Atelopus_spurrelli.png",
             scale=.07,x=.065,y=-0.121) +
  draw_image("../Grupos/Siluetas/LycophytaBrazil_Huperzia_biformis2.png",
             scale=0.1,x=-.044,y=-0.126) + 
  draw_image("../Grupos/Siluetas/MonilophytaMadagascar_Alsophila_hyacinthei.png",
             scale=.08,x=-0.11,y=-0.05) + 
  draw_image("../Grupos/Siluetas/PinophytaChina_Ginkgo_biloba.png",
             scale=.06,x=-0.125,y=0.055) + 
  draw_image("../Grupos/Siluetas/MagnoliophytaPeru_Nasa_ranunculifolia.png",
             scale=.06,x=-0.06,y=0.138) + 
  draw_plot(p) +
  plot_annotation(title = "Vascular Plants and Vertebrates",caption="Santiago Ramírez Barahona") &
  theme(plot.title = element_text(hjust = 0.5,vjust=-20),title=element_text(size=12)) &
  theme(plot.caption = element_text(hjust = 0.5,vjust=40))
ggsave(filename = "output/SpeciesbyGroup.pdf")




##### Format data to create circular barplot by country ####
data2 -> data
x <- unlist(lapply(strsplit(data$code_2,"_"),"[[",1))
y <- unlist(lapply(strsplit(data$code_2,"_"),"[[",2))
x
y
data$group <- as.character(data$group)
data[which(is.na(data$group)),2] <- y[which(is.na(data$group))]
data[which(is.na(data$country)),1] <- x[which(is.na(data$country))]
which(is.na(data$country))
filling = as.factor(data$country)
data$country <- factor(data$country)
empty_bar <- 3
data <- data.frame( matrix(NA,empty_bar*nlevels(filling), ncol(data)) ) %>% 
  'colnames<-' (colnames(data)) %>% mutate(.,country=rep(levels(filling), each=empty_bar)) %>% 
  rbind(data, .) %>% arrange(country)  %>% mutate(.,id=seq(1, nrow(.))) %>% 
  mutate(spp=log(spp)+.1)
data[which(!is.finite(data$spp)&!is.na(data$spp)),"spp"] <- NA
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
base_data <- data %>% 
  group_by(country) %>%  ##### here!!
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]
data
grouping = data$group
y.limits <- floor(range(data$spp,na.rm = T)[2] + 3)
ys <- seq(2,12,2)
p <- ggplot(data, aes(x=as.factor(id), y=spp, fill=grouping)) +  ### here!!
  geom_bar(aes(x=as.factor(id), y=spp, fill=grouping), ### here!!
           stat="identity", alpha=0.5) +
  scale_fill_viridis_d(option="D") +
  geom_segment(data=grid_data, aes(x = end, y = ys[length(ys)-3], xend = start, yend = ys[length(ys)-3]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = ys[length(ys)-4], xend = start, yend = ys[length(ys)-4]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = ys[length(ys)-5], xend = start, yend = ys[length(ys)-5]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  annotate("text", x = rep(max(data$id),length(ys)-2), y = ys[1:4], 
           label = as.character(ys[1:4]), color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  geom_bar(aes(x=as.factor(id), y=spp, fill=grouping), stat="identity", alpha=0.5) + ### here!!
  ylim(-10,12) + theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")) + coord_polar() + 
  geom_text(data=label_data, aes(x=id, y = spp + 0.2, label=group, hjust=hjust), 
            color="black", fontface="bold",alpha=0.6, size=2, 
            angle= label_data$angle, inherit.aes = FALSE ) +
  geom_segment(data=base_data, aes(x = start, y = -0.5, xend = end, yend = -.5), 
               colour = "black", alpha=0.8, size=0.6, inherit.aes = FALSE )
ggdraw() + 
  draw_image("../Grupos/Siluetas/MammaliaAustralia_Tarsipes_rostratus.png",
             scale=.06,x=.06,y=0.135) + 
  draw_image("../Grupos/Siluetas/SquamataVenezuela_Gonatodes_rozei.png",
             scale=.07,x=.135,y=0.072) +
  draw_image("../Grupos/Siluetas/AvesMexico_Phaetornis_mexicanus.png",
             scale=.06,x=.125,y=-0.05) +
  draw_image("../Grupos/Siluetas/AmphibiaColombia_Atelopus_spurrelli.png",
             scale=.07,x=.065,y=-0.121) +
  draw_image("../Grupos/Siluetas/LycophytaBrazil_Huperzia_biformis2.png",
             scale=0.1,x=-.044,y=-0.126) + 
  draw_image("../Grupos/Siluetas/MonilophytaMadagascar_Alsophila_hyacinthei.png",
             scale=.08,x=-0.11,y=-0.05) + 
  draw_image("../Grupos/Siluetas/PinophytaChina_Ginkgo_biloba.png",
             scale=.06,x=-0.125,y=0.055) + 
  draw_image("../Grupos/Siluetas/MagnoliophytaPeru_Nasa_ranunculifolia.png",
             scale=.06,x=-0.06,y=0.138) + 
  draw_plot(p) +
  plot_annotation(title = "Vascular Plants and Vertebrates",caption="Santiago Ramírez Barahona") &
  theme(plot.title = element_text(hjust = 0.5,vjust=-20),title=element_text(size=12)) &
  theme(plot.caption = element_text(hjust = 0.5,vjust=40))
ggsave(filename = "output/SpeciesbyGroup.pdf")



#### List species ####
setwd("~/Documents/1.PROYECTOS/1.Termohalina/REGISTROS_TODAS_v1/")
lista <- list.files(pattern=".csv",recursive = T)
length(lista)
m <- strsplit(lista,"/")
group <- unlist(lapply(m,"[",1))
country <- unlist(lapply(m,"[",2))
species <- unlist(lapply(m,"[",3)) %>% substr(.,1,nchar(.)-4)
data <- cbind(group,country,species) %>% as_data_frame(.)
head(data);dim(data)
data$species[which(duplicated(data$species))]
data %>% as_data_frame(.) %>% dplyr::select(.,species) %>% unique(.) %>% dim(.)
data %>% as_data_frame(.) %>% dplyr::select(.,group) %>% unique(.) %>% dim(.)
pp <- data %>% group_by(.,group,country) %>% summarise(spp = n())
tabla_riqueza <- read.table("../output/tabla_spp_con_nuevas.csv",header=F,sep=",") %>% as_tibble(.)
mods <- tabla_riqueza %>%  group_by(V2) %>% summarise(.,tot=sum(V3,na.rm=T))
tots <- tabla_riqueza %>%  group_by(V2) %>% summarise(.,tot=sum(V5,na.rm=T))
mods <- mods %>% mutate(.,prop= 100*(tot / tots$tot)) %>% mutate(.,code=(substring(mods$V2,1,4))) %>% 
  arrange(.,prop) %>% filter(.,!is.na(V2))
p <- mods %>%  
ggplot(.,aes(x=as_factor(V2),y=prop,fill=fct_reorder(V2,prop))) + 
  geom_bar(stat="identity",col="black") +
  scale_fill_viridis_d(option="D") +
  ylim(c(0,30)) +
  theme_bw() +
  theme(legend.position = "none", legend.key.width = unit(1, "cm"),
        panel.grid = element_blank(),
        plot.background = element_rect(fill=NULL,inherit.blank = T),
              axis.title.x=element_blank(),
              axis.text.x=element_text(),
              legend.title=element_blank()) + 
  geom_text(aes(label=formatC(tot,big.mark = ","),fontface="bold"), vjust=1.6, color=c(rep("white",3),rep("black",5)), size=3.5) + 
  scale_x_discrete(labels = mods$code) + 
  labs(
    title="Number of species: vascular plants and vertebrates",
    caption = paste("Total number of species: ", formatC(sum(mods$tot,na.rm=T),big.mark = ",") ),
    subtitle = paste("Proportion of total sum across 12 mega-diverse countries")) + 
  ylab("Proportion (raw number in bars)")
p

ggdraw(p) +
  draw_image("../../Grupos/Siluetas/MonilophytaMadagascar_Alsophila_hyacinthei.png",
             scale=.09,x=-0.32,y=-.25) +
  draw_image("../../Grupos/Siluetas/LycophytaBrazil_Huperzia_biformis2.png",
             scale=0.15,x=-.22,y=-.21) + 
  draw_image("../../Grupos/Siluetas/MagnoliophytaPeru_Nasa_ranunculifolia.png",
             scale=.07,x=-0.12,y=-0.17) +
  draw_image("../../Grupos/Siluetas/AvesMexico_Phaetornis_mexicanus.png",
             scale=.06,x=-.02,y=-0.09) +
  draw_image("../../Grupos/Siluetas/PinophytaChina_Ginkgo_biloba.png",
             scale=.07,x=0.09,y=-0.01) + 
  draw_image("../../Grupos/Siluetas/AmphibiaColombia_Atelopus_spurrelli.png",
             scale=.07,x=.19,y=0.09) +
draw_image("../../Grupos/Siluetas/MammaliaAustralia_Tarsipes_rostratus.png",
           scale=.07,x=.32,y=0.112) + 
  draw_image("../../Grupos/Siluetas/SquamataVenezuela_Gonatodes_rozei.png",
             scale=.08,x=.42,y=0.246)



 




ggsave("../../TEST_v2/Numbers_prop.png")



#### PROCESS RASTERS (version 5 of PAMA_cluster)  ####
#### This reads all SDMs (rasters) and produces a tibble of presences ####
#### It produces three types of files:
#### - MasterTib* : presence data by cell (CellID)
#### - MasterList* : list of species modeled
#### - MasterModels : list of models processed
setwd("~/Desktop/3.Termohalina/TEST_v2/") #### This should be the master folder
#### The structure of the folder system in CONABIO is: /Group/Country[i]
#### The structure of the folder system in CCA is: /Group/GroupCountry[i]
if(!dir.exists("PAM")) dir.create("PAM")
lista <- list.files(pattern=".tif",recursive = T) %>% .[grep("EMwmeanByROC_",.)] %>% .[grep("TSSbin_",.)]
lista2 <- list.files(pattern=".tif",recursive = T) %>% .[grep("EMcvByROC_",.)] %>% .[grep("/c_",.)]
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
plan(sequential)


###### READ master_tibs #####
rm(list=ls())
load("DuplicatedSpecies.Rdata") ### duplicate species in the original data (source uncertain)
names(dobles)[c(6,10)] <- c("Brazil","Brazil_IUCN")
pais_target <- "India"
## Condition to remove duplicate or not (some are very specific)
remove_duplicates = FALSE
doube_china = FALSE
doube_brazil = FALSE
remove_within = FALSE
######
list_arch <- list.files(pattern="MasterTib",recursive = T)
list_arch <- list_arch[grep(pais_target,list_arch)]
#list_arch <- list_arch[1] ### for partial Australia
cat("MasterTibs to read:",list_arch,"\n")
list_tibs <- lapply(list_arch, function(x) mget(load(x)))
names(list_tibs) <- list_arch
for (i in 1:length(list_tibs)){
  list_tibs[[i]] <- list_tibs[[i]]$master_tib %>% #mutate(.data = .,Group=group) %>% 
    mutate(.data = .,Country=all_of(pais_target)) %>% bind_rows(.)
}
list_tibs

list_tibs[[2]] %>% select(name) %>% as_vector() -> spp 
sub("_.*","",spp) %>% unique(.)

#### for duplicates ####
who = 3
if(remove_within){  dobles %>% dplyr::select(1,grep(pais_target,names(dobles))) %>% 
    mutate(Sum=rowSums(.[,-1],na.rm=T)) %>% filter(Sum!=0) -> dobles
  dobles %>% select(Sum) %>% unique(); dobles
  dupli_spp <- dobles %>% select(especie) %>% unique() %>% as_vector()
  list_tibs[[who]] %>% select(name) -> target_names
  target_names %>% as_vector() %>% stringi::stri_split_fixed(.,pattern="_",simplify = F) %>% 
    lapply(.,"[",1) %>% unlist() -> target_names
  list_tibs[[who]] %>% distinct(name,CellID,.keep_all = T) -> list_tibs[[who]] 
  
}
if(remove_duplicates){ 
  dobles %>% dplyr::select(1,grep(pais_target,names(dobles))) %>% 
    mutate(Sum=rowSums(.[,-1],na.rm=T)) %>% filter(Sum!=0) -> dobles
  dobles %>% select(Sum) %>% unique(); dobles
  
  dupli_spp <- dobles %>% select(especie) %>% unique() %>% as_vector()
  list_tibs[[grep("IUCN",names(list_tibs))]] %>% select(name) -> target_names
  target_names %>% as_vector() %>% stringi::stri_split_fixed(.,pattern="_",simplify = F) %>% 
    lapply(.,"[",1) %>% unlist() -> target_names
  list_tibs[[grep("IUCN",names(list_tibs))]] %>% .[which(!target_names %in% dupli_spp),] -> list_tibs[[grep("IUCN",names(list_tibs))]]
  if(doube_brazil){list_tibs[[3]] %>% filter(!grepl("Anemia.organensis",name)) -> list_tibs[[3]]}
  
  if(doube_china){list_tibs[[4]] %>% filter(!grepl("China",name)) -> list_tibs[[4]]}
}
rm(target_names,dobles)
#### write final FINAL TIB : join data from mulitple sources
list_tibs <- bind_rows(list_tibs)
lista_mods <- list_tibs %>% select(name) %>% as_vector() %>% unique()
list_tibs <- bind_cols(as.data.table(stringi::stri_split_fixed(list_tibs$name,pattern="_",simplify = T)),list_tibs[,2:4]) %>% as_tibble(.)
cha <- c("Species"="V1","Scenario"="V2","Time"="V3")
list_tibs <- rename(list_tibs,all_of(cha))
cat(bold(" - Number of models:",length(lista_mods),"\n")) #### Este es el número de modelos que tenemos. 
spp <- sub("_.*","",lista_mods) %>% unique(.)
cat(green(bold("  -- Number of species:",length(spp),"\n"))) #### Este es el número de especies que tenemos. 
cat("CORRECT:",length(spp)*16 == length(lista_mods),"\n") #### Checando que todas las especies estén completamente modeladas
list_tibs <- list_tibs %>% .[which(!is.na(list_tibs$CellID)),]
list_tibs %>% mutate(Group=NA) -> list_tibs

list_arch <- list.files(pattern="MasterList",recursive = F)
list_arch <- list_arch[grep(pais_target,list_arch)]
lista <- lapply(list_arch, function(x) mget(load(x))) %>% bind_rows(.)
lista <- lapply(lista,strsplit,split="/")
grupo <- sapply(lista[[1]],"[",1)
pais <- sapply(lista[[1]],"[",2) #%>% mapply(sub,grupo,"",.,SIMPLIFY = T,USE.NAMES = F)
especie <- sapply(lista[[1]],"[",3)
modelo <- sapply(lista[[1]],"[",4)
raster <- sapply(lista[[1]],"[",5)
lista2 <- tibble(group=as.factor(grupo),pais=as.factor(pais),especie,modelo,raster)
lista2 %>% distinct(group,pais,especie,.keep_all = T) -> lista2
lista2 %>% filter(group %in% c("IUCN","UICN_Plantas")) -> IUCN
cat("Number of records to match:",dim(IUCN)[1],"\n")
if(dim(IUCN)[1]!=0){
  list.files(path="../../IUCN/", pattern = "taxonomy.csv",recursive = T) %>% 
    paste("../../IUCN/", .,sep="")  %>% lapply(.,FUN = read.csv) %>% rbindlist() %>% as_tibble() -> taxonomies
  taxonomies$className <- as.factor(taxonomies$className)
  levels(taxonomies$className) <- c("Pinophyta","Pinophyta","Pinophyta","Magnoliophyta","Lycophyta","Magnoliophyta","Pinophyta","Monilophyta")
  taxonomies$className <- as.character(taxonomies$className)
  taxonomies$scientificName <- sub(" ",".",taxonomies$scientificName)
  IUCN$especie %>% unique() -> especies
  taxonomies %>% filter(scientificName %in% especies) -> taxonomies
  IUCN$taxonomies <- taxonomies$className[match(IUCN$especie, taxonomies$scientificName)]
  IUCN$taxonomies %>% as_factor() -> IUCN$taxonomies
  levels(IUCN$taxonomies)
  IUCN %>% dplyr::select(especie,taxonomies) %>% filter(!duplicated(especie)) -> IUCN
  IUCN$taxonomies %>% as.character() -> IUCN$taxonomies
  length(IUCN$taxonomies)
  length(which(lista2$group %in% c("IUCN","UICN_Plantas")))
  lista2$group <- as.character(lista2$group)
  lista2[which(lista2$group %in% c("IUCN","UICN_Plantas")),]$group <- IUCN$taxonomies[match(lista2[which(lista2$group %in% c("IUCN","UICN_Plantas")),"especie"]$especie,IUCN$especie)]
  lista2 %>% select(group) %>% unique()
  list_tibs$Group <- lista2$group[match(list_tibs$Species,lista2$especie)]
}
#save(list_tibs,file=paste("../PROCESSED/",pais_target,"FinalTib.Rdata",sep=""))
save(list_tibs,file=paste("../PROCESSED/",pais_target,"sixPartial_","FinalTib.Rdata",sep=""))

#### Estimate species richness across scenarios/times per cell : SumTib*        ####
#### Also estimate range changes per species across scenarios/times : RangeTib* ####
#### Australia had to be divided into partial files to process.                 ####
list_tibs %>% select(Group) %>% unique() %>% print()
list_tibs %>% group_by(Group,CellID,Scenario,Time) %>% summarise(SR=n(),mean_cv=mean(CV,na.rm=T),sd_cv=sd(CV,na.rm=T))  -> summary_tib
g <- raster(nrows=180*12,ncols=360*12,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1,crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% as(., 'SpatialPixels')
g@coords[match(summary_tib$CellID,g@grid.index),] %>% as_tibble() %>% rename_with(~c("Longitud","Latitude")) -> coordas
summary_tib %>% bind_cols(.,coordas) -> summary_tib
summary_tib %>% mutate(Country=all_of(pais_target)) -> summary_tib
summary_tib %>% select(Group,Country,Scenario:sd_cv,CellID,Longitud,Latitude) -> summary_tib
summary_tib
#save(summary_tib,file=paste("../PROCESSED/",pais_target,"SumTib.Rdata",sep=""))
save(summary_tib,file=paste("../PROCESSED/",pais_target,"sixPartial_","SumTib.Rdata",sep=""))
list_tibs %>% group_by(Group,Species,Scenario,Time) %>% 
  summarise(Range=n()) %>% pivot_wider(names_from=c(Scenario,Time),values_from=Range) -> range_tib
range_tib %>% mutate(Country=all_of(pais_target)) -> range_tib
range_tib %>% select(Group,Country,Species,Present_, GreenControl_T1,GreenControl_T2,GreenControl_T3,
                     Green0.5_T1,Green0.5_T2,Green0.5_T3,Green1_T1,Green1_T2,Green1_T3,Green1.5_T1,Green1.5_T2,Green1.5_T3,
                     Green3_T1,Green3_T2,Green3_T3) -> range_tib
range_tib
save(range_tib,file=paste("../PROCESSED/",pais_target,"RangeTib.Rdata",sep=""))
#save(range_tib,file=paste("../PROCESSED/",pais_target,"sixPartial_","RangeTib.Rdata",sep=""))


##### Australia only : join partial files #####
pais_target="Australia"
lista_ranges <- list.files("../PROCESSED/",pattern = "Partial_SumTib")
lista_ranges <- paste0("../PROCESSED/",lista_ranges)
sum_tibs <- lapply(lista_ranges, function(x) mget(load(x)))
sum_tibs
for (i in 1:length(sum_tibs)){
  sum_tibs[[i]] <- sum_tibs[[i]]$summary_tib  %>% ungroup()
}
sum_tibs %>% bind_rows() -> sum_tibs

sum_tibs %>% group_by(Group,Country,Scenario,Time,CellID) %>% summarise(SR=sum(SR),mean_cv=mean(mean_cv),sd_cv=max(sd_cv),
                                                                        Longitud=first(Longitud),Latitude=first(Latitude)) %>% 
  select(Group,Country,Scenario,Time,SR,mean_cv,sd_cv,CellID,Longitud,Latitude) -> sum_tibs
sum_tibs
save(sum_tibs,file=paste("../PROCESSED/",pais_target,"SumTib.Rdata",sep=""))

pais_target="Australia"
lista_ranges <- list.files("../PROCESSED/",pattern = "Partial_RangeTib")
lista_ranges <- paste0("../PROCESSED/",lista_ranges)
sum_tibs <- lapply(lista_ranges, function(x) mget(load(x)))
sum_tibs
for (i in 1:length(sum_tibs)){
  sum_tibs[[i]] <- sum_tibs[[i]]$range_tib
}
sum_tibs %>% bind_rows() -> sum_tibs
sum_tibs
sum_tibs %>% distinct()
save(sum_tibs,file=paste("../PROCESSED/",pais_target,"RangeTib.Rdata",sep=""))



########

##### Files generated #####
##### MasterTib.R: data obtained from processing raster files.
##### FinalTib.R: filtered and processed MasterTib data (these are the data to upload).
##### RangeTib.R: per species range size estimates.
##### SummaryTib: per grid-cell richness estimates.


#### FROM HERE ON IS CODE TO PROCESS THE MASTER TIBS, SOME OF IT MIGHT BE REPETITIVE (Santiago FEB-2021) #####

#### Estimate species range loss across scenarios/times ####
library(GGally)
library(viridis)
library(tidyverse)
setwd("~/Documents/1.PROYECTOS/1.Termohalina/Final_v4/TIBS/")
list_arch <- list.files(pattern="RangeTib",recursive = F)
list_tibs <- lapply(list_arch, function(x) mget(load(x))) %>% lapply(function(x)x[[1]]) %>% bind_rows() %>% 
mutate_if(is.numeric , replace_na, replace = 0) %>% rename("Present"=Present_) -> range_tib
range_tib %>% mutate_if(is.numeric , replace_na, replace = 0) %>% mutate_at(5:ncol(.), ~((.x - Present)/Present)) %>% 
  mutate(Present = 0) -> STDrange_tib

#### Species range loss by group
STDrange_tib %>% ungroup() %>% group_by(Group) %>% summarise_at(3:18, quantile,probs=0.5) -> df_group
STDrange_tib %>% ungroup() %>% group_by(Country) %>% summarise_at(3:18, quantile,probs=0.5) -> df_country

STDrange_tib %>% group_by(Group,Country) %>% summarise(N=n()) %>% pull(N) -> N
STDrange_tib %>% ungroup() %>% mutate_if(is.numeric, ~ifelse(. > -1, NA, .)) %>% group_by(Group,Country) %>% summarise(across(.cols = -c(1),.fns = sum,na.rm=T)) %>% ungroup() %>% mutate_if(is.numeric, ~ . / all_of(N)) %>% group_by(Group,Country) %>% select(ends_with("_T1")) %>% select(-GreenControl_T1) %>% view()
una$Range

#



STDrange_tib %>% ungroup() %>% group_by(Group) %>% summarise(M1=median(GreenControl_T1),M2=median(Green0.5_T1),M3=median(Green1_T1),M4=median(Green1.5_T1),M5=median(Green3_T1)) %>% 
  dplyr::select(-1) %>% mutate(MT=apply(.,1,mean)) %>% dplyr::select(MT) %>% .[-c(1,2,5,8),] %>% colMeans()
nms=df_group$Group
sce=names(df_group)[-1]
df_group %>% dplyr::select(-1) %>% t(.) %>% as_tibble() %>% bind_cols(Scenario=sce,.) %>% rename_at(-1,~ all_of(nms)) %>% mutate_at(-1,~format(floor(.x * -1000) /-1000, nsmall = 2)) %>% .[c(1,2,5,8,11,14,3,6,9,12,15,4,7,10,13,16),] %>% filter(Scenario!="Present") %>% .[,c(1,4,5,7,8,2,3,6,9)] %>% flextable::flextable() %>% flextable::save_as_docx(.,path="Table_by_Group.docx")

STDrange_tib %>% ungroup() %>% distinct(Group) %>% pull() -> grupos
STDrange_tib %>% ungroup() %>% dplyr::select(-Present) %>% filter(Group%in%all_of(grupos[-c(1:4)])) %>% summarise_at(-c(1:3), quantile,probs=0.5) 

STDrange_tib %>% ungroup()  %>% group_by(Country)  %>% dplyr::select(-Present,-starts_with("GreenControl")) %>% summarise_at(-c(1:2), quantile,probs=0.5)  %>% select(ends_with("T3")) %>% summarise_all(list(max=min,min=max)) 



nms=df_country$Country
sce=names(df_country)[-1]
df_country %>% dplyr::select(-1) %>% t(.) %>% as_tibble() %>% bind_cols(Scenario=sce,.) %>% rename_at(-1,~ all_of(nms)) %>% mutate_at(-1,~format(floor(.x * -1000) /-1000, nsmall = 2))%>% .[c(1,2,5,8,11,14,3,6,9,12,15,4,7,10,13,16),] %>% filter(Scenario!="Present") %>% flextable::flextable() %>% flextable::save_as_docx(.,path="Table_by_Country.docx")

STDrange_tib %>% dplyr::select(-Present) %>% group_by(Country) %>% tally() %>% pull(n) -> tot_spp
case_complete <- function(x) {case_when(x == -1 ~ 1, x!= -1 ~ 0)}
STDrange_tib %>% dplyr::select(-Present) %>% ungroup() %>% mutate_if(is.numeric , case_complete) %>% group_by(Country) %>% summarise_at(-c(1:2),sum) %>% mutate_if(is.numeric, ~ .x / all_of(tot_spp)) %>% pivot_longer(cols=-1,names_to = "Scenario",values_to = "Lost_spp") %>% pivot_wider(names_from = Country,values_from = Lost_spp) %>% mutate_at(-1,~format(floor(.x * -1000) /-1000, nsmall = 2)) %>% .[c(1,4,7,10,13,2,5,8,11,14,3,6,8,12,15),] %>% flextable::flextable() %>% flextable::save_as_docx(.,path="../../TERMOS_MS/EXTRAS/Table_LostSpp.docx")

STDrange_tib %>% dplyr::select(-Present) %>% group_by(Group) %>% tally() %>% pull(n) -> tot_spp
STDrange_tib %>% dplyr::select(-Present) %>% ungroup() %>% mutate_if(is.numeric , case_complete) %>% group_by(Group) %>% summarise_at(-c(1:2),sum) %>% mutate_if(is.numeric, ~ .x / all_of(tot_spp)) %>% pivot_longer(cols=-1,names_to = "Scenario",values_to = "Lost_spp") %>% pivot_wider(names_from = Group,values_from = Lost_spp) %>% mutate_at(-1,~format(floor(.x * -1000) /-1000, nsmall = 2)) %>% .[c(1,4,7,10,13,2,5,8,11,14,3,6,8,12,15),c(1,2,3,6,9,4,5,7,8)] %>% flextable::flextable() %>% flextable::save_as_docx(.,path="../../TERMOS_MS/EXTRAS/Table_LostSpp_group.docx")


STDrange_tib %>% filter(Green0.5_T1==-1) -> a
a$Country %>% table(.) -> loss
a$Country %>% table(.) %>% sum() -> loss_sum
(loss/loss_sum * 100) [c(1,2,3,8,9)] %>% sum()
STDrange_tib$Country %>% table(.) -> total
plot(y=as.vector(loss)/as.vector(total),x=as.vector(total),xlab="Species modelled",
ylab="Proportion of complete loss"     )
id <- identify(y=as.vector(loss)/as.vector(total),x=as.vector(total))
sum(total[id])/sum(total) *100
sum(loss[id])/sum(loss) *100
(loss/total * 100)

STDrange_tib$Country %>% table(.) %>% sum() -> zz
z/zz * 100 -> z
z[order(z,decreasing = T)][1:5] %>% sum()
STDrange_tib$Group %>% table() -> z
z <- names(z)
STDrange_tib %>% ungroup() %>% 
  dplyr::select(Group,GreenControl_T1) %>% group_by(Group) %>% summarise(Median=median(GreenControl_T1,na.rm=T),N=n()) %>% ungroup() %>% as.data.frame() -> z

STDrange_tib %>% ungroup() %>% 
  dplyr::select(Group,Green0.5_T1) %>% group_by(Group) %>% summarise(Median=median(Green0.5_T1,na.rm=T),N=n()) %>% ungroup() %>% as.data.frame() %>% bind_cols(z,.) -> z
z[order(z[,2]),]
z[order(z[,5]),]

df_country$Country
df_country$Green0.5_T1
df_country$GreenControl_T3 %>% median()


time <- c("T1","T2","T3")
t_title <- c("2020-2060","2050-2080","2070-2100")
ylims <- range(df_group[,-1])
plotas <- list()
for (i in 1:3){ 
  if(i!=1) {
    plotas[[i]] <-  df_group %>% ungroup() %>% mutate(Group=as_factor(Group)) %>% dplyr::select(1,contains(time[i])) %>% 
    rename_with(~all_of(c("Group","RCP_8.5","Green_0.5","Green_1","Green_1.5","Green_3"))) %>% 
    ggparcoord(columns = c(2:ncol(.)), groupColumn = 1,scale="globalminmax",showPoints = TRUE, boxplot=F,
               title=t_title[i], alphaLines = 1) + 
    ylim(ylims) + ylab("Range Loss") +
      #scale_color_viridis(discrete=TRUE) +
      scale_color_brewer(palette = "Set1")+
   # theme_minimal() +
    theme(
      legend.position="bottom",
      plot.title = element_text(size=13),
      panel.background=element_blank(),
      panel.border=element_rect(colour="black",fill=NA),
      #axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    ) +
    xlab("")
  } 
  else { 
plotas[[i]] <- df_group %>% ungroup() %>% mutate(Group=as_factor(Group)) %>% dplyr::select(1,contains(time[i])) %>% 
  rename_with(~all_of(c("Group","RCP_8.5","Green_0.5","Green_1","Green_1.5","Green_3"))) %>% 
ggparcoord(columns = c(2:ncol(.)), groupColumn = 1,
           scale="globalminmax",
           showPoints = TRUE, boxplot=F,
           title=t_title[i],
           alphaLines = 1) + 
  ylim(ylims) + ylab("Range Loss") +
  #scale_color_viridis(discrete=TRUE) +
  scale_color_brewer(palette = "Set1") +
  #theme_minimal() +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=13),
    panel.background=element_blank(),
    panel.border=element_rect(colour="black",fill=NA),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  xlab("")
}
}
pdf("../../TERMOS_MS/EXTRAS/Loss_by_Group.pdf")
wrap_plots(plotas) +
  plot_layout(ncol=3,byrow=T,guides = "collect",tag_level = 'new') + 
  plot_annotation(title = paste("Median range loss by Taxonomic group"),
                  caption='Ureta et al. 2020') & theme(legend.position = 'bottom')
dev.off()

#### Species range loss by country
ylims <- range(df_country[,-1])
plotas <- list()
for (i in 1:3){ 
  if(i!=1) {
    plotas[[i]] <-  df_country %>% ungroup() %>% mutate(Country=as_factor(Country)) %>% dplyr::select(1,contains(time[i])) %>% 
      rename_with(~all_of(c("Country","RCP_8.5","Green_0.5","Green_1","Green_1.5","Green_3"))) %>% 
      ggparcoord(columns = c(2:ncol(.)), groupColumn = 1,scale="globalminmax",showPoints = TRUE, boxplot=F,
                 title=t_title[i], alphaLines = 1) + 
      ylim(ylims) + ylab("Range Loss") +
      #scale_color_viridis(discrete=TRUE) +
      scale_color_brewer(palette = "Set3")+
    #theme_minimal()+
      theme(
        legend.position="bottom",
        plot.title = element_text(size=13),
        panel.background=element_blank(),
        panel.border=element_rect(colour="black",fill=NA),
        #axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
      ) +
      xlab("")
  } 
  else { 
    plotas[[i]] <- df_country %>% ungroup() %>% mutate(Country=as_factor(Country)) %>% dplyr::select(1,contains(time[i])) %>% 
      rename_with(~all_of(c("Country","RCP_8.5","Green_0.5","Green_1","Green_1.5","Green_3"))) %>% 
      ggparcoord(columns = c(2:ncol(.)), groupColumn = 1,
                 scale="globalminmax",
                 showPoints = TRUE, boxplot=F,
                 title=t_title[i],
                 alphaLines = 1) + 
      ylim(ylims) + ylab("Range Loss") +
      #scale_color_viridis(discrete=TRUE) +
      scale_color_brewer(palette = "Set3")+
   # theme_minimal()+
      theme(
        legend.position="bottom",
        plot.title = element_text(size=13),
        panel.background=element_blank(),
        panel.border=element_rect(colour="black",fill=NA),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
      ) +
      xlab("")
  }
}
pdf("../../TERMOS_MS/EXTRAS/Loss_by_Country.pdf")
wrap_plots(plotas) +
  plot_layout(ncol=3,byrow=T,guides = "collect",tag_level = 'new') + 
  plot_annotation(title = paste("Median range loss by Country"),
                  caption='Ureta et al. 2020') & theme(legend.position = 'bottom')
dev.off()




setwd("~/Documents/1.PROYECTOS/1.Termohalina/Final_v4/TIBS/")
list_arch <- list.files(pattern="SumTib",recursive = F)
list_tibs <- lapply(list_arch, function(x) mget(load(x))) %>% lapply(function(x)x[[1]]) %>% bind_rows() 
list_tibs %>% ungroup() %>% group_by(Country,Scenario,Time,CellID) -> sum_tib
sum_tib %>% dplyr::select(Longitud,Latitude,SR) %>% 
  summarise(SR=sum(SR),Longitud=first(Longitud),Latitude=first(Latitude)) %>% 
  ungroup() %>% group_by(Country) %>% mutate(SR_std = (SR-min(SR))/(max(SR)-min(SR))) -> df_tib

lista_ras <- list.files(path="~/Documents/7.MAPOTECA/WorldClim_2/mawc2.0_10m_bio/",pattern=".tif$",full.names = T)
stacka <- raster::stack(lista_ras)
names(stacka)
g <- raster(nrows=180*12,ncols=360*12,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1,crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% as(., 'SpatialPixels')
cordas <- g@coords[df_tib$CellID,]
raster::extract(stacka[[c(1,4,12,15)]],cordas) -> bio_clims
raster::extract(raster("~/Documents/7.MAPOTECA/WorldClim_2/wc2.1_30s_elev.tif"),cordas) -> altas
df_tib %>% ungroup() %>% mutate(AnnTemp = bio_clims[,1], AnnPrec = bio_clims[,3],Alt = altas) -> df_tib

df_tib %>% filter(CellID==4811874) %>% view
df_tib %>% pivot_wider(names_from = c(Scenario,Time), values_from = SR)%>% filter(CellID==4811874)


df %>% ggplot(aes(x=AnnTemp,y=Scenario,fill=Scenario)) + 
  geom_density_ridges_gradient(scale = 2.5, rel_min_height = 0.001,stat="binline",binwidth=.5) +
  NULL
  #scale_fill_gradientn(colours = MaizePal::maize_pal("HighlandMAGIC")) +
  #labs(title =paste('Gene flow into ',meta_target$SourcePop," (",round(meta_target$Score_source,2),")",sep="")) +
  ylab("") +
  geom_vline(xintercept = 1,linetype="dashed",size=1) +
  theme(
    panel.grid=element_blank(),
    panel.grid.major.y=element_line(colour="grey80"),
    panel.border = element_rect(colour="grey50",fill=NA),
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )


  
  
  scenario="Present"
df_tib %>% filter(Scenario==scenario) -> df
ylims <- sum_tib %>% ungroup() %>%  dplyr::select(Latitude) %>% range(.,na.rm=T) %>% + c(-2,2)
xlims <- sum_tib %>% ungroup() %>% dplyr::select(Longitud) %>% range(.,na.rm=T) %>% + c(-2,2)
val_lims <- range(df$SR)
nonstd <- ggplot() + xlim(xlims) + ylim(ylims) +
  geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="grey80",fill="grey80") +
  geom_tile(data=df,aes(y=Latitude, x=Longitud,fill = SR),colour = "white",size=0.00001,height=0.1,width=0.1) + 
  scale_fill_viridis(option = "inferno",direction=-1,begin=0,end=1,limits=val_lims) +
  # geom_sf(data = world) +
  theme(legend.position = "left", legend.key.width = unit(2, "cm"),
        panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),
        panel.background=element_rect(colour="black",fill="azure2"),
        panel.border=element_rect(colour="black",fill=NA)) + 
  labs(subtitle="non-standarized") +
  NULL

val_lims <- range(df$SR_log)
std <- ggplot() + xlim(xlims) + ylim(ylims) +
  geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="grey80",fill="grey80") +
  geom_tile(data=df,aes(y=Latitude, x=Longitud,fill = SR_log),colour = "white",size=0.00001,height=0.1,width=0.1) + 
  scale_fill_viridis(option = "inferno",direction=-1,begin=0,end=1,limits=val_lims) +
  # geom_sf(data = world) +
  theme(legend.position = "left", legend.key.width = unit(2, "cm"),
        panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),
        panel.background=element_rect(colour="black",fill="azure2"),
        panel.border=element_rect(colour="black",fill=NA)) + 
  labs(subtitle="standardized (0-1)") +
  NULL
wrap_plots(list(nonstd,std)) +
  plot_layout(ncol=1,nrow = 2,byrow=F,tag_level = 'new') + 
  plot_annotation(title = paste("Species Richness"),
                  caption='Ureta et al. 2020') & theme(legend.position = 'bottom')


scenario="GreenControl"
time="T1"
df_tib %>% filter(Scenario==scenario,Time==time) -> df
ylims <- sum_tib %>% ungroup() %>%  dplyr::select(Latitude) %>% range(.,na.rm=T) %>% + c(-2,2)
xlims <- sum_tib %>% ungroup() %>% dplyr::select(Longitud) %>% range(.,na.rm=T) %>% + c(-2,2)
val_lims <- range(df$SR)
nonstd_control <- ggplot() + xlim(xlims) + ylim(ylims) +
  geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="grey80",fill="grey80") +
  geom_tile(data=df,aes(y=Latitude, x=Longitud,fill = SR),colour = "white",size=0.00001,height=0.1,width=0.1) + 
  scale_fill_viridis(option = "inferno",direction=-1,begin=0,end=1,limits=val_lims) +
  # geom_sf(data = world) +
  theme(legend.position = "left", legend.key.width = unit(2, "cm"),
        panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),
        panel.background=element_rect(colour="black",fill="azure2"),
        panel.border=element_rect(colour="black",fill=NA)) + 
  labs(subtitle="non-standarized") +
  NULL
val_lims <- range(df$SR_std)
std_control <- ggplot() + xlim(xlims) + ylim(ylims) +
  geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="grey80",fill="grey80") +
  geom_tile(data=df,aes(y=Latitude, x=Longitud,fill = SR_std),colour = "white",size=0.00001,height=0.1,width=0.1) + 
  scale_fill_viridis(option = "inferno",direction=-1,begin=0,end=1,limits=val_lims) +
  # geom_sf(data = world) +
  theme(legend.position = "left", legend.key.width = unit(2, "cm"),
        panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),
        panel.background=element_rect(colour="black",fill="azure2"),
        panel.border=element_rect(colour="black",fill=NA)) + 
  labs(subtitle="standardized (0-1)") +
  NULL

wrap_plots(list(nonstd_control,std_control)) +
  plot_layout(ncol=1,nrow = 2,byrow=F,tag_level = 'new') + 
  plot_annotation(title = paste("Species Richness",scenario,time),
                  caption='Ureta et al. 2020') & theme(legend.position = 'bottom')



scenario="Green0.5"
time="T1"
df_tib %>% filter(Scenario==scenario,Time==time) -> df
val_lims <- range(df$SR)
nonstd_green0.5 <- ggplot() + xlim(xlims) + ylim(ylims) +
  geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="grey80",fill="grey80") +
  geom_tile(data=df,aes(y=Latitude, x=Longitud,fill = SR),colour = "white",size=0.00001,height=0.1,width=0.1) + 
  scale_fill_viridis(option = "inferno",direction=-1,begin=0,end=1,limits=val_lims) +
  # geom_sf(data = world) +
  theme(legend.position = "left", legend.key.width = unit(2, "cm"),
        panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),
        panel.background=element_rect(colour="black",fill="azure2"),
        panel.border=element_rect(colour="black",fill=NA)) + 
  labs(subtitle="non-standarized") +
  NULL

val_lims <- range(df$SR_std)
std_green0.5 <- ggplot() + xlim(xlims) + ylim(ylims) +
  geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="grey80",fill="grey80") +
  geom_tile(data=df,aes(y=Latitude, x=Longitud,fill = SR_std),colour = "white",size=0.00001,height=0.1,width=0.1) + 
  scale_fill_viridis(option = "inferno",direction=-1,begin=0,end=1,limits=val_lims) +
  # geom_sf(data = world) +
  theme(legend.position = "left", legend.key.width = unit(2, "cm"),
        panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),
        panel.background=element_rect(colour="black",fill="azure2"),
        panel.border=element_rect(colour="black",fill=NA)) + 
  labs(subtitle="standardized (0-1)") +
  NULL

wrap_plots(list(nonstd_green0.5,std_green0.5)) +
  plot_layout(ncol=1,nrow = 2,byrow=F,tag_level = 'new') + 
  plot_annotation(title = paste("Species Richness",scenario,time),
      caption='Ureta et al. 2020') & theme(legend.position = 'bottom')



setwd("~/Documents/1.PROYECTOS/1.Termohalina/Final_v4/TIBS/")
list_arch <- list.files(pattern="SumTib",recursive = F)
list_tibs <- lapply(list_arch, function(x) mget(load(x))) %>% lapply(function(x)x[[1]]) %>% bind_rows() 
list_tibs %>% ungroup() %>% group_by(Country,Scenario,Time,CellID) -> sum_tib
sum_tib %>% dplyr::select(Longitud,Latitude,SR) %>% 
  summarise(SR=sum(SR),Longitud=first(Longitud),Latitude=first(Latitude)) %>% 
  ungroup() %>% group_by(Country) %>% mutate(SR_std = (SR-min(SR))/(max(SR)-min(SR))) -> df_tib

df_tib %>% filter(Scenario=="Present") %>% summarize(MED=ceiling((max(SR)*0.6))) -> umbral

plotas <- list()
result_tib <- list()
for (i in 1:12){ 
   df_tib %>% filter(Country == umbral$Country[i], SR > umbral$MED[i]) %>% group_by(Scenario,Time) %>% summarize(PSH=n()) %>% mutate(Scenario=factor(Scenario,levels=c("Present","GreenControl","Green0.5","Green1","Green1.5","Green3"))) -> temp 
  present <- temp %>% filter(Time=="") %>% ungroup() %>%  dplyr::select(PSH) %>% as_vector()
  temp %>% mutate(PSH = PSH/ present) -> result_tib[[i]]
 # plotas[[i]] <- temp %>% mutate(PSH = PSH/ present) %>% 
 #  ggplot(aes(#x=Scenario,
 #             x=Time,
 #             y=PSH)) + geom_bar(stat="identity") + ggtitle(umbral$Country[i])
}
names(result_tib) <- umbral$Country
for (i in 1:12) result_tib[[i]]<-result_tib[[i]] %>% mutate(Country=umbral$Country[i],.before=Scenario)
result_tib %>% do.call(rbind,.) %>% pivot_wider(names_from = c(Scenario,Time),values_from = PSH) %>% select(-starts_with("Present")) %>% mutate_if(is.numeric,replace_na,replace=0) %>% mutate_at(-1,~format(floor(.x * -1000) /-1000, nsmall = 2)) %>% pivot_longer(cols=-1,names_to = "Scenario",values_to = "PSH") %>% pivot_wider(names_from = Country,values_from = PSH) %>% .[c(1,4,7,10,13,2,5,8,11,14,3,6,9,12,15),] %>% flextable::flextable() %>% flextable::save_as_docx(.,path="../../TERMOS_MS/EXTRAS/Table_PSH.docx")

result_tib %>% do.call(rbind,.) %>% pivot_wider(names_from = c(Scenario,Time),values_from = PSH) %>% select(-starts_with("Present")) %>% mutate_if(is.numeric,replace_na,replace=0) %>% mutate_at(-1,~format(floor(.x * -1000) /-1000, nsmall = 2)) %>% pivot_longer(cols=-1,names_to = "Scenario",values_to = "PSH") %>% pivot_wider(names_from = Country,values_from = PSH) %>% .[c(1,4,7,10,13,2,5,8,11,14,3,6,9,12,15),] %>% mutate_at(-1,as.numeric) %>% rowwise(Scenario) %>% mutate(Median= 1 - median(c_across(-1))) %>% arrange(Median)


pdf("~/Desktop/FIGURES/PSH_green0.5.pdf")
wrap_plots(plotas) +
  plot_layout(ncol=4,nrow = 3,byrow=F,tag_level = 'new') + 
  plot_annotation(title = "Species Richness (Green 0.5)",
                  caption='Ureta et al. 2020') & theme(legend.position = 'bottom')
dev.off()


### ESTIMATE CHANGE IN RICHNESS
sum_tib %>% summarise(Tot_SR=sum(SR)) -> una
una
una %>% arrange(Time) %>% ungroup() %>% pivot_wider(names_from=c(Scenario,Time),values_from=Tot_SR) %>% rename_with(~"Present",Present_) %>% ungroup() -> una_dos

una_dos %>%  mutate_at(-c(1:2),replace_na,replace=0) %>% mutate_at(-c(1:3), ~.x - Present) %>% group_by(Country) %>% summarise_at(-c(1:2),quantile,prob=0.25) %>% mutate_if(is.numeric,~case_when(.x > 0 ~ 0,.x <= 0 ~.x)) -> xx_neg

una_dos %>%  mutate_at(-c(1:2),replace_na,replace=0) %>% mutate_at(-c(1:3), ~.x - Present) %>% group_by(Country) %>% summarise_at(-c(1:2),quantile,prob=0.75) %>% mutate_if(is.numeric,~case_when(.x < 0 ~ 0,.x >= 0 ~.x)) -> xx_pos

una_dos %>% mutate_at(-c(1:2),replace_na,replace=0) %>% mutate_at(-c(1:3), ~.x - Present) %>% group_by(Country) %>% filter(Country=="Australia") %>% summarise_at(-c(1:2), .funs = c(Positive = ~ sum(.x > xx_pos[1,1])))

xx_pos %>% filter(Country=="Australia")

pivot_longer(cols=-1,names_to = "Scenario",values_to = "Extent") -> una_tres

strsplit(una$Scenario,"_") -> dos
una %>% mutate(Scenario=lapply(dos,"[",1) %>% unlist,Time=lapply(dos,"[",2) %>% unlist,Change=lapply(dos,"[",3) %>% unlist) %>% filter(Change=="Positive") %>% pivot_wider(names_from = Country,values_from = Extent) %>% .[c(5,10,15,1,6,11,2,7,12,3,8,13,4,9,14),]



#### PROCESS MASTERTIBS #####
list_arch <- list.files(pattern="MasterList",recursive = T)
lista <- lapply(list_arch, function(x) mget(load(x))) %>% bind_rows(.)
long <- dim(lista)[1]
cat("Number of models in tibble:", format(dim(lista)[1],big.mark   =","),"\n") #### Este es el número de modelos que tenemos. 
lista <- lapply(lista,strsplit,split="/")
grupo <- sapply(lista[[1]],"[",1)
pais <- sapply(lista[[1]],"[",2) #%>% mapply(sub,grupo,"",.,SIMPLIFY = T,USE.NAMES = F)
especie <- sapply(lista[[1]],"[",3)
modelo <- sapply(lista[[1]],"[",4)
raster <- sapply(lista[[1]],"[",5)
lista2 <- tibble(grupo,pais,especie,modelo,raster)
dim(lista2)[1] == length(lista[[1]]) ### sólo checando
cat("Number of species in tibble:", format(length(unique(lista2$especie)) ,big.mark   =","),"\n")
cat("All species were modelled completely (16 models):",length(unique(lista2$especie)) * 16 == long,"\n") # Checking if all species were modelled completely

list_arch <- list.files(pattern="MasterTib",recursive = T)
list_tibs <- lapply(list_arch, function(x) mget(load(x)))
tabla <- lista2 %>% filter(!duplicated(.[,1:3]))
pais <- sub("MasterTib","",list_arch) %>% sub(".Rdata","",.)

list_info <- list()
for (i in 1:length(list_tibs)){ 
 list_info[[i]] <- list_tibs[[i]]$master_tib$name %>% strsplit(.,"_") %>% lapply(.,"[",1)  %>% 
   unlist(.) %>% match(.,tabla$especie) %>% tabla[.,]  %>% select(1:2)
 list_tibs[[i]] <- list_tibs[[i]]$master_tib
}
list_tibs <- bind_rows(list_tibs);list_tibs
list_tibs
list_info
groups = unique(list_info[[1]]$grupo)
countries = unique(list_info[[1]]$pais)
country_choice <- countries[1]
group_choice <- groups

sub_tibs <- list_tibs[which(list_info[[1]]$pais == country_choice & list_info[[1]]$grupo %in% group_choice),]
lista_mods <- unique(sub_tibs$name)
cat(bold(" - Number of models:",length(lista_mods),"\n")) #### Este es el número de modelos que tenemos. 
spp <- sub("_.*","",lista_mods) %>% unique(.)
cat(green(bold("  -- Number of species:",length(spp),"\n")))
sub_tibs <- sub_tibs %>% .[which(!is.na(sub_tibs$CellID)),]
sub_tibs <- bind_cols(as.data.table(stringi::stri_split_fixed(sub_tibs$name,pattern="_",simplify = T)),sub_tibs[,2:3]) %>% as_tibble(.)
cha <- c("Species"="V1","Scenario"="V2","Time"="V3")
sub_tibs <- rename(sub_tibs,all_of(cha))
sub_tibs <- sub_tibs %>% mutate(.,Country=pais)
sub_tibs <- sub_tibs %>% bind_cols(.,tabla[match(sub_tibs$Species,tabla$especie),1])
sub_tibs


#### ESTIMATE RICHNESS ####
#plan(sequential)
system.time({
  rich_MASTER <- 1:length(group_choice) %>% map(function(i){
    cat("Processing",group_choice[i],"\n")
    r1 <- sub_tibs %>% filter(.,grupo == group_choice[i]) %>% group_by(Scenario,Time,CellID) %>% 
      summarize(SR=n(),mean_cv=mean(CV,na.rm=T),sd_cv=sd(CV,na.rm=T))
    temp <- as.data.frame(g@coords[match(r1$CellID,g@grid.index),])
    r1 <- bind_cols(temp,r1) %>% as_tibble(.) %>% mutate(.,Country=countries[i],Group=group_choice[i])
    return(r1)
  })
})



list_tibs

#rich_MASTER <- rbindlist(rich_MASTER)
poly<-rgdal::readOGR("~/Documents/7.MAPOTECA/official/wwf_terr_ecos.shp")
for (i in 1:length(rich_MASTER)){
pointos <- SpatialPoints(rich_MASTER[[i]][,1:2],proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
points_eco <- sp::over(pointos,poly)
points_eco$BIOME <- as.factor(points_eco$BIOME)
rich_MASTER[[i]] <-  mutate(rich_MASTER[[i]],Biome=points_eco$BIOME)
rich_MASTER[[i]]$Biome <- as.factor(rich_MASTER[[i]]$Biome)
levels(rich_MASTER[[i]]$Biome) <- c("TropMoistForest","TropDryForests","TropConiferForests","TempBroadleafForests",
                               "TempConiferForests","TropGrasslands","TempGrasslands","FloodedGrasslands",
                               "MontaneGrasslands","MediterraneanForests","Deserts","Mangroves",NA)
rich_MASTER[[i]]$Biome <- forcats::fct_explicit_na(rich_MASTER[[i]]$Biome,"Undet. Biome")
}
i=2
xlims <- rich_MASTER[[i]] %>% #filter(.,Country==country) %>% 
  dplyr::select(.,x)  %>% range(.+2)
ylims <- rich_MASTER[[i]] %>% #filter(.,Country==country) %>% 
  dplyr::select(.,y)  %>% range(.+2)
val_lims <- rich_MASTER[[i]] %>% #filter(.,Country==country) %>% 
  dplyr::select(.,SR)  %>% range(.)


###### Richness: Same time, different scenarios ####
setwd("~/Documents/1.PROYECTOS/1.Termohalina/")
load("~/Downloads/MadagascarSumTib.Rdata")
plotas <- list()
summary_tib
ord <- unique(rich_MASTER[[i]]$Scenario)
ord <- ord[c(6,5,1,2,3,4)]
world <- ne_coastline(scale=110,returnclass = "sf")
times = "T3"
years = 2070
for (sc in 1:length(ord)) {
if(ord[sc] == "Present") temp <- rich_MASTER[[i]] %>% filter(.,Time=="") %>% filter(.,Scenario==ord[sc])
if(ord[sc] != "Present") temp <- rich_MASTER[[i]] %>% filter(.,Time==times) %>% filter(.,Scenario==ord[sc])
plotas[[sc]] <- ggplot() + xlim(xlims) + ylim(ylims) +
  geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="grey80",fill="grey80") +
  geom_tile(data=temp,aes(y=y, x=x,fill = SR),colour = "white",size=0.00001,height=0.1,width=0.1) + 
  scale_fill_viridis(option = "inferno",direction=-1,begin=0,end=1,limits=val_lims) +
  geom_sf(data = world) +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),
        panel.background=element_rect(colour="black",fill="azure2"),
        panel.border=element_rect(colour="black",fill=NA)) + 
  labs(subtitle=paste(LETTERS[1:6][sc]))
}
wrap_plots(plotas) +
  plot_layout(byrow=T,guides = "collect",tag_level = 'new') + 
    plot_annotation(title = paste("Species richness for",country_choice),subtitle = paste(years,sep="-"),
                  caption='Ureta et al. 2020') & theme(legend.position = 'bottom')

myList <- setNames(lapply(vector("list", ncol(range_tib)), function(x) x <- 0), names(range_tib))
range_tib %>% ungroup() %>% replace_na(myList) %>% mutate(GC_1= (GreenControl_T1 - Present_) / Present_) %>% select(GC_1) %>% 
  as_vector() %>% max()
range_tib %>% ungroup() %>% 
  mutate(Group=dplyr::recode(Group,"Anfibios"="Amphibia", "Reptiles"="Squamata")) %>% 
  select(Group) %>% unique()


#### Coef.var: NEED to loop over times and standardise mean_cv ####
normalized <- function(x) (x- min(x))/(max(x) - min(x))
rich_MASTER[[i]] <- rich_MASTER[[i]] %>% 
  mutate(Norm_mean_cv = scale(mean_cv,scale=F,center=F))
rich_MASTER[[i]] 
xlims <- rich_MASTER[[i]] %>% #filter(.,Country==country) %>% 
  dplyr::select(.,x)  %>% range(.+2)
ylims <- rich_MASTER[[i]] %>% #filter(.,Country==country) %>% 
  dplyr::select(.,y)  %>% range(.+2)
val_lims <- rich_MASTER[[i]] %>% #filter(.,Country==country) %>% 
  dplyr::select(.,Norm_mean_cv)  %>% range(.)
plotas <- list()
ord <- unique(rich_MASTER[[i]]$Scenario)
ord <- ord[c(6,5,1,2,3,4)]
world <- ne_coastline(scale=110,returnclass = "sf")
for (sc in 1:length(ord)) {
  if(ord[sc] == "Present") temp <- rich_MASTER[[i]] %>% filter(.,Time=="") %>% filter(.,Scenario==ord[sc])
  if(ord[sc] != "Present") temp <- rich_MASTER[[i]]  %>% filter(.,Time=="T1") %>% filter(.,Scenario==ord[sc])
  plotas[[sc]] <- ggplot() + xlim(xlims) + ylim(ylims) +
    geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="grey80",fill="grey80") +
    geom_tile(data=temp,aes(y=y, x=x,fill = Norm_mean_cv),colour = "white",size=0.00001,height=0.1,width=0.1) + 
    scale_fill_viridis(option = "inferno",direction=1,begin=0.2,end=0.95,limits=val_lims) +
    geom_sf(data = world) +
    theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
          panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
          axis.text.x = element_blank(),axis.text.y = element_blank(),
          panel.background=element_rect(colour="black",fill="azure2"),panel.border=element_rect(colour="black",fill=NA)) + 
    labs(subtitle=paste(ord[sc]))
}
wrap_plots(plotas) +
  plot_layout(byrow=T,guides = "collect",tag_level = 'new') + 
  plot_annotation(title = paste("Species richness for",country_choice),subtitle = paste(years,sep="-"),
                  caption='Ureta et al. 2020') & theme(legend.position = 'bottom')

######## ESTIMATE CHANGE IN SPECIES RIHCNESS: UNLIMITED DISPERSAL (this is deprecated) ####
mm <- rich_MASTER[[i]] %>% #filter(.,Country==country) %>% 
  as_tibble(.)
scenario <- c("GreenControl","Green0.5","Green1","Green1.5","Green3")
sc=2
mm_fut <- mm %>% filter(.,Scenario==scenario[sc])
times <- c("Pre","T1","T2","T3")
mm_temp <- list()
mm_temp[[1]] <- mm %>% filter(.,Scenario=="Present") %>% dplyr::select(.,CellID,SR,x,y,Biome)
for (ti in 2:length(times)){ 
mm_temp[[ti]] <- mm_fut %>% filter(.,Time==times[ti]) %>% dplyr::select(.,CellID,SR,x,y)
}
RT <- Reduce(function(...) merge(...,by="CellID",all.x=TRUE,all.y=TRUE),mm_temp)
head(RT)
colnames(RT)[2:4] <- c("SR_pre","x","y")
colnames(RT)[6:8] <- c("SR_T1","x_1","y_1")
colnames(RT)[9:11] <- c("SR_T2","x_2","y_2")
colnames(RT)[12:14] <- c("SR_T3","x_3","y_3")
RT$x[is.na(RT$x)] = RT$x_1[is.na(RT$x)];RT$x[is.na(RT$x)] = RT$x_2[is.na(RT$x)];RT$x[is.na(RT$x)] = RT$x_3[is.na(RT$x)]
RT$y[is.na(RT$y)] = RT$y_1[is.na(RT$y)];RT$y[is.na(RT$y)] = RT$y_2[is.na(RT$y)];RT$y[is.na(RT$y)] = RT$y_3[is.na(RT$y)]
RT <- RT[,-c(7,8,10,11,13,14)]
RT$SR_pre[is.na(RT$SR_pre)] = 0
RT$SR_T1[is.na(RT$SR_T1)] = 0
RT$SR_T2[is.na(RT$SR_T2)] = 0; RT$SR_T3[is.na(RT$SR_T3)] = 0
RT <- RT %>%mutate(.,Change_T1=SR_T1-SR_pre) %>% mutate(.,Change_T2=SR_T2-SR_pre) %>% mutate(.,Change_T3=SR_T3-SR_pre)
head(RT);dim(RT);length(unique(RT$CellID))

df2 <- as_tibble(RT)
xlims <- df2 %>% dplyr::select(.,x)  %>% range(.+2)
ylims <- df2 %>% dplyr::select(.,y)  %>% range(.+2)
val_lims <- range(df2[,8:10],na.rm=T)
plotas <- list()
ord <- unique(rich_MASTER[[i]]$Time)
ord <- ord[c(1:3)]
world <- ne_coastline(scale=110,returnclass = "sf")
this <- c(9:11)

for(ti in 1:length(ord)){
  df3 <- df2 %>% dplyr::select(.,x,y,this[ti])
  colnames(df3)[3] <- "Change"
plotas[[ti]] <- ggplot() + xlim(xlims) + ylim(ylims) +
  geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="grey80",fill="grey80") +
  geom_tile(data=df3,aes(y=y, x=x,fill = Change),colour = "white",size=0.00001,height=0.1,width=0.1) + 
  scale_fill_viridis(option = "inferno",direction=1,begin=0.0,end=1,limits=val_lims) +
  geom_sf(data = world) +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),
        panel.background=element_rect(colour="black",fill="azure2"),panel.border=element_rect(colour="black",fill=NA)) + 
  labs(subtitle=paste(ord[ti]))
}

titulo = "Projected change in species richness (unlimited dispersal)"
titulo = "Projected change in species richness (limited dispersal)"
wrap_plots(plotas) +
  plot_layout(byrow=T,guides = "collect",tag_level = 'new') + 
  plot_annotation(title = paste(titulo,scenario[sc]),subtitle = paste(country_choice,sep="-"),
                  caption='Ureta et al. 2020') & theme(legend.position = 'bottom')

######## ESTIMATE CHANGE IN SPECIES RIHCNESS: LIMITED & UNLIMITED DISPERSAL ####
setwd("~/Documents/1.PROYECTOS/1.Termohalina/TEST_v4/")
load("PeruFinalTib.Rdata")
list_tibs$Group %>% unique()
g <- raster(nrows=180*12,ncols=360*12,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1,crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% as(., 'SpatialPixels')
### Assign coordinates to all presences
temp <- as.data.frame(g@coords[match(list_tibs$CellID,g@grid.index),])
list_tibs <- bind_cols(temp,list_tibs) %>% as_tibble(.)
rich_MASTER <- list_tibs
rich_MASTER <- rich_MASTER %>% as_tibble(.) %>% mutate(SR=1)
rich_MASTER
species_list <- rich_MASTER %>% dplyr::select(.,Species) %>% unique(.) %>% as_vector(.)
scenario <- unique(rich_MASTER$Scenario); scenario <- scenario[which(scenario!="Present")]
rich_MASTER$Time[which(rich_MASTER$Time=="")] <- "Pre"
times <- c("Pre","T1","T2","T3")
#### Global analysis..... for all species
weights = c(0) ### These are the weight to be applied: c(0) for limited dispersal, c(1) for unlimited dispersal
df3 <- list()
counter <- 0
for (spp in 1:length(species_list)){
  counter <- counter + 1
  if(counter==1) tictoc::tic()
  cat(yellow(bold("Processing",species_list[spp],"\n")))
mm <- rich_MASTER %>% #filter(.,Country==country) %>% 
  filter(.,Species==species_list[spp]) %>% as_tibble(.)
df2 <- list()
for (sc in 1:length(scenario)){
  cat(green(bold("Processing",scenario[sc],"\n")))
mm_fut <- mm %>% filter(.,Scenario==scenario[sc])
if(dim(mm_fut)[1]==0) {
  df2[[sc]] <- mm %>% filter(.,Scenario=="Present") %>% dplyr::select(.,Species,Scenario,CellID,SR,x,y) %>% 
    mutate(.,ChangeT1=-1,Change_T2=-1,Change_T3=-1) %>% mutate(.,Scenario=scenario[sc])
} #### This creates -1 changes  for all times when the scenario is missing
mm_temp <- list()
mm_temp[[1]] <- mm %>% filter(.,Scenario=="Present") %>% dplyr::select(.,Species,CellID,SR,x,y)
for (ti in 2:length(times)){ 
  mm_temp[[ti]] <- mm_fut %>% filter(.,Time==times[ti]) %>% dplyr::select(.,Species,Scenario,CellID,SR,x,y)
}
mm_temp #### The order is as in the vecrtor times: present, T1, T2, T3.
RT <- Reduce(function(...) merge(...,by="CellID",all.x=TRUE,all.y=TRUE),mm_temp)
colnames(RT)[2:5] <- c("Species","SR_pre","x","y")
colnames(RT)[6:10] <- c("Species_T1","Scenario_T1","SR_T1","x_1","y_1")
colnames(RT)[11:15] <- c("Species_T2","Scenario_T2","SR_T2","x_2","y_2")
colnames(RT)[16:20] <- c("Species_T3","Scenario_T3","SR_T3","x_3","y_3")
RT$Species <- species_list[spp]
RT$x[is.na(RT$x)] = RT$x_1[is.na(RT$x)];RT$x[is.na(RT$x)] = RT$x_2[is.na(RT$x)];RT$x[is.na(RT$x)] = RT$x_3[is.na(RT$x)]
RT$y[is.na(RT$y)] = RT$y_1[is.na(RT$y)];RT$y[is.na(RT$y)] = RT$y_2[is.na(RT$y)];RT$y[is.na(RT$y)] = RT$y_3[is.na(RT$y)]
if(dim(mm_temp[[4]])[1]==0) RT <- RT[,-c(16:20)]
if(dim(mm_temp[[3]])[1]==0) RT <- RT[,-c(11:15)]
if(dim(mm_temp[[2]])[1]==0) RT <- RT[,-c(6:10)]
head(RT)
RT$SR_pre[is.na(RT$SR_pre)] <- 0

if(dim(mm_temp[[2]])[1]!=0) {RT$SR_T1[is.na(RT$SR_T1)] = 0
    RT <- RT %>% mutate(.,Change_T1=SR_T1 - SR_pre)} else  RT <- RT %>% mutate(.,Change_T1=-1)

if(dim(mm_temp[[3]])[1]!=0) {RT$SR_T2[is.na(RT$SR_T2)] = 0
    RT <- RT %>% mutate(.,Change_T2=SR_T2 - SR_pre)} else  RT <- RT %>% mutate(.,Change_T2=-1)

if(dim(mm_temp[[4]])[1]!=0) {RT$SR_T3[is.na(RT$SR_T3)] = 0
    RT <- RT %>% mutate(.,Change_T3=SR_T3 - SR_pre)} else  RT <- RT %>% mutate(.,Change_T3=-1)
df2[[sc]] <- as_tibble(RT) %>% dplyr::select(.,CellID,Species,x,y,Change_T1,Change_T2,Change_T3) %>% 
  mutate(.,Scenario=scenario[sc])
}
df3[[spp]] <- bind_rows(df2)

ct1 <- df3[[spp]] %>% dplyr::select(.,Change_T1) %>% as_vector(.) %>% dplyr::recode(.,`1` = 1L*weights)

ct2 <- df3[[spp]] %>% dplyr::select(.,Change_T2) %>% as_vector(.) %>% dplyr::recode(.,`1` = 1L*weights)

ct3 <- df3[[spp]] %>% dplyr::select(.,Change_T3) %>% as_vector(.) %>% dplyr::recode(.,`1` = 1L*weights)

df3[[spp]] <- df3[[spp]] %>% mutate(.,Change_T1=ct1) %>% mutate(.,Change_T2=ct2) %>% mutate(.,Change_T3=ct3)

if(counter==20) {clocka <- tictoc::toc()
  counter <- 0
  remaining <- (clocka$toc - clocka$tic) * (length(species_list)/20)
  cat(red(bold("ETA (hrs)",remaining/3600,"\n")))}
}
df4 <- bind_rows(df3)

#### Plot the desired scenario.....
scenario_choice = "GreenControl"
df5 <- df4 %>%  filter(.,Scenario==scenario_choice) %>% group_by(.,CellID) %>% 
  summarise(.,Global_changeT1=sum(Change_T1,na.rm=T),Global_changeT2=sum(Change_T2,na.rm=T),Global_changeT3=sum(Change_T3,na.rm=T))
df5 <- df5 %>% mutate(.,x=df4$x[match(.$CellID,df4$CellID)])
df5 <- df5 %>% mutate(.,y=df4$y[match(.$CellID,df4$CellID)])
xlims <- rich_MASTER %>% #filter(.,Country %in% country) %>% 
  dplyr::select(.,x)  %>% range(.+2)
ylims <- rich_MASTER %>% #filter(.,Country %in% country) %>% 
  dplyr::select(.,y)  %>% range(.+2)
val_lims = range(c(df5$Global_changeT1,df5$Global_changeT2,df5$Global_changeT2))
plotas <- list()
ord <- unique(rich_MASTER$Time)
ord <- ord[c(1:3)]
years <- c(2030,2050,2070)
world <- ne_coastline(scale=110,returnclass = "sf")
this <- c(2:4)
df5
for(ti in 1:length(ord)){
  df6 <- df5 %>% dplyr::select(.,x,y,this[ti])
  colnames(df6)[3] <- "Change"
  plotas[[ti]] <- ggplot() + xlim(xlims) + ylim(ylims) +
  geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="grey80",fill="grey80") +
  geom_tile(data=df6,aes(y=y, x=x,fill = Change),colour = "white",size=0.00001,height=0.1,width=0.1) + 
  scale_fill_viridis(option = "inferno",direction=1,begin=0.0,end=1,limits=val_lims) +
  geom_sf(data = world) +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),
        panel.background=element_rect(colour="black",fill="azure2"),panel.border=element_rect(colour="black",fill=NA)) + 
  labs(subtitle=paste(ord[ti],":",years[ti]))
}
titulo = "Projected change in species richness (unlimited dispersal)"
wrap_plots(plotas) +
  plot_layout(byrow=T,guides = "collect",tag_level = 'new') + 
  plot_annotation(title = paste(titulo,scenario[sc]),
                  caption='Ureta et al. 2020') & theme(legend.position = 'bottom')

unique(df5$CellID==unlim$CellID)
plot(unlim$Global_changeT1, df5$Global_changeT1)
abline(h=0);abline(v=0)

limired <- df5




####### ESTIMATE CHANGE IN SPECIES RIHCNESS: CONTINUOS DISPERSAL (takes a while) ####
multiplier = 50
counter = 0
df3 <- list()
plan(multisession)
for (spp in 1:length(species_list)){
  counter <- counter + 1
  if(counter==1) tictoc::tic()
  mm <- rich_MASTER %>% #filter(.,Country==country) %>% 
    filter(.,Species==species_list[spp]) %>% as_tibble(.)
group_spp <- unique(mm$Group)
cat(yellow(bold(species_list[spp],"\n")))
  list1.sf <- mm %>% filter(.,Scenario == "Present") %>% 
    st_as_sf( coords = c("x", "y"), crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") )
  list2.sf <- mm %>% filter(.,Scenario != "Present") %>%
    st_as_sf( coords = c("x", "y"), crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") )
  
  # if(dim(list2.sf)[1]!=0) {
  #   present_dists <- list1.sf %>% dplyr::mutate(., target.Cell = sf::st_nearest_feature( geometry, list1.sf )) 
  # present_dists <- sf::st_distance(present_dists, present_dists[present_dists$target.Cell,] ,"Great Circle",by_element = F)
  # diag(present_dists) = NA
  # present_dists <- ceiling(max(apply(present_dists,1,max,na.rm=T))/1000)
  # tt <- list2.sf %>% dplyr::mutate(., target.Cell = sf::st_nearest_feature( geometry, list1.sf )) %>% 
  #   mutate(., Distance = sf::st_distance(geometry, list1.sf[target.Cell,] ,"Great Circle", by_element = TRUE)) %>% as_tibble(.) %>% 
  #   mutate(.,Distance = as.numeric(Distance/1000))
  # tt$w <- 1-(tt$Distance/(present_dists/multiplier))
  # }
  
  ##### TEST.....
  if(dim(list2.sf)[1]!=0) {
  system.time({
  res <- split(list2.sf,cut(1:dim(list2.sf)[1],10)) %>% future_map(function(i){
  tt <- i %>% dplyr::mutate(., target.Cell = sf::st_nearest_feature( geometry, list1.sf )) %>% 
    mutate(., Distance = sf::st_distance(geometry, list1.sf[target.Cell,] ,"Great Circle", by_element = TRUE)) %>% as_tibble(.) %>% 
    mutate(.,Distance = as.numeric(Distance/1000))
  tt$w <- 1-(tt$Distance/(present_dists/multiplier))
  return(tt)
  },.progress = TRUE)
  })
}
  
  df2 <- list()
  for (sc in 1:length(scenario)){
    cat(green(bold("------ Processing",scenario[sc],"\n")))
    mm_fut <- mm %>% filter(.,Scenario==scenario[sc])
    if(dim(mm_fut)[1]==0) {
      df2[[sc]] <- mm %>% filter(.,Scenario=="Present") %>% dplyr::select(.,Species,Scenario,CellID,SR,x,y) %>% 
        mutate(.,ChangeT1=-1,Change_T2=-1,Change_T3=-1) %>% mutate(.,Scenario=scenario[sc])
    } #### This creates -1 changes  for all times when the scenario is missing
    mm_temp <- list()
    mm_temp[[1]] <- mm %>% filter(.,Scenario=="Present") %>% dplyr::select(.,Species,CellID,SR,x,y)
    for (ti in 2:length(times)){ 
      mm_temp[[ti]] <- mm_fut %>% filter(.,Time==times[ti]) %>% dplyr::select(.,Species,Scenario,CellID,SR,x,y)
    }
    mm_temp #### The order is as in the vecrtor times: present, T1, T2, T3.
    RT <- Reduce(function(...) merge(...,by="CellID",all.x=TRUE,all.y=TRUE),mm_temp)
    colnames(RT)[2:5] <- c("Species","SR_pre","x","y")
    colnames(RT)[6:10] <- c("Species_T1","Scenario_T1","SR_T1","x_1","y_1")
    colnames(RT)[11:15] <- c("Species_T2","Scenario_T2","SR_T2","x_2","y_2")
    colnames(RT)[16:20] <- c("Species_T3","Scenario_T3","SR_T3","x_3","y_3")
    RT$Species <- species_list[spp]
    RT$x[is.na(RT$x)] = RT$x_1[is.na(RT$x)];RT$x[is.na(RT$x)] = RT$x_2[is.na(RT$x)];RT$x[is.na(RT$x)] = RT$x_3[is.na(RT$x)]
    RT$y[is.na(RT$y)] = RT$y_1[is.na(RT$y)];RT$y[is.na(RT$y)] = RT$y_2[is.na(RT$y)];RT$y[is.na(RT$y)] = RT$y_3[is.na(RT$y)]
    if(dim(mm_temp[[4]])[1]==0) RT <- RT[,-c(16:20)]
    if(dim(mm_temp[[3]])[1]==0) RT <- RT[,-c(11:15)]
    if(dim(mm_temp[[2]])[1]==0) RT <- RT[,-c(6:10)]
    head(RT)
    RT$SR_pre[is.na(RT$SR_pre)] <- 0
    
    if(dim(mm_temp[[2]])[1]!=0) {RT$SR_T1[is.na(RT$SR_T1)] = 0
    RT <- RT %>% mutate(.,Change_T1=SR_T1 - SR_pre)} else  RT <- RT %>% mutate(.,Change_T1=-1)
    
    if(dim(mm_temp[[3]])[1]!=0) {RT$SR_T2[is.na(RT$SR_T2)] = 0
    RT <- RT %>% mutate(.,Change_T2=SR_T2 - SR_pre)} else  RT <- RT %>% mutate(.,Change_T2=-1)
    
    if(dim(mm_temp[[4]])[1]!=0) {RT$SR_T3[is.na(RT$SR_T3)] = 0
    RT <- RT %>% mutate(.,Change_T3=SR_T3 - SR_pre)} else  RT <- RT %>% mutate(.,Change_T3=-1)
    df2[[sc]] <- as_tibble(RT) %>% select(.,CellID,Species,x,y,Change_T1,Change_T2,Change_T3) %>% 
      mutate(.,Scenario=scenario[sc])
  }
  
  df3[[spp]] <- bind_rows(df2)
  df3[[spp]] <- df3[[spp]] %>% mutate(.,Country=country_spp) %>% mutate(.,Group=group_spp)
  
  if(dim(list2.sf)[1]==0) {
    df3[[spp]] <- df3[[spp]] %>% mutate(.,w = 1)
  }
  
  df3[[spp]] <- df3[[spp]] %>% mutate(.,w=tt$w[match(df3[[spp]]$CellID,tt$CellID)])
  df3[[spp]] <- df3[[spp]] %>% 
    mutate(.,w= case_when(.$w < 0  ~  0,
                                  TRUE ~ as.numeric(.$w))) %>%
        mutate(.,Change_T1 = case_when(.$Change_T1 == 1  ~  .$Change_T1*.$w,
                                                                TRUE ~ as.numeric(.$Change_T1))) %>% 
    mutate(.,Change_T2 = case_when(.$Change_T2 == 1  ~  .$Change_T2*.$w,
                                   TRUE ~ as.numeric(.$Change_T2))) %>% 
    mutate(.,Change_T3 = case_when(.$Change_T3 == 1  ~  .$Change_T3*.$w,
                                   TRUE ~ as.numeric(.$Change_T3)))
  if(counter==20) {clocka <- tictoc::toc()
  counter <- 0
  remaining <- (clocka$toc - clocka$tic) * (length(species_list)/20)
  cat(red(bold("ETA (hrs)",remaining/3600,"\n")))}
}
df4 <- bind_rows(df3)
plan(sequential)
# AQUI SE ESPECIFICA EL ESCENARIO PARA HACER LAS GRAFICAS
scenario
sc=1
df5 <- df4 %>%  filter(.,Scenario==scenario[sc]) %>% group_by(.,CellID) %>% 
  summarise(.,Global_changeT1=sum(Change_T1,na.rm=T),
            Global_changeT2=sum(Change_T2,na.rm=T),
            Global_changeT3=sum(Change_T3,na.rm=T))

df5 <- df5 %>% mutate(.,x=df4$x[match(.$CellID,df4$CellID)]) %>% 
  mutate(.,y=df4$y[match(.$CellID,df4$CellID)])
xlims <- rich_MASTER %>% filter(.,Country %in% country) %>% dplyr::select(.,x)  %>% range(.+2)
ylims <- rich_MASTER %>% filter(.,Country %in% country) %>% dplyr::select(.,y)  %>% range(.+2)

val_lims_real = round(range(c(df5$Global_changeT1,df5$Global_changeT2,df5$Global_changeT2)),1)
val_lims = c(-8,3)
plotas <- list()
ord <- unique(rich_MASTER$Time)
ord <- ord[c(2:4)]
years <- c(2030,2050,2070)
world <- ne_coastline(scale=110,returnclass = "sf")
this <- c(2:4)
df5
breaks = diff(val_lims)
colores = inferno(breaks+1,direction = 1)

for(ti in 1:length(ord)){
  df6 <- df5 %>% dplyr::select(.,x,y,this[ti])
  colnames(df6)[3] <- "Change"
  plotas[[ti]] <- ggplot() + xlim(xlims) + ylim(ylims) +
    geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="grey80",fill="grey80") +
    geom_tile(data=df6,aes(y=y, x=x,fill = Change),colour = "white",size=0.00001,height=0.1,width=0.1) + 
    #scale_fill_viridis(option = "inferno",direction=1,begin=0.0,end=1,limits=val_lims) +
    scale_fill_stepsn(colours = colores,limits=val_lims,n.breaks=breaks,guide="legend") +
    geom_sf(data = world) +
    theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
          panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
          axis.text.x = element_blank(),axis.text.y = element_blank(),
          panel.background=element_rect(colour="black",fill="azure2"),panel.border=element_rect(colour="black",fill=NA)) + 
    labs(subtitle=paste(ord[ti],":",years[ti]))
}
#pdf("~/Desktop/Change_CD.5.pdf")
titulo = "Projected change in species richness (max_dist mpl=100)"
wrap_plots(plotas) +
  plot_layout(byrow=T,guides = "collect",tag_level = 'new') + 
  plot_annotation(title = paste(titulo,scenario[sc]),subtitle = paste(group,country,sep="-"),
      caption=paste('Range of change:',val_lims_real[1],"-",val_lims_real[2])) & theme(legend.position = 'bottom')
#dev.off()





####### ESTIMATE CHANGE IN SPECIES RIHCNESS: CONTINUOS DISPERSAL (using futures!!) ####
rich_MASTER <- 1:length(group_choice) %>% map(function(i){
    cat("Processing",group_choice[i],"\n")
    r1 <- sub_tibs %>% filter(.,grupo == group_choice[i])
    temp <- as.data.frame(g@coords[match(r1$CellID,g@grid.index),])
    r1 <- bind_cols(temp,r1) %>% as_tibble(.) %>% mutate(.,Country=countries[i],Group=group_choice[i])
    return(r1)
  })

rich_MASTER <- rbindlist(rich_MASTER)
rich_MASTER <- rich_MASTER %>% as_tibble(.) %>% mutate(Country=country_choice) %>% 
  mutate(SR=1)
rich_MASTER

species_list <- rich_MASTER %>% #filter(.,Country %in% country) %>% 
  dplyr::select(.,Species) %>% unique(.) %>% as_vector(.)
scenario <- unique(rich_MASTER$Scenario) %>% .[which(. !="Present")]
times <- c("Pre","T1","T2","T3")
plan(multisession)
multiplier = 50
counter = 0
df3 <- list()
for (spp in 1:length(species_list)){
  counter <- counter + 1
  if(counter==1) tictoc::tic()
  mm <- rich_MASTER %>% #filter(.,Country==country) %>% 
    filter(.,Species==species_list[spp]) %>% as_tibble(.)
  group_spp <- unique(mm$Group)
  cat(yellow(bold(species_list[spp],"\n")))
  list1.sf <- mm %>% filter(.,Scenario == "Present") %>% 
    st_as_sf( coords = c("x", "y"), crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") )
  list2.sf <- mm %>% filter(.,Scenario != "Present") %>%
    st_as_sf( coords = c("x", "y"), crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") )
  
  if(dim(list2.sf)[1]!=0) {
    present_dists <- list1.sf %>% dplyr::mutate(., target.Cell = sf::st_nearest_feature( geometry, list1.sf )) 
    present_dists <- sf::st_distance(present_dists, present_dists[present_dists$target.Cell,] ,"Great Circle",by_element = F)
    diag(present_dists) = NA
    present_dists <- ceiling(max(apply(present_dists,1,max,na.rm=T))/1000)
    #system.time({
    # tt <- list2.sf %>% dplyr::mutate(., target.Cell = sf::st_nearest_feature( geometry, list1.sf )) %>% 
    # mutate(., Distance = sf::st_distance(geometry, list1.sf[target.Cell,] ,"Great Circle", by_element = TRUE)) %>% as_tibble(.) %>% 
    # mutate(.,Distance = as.numeric(Distance/1000))
    # tt$w <- 1-(tt$Distance/(present_dists/multiplier))
    # }
    # })
  ##### TEST.....
  #system.time({
    if(dim(list2.sf)[1] < 100) {
      tt <- list2.sf %>% dplyr::mutate(., target.Cell = sf::st_nearest_feature( geometry, list1.sf )) %>% 
       mutate(., Distance = sf::st_distance(geometry, list1.sf[target.Cell,] ,"Great Circle", by_element = TRUE)) %>% as_tibble(.) %>% 
       mutate(.,Distance = as.numeric(Distance/1000))
       tt$w <- 1-(tt$Distance/(present_dists/multiplier))
    }
    if(dim(list2.sf)[1] >= 100) {
      tt <- split(list2.sf,cut(1:dim(list2.sf)[1],10)) %>% future_map_dfr(function(i){
      res <- i %>% dplyr::mutate(., target.Cell = sf::st_nearest_feature( geometry, list1.sf )) %>% 
        mutate(., Distance = sf::st_distance(geometry, list1.sf[target.Cell,] ,"Great Circle", by_element = TRUE)) %>% as_tibble(.) %>% 
        mutate(.,Distance = as.numeric(Distance/1000))
      res$w <- 1-(res$Distance/(present_dists/multiplier))
      return(res)
    },.progress = TRUE)}
  #})
  } 
  
  df2 <- list()
  for (sc in 1:length(scenario)){
    cat(green(bold("------ Processing",scenario[sc],"\n")))
    mm_fut <- mm %>% filter(.,Scenario==scenario[sc])
    if(dim(mm_fut)[1]==0) {
      df2[[sc]] <- mm %>% filter(.,Scenario=="Present") %>% dplyr::select(.,Species,Scenario,CellID,SR,x,y) %>% 
        mutate(.,ChangeT1=-1,Change_T2=-1,Change_T3=-1) %>% mutate(.,Scenario=scenario[sc])
    } #### This creates -1 changes  for all times when the scenario is missing
    mm_temp <- list()
    mm_temp[[1]] <- mm %>% filter(.,Scenario=="Present") %>% dplyr::select(.,Species,CellID,SR,x,y)
    for (ti in 2:length(times)){ 
      mm_temp[[ti]] <- mm_fut %>% filter(.,Time==times[ti]) %>% dplyr::select(.,Species,Scenario,CellID,SR,x,y)
    }
    mm_temp #### The order is as in the vecrtor times: present, T1, T2, T3.
    RT <- Reduce(function(...) merge(...,by="CellID",all.x=TRUE,all.y=TRUE),mm_temp)
    colnames(RT)[2:5] <- c("Species","SR_pre","x","y")
    colnames(RT)[6:10] <- c("Species_T1","Scenario_T1","SR_T1","x_1","y_1")
    colnames(RT)[11:15] <- c("Species_T2","Scenario_T2","SR_T2","x_2","y_2")
    colnames(RT)[16:20] <- c("Species_T3","Scenario_T3","SR_T3","x_3","y_3")
    RT$Species <- species_list[spp]
    RT$x[is.na(RT$x)] = RT$x_1[is.na(RT$x)];RT$x[is.na(RT$x)] = RT$x_2[is.na(RT$x)];RT$x[is.na(RT$x)] = RT$x_3[is.na(RT$x)]
    RT$y[is.na(RT$y)] = RT$y_1[is.na(RT$y)];RT$y[is.na(RT$y)] = RT$y_2[is.na(RT$y)];RT$y[is.na(RT$y)] = RT$y_3[is.na(RT$y)]
    if(dim(mm_temp[[4]])[1]==0) RT <- RT[,-c(16:20)]
    if(dim(mm_temp[[3]])[1]==0) RT <- RT[,-c(11:15)]
    if(dim(mm_temp[[2]])[1]==0) RT <- RT[,-c(6:10)]
    head(RT)
    RT$SR_pre[is.na(RT$SR_pre)] <- 0
    
    if(dim(mm_temp[[2]])[1]!=0) {RT$SR_T1[is.na(RT$SR_T1)] = 0
    RT <- RT %>% mutate(.,Change_T1=SR_T1 - SR_pre)} else  RT <- RT %>% mutate(.,Change_T1=-1)
    
    if(dim(mm_temp[[3]])[1]!=0) {RT$SR_T2[is.na(RT$SR_T2)] = 0
    RT <- RT %>% mutate(.,Change_T2=SR_T2 - SR_pre)} else  RT <- RT %>% mutate(.,Change_T2=-1)
    
    if(dim(mm_temp[[4]])[1]!=0) {RT$SR_T3[is.na(RT$SR_T3)] = 0
    RT <- RT %>% mutate(.,Change_T3=SR_T3 - SR_pre)} else  RT <- RT %>% mutate(.,Change_T3=-1)
    df2[[sc]] <- as_tibble(RT) %>% dplyr::select(.,CellID,Species,x,y,Change_T1,Change_T2,Change_T3) %>% 
      mutate(.,Scenario=scenario[sc])
  }
  
  df3[[spp]] <- bind_rows(df2)
  df3[[spp]] <- df3[[spp]] %>% #mutate(.,Country=country_spp) %>% 
    mutate(.,Group=group_spp)
  
  if(dim(list2.sf)[1]==0) {
    df3[[spp]] <- df3[[spp]] %>% mutate(.,w = 1)
  }
  
  df3[[spp]] <- df3[[spp]] %>% mutate(.,w=tt$w[match(df3[[spp]]$CellID,tt$CellID)])
  df3[[spp]] <- df3[[spp]] %>% 
    mutate(.,w= case_when(.$w < 0  ~  0,
                          TRUE ~ as.numeric(.$w))) %>%
    mutate(.,Change_T1 = case_when(.$Change_T1 == 1  ~  .$Change_T1*.$w,
                                   TRUE ~ as.numeric(.$Change_T1))) %>% 
    mutate(.,Change_T2 = case_when(.$Change_T2 == 1  ~  .$Change_T2*.$w,
                                   TRUE ~ as.numeric(.$Change_T2))) %>% 
    mutate(.,Change_T3 = case_when(.$Change_T3 == 1  ~  .$Change_T3*.$w,
                                   TRUE ~ as.numeric(.$Change_T3)))
  if(counter==20) {clocka <- tictoc::toc()
  counter <- 0
  remaining <- (clocka$toc - clocka$tic) * (length(species_list)/20)
  cat(red(bold("ETA (hrs)",remaining/3600,"\n")))}
}
plan(sequential)
df4 <- bind_rows(df3)

# AQUI SE ESPECIFICA EL ESCENARIO PARA HACER LAS GRAFICAS
scenario_choice = "GreenControl"
df5 <- df4 %>%  filter(.,Scenario==scenario_choice) %>% group_by(.,CellID) %>% 
  summarise(.,Global_changeT1=sum(Change_T1,na.rm=T),Global_changeT2=sum(Change_T2,na.rm=T),Global_changeT3=sum(Change_T3,na.rm=T))
df5 <- df5 %>% mutate(.,x=df4$x[match(.$CellID,df4$CellID)])
df5 <- df5 %>% mutate(.,y=df4$y[match(.$CellID,df4$CellID)])
xlims <- rich_MASTER %>% #filter(.,Country %in% country) %>% 
  dplyr::select(.,x)  %>% range(.+2)
ylims <- rich_MASTER %>% #filter(.,Country %in% country) %>% 
  dplyr::select(.,y)  %>% range(.+2)
val_lims = range(c(df5$Global_changeT1,df5$Global_changeT2,df5$Global_changeT2))
plotas <- list()
ord <- unique(rich_MASTER$Time)
ord <- ord[c(1:3)]
years <- c(2030,2050,2070)
world <- ne_coastline(scale=110,returnclass = "sf")
this <- c(2:4)
df5
for(ti in 1:length(ord)){
  df6 <- df5 %>% dplyr::select(.,x,y,this[ti])
  colnames(df6)[3] <- "Change"
  plotas[[ti]] <- ggplot() + xlim(xlims) + ylim(ylims) +
    geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="grey80",fill="grey80") +
    geom_tile(data=df6,aes(y=y, x=x,fill = Change),colour = "white",size=0.00001,height=0.1,width=0.1) + 
    scale_fill_viridis(option = "inferno",direction=1,begin=0.0,end=1,limits=val_lims) +
    geom_sf(data = world) +
    theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
          panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
          axis.text.x = element_blank(),axis.text.y = element_blank(),
          panel.background=element_rect(colour="black",fill="azure2"),panel.border=element_rect(colour="black",fill=NA)) + 
    labs(subtitle=paste(ord[ti],":",years[ti]))
}
titulo = "Projected change in species richness (unlimited dispersal)"
wrap_plots(plotas) +
  plot_layout(byrow=T,guides = "collect",tag_level = 'new') + 
  plot_annotation(title = paste(titulo,scenario[sc]),
                  caption='Ureta et al. 2020') & theme(legend.position = 'bottom')







#####
histas <- list()
y_lims <- c(0,max(signif(ceiling(max(table(c(df5$Global_changeT1)))),1),
signif(ceiling(max(table(c(df5$Global_changeT2)))),1),
signif(ceiling(max(table(c(df5$Global_changeT3)))),1)))
val_lims = range(c(df5$Global_changeT1,df5$Global_changeT2,df5$Global_changeT2))

for(ti in 1:length(ord)){
df6 <- df5 %>% dplyr::select(.,x,y,this[ti])
colnames(df6)[3] <- "Change"
histas[[ti]] <- ggplot(df6, aes(Change,fill=..x..)) + ylim(y_lims) + geom_histogram(stat="count")  +
  scale_fill_viridis(option = "inferno",direction=1,begin=0.0,end=1,limits=val_lims) +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
              panel.background=element_rect(colour="black"),panel.border=element_rect(colour="black",fill=NA)) +
  labs(subtitle=paste(ord[ti],":",years[ti]))
}
titulo = "Projected change in species richness (limited dispersal)"
wrap_plots(histas) +
  plot_layout(byrow=T,guides = "collect",tag_level = 'new') + 
  plot_annotation(title = paste(titulo,scenario[sc]),subtitle = paste(group,country,sep="-"),
                  caption='Ureta et al. 2020') & theme(legend.position = 'bottom')





#### RANGE
system.time({
  range_MASTER <- 1:length(countries) %>% future_map(function(i){
    r1 <- list_tibs %>% filter(.,Country == countries[i]) %>% group_by(Scenario,Time,Species) %>% 
      summarize(Range=n(),mean_cv=mean(CV,na.rm=T),sd_cv=sd(CV,na.rm=T)) %>% as_tibble(.) %>% 
      mutate(.,Country=countries[i],Group="Pinophyta")
     return(r1)
  },.progress = TRUE)
})
range_MASTER <- rbindlist(range_MASTER)
range_MASTER

country = "Mexico"
df2 <- range_MASTER %>% group_by(Species,Scenario,Time) %>% 
  summarize(.,n=Range)
df_pre <- df2[which(df2$Scenario=="Present"),]
df2 <- df2[which(df2$Scenario!="Present"),]
df2
df2 <- df2 %>% ungroup(.) %>%  mutate(.,Prop_range = 100*(.$n/df_pre$n[match(.$Species,df_pre$Species)])) %>% 
      group_by(Species,Scenario,Time) %>% 
      summarise(n = Prop_range)
df2 <- ungroup(df2) %>% mutate(.,Country=range_MASTER$Country[match(.$Species,range_MASTER$Species)]) %>% 
  mutate(.,Group=range_MASTER$Group[match(.$Species,range_MASTER$Species)])
df2$Scenario <- as.factor(df2$Scenario)
df2$Scenario <- factor(df2$Scenario, levels = levels(df2$Scenario)[c(5,1,2,3,4)])



#pdf("~/Desktop/Test2070_v1.pdf",height = 7,width=11)
ggplot(df2, aes(y=Species,x=Time)) + 
  geom_tile(aes(fill = n), colour = "white") + 
  scale_fill_viridis(option = "inferno",direction=1,begin=0.2,end = 0.9,na.value="grey80") + 
  theme_minimal() +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        panel.grid = element_blank(),panel.background=element_rect(fill="grey80"),
        strip.text.y=element_text(angle=0),axis.text.y = element_text(size=4)) + 
  labs(
    title="Forecast of Range size (proportional to present)",
    subtitle=paste("Endemic",group,"by country"),
    caption='Ureta et al. 2020') +
  facet_grid(Country~Scenario,shrink=FALSE,space = "free_y",scales = "free_y",labeller = label_context)
#dev.off()


#### RICHNESS











# "~/Documents/1.PROYECTOS/3.Termohalina/PAM17.6.2020/" This are the master_tibs at 1X1º resolution.
# "~/Documents/1.PROYECTOS/3.Termohalina/PAM_4/" This are the master_tibs at full resolution.
poly<-rgdal::readOGR("~/Documents/7.MAPOTECA/official/wwf_terr_ecos.shp")
g <- raster(nrows=180*12,ncols=360*12,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1,crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% as(., 'SpatialPixels')
rich_MASTER <- list()
#country=c("Australia","China","India","Colombia","Ecuador","Indonesia","Madagascar","Peru","Phillipines","Venezuela","Mexico","Brasil")
country=c("Brasil","Mexico")
for (cy in 1:length(country)){
cat(red(bold("-----","PROCESSING:",country[cy],"-----","\n")))
rich_MASTER[[country[cy]]] <- list()
scenario = c("GreenControl","Green0.5","Green1","Green1.5","Green3")
for (sc in 1:length(scenario)){
cat(yellow(bold("-+- Processing:",scenario[sc],"-+-","\n")))
rich_MASTER[[country[cy]]][[scenario[sc]]] <- list()
group =c("Mag","Pino","Moni","Lyco","Aves","Anfibios","Mamiferos","Reptiles")
#group =c("Moni","Lyco","Anfibios","Aves")
for (gr in 1:length(group)) {
time_s =c("Present",paste(scenario[sc],c("T1","T2","T3"),sep="_"))
list_arch <- list.files(pattern="MasterTib") %>%  .[grep(country[cy],.)] %>% .[grep(group[gr],.)]
if(length(list_arch)==0) {cat("No models for:",group[gr],"\n"); next}
name <- list_arch %>%  sub("MasterTib","",.) %>% sub(".Rdata","",.)
cat("Processing data for:", name,"\n")
#list_tibs <- lapply(list_arch, function(x) mget(load(x)))
list_tibs <- mget(load(list_arch))
lista <- unique(list_tibs$master_tib$name)
cat(bold(" - Number of models:",length(lista),"\n")) #### Este es el número de modelos que tenemos. 
spp <- sub("_.*","",lista) %>% unique(.)
cat(green(bold("  -- Number of species:",length(spp),"\n")))
richness <- list()
pre <- list_tibs$master_tib %>%  .[which(.$name %in% lista[grep(time_s[1],lista)]),]
length(unique(pre$name))
master_PAM <- table(pre$value)
master_PAM <- cbind(g@coords[match(names(master_PAM),g@grid.index),],master_PAM)
richness[[1]] <- as.data.frame(master_PAM[,1:3])
for (j in 2:4){
  cat(blue("  -- processing",time_s[j],"\n"))
  time_x <- list_tibs$master_tib %>%  .[which(.$name %in% lista[grep(time_s[j],lista)]),]
  if(dim(time_x)[1]==0) {cat("  ------- no data","\n");
    master_PAM <- data.frame(x=NA,y=NA,master_PAM=NA)
    rownames(master_PAM) <- "999999"
    richness[[j]] <- master_PAM[,1:3]
    next} else { 
  if(length(table(time_x$value))==0) {cat("  ------- no data","\n");
    master_PAM <- data.frame(x=NA,y=NA,master_PAM=NA)
    rownames(master_PAM) <- "999999"
    richness[[j]] <- master_PAM[,1:3]
    next} else { 
  master_PAM <- table(time_x$value)
  if(length(master_PAM)==1){
  master_PAM <- c(g@coords[match(names(master_PAM),g@grid.index),],master_PAM)
  richness[[j]] <- as.data.frame(t(master_PAM[1:3]));cat("  ------- one datum!","\n");
  rownames(richness[[j]]) <- colnames( richness[[j]]) [3]
  colnames( richness[[j]]) [3] <- "master_PAM"
  next} else { 
  if(length(master_PAM) > 1) {master_PAM <- cbind(g@coords[match(names(master_PAM),g@grid.index),],master_PAM)
  richness[[j]] <- as.data.frame(master_PAM[,1:3])}
    }}}
  
  }
for (ll in 1:length(richness)){
  richness[[ll]] <- cbind(rownames(richness[[ll]]),richness[[ll]])
  colnames(richness[[ll]])[1] <- "CellID"
}

cat(cyan("  -- merging data",scenario[sc],"\n"))
lapply(richness,dim)
richness.total <- Reduce(function(...) merge(...,by="CellID",all.x=TRUE,all.y=TRUE),richness)
names(richness.total) <- c("CellID","x","y","richness_pre","x.y","y.y","richness_t1","x.x","y.x","richness_t2","x.y","y.y","richness_t3")
richness.total[which(is.na(richness.total[,2])),2] <- richness.total[which(is.na(richness.total[,2])),5]
richness.total[which(is.na(richness.total[,2])),2] <- richness.total[which(is.na(richness.total[,2])),8]
richness.total[which(is.na(richness.total[,2])),2] <- richness.total[which(is.na(richness.total[,2])),11]
richness.total[which(is.na(richness.total[,3])),3] <- richness.total[which(is.na(richness.total[,3])),6]
richness.total[which(is.na(richness.total[,3])),3] <- richness.total[which(is.na(richness.total[,3])),9]
richness.total[which(is.na(richness.total[,3])),3] <- richness.total[which(is.na(richness.total[,3])),12]

richness.total <-  richness.total[,-c(5,6,8,9,11,12)]
richness.total$richness_pre[which(is.na(richness.total$richness_pre))] <- 0
richness.total$richness_t1[which(is.na(richness.total$richness_t1))] <- 0
richness.total$richness_t2[which(is.na(richness.total$richness_t2))] <- 0
richness.total$richness_t3[which(is.na(richness.total$richness_t3))] <- 0
iden <- which(richness.total$CellID=="999999")
if(length(iden)!=0) richness.total <- richness.total[-iden,]

pointos <- SpatialPoints(richness.total[,2:3])
proj4string(pointos)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
points_eco <- sp::over(pointos,poly)
points_eco$BIOME <- as.factor(points_eco$BIOME)
richness.total <- cbind(richness.total,points_eco$BIOME)
colnames(richness.total)[8] <- "Biome"
richness.total$Biome <- as.character(richness.total$Biome)
richness.total$change_t1 <- richness.total$richness_t1 - richness.total$richness_pre
richness.total$change_t2 <- richness.total$richness_t2 - richness.total$richness_pre
richness.total$change_t3 <- richness.total$richness_t3 - richness.total$richness_pre
richness.total <- cbind(Scenario=rep(scenario[sc],dim(richness.total)[1]), Country=rep(country[cy],dim(richness.total)[1]),Group=rep(group[gr],dim(richness.total)[1]),richness.total)
rich_MASTER[[country[cy]]][[scenario[sc]]][[group[gr]]] <- richness.total
}
rich_MASTER[[country[cy]]][[scenario[sc]]] <- data.table::rbindlist(rich_MASTER[[country[cy]]][[scenario[sc]]])
cat(yellow("-+- Done:","-+-","\n"))
}
rich_MASTER[[country[cy]]] <- data.table::rbindlist(rich_MASTER[[country[cy]]])
cat(red $ bgWhite $ bold("SUPER DONE !!!","\n"))
}
str(rich_MASTER, give.attr=FALSE, give.length=FALSE, give.head=FALSE, vec.len=0, max.level = 2,
    indent.str="|", comp.str="----")
#save(rich_MASTER,file="richness_MASTER.Rdata")
#load("richness_MASTER.Rdata")
r_MASTER <- as.data.frame(data.table::rbindlist(rich_MASTER))
table(r_MASTER$Group)
table(r_MASTER$Scenario)
dim(r_MASTER)
target_group = "all"
target_scenario = "all"
if(target_scenario!="all") r_MASTER <- (r_MASTER[which(r_MASTER$Scenario==target_scenario & r_MASTER$Group == target_group),])
if(target_scenario=="all") {
if(target_group!="all") r_MASTER <- (r_MASTER[which(r_MASTER$Group == target_group),])
if(target_group=="all") 
r_MASTER <- r_MASTER
g <- raster(nrows=180*12,ncols=360*12,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1,crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
lon.lim <- range(r_MASTER$x) + c(-5,5)
lat.lim <- range(r_MASTER$y) + c(-5,5)
gr <- crop(g,c(lon.lim,lat.lim))
r_test <- r_MASTER
r_test$Biome <- as.factor(r_test$Biome)
levels(r_test$Biome) <- c("TropMoistForest","MediterraneanForests","Deserts","Mangroves","TropDryForests",
                          "TropConiferForests","TempConiferForests","TropGrasslands","FloodedGrasslands")
r_test$Biome <- forcats::fct_explicit_na(r_test$Biome,"Undet. Biome")
head(r_test)
df2 <- r_test %>% .[which(.$change_t1<0 & .$Country=="Brasil"),]
df2 <- df2 %>% 
  group_by(Group,Scenario,Biome) %>% 
  summarise(n= abs(mean(change_t1 / richness_pre)*100))
#pdf("~/Desktop/Test2070_v1.pdf",height = 7,width=11)
ggplot(df2, aes(y=Scenario,x=Group)) + 
  geom_tile(aes(fill = n), colour = "white") + 
  scale_fill_viridis(option = "viridis",direction=1,begin=0.2,end=0.95) + 
  theme_minimal() +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        panel.grid = element_blank()) + 
  labs(
    title="Forecast of Species Loss (mean %)",
    subtitle=paste(country,"-2070"),
    caption='Ureta et al. 2020') + 
facet_wrap(~ Biome,ncol = 3,shrink=TRUE)
#dev.off()


}


as.character(unique(r_MASTER$Country))
unique(r_MASTER$Scenario)
table(r_MASTER$Group)
head(r_MASTER)
which(is.na(r_MASTER[,7]))
g <- raster(nrows=180*12,ncols=360*12,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1,crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
lon.lim <- range(r_MASTER$x) + c(-5,5)
lat.lim <- range(r_MASTER$y) + c(-5,5)
gr <- crop(g,c(lon.lim,lat.lim))
r <- rasterize(r_MASTER[,5:6],gr,field=r_MASTER$richness_t1)
r@data@values[which(r@data@values==0)] <- NA
wld <- maps::map('world', xlim=lon.lim,
           ylim=lat.lim,plot=FALSE)
wld <- data.frame(lon=wld$x, lat=wld$y)
col.l <- rev(magma(30))

rasterVis::levelplot(r,col.regions=col.l,main=paste("Species richness for",target_group)) +
  xyplot(lat ~ lon, wld, type='l', lty=1, lwd=0.5, col='black')


r <- rasterize(r_MASTER[,5:6],gr,field=r_MASTER$change_t1)
col.l <- (viridis(30))
rasterVis::levelplot(r,zscaleLog=FALSE,col.regions=col.l,margin=list(draw=T,FUN=mean,
              axis=T),pretty=T,
                main=paste("Mean change in richness for",target_group,"in",target_scenario,"(time t1)")) +
  xyplot(lat ~ lon, wld, type='l', lty=1, lwd=0.5, col='black')

dim(r_MASTER)


richness <- list_tibs$value %>% table(.) %>% enframe(.)




# plotas[[sc]] <-  ggplot(data = world) + 
#   xlim(xlims) + ylim(ylims) +
#   geom_sf() +
#   #borders("world",xlim = xlims,ylim=ylims,colour = "gray80",fill="gray80") +
#   geom_tile(data=temp,aes(y=y, x=x,fill = SR), colour = "white",size=0.00001,height=0.1,width=0.1) + 
#   scale_fill_viridis(option = "inferno",direction=1,begin=0.2,end=0.95,limits=val_lims) +
#   theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
#         panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
#         axis.text.x = element_blank(),axis.text.y = element_blank(),
#         panel.background=element_rect(colour="black")) + 
#   labs(
#     #title=paste("Species Richness (",ord[1],")",sep=""),
#     subtitle=paste(ord[sc])




##### ##### 



load("~/Downloads/TRAQUEOS_Corrected_overlay.r")
levels(as_factor(p_clean$PAÍS)) -> paises
p_clean$class <- (as_factor(p_clean$class))
levels(p_clean$class) <- c("Pinophyta","Magnoliophyta","Monilophyta",
                           "Magnoliophyta","Monilophyta","Pinophyta",
                           "Lycophyta","Monilophyta","Pinophyta",
                           "Monilophyta","Pinophyta")

master <- list()
for (y in 1:nlevels(p_clean$class)){
  p_filter <- p_clean %>% filter(.,class==levels(p_clean$class)[y])
  spp <- paises %>% map_int(function(x){ 
    r1 <- p_filter %>% filter(.,PAÍS==x) %>% dplyr::select(., Resolved_ACCEPTED) %>% unique(.) %>% dim(.) %>% .[1]
    return(r1)})
  master[[y]] <- spp
}
master
names(master) <- levels(p_clean$class)
spp <- do.call(cbind,master) %>% as_tibble(.)
spp <- cbind(paises,spp)
spp
spp[,1] <- c("Australia","Brazil","China","Colombia","Ecuador","Phillipines",
             "India","Indonesia","Madagascar","Mexico","Peru","Venezuela" )
che <- read.table("../Checklist.txt",sep="\t",header=T) %>% as_tibble()
che
spp <- spp %>% .[match(che$Country,.$paises),] %>% dplyr::select(.,-paises) %>% bind_cols(che,.) %>% dplyr::select(.,-Plants)
spp
write.table(spp,file="SpeciesbyGroup.txt",sep=",",col.names = T,row.names = F,quote = F)


false.pos <- c()
sample.size <- c()
tt=10
ns <- c(.2,.15,.1,0.05)
for (j in 1:length(ns)){ 
pis <- c()
r.sq <- c()
t=seq(0.1,1,by=ns[j])*tt
n=length(t)
for (i in 1:1000){ 
  cat(i,"\r")
x <- sample(1:1000,n)
y <- sample(1:1000,n)
xx <- x/(t)
yy <- y/(t)
yy <- log(yy)
summ.s <- summary(lm(yy~xx))
pis <- c(pis,summ.s$coefficients[2,4])
r.sq <- c(r.sq,summ.s$r.squared)
}
#hist(pis,breaks=seq(0,1,0.05))
#hist(r.sq,breaks=seq(0,1,0.05))
false.pos <- c(false.pos,length(which(pis<0.05)))
sample.size <- c(sample.size,length(t))
}

plot(sample.size,false.pos)
title(paste("Time = ", tt))


r %>% levelplot(contour=TRUE, par.settings=theonlytwocolors)






#####################################################################
# R script use for the GLMM mini-Workshop on 11th March in Freising #
#       Author: Lionel Hertzog and Franzi Korner-Nievergelt         #
#####################################################################


#load the libraries
library(lme4)
library(nlme)
library(arm)
library(RCurl) #to directly download the rikz dataset

########################
# part 0 fitting GLMM #
#  #  #  #  #  #  # 

load("~/Downloads/MexicoSumTib.Rdata")
summary_tib
summary_tib %>% ungroup() %>% select(Scenario) %>% unique()
summary_tib %>% ungroup() %>% filter(Time=="") %>% 
  mutate(ID_complete = paste(CellID,Group,sep="_")) -> present
  
summary_tib %>% ungroup() %>% filter(Scenario=="GreenControl") %>% filter(Time=="T1") %>% 
  mutate(ID_complete = paste(CellID,Group,sep="_")) -> GreenCtr_T1

present %>% group_by(CellID) %>% summarize(Total_SR=sum(SR)) -> present
GreenCtr_T1 %>% group_by(CellID) %>% summarize(Total_SR=sum(SR)) -> GreenCtr_T1

merge(present,GreenCtr_T1,by="CellID",all=T) %>% as_tibble() %>% 
  mutate_if(is.numeric, funs(replace_na(., 0))) -> merged

#merged$Group.x[which(is.na(merged$Group.x))] <- merged$Group.y[which(is.na(merged$Group.x))]
#colores <- viridis::inferno(8)
#merged %>% mutate(Color=as_factor(Group.x)) -> merged
#levels(merged$Color) <- colores
#merged$Color <- as.character(merged$Color)
#id <- sample(dim(merged)[1],50000)
plot(merged$Total_SR.x,merged$Total_SR.y,xlab="Species Richness (present)",
     ylab="Species Richness (Green Control: T1)",lwd=0.2,bg=rgb (0,0,0,0.3),pch=21)
plot(merged$Total_SR.y - merged$Total_SR.x,lwd=0.2,bg=rgb (0,0,0,0.3),pch=21)
plot(density(merged$Total_SR.y,bw=10,window="cosine"))

#first a random intercept model
#plot the fit of the model
pal<-colores
f <- as.factor(merged$Group.x[id])
x=merged$SR.x[id]
y=merged$SR.y[id]
 idd <- which(y==0 | x==0)
 x=x[-idd]
 y=y[-idd]
 f=f[-idd]
x <- tapply(x,f,range_01) %>% unlist()
y <- tapply(y,f,range_01) %>% unlist()
mod_lme1<-lme(y~x,random=~1|f,control = lmeControl(opt = "optim"))
plot(x,y,col=pal[f],pch=16,main="Linear Model")
for(i in 1:length(levels(f))){
  lines(x,fixef(mod_lme1)[1]+fixef(mod_lme1)[2]*x+ranef(mod_lme1)[i,1],col=pal[i],lwd=1.5) 
}

mod_lme2<-lme(y~x,random = ~ x|f,control = lmeControl(opt = "optim"))
plot(y~x,col=pal[f],pch=16,main="Linear Mixed Model")
for(i in 1:length(levels(f))){
  lines(x,fixef(mod_lme2)[1]+(fixef(mod_lme2)[2]+ranef(mod_lme2)[i,2])*x+ranef(mod_lme2)[i,1],col=pal[i],lwd=1.5) 
}



#then a random slope plus intercept model
mod_lme2<-lme(merged$SR.y~merged$SR.x,random=merged$SR.x|merged$Group.x)
mod_lmer2<-lmer(SR.y~SR.x+(SR.x|Group.x),data=merged)

mod_lme2<-lme(data=data,Richness~NAP,random=NAP|Beach)
mod_lmer2<-lmer(Richness~NAP+(NAP|Beach),data=data)

#Poisson model
mod_glmer1<-glmer(Richness~NAP+(1|Beach),data=data,family="poisson")
#nested and crossed random effect??

##################################
#   part 1 mixed vs fixed effect #
#   #   #   #   #   #   #   #
#factor variable with intercept only effect
#simulate data in a fixed effect way
x<-rnorm(120,0,1)
f<-gl(n=6,k=20)
modmat<-model.matrix(~x+f,data.frame(x=x,f=f))
betas<-c(1,2,0.3,-3,0.5,1.2,0.8)
y<-modmat%*%betas+rnorm(120,0,1)

#the fixed effect model
m_lm<-lm(y~x+f)
#the mixed effect model
m_lme<-lme(y~x,random=~1|f)

#checking model assumptions in both case
par(mfrow=c(2,2))
plot(m_lm)

plot(fitted(m_lme),resid(m_lme))
qqnorm(resid(m_lme))
qqnorm(ranef(m_lme)[,1])
plot(x,resid(m_lme))

#looking at the fitted parameters
summary(m_lm)
summary(m_lme)

#plot the fit of the model
par(mfrow=c(1,1))
library(RColorBrewer)
pal<-brewer.pal(6,"Set1")
plot(y~x,col=pal[f],pch=16,main="Linear Model")
for(i in 1:length(levels(f))){
  if(i==1){
    lines(x,coef(m_lm)[1]+coef(m_lm)[2]*x,col=pal[i],lwd=1.5)
  }
  else{
    lines(x,coef(m_lm)[1]+coef(m_lm)[2]*x+coef(m_lm)[i+1],col=pal[i],lwd=1.5)
  }
}

plot(y~x,col=pal[f],pch=16,main="Linear Mixed Model")
for(i in 1:length(levels(f))){
  lines(x,fixef(m_lme)[1]+fixef(m_lme)[2]*x+ranef(m_lme)[i,1],col=pal[i],lwd=1.5) 
}
#no clear difference visible

#now generqte a random slope/intercept model through the mixed effect technique
rnd.eff<-rnorm(5,0,2)
slp.eff<-rnorm(5,0,1)
all.eff<-c(1,2,rnd.eff,slp.eff)
modmat<-model.matrix(~x*f,data.frame(x=x,f=f))
y<-modmat%*%all.eff+rnorm(120,0,1)

#build the two model
m_lm2<-lm(y~x*f)
m_lme2<-lme(y~x,random=~x|f)

#checking model assumptions
par(mfrow=c(2,2))
plot(m_lm2)
plot(fitted(m_lme2),resid(m_lme2))
abline(h=0,lty=2,col="red")
qqnorm(resid(m_lme2))
qqnorm(ranef(m_lme2)[,1])
qqnorm(ranef(m_lme2)[,2])

#summary of the models
summary(m_lm2)
summary(m_lme2)

#plot the model fitted values
par(mfrow=c(1,2))
plot(y~x,col=pal[f],pch=16,main="Linear Model")
for(i in 1:length(levels(f))){
  if(i==1){
    lines(x,coef(m_lm2)[1]+coef(m_lm2)[2]*x,col=pal[i],lwd=1.5)
  }
  else{
    lines(x,coef(m_lm2)[1]+(coef(m_lm2)[2]+coef(m_lm2)[i+6])*x+coef(m_lm2)[i+1],col=pal[i],lwd=1.5)
  }
}

plot(y~x,col=pal[f],pch=16,main="Linear Mixed Model")
for(i in 1:length(levels(f))){
  lines(x,fixef(m_lme2)[1]+(fixef(m_lme2)[2]+ranef(m_lme2)[i,2])*x+ranef(m_lme2)[i,1],col=pal[i],lwd=1.5) 
}

#again no clear difference can be seen ...

#conclusion
#end of Practical 1

#######################
#   Practical 2      #
#   #   #   #   #  

#checking model assumptions
par(mfrow=c(2,2))
plot(fitted(m_lme2),resid(m_lme2))
abline(h=0,lty=2,col="red")
qqnorm(resid(m_lme2))
qqline(resid(m_lme2))
qqnorm(ranef(m_lme2)[,1])
qqline(ranef(m_lme2)[,1])
qqnorm(ranef(m_lme2)[,2])
qqline(ranef(m_lme2)[,2])
scatter.smooth(fitted(m_lme2),sqrt(abs(resid(m_lme2))))

#wrong data
modmat[,2]<-log(modmat[,2]+10)
y<-modmat%*%all.eff+runif(120,0,5)
m_wrg<-lme(y~x,random=~x|f)

plot(fitted(m_wrg),resid(m_wrg))
abline(h=0,lty=2,col="red")
qqnorm(resid(m_wrg))
qqline(resid(m_wrg))
qqnorm(ranef(m_wrg)[,1])
qqline(ranef(m_wrg)[,1])
qqnorm(ranef(m_wrg)[,2])
qqline(ranef(m_wrg)[,2])
scatter.smooth(fitted(m_wrg),sqrt(abs(resid(m_wrg))))

#plot fitted values vs resid, qqnorm the residuals and all random effects

#end of practical 2

###################
#  Practical 3   #
#  #  #  #  #  #

#Model selection
#work with the RIKZ dataset from Zuur et al

data<-read.table("/home/lionel/Documents/PhD/GLMM_WS/data/rikz.txt",sep=" ",head=TRUE)

#testing the random effect
#a first model
mod1<-lme(Richness~NAP+Exposure,data=data,random=~1|Beach,method="REML")
#a second model without the random term, gls is used because it also fit the model through REML
mod2<-gls(Richness~NAP+Exposure,data=data,method="REML")
#likelihood ratio test, not very precise for low sample size
anova(mod1,mod2)

# parameteric bootstrap
lrt.obs <- anova(mod1, mod2)$L.Ratio[2] # save the observed likelihood ratio test statistic
n.sim <- 1000  # use 1000 for a real data analysis
lrt.sim <- numeric(n.sim)
dattemp <- data
for(i in 1:n.sim){
  dattemp$ysim <- simulate(lm(Richness ~ NAP+Exposure, data=dattemp))$sim_1 # simulate new observations from the null-model
  modnullsim <- gls(ysim ~ NAP+Exposure, data=dattemp)   # fit the null-model
  modaltsim <-lme(ysim ~ NAP+Exposure, random=~1|Beach, data=dattemp)  # fit the alternative model
  lrt.sim[i] <- anova(modnullsim, modaltsim)$L.Ratio[2] # save the likelihood ratio test statistic
}

(sum(lrt.sim>=lrt.obs)+1)/(n.sim+1)  # p-value

#plot
par(mfrow=c(1,1))
hist(lrt.sim, xlim=c(0, max(c(lrt.sim, lrt.obs))), col="blue", xlab="likelihood ratio test statistic", ylab="density", cex.lab=1.5, cex.axis=1.2)
abline(v=lrt.obs, col="orange", lwd=3)

#model selection for the fixed effect part, to use the simulate function we need MER object
mod1_ML<-lme(Richness~NAP+Exposure,data,random=~1|Beach,method="ML")
mod3<-lme(Richness~NAP,data,random=~1|Beach,method="ML")
mod1_lmer<-lmer(Richness~NAP+Exposure+(1|Beach),data=data,REML=FALSE)
mod3_lmer<-lmer(Richness~NAP+(1|Beach),data=data,REML=FALSE)
#compare with lme results
summary(mod1_lmer)
summary(mod1_ML)
#anova
anova(mod1_lmer,mod3_lmer)

#again parametric boostrapping of the LRT
lrt.obs<-anova(mod1_lmer, mod3_lmer)$Chisq[2]
n.sim <- 1000  # use 1000 for a real data analysis
lrt.sim <- numeric(n.sim)
dattemp <- data
for(i in 1:n.sim){
  dattemp$ysim <-  unlist(simulate(mod3_lmer)) # simulate new observations from the null-model
  modnullsim <- lmer(ysim ~ NAP+(1|Beach), data=dattemp,REML=FALSE)   # fit the null-model
  modaltsim <-lmer(ysim ~ NAP+Exposure+(1|Beach), data=dattemp,REML=FALSE)  # fit the alternative model
  lrt.sim[i] <- anova(modnullsim, modaltsim)$Chisq[2] # save the likelihood ratio test statistic
}

(sum(lrt.sim>=lrt.obs)+1)/(n.sim+1)  # p-value

#plot
hist(lrt.sim, xlim=c(0, max(c(lrt.sim, lrt.obs))), col="blue", xlab="likelihood ratio test statistic", ylab="density", cex.lab=1.5, cex.axis=1.2)
abline(v=lrt.obs, col="orange", lwd=3)

#the next step would be to drop NAP first and then see if the likelihood ratio test is significant and if dropping Exposure first always
#lead to higher LRT statistic
#other methods, AIC..
#R square computation for GLMM, see supplementary material from Nakagawa 2013 MEE
VarF <- var(as.vector(fixef(mod1_lmer) %*% t(mod1_lmer@pp$X)))
# VarCorr() extracts variance components
# attr(VarCorr(lmer.model),’sc’)^2 extracts the residual variance, VarCorr()$plot extract the variance of the plot effect
VarF/(VarF + VarCorr(mod1_lmer)$Beach[1] + attr(VarCorr(mod1_lmer), "sc")^2 )

#compute the conditionnal R-square
(VarF + VarCorr(mod1_lmer)$Beach[1])/(VarF + VarCorr(mod1_lmer)$Beach[1] + (attr(VarCorr(mod1_lmer), "sc")^2))

#end of practical 3


######################
#    Practical 4    #
#  #  #  #  #  #  #

#drawing inference from a model
#p-values can be retrieved from lme and glmer but not from lmer call
summary(mod1)
summary(mod1_lmer)

mod1_glmer<-glmer(Richness~NAP+Exposure+(1|Beach),data=data,family="poisson")
summary(mod1_glmer)

#using sim from the arm package
n.sim<-1000
simu<-sim(mod1_glmer,n.sims=n.sim)
head(simu@fixef)
#95% credible interval
apply(simu@fixef,2,quantile,prob=c(0.025,0.5,0.975))
#plotting the effect of NAP on the richness
nsim <- 1000
bsim <- sim(mod1_glmer, n.sim=nsim)
newdat <- data.frame(NAP=seq(-1.5, 2.5, length=100),Exposure=mean(data$Exposure))
Xmat <- model.matrix(~NAP+Exposure, data=newdat)
predmat <- matrix(ncol=nsim, nrow=nrow(newdat))
predmat<-apply(bsim@fixef,1,function(x) exp(Xmat%*%x))
newdat$lower <- apply(predmat, 1, quantile, prob=0.025)
newdat$upper <- apply(predmat, 1, quantile, prob=0.975)
newdat$med<-apply(predmat, 1, quantile, prob=0.5)

plot(Richness~NAP, data=data, pch=16, las=1, cex.lab=1.4, cex.axis=1.2)
lines(newdat$NAP,newdat$med,col="blue",lty=1,lwd=1.5)
lines(newdat$NAP,newdat$upper,col="red",lty=2,lwd=1.5)
lines(newdat$NAP,newdat$lower,col="red",lty=2,lwd=1.5)

#to compare the coefficient between the different terms standardize the variable
data$stdNAP<-scale(data$NAP)
data$stdExposure<-scale(data$Exposure)
mod2_glmer<-glmer(Richness~stdNAP+stdExposure+(1|Beach),data=data,family="poisson")

#simulate to draw the posterior distribution of the coefficients
n.sim<-1000
simu<-sim(mod2_glmer,n.sims=n.sim)
head(simu@fixef)
#95% credible interval
coeff<-t(apply(simu@fixef,2,quantile,prob=c(0.025,0.5,0.975)))
#plot
plot(1:3,coeff[,2],ylim=c(-0.8,2),xaxt="n",xlab="",ylab="Estimated values")
axis(side=1,at=1:3,labels=attr(fixef(mod2_glmer),"names"))
segments(x0=1:3,y0=coeff[,1],x1=1:3,y1=coeff[,3],lend=1)
abline(h=0,lty=2,col="red")

#end of practical 4






##### Format data to create circular barplot by group ####
setwd("~/Documents/1.PROYECTOS/1.Termohalina/Final_v4/TIBS/")
list_arch <- list.files(pattern="RangeTib",recursive = F)
list_tibs <- lapply(list_arch, function(x) mget(load(x))) %>% lapply(function(x)x[[1]]) %>% bind_rows() %>% 
  mutate_if(is.numeric , replace_na, replace = 0) %>% rename("Present"=Present_) -> range_tib
range_tib %>% mutate_if(is.numeric , replace_na, replace = 0) %>% mutate_at(5:ncol(.), ~((.x - Present)/Present)) %>%  mutate(Present = 0) -> STDrange_tib
STDrange_tib %>% group_by(Country,Group) %>% summarise(Species=n()) %>% arrange(Group) -> dataTib
data2
x <- unlist(lapply(strsplit(data$code_2,"_"),"[[",1))
y <- unlist(lapply(strsplit(data$code_2,"_"),"[[",2))
data$group <- as.character(data$group)
data[which(is.na(data$group)),2] <- y[which(is.na(data$group))]
data[which(is.na(data$country)),1] <- x[which(is.na(data$country))]
which(is.na(data$country))
data.sum <- data %>% group_by(group) %>% summarise(Modeladas=sum(spp,na.rm=T)) %>% 
  mutate(Linaje=c("animal","animal","planta","planta","animal","planta","planta","animal")) %>% 
  mutate(Total=c(8181,11147,1340,291000,6495,11916,1058,11136)) %>% 
  mutate(Proporcion=Modeladas/Total * 100)

data.sum %>% group_by(Linaje) %>% summarise(Total_Mod=sum(Modeladas,na.rm=T),Total_Tot=sum(Total,na.rm=T)) %>% 
  mutate(Total_Mod/Total_Tot * 100)


dataTib

filling = as.factor(dataTib$Group)
empty_bar <- 3
data <- data.frame( matrix(NA,empty_bar*nlevels(filling), ncol(dataTib)) ) %>% 
  'colnames<-' (colnames(dataTib)) %>% mutate(.,Group=rep(levels(filling), each=empty_bar)) %>% 
  rbind(dataTib, .) %>% arrange(Group) %>% ungroup() %>%  mutate(.,id=seq(1, nrow(.))) %>% 
  mutate(spp=log(Species)+.1)
data[which(!is.finite(data$Species)&!is.na(data$Species)),"Species"] <- NA
factor(data$Group,levels = c("Anfibios", "Aves", "Mamiferos", "Reptiles", "Lycophyta", "Magnoliophyta", "Monilophyta", "Pinophyta")) -> data$Group
data %>% arrange(Group) %>%  mutate(.,id=seq(1, nrow(.))) -> data

#mutate(spp=(spp/total)*100)
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
base_data <- data %>% 
  group_by(Group) %>%  ##### here!!
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]
grouping = data$Country
y.limits <- floor(range(data$spp,na.rm = T)[2] + 3)
ys <- seq(2,12,2)
p <- ggplot(data, aes(x=as.factor(id), y=spp, fill=grouping)) +  ### here!!
  geom_bar(aes(x=as.factor(id), y=spp, fill=grouping), ### here!!
           stat="identity", alpha=0.5) +
  scale_fill_viridis_d(option="D") +
  geom_segment(data=grid_data, aes(x = end, y = ys[length(ys)-3], xend = start, yend = ys[length(ys)-3]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = ys[length(ys)-4], xend = start, yend = ys[length(ys)-4]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = ys[length(ys)-5], xend = start, yend = ys[length(ys)-5]), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  annotate("text", x = rep(max(data$id),length(ys)-2), y = ys[1:4], 
           label = as.character(ys[1:4]), color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  geom_bar(aes(x=as.factor(id), y=spp, fill=grouping), stat="identity", alpha=0.5) + ### here!!
  ylim(-10,12) + theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")) + coord_polar() +
  geom_text(data=label_data, aes(x=id, y = spp + 0.2, label=Country, hjust=hjust), 
            color="black", fontface="bold",alpha=0.6, size=2, 
            angle= label_data$angle, inherit.aes = FALSE ) +
  geom_segment(data=base_data, aes(x = start, y = -0.5, xend = end, yend = -.5), 
               colour = "black", alpha=0.8, size=0.6, inherit.aes = FALSE )
pdf("Figure_1.pdf",useDingbats = F)
p
dev.off()


#### Get validation statistics
library(tidyverse)
list_arch <- list.files(pattern="myBiomodEMEval.csv",full.names = F,recursive = T)
list_tss <- lapply(list_arch, function(x) read.csv(x) %>% .[,17]) %>% bind_cols() %>% t(.) %>% as_tibble() %>% rename_with(~c("Kappa","TSS","ROC"))
list_arch %>% strsplit(.,split="/") -> nombres
list_tss <- list_tss %>% mutate(Group=nombres %>% lapply(.,"[",1) %>% unlist(),.before=Kappa) %>% mutate(Country=nombres %>% lapply(.,"[",2) %>% unlist(),.before=Kappa) %>% mutate(Species= nombres %>% lapply(.,"[",3) %>% unlist(),.before=Kappa)
save(list_tss,file="ValidationSummary.Rdata")

list_arch <- list.files(pattern="myVarImportEM",full.names = F,recursive = T)
list_vars <- lapply(list_arch, function(x) read.csv(x) %>% as_tibble() %>% rename_with(~c("Variable","Importance")))
list_arch %>% strsplit(.,split="/") -> nombres
tibble(Group=nombres %>% lapply(.,"[",1) %>% unlist(), Country=nombres %>% lapply(.,"[",2) %>% unlist(),Species=nombres %>% lapply(.,"[",3) %>% unlist()) %>% mutate(Bios=list_vars) -> list_vars
save(list_vars,file="VariablesSummary.Rdata")



##### Process validation statistics
load("~/Desktop/TERMOS_MS/Validation/ValidationSummary_Animals.Rdata")
list_tss -> animals
list_arch <-  list.files(path="~/Desktop/TERMOS_MS/Validation/Plants",full.names = T)
list_tss <- lapply(list_arch, function(x) mget(load(x)))
for (i in 1:length(list_tss)){ 
list_tss[[i]]$list_tss %>% select(-Species) %>% rename_with(.cols=1:2,~c("Country","Species")) %>% mutate(Group=NA) -> list_tss[[i]]}
list_tss %>% do.call(rbind,.) -> plantas

load("~/Desktop/TERMOS_MS/Validation/ValidationSummary.Rdata")
list_tss %>% bind_rows(.,plantas,animals) %>% distinct(Species,.keep_all = T) -> list_tss

list_tss %>% select(-Group,-Country) -> list_tss
STDrange_tib %>% distinct(Species,.keep_all = T) %>% left_join(.,list_tss,by="Species") %>% select(1:3,Kappa,ROC,TSS) %>% group_by(Group) %>% summarize(N=n())
STDrange_tib %>% ungroup() %>% distinct(Species,.keep_all = T) %>% left_join(.,list_tss,by="Species") %>% select(1:3,Kappa,ROC,TSS) %>% mutate(Group=ifelse(Group %in% c("Anfibios","Aves","Mamiferos","Reptiles"),"Tetrapods", "Vascular Plants")) -> stats
stats %>% group_by(Group) %>% filter(!is.na(TSS)) %>%
  ggplot(aes(x=TSS,fill=Group)) +
  geom_histogram(color = "black",position = "identity")+ facet_wrap(facets = "Group")+
  geom_hline(yintercept = 0)+theme(panel.grid = element_blank(),legend.position="none",panel.background = element_blank(),axis.line = element_line()) + scale_x_continuous(breaks = pretty(stats$TSS),labels = abs(pretty(stats$TSS)))  + geom_vline(xintercept = 0.7) +
    NULL
  
  

load("~/Desktop/TERMOS_MS/Validation/VariablesSummary.Rdata")
list_vars -> plantas
load("~/Desktop/TERMOS_MS/Validation/VariablesSummary_Animals.Rdata")
list_vars %>% bind_rows(.,plantas) -> list_vars
MostImp <- lapply(list_vars$Bios, function(x) x %>% arrange(desc(Importance)) %>% .[1,1] %>%  mutate_if(is.factor,as.character) %>% pull())

list_vars %>% mutate(MostImp = (MostImp)) %>% select(-Bios) %>% unnest(MostImp) %>% select(Species,MostImp) %>% inner_join(STDrange_tib %>% ungroup(),.,by="Species") %>% select(1:3,MostImp,ends_with("_T2")) %>%  pivot_longer(cols=ends_with("T2"),names_to = "Scenario",values_to = "Range") %>% filter(Scenario%in%c("GreenControl_T2","Green0.5_T2"),Range <= - 0.95) %>% mutate(Group=case_when(Group == "Anfibios"|Group =="Aves"|Group =="Mamiferos"|Group =="Reptiles" ~ "Tetrapods", Group=="Lycophyta"|Group =="Magnoliophyta"|Group =="Monilophyta"|Group =="Pinophyta" ~ "Vascular Plants")) %>% group_by(MostImp,Scenario) %>% summarise(Freq=n()) %>% ungroup()  %>% pivot_wider(names_from = Scenario,values_from = Freq) %>% mutate(Green0.5_T2=Green0.5_T2/sum(Green0.5_T2),GreenControl_T2=GreenControl_T2/sum(GreenControl_T2)) %>% pivot_longer(cols=2:3,names_to = "Group",values_to = "Freq") %>% mutate(Freq= ifelse(Group =="Green0.5_T2",Freq,Freq*-1)) -> variables 
variables %>% ggplot(aes(x=MostImp,y=Freq,fill=Group)) + 
  geom_bar(stat="identity",position="identity") + 
  coord_flip()+
  theme(legend.position="bottom",axis.text.x = element_text(angle = 90),panel.background = element_blank(),axis.line = element_line(),axis.text.x.bottom  = element_text(angle=0),panel.grid.major.x = element_line(colour = "grey80")) + scale_y_continuous(limits = c(-0.35,0.35),breaks = pretty(variables$Freq),labels = abs(pretty(variables$Freq))) + labs(title="Most important variables in SDMs",subtitle="for species with range loss > 95%") +
  NULL











#### Alluvial plots
library(ggalluvial)
list_vars %>% mutate(MostImp = unlist(MostImp)) %>% select(Species,Bios,MostImp) %>% inner_join(STDrange_tib %>% ungroup(),.,by="Species") %>% select(1:3,MostImp,ends_with("_T1")) %>%  pivot_longer(cols=ends_with("T1"),names_to = "Scenario",values_to = "Range") %>% group_by(MostImp,Country) %>% summarise(Freq=n()) %>% 
  ggplot(aes( axis1=Country, axis2=MostImp, y = Freq)) +
  scale_x_discrete(limits = c("Country","MostImp"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = MostImp)) +
  geom_stratum() +
  scale_fill_viridis_d(option="B",direction=-1,begin = 0.4,end=0.9) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void() + theme(legend.position = "none") + labs(title="Range Loss at 2070")
NULL


pdf("~/Documents/1.PROYECTOS/1.Termohalina/TERMOS_MS/FINAL_v2/EDIT/Figure_SAlluvials.pdf")
STDrange_tib %>% ungroup() %>% #distinct(Species,.keep_all = T) %>% 
pivot_longer(cols=4:19,names_to = "Scenario",values_to = "Loss") %>% mutate(LossCat=case_when(Loss >= 0 ~ "Gain",Loss < 0 & Loss > -0.45 ~"Moderate Loss",Loss <= -0.45 & Loss > -0.8 ~ "Severe Loss",Loss <= -0.8 & Loss > -0.99 ~"Extreme Loss",Loss == -1 ~ "Complete Loss")) %>% filter(Scenario!="Present") %>% separate(Scenario,sep="_",into = c("Scenario","Time")) %>% mutate(Scenario=sub("GreenControl","RCP8.5",Scenario)) %>% mutate(Scenario=sub("Green","Melting",Scenario)) %>% filter(Time=="T3")%>% group_by(LossCat,Scenario) %>% summarise(Freq=n()) %>% mutate(LossCat=factor(LossCat,levels=c("Gain","Moderate Loss","Severe Loss","Extreme Loss","Complete Loss"))) %>% mutate(Scenario=factor(Scenario,levels=c("RCP8.5","Melting0.5","Melting1","Melting1.5","Melting3"))) %>% 
  ggplot(aes( axis1=Scenario, axis2=LossCat, y = Freq)) +
  #scale_x_discrete(limits = c("LossCat","Scenario"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = LossCat)) +
  geom_stratum() +
  scale_fill_viridis_d(option="B",direction=-1,begin = 0.4,end=0.9) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void() + theme(legend.position = "none") +
NULL
dev.off()



