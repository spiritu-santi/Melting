library(tidyverse)
library(viridis)
library(patchwork)
library("rnaturalearth")
library("rnaturalearthdata")
#### PROCESS FINAL TIBS #####
setwd("~/Documents/1.PROYECTOS/1.Termohalina/TERMOS_MS/EXTRAS/BETA")
pais="Australia"
g <- raster::raster(nrows=180*12,ncols=360*12,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1,crs = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% as(., 'SpatialPixels')

if(pais=="Australia") {
  lista_ranges <- list.files("Partials",pattern = "Partial_FinalTib",full.names = T)
  testGlobal_tib <- list()

  lista_subs <- list.files("Partials",pattern = "Partial_SumTib",full.names = T)
  sum_tibs <- lapply(lista_subs, function(x) mget(load(x)))
  sum_tibs
  for (i in 1:length(sum_tibs)){
    sum_tibs[[i]] <- sum_tibs[[i]]$summary_tib  %>% ungroup()
  }
  sum_tibs %>% bind_rows() -> sum_tibs
  sum_tibs %>% group_by(Group,Country,Scenario,Time,CellID) %>% summarise(SR=sum(SR),mean_cv=mean(mean_cv),sd_cv=max(sd_cv),Longitud=first(Longitud),Latitude=first(Latitude)) %>%  select(Group,Country,Scenario,Time,SR,mean_cv,sd_cv,CellID,Longitud,Latitude) -> sum_tibs
 sum_tibs %>%  filter(Scenario=="Present") %>% group_by(CellID) %>% summarise(SR=sum(SR)) %>% select(SR) %>% summarize(MED=ceiling((max(SR) * 0.6))) -> umbral
  xlims=range(sum_tibs$Longitud)
  ylims=range(sum_tibs$Latitude)
  sum_tibs %>% ungroup() %>% group_by(CellID,Scenario) %>% summarise(SR=sum(SR)) -> SRs
  SRs %>% pull(CellID) %>% match(.,g@grid.index) -> index
  SRs %>% bind_cols(x=g@coords[index,"x"],y=g@coords[index,"y"]) -> SRs
  SRs %>% filter(Scenario=="Present") %>% filter(SR > umbral$MED[1]) %>% pull(CellID) -> celdas
  # val_lims = SRs %>% filter(Scenario=="Present") %>% ungroup() %>% select(SR) %>% range()
  #  SRs %>% filter(Scenario=="Present") %>% ggplot() + xlim(xlims) + ylim(ylims) +
  #    geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="black",fill="grey80") +
  #    geom_tile(data=SRs,aes(y=y, x=x,fill = SR),colour = "white",size=0.00001,height=0.1,width=0.1) + 
  #    scale_fill_viridis(option = "inferno",direction=-1,begin=0,end=1,limits=val_lims) +
  #    theme(legend.position = "bottom", legend.key.width = unit(0.6, "cm"),
  #          panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
  #          axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),
  #          axis.ticks.x=element_blank(),
  #          panel.background=element_rect(colour="black",fill="white"),
  #          panel.border=element_rect(colour="black",fill=NA)) +
  #    NULL
  rm(sum_tibs)
  for(k in 1:length(lista_ranges)){
 cat("Loading partial tib",k,"\n")
  load(lista_ranges[k])
  cat("--- Processing partial tib",k,"\n")
  list_tibs %>% group_by(CellID) -> sum_tib
  sum_tib %>%  mutate(SceTime = paste(Scenario,Time,sep="")) %>% group_by(CellID,SceTime) %>% filter(CellID %in% all_of(celdas)) -> testGlobal_tib[[k]]
  }
  testGlobal_tib %>% bind_rows() -> test_tib
  test_tib %>% group_by(CellID,SceTime) %>% nest() -> test_tib
  test_tib %>% mutate(SpeciesList = map2(data, 1, select)) -> test_tib
  test_tib %>% select(CellID,SceTime,SpeciesList) %>% #filter(SceTime %in% all_of(times)) %>% 
    pivot_wider(names_from = SceTime,values_from = SpeciesList) -> test_tib
  save(test_tib,file=paste(pais,"_BetaTib.R",sep=""))
}
if(pais!="Australia") {
input=paste(pais,"FinalTib.Rdata",sep="")
load(input)
list_tibs %>% group_by(CellID) -> sum_tib
sum_tib %>% group_by(CellID,Scenario) %>% summarise(SR=n()) -> SRs
SRs %>% pull(CellID) %>% match(.,g@grid.index) -> index
SRs %>% bind_cols(x=g@coords[index,"x"],y=g@coords[index,"y"]) -> SRs
xlims=range(SRs$x)
ylims=range(SRs$y)
val_lims = SRs %>% filter(Scenario=="Present") %>% ungroup() %>% select(SR) %>% range()
# SRs %>% filter(Scenario=="Present") %>% ggplot() + xlim(xlims) + ylim(ylims) +
#   geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="black",fill="grey80") +
#   geom_tile(data=SRs,aes(y=y, x=x,fill = SR),colour = "white",size=0.00001,height=0.1,width=0.1) + 
#   scale_fill_viridis(option = "inferno",direction=-1,begin=0,end=1,limits=val_lims) +
#   theme(legend.position = "bottom", legend.key.width = unit(0.6, "cm"),
#         panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
#         axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),
#         axis.ticks.x=element_blank(),
#         panel.background=element_rect(colour="black",fill="white"),
#         panel.border=element_rect(colour="black",fill=NA)) +
#   NULL

sum_tib %>% filter(Scenario=="Present") %>% group_by(CellID) %>% summarise(SR=n()) %>% select(SR) %>% summarize(MED=ceiling((max(SR) * 0.6))) -> umbral
SRs %>% filter(Scenario=="Present") %>% filter(SR > umbral$MED[1]) %>% pull(CellID) -> celdas
sum_tib %>%  mutate(SceTime = paste(Scenario,Time,sep="")) %>% group_by(CellID,SceTime) %>% 
  filter(CellID %in% all_of(celdas)) %>% nest() -> test_tib
test_tib %>% mutate(SpeciesList = map2(data, 1, select)) -> test_tib
test_tib %>% select(CellID,SceTime,SpeciesList) %>% #filter(SceTime %in% all_of(times)) %>% 
  pivot_wider(names_from = SceTime,values_from = SpeciesList) -> test_tib
save(test_tib,file=paste(pais,"_BetaTib.R",sep=""))
}
#####

#### ESTIMATE BETAS ######
# "Present"        "Green0.5T1"     "Green0.5T2"     "Green0.5T3"     "Green1T1"      "Green1T2"       "Green1T3"       "Green1.5T1"     "Green1.5T2"     "Green1.5T3"     "Green3T1"      "Green3T2"       "Green3T3"       "GreenControlT1" "GreenControlT2" "GreenControlT3"
rm(list=ls())
betas <- function(x,y){ 
  a <- intersect(x,y) %>% length()
  b <- setdiff(x,y) %>% length()
  c <-  setdiff(y,x) %>% length()
  beta.sor <-  (b+c)/(2*a+b+c)
  beta.sim <- min(b,c)/(a+min(b,c))
  beta.nes <- beta.sor-beta.sim
  return(list(Bsor=beta.sor,Bsim=beta.sim,Bnes=beta.nes))
}
times = c("Present","GreenControlT1","GreenControlT2","GreenControlT3","Green0.5T1","Green0.5T2","Green0.5T3",
          "Green1T1","Green1T2","Green1T3","Green1.5T1","Green1.5T2","Green1.5T3","Green3T1","Green3T2","Green3T3")
g <- raster::raster(nrows=180*12,ncols=360*12,xmn=-180,xmx=180,ymn=-90,ymx=90,vals=1,crs = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) %>% as(., 'SpatialPixels')
paises <- c("Australia","Brazil","China","Colombia","Ecuador","India",
            "Indonesia","Madagascar","Mexico","Peru","Phillipines","Venezuela")
beta_gral <- list()
for (kc in 2:length(times)){ 
  cat("SCENARIO", times[kc],"\n")
#for (p in 1:length(paises)) {
res_beta <- list()
for(p in 1:12){ 
pais=paises[p]
cat("---- Estimating betas and mapping for", pais,"\n")
if(pais=="Phillipines") { m <- maps::map("world","Philippines",plot=F)}
if(pais!="Phillipines"){ m <- maps::map("world",pais,plot=F)}
xlims=m$range[1:2]
ylims=m$range[3:4]
load(paste(pais,"_BetaTib.R",sep=""))
time0=times[1]
time1=times[kc]
lapply(test_tib[[time0]],is.null) %>% unlist() %>% which(.==TRUE) -> uno
c(uno,lapply(test_tib[[time1]],is.null) %>% unlist() %>% which(.==TRUE)) -> uno
if(length(uno) > 0){
test_tib[[time0]][-uno] %>% map(~.x %>% pull(Species)) -> sites1
test_tib[[time1]][-uno] %>% map(~.x %>% pull(Species)) -> sites2
}
if(length(uno) == 0) {
  test_tib[[time0]] %>% map(~.x %>% pull(Species)) -> sites1
  test_tib[[time1]] %>% map(~.x %>% pull(Species)) -> sites2
}
res <- list()
for (i in 1:length(sites1)) {
  betas(sites1[[i]],sites2[[i]])-> res[[i]]
}
Bsor <- unlist(lapply(res,"[[",'Bsor'))
Bsim <- unlist(lapply(res,"[[",'Bsim'))
Bnes <- unlist(lapply(res,"[[",'Bnes'))
if(length(uno) > 0){test_tib %>% .[-uno,] -> inter}
if(length(uno) == 0){test_tib -> inter}
inter %>% dplyr::select(CellID) %>% bind_cols(Bsor=Bsor,Bsim=Bsim,Bnes=Bnes) -> beta_tib
beta_tib %>% ungroup() %>% summarise(M_Bsor=median(Bsor),M_Bsim=median(Bsim),M_Bnes=median(Bnes)) %>% mutate(Country=pais,Time=time1,.before=1) -> res_beta[[p]]
}
res_beta %>% do.call(rbind,.) -> beta_gral[[kc]]
}
beta_gral %>% do.call(rbind,.) %>% arrange(Country) %>% dplyr::select(Country,Time,M_Bsor) %>% 
  pivot_wider(names_from = Country, values_from = starts_with("M_")) %>% .[c(1,4,7,10,13,2,5,8,11,14,3,6,9,12,15),] %>% mutate_at(-1,~format(floor(.x * -1000) /-1000, nsmall = 2)) %>% flextable::flextable() %>% flextable::save_as_docx(.,path="../Table_Betas_Country.docx")

beta_gral %>% do.call(rbind,.) %>% arrange(Country) %>% dplyr::select(Country,Time,M_Bsor) %>% 
  pivot_wider(names_from = Country, values_from = starts_with("M_")) %>% .[c(1,4,7,10,13,2,5,8,11,14,3,6,9,12,15),] %>% rowwise(Time) %>% mutate(Max=max(c_across(where(is.numeric))),Min=min(c_across(where(is.numeric)))) %>% dplyr::select(Time,Max,Min) %>% ungroup() %>% filter(!grepl("Control",Time)) %>% filter(grepl("T1",Time)) %>% 
  summarise(mean(Max),mean(Min))

beta_tib %>% pull(CellID) %>% match(.,g@grid.index) -> index
beta_tib %>% bind_cols(x=g@coords[index,"x"],y=g@coords[index,"y"]) -> beta_tib
val_lims= c(0,1)
Bsor <- ggplot() + xlim(xlims) + ylim(ylims) +
  geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="black",fill="grey80") +
  geom_tile(data=beta_tib,aes(y=y, x=x,fill = Bsor),colour = "white",size=0.00001,height=0.1,width=0.1) + 
  scale_fill_viridis(option = "inferno",direction=-1,begin=0,end=1,limits=val_lims) +
  labs(title="Beta diversity (Bnes)") +
  theme(legend.position = "bottom", legend.key.width = unit(0.6, "cm"),
        panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background=element_rect(colour="black",fill="white"),
        panel.border=element_rect(colour="black",fill=NA)) +
  NULL
Bsim <- ggplot() + xlim(xlims) + ylim(ylims) +
  geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="black",fill="grey80") +
  geom_tile(data=beta_tib,aes(y=y, x=x,fill = Bsim),colour = "white",size=0.00001,height=0.1,width=0.1) + 
  scale_fill_viridis(option = "inferno",direction=-1,begin=0,end=1,limits=val_lims) +
  labs(title="Turnover (Bsim)") +
  theme(legend.position = "bottom", legend.key.width = unit(0.6, "cm"),
        panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background=element_rect(colour="black",fill="white"),
        panel.border=element_rect(colour="black",fill=NA)) +
  NULL
Bnes <- ggplot() + xlim(xlims) + ylim(ylims) +
  geom_sf(data=ne_countries(scale=110,type="countries",returnclass = "sf"),colour="black",fill="grey80") +
  geom_tile(data=beta_tib,aes(y=y, x=x,fill = Bnes),colour = "white",size=0.00001,height=0.1,width=0.1) + 
  scale_fill_viridis(option = "inferno",direction=-1,begin=0,end=1,limits=val_lims) +
  labs(title="Nestedness (Bnes)") +
  theme(legend.position = "bottom", legend.key.width = unit(0.6, "cm"),
        panel.grid = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background=element_rect(colour="black",fill="white"),
        panel.border=element_rect(colour="black",fill=NA)) +
  NULL
hists <- beta_tib %>% pivot_longer(cols=starts_with("B"),names_to = "Index",values_to = "Value") %>% ggplot(aes(y=Value,x=Index)) + ylim(0,1) +
  geom_jitter(aes(color = Value)) + scale_color_viridis_c(option="B",direction = -1,limits=val_lims) + ylab("") + xlab("") +  geom_violin(fill=NA,size=1,color="white",bw=.01) + 
  geom_violin(fill=NA,size=0.7,color="black",bw=.01) + theme(legend.position = "none",
               panel.grid = element_blank(),
               plot.background = element_blank(),
               panel.background= element_blank()) + #coord_flip() +
  NULL
pdf(paste("Beta_",time1,"/",pais,"_betaTib.pdf",sep=""))
  ((Bsim / Bnes /Bsor)|hists) +
  plot_layout(byrow=F,guides = "collect",tag_level = 'new') + 
  plot_annotation(title = "Temporal Beta Diversity",subtitle = "Present - 2030",
                  caption= "Ureta et al.") & theme(legend.position = '')
dev.off()
#}


