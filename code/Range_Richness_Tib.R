library(GGally)
library(viridis)
library(tidyverse)
library(sf)
library(raster)
library(patchwork)

setwd("/TIBS")
###### ESTIMATE PROPORTIONAL RANGE CHANGE #######
list_arch <- list.files(pattern="RangeTib",recursive = F)
list_tibs <- lapply(list_arch, function(x) mget(load(x))) %>% lapply(function(x)x[[1]]) %>% bind_rows() %>% 
mutate_if(is.numeric , replace_na, replace = 0) %>% rename("Present"=Present_) -> range_tib
range_tib

range_tib %>% mutate_if(is.numeric , replace_na, replace = 0) %>% mutate_at(5:ncol(.), ~((.x - Present)/Present)) %>% mutate(Present = 0) -> STDrange_tib
STDrange_tib
range_tib

#### STANDARDIZE RICHNESS VALUES BY COUNTRY ####
#### The procedure considers the maximum and minimum richness values observed within a country #### 
#### across all scenarios and times (including present).                                       ####

list_arch <- list.files(pattern="SumTib",recursive = F)
list_tibs <- lapply(list_arch, function(x) mget(load(x))) %>% lapply(function(x)x[[1]]) %>% bind_rows() 
list_tibs %>% ungroup() %>% group_by(Country,Scenario,Time,CellID) -> sum_tib

#Group
sum_tib %>% dplyr::select(Longitud,Latitude,SR) %>% 
  summarise(SR=sum(SR),Longitud=first(Longitud),Latitude=first(Latitude)) %>% 
  ungroup() %>% group_by(Country) %>% mutate(SR_std = (SR-min(SR))/(max(SR)-min(SR))) -> df_tib


#plants<-c("Lycophyta","Magnoliophyta","Monilophyta","Pinophyta")
sum_tib %>% dplyr::filter(Group %in%  plants) %>% dplyr::select(Longitud,Latitude,SR) %>% 
  summarise(SR=sum(SR),Longitud=first(Longitud),Latitude=first(Latitude)) %>% 
  ungroup() %>% group_by(Country) %>% mutate(SR_std = (SR-min(SR))/(max(SR)-min(SR))) -> df_tib


#vertebrates <-c("Anfibios","Aves","Mamiferos","Reptiles")
sum_tib %>% dplyr::filter(Group %in%  vertebrates) %>% dplyr::select(Longitud,Latitude,SR) %>% 
  summarise(SR=sum(SR),Longitud=first(Longitud),Latitude=first(Latitude)) %>% 
  ungroup() %>% group_by(Country) %>% mutate(SR_std = (SR-min(SR))/(max(SR)-min(SR))) -> df_tib



