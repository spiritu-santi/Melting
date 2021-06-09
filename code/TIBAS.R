library(GGally)
library(viridis)
library(tidyverse)
library(patchwork)
setwd("Documents/1.PROYECTOS/1.Termohalina/Final_v4/TIBS")
###### ESTIMATE PROPORTIONAL RANGE CHANGE #######
list_arch <- list.files(pattern="RangeTib",recursive = F)
list_tibs <- lapply(list_arch, function(x) mget(load(x))) %>% lapply(function(x)x[[1]]) %>% bind_rows() %>% 
mutate_if(is.numeric , replace_na, replace = 0) %>% dplyr::rename("Present"=Present_) -> range_tib
range_tib
range_tib %>% mutate_at(5:ncol(.), ~((.x - Present)/Present)) %>% 
  mutate(Present = 0) -> STDrange_tib
STDrange_tib
range_tib

###### PLOT CHANGE BY GROUP AND COUNTRY ########
STDrange_tib %>% ungroup() %>% group_by(Group) %>% summarise_at(3:18, quantile, probs=0.5) -> df_group
STDrange_tib %>% ungroup() %>% group_by(Country) %>% summarise_at(3:18, quantile, probs=0.5) -> df_country

pdf("~/Documents/1.PROYECTOS/1.Termohalina/TERMOS_MS/FINAL_v2/EDIT/Figure_ED2.pdf")
scen <- c("GreenControl_","Green0.5_","Green1_","Green1.5_","Green3_")
t_title <- c("RCP_8.5","Scenario_0.5","Scenario_1","Scenario_1.5","Scenario_3")
limites = range(df_country[,-1])
plotas <- list()
for (i in 1:length(scen)){ 
  if(i!=1) {
    plotas[[i]] <-  df_group %>% dplyr::rename("Present_T0"=Present) %>% #pivot_longer(cols=-1,names_to = c("Scenario","Time"),names_pattern = "(.*)_(.*)",values_to = "Estimate") -> df_group
ungroup() %>% mutate(Group=as_factor(Group)) %>% dplyr::select(1,2,contains(scen[i])) %>% 
  rename_with(~all_of(c("Group","T0","T1","T2","T3"))) %>% 
  ggparcoord(columns = c(3:ncol(.)), groupColumn = 1,scale="globalminmax",showPoints = TRUE, boxplot=F,
             title=t_title[i], alphaLines = 1) + 
 ylab("Range Loss") +
      scale_y_continuous(limits=limites) +
  scale_color_brewer(palette = "Set1") +
  theme(
    legend.position="bottom",
    plot.title = element_text(size=13),
    panel.background=element_blank(),
    panel.border=element_rect(colour="black",fill=NA),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  xlab("")

  }
  else { 
    plotas[[i]] <- df_group %>% dplyr::rename("Present_T0"=Present) %>% #pivot_longer(cols=-1,names_to = c("Scenario","Time"),names_pattern = "(.*)_(.*)",values_to = "Estimate") -> df_group
      ungroup() %>% mutate(Group=as_factor(Group)) %>% dplyr::select(1,2,contains(scen[i])) %>% 
      rename_with(~all_of(c("Group","T0","T1","T2","T3"))) %>% 
      ggparcoord(columns = c(3:ncol(.)), groupColumn = 1,scale="globalminmax",showPoints = TRUE, boxplot=F,title=t_title[i], alphaLines = 1) + 
      ylim(c(1,0)) + ylab("Range Loss") +
      scale_y_continuous(limits=limites) +
      scale_color_brewer(palette = "Set1")+
      theme(
        legend.position="bottom",
        plot.title = element_text(size=13),
        panel.background=element_blank(),
        panel.border=element_rect(colour="black",fill=NA)
  ) +
      xlab("")
  
}
}
wrap_plots(plotas) +
  plot_layout(ncol=3,byrow=T,guides = "collect",tag_level = 'new') + 
  plot_annotation(title = paste("Median range loss by Taxonomic group"),
                  caption='Ureta et al. 2020') & theme(legend.position = 'bottom')

scen <- c("GreenControl_","Green0.5_","Green1_","Green1.5_","Green3_")
t_title <- c("RCP_8.5","Scenario_0.5","Scenario_1","Scenario_1.5","Scenario_3")
limites = range(df_country[,-1])
plotas <- list()
for (i in 1:length(scen)){ 
  if(i!=1) {
    plotas[[i]] <-  df_country %>% dplyr::rename("Present_T0"=Present) %>%
      ungroup() %>% mutate(Group=as_factor(Country)) %>% dplyr::select(1,2,contains(scen[i])) %>% 
      rename_with(~all_of(c("Country","T0","T1","T2","T3"))) %>% 
      ggparcoord(columns = c(3:ncol(.)), groupColumn = 1,scale="globalminmax",showPoints = TRUE, boxplot=F,
                 title=t_title[i], alphaLines = 1) + 
      ylab("Range Loss") +
      scale_y_continuous(limits=limites) +
  theme(
        legend.position="bottom",
        plot.title = element_text(size=13),
        panel.background=element_blank(),
        panel.border=element_rect(colour="black",fill=NA),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
      ) +
      xlab("")
    
  }
  else { 
    plotas[[i]] <- df_country %>% dplyr::rename("Present_T0"=Present) %>%
      ungroup() %>% mutate(Group=as_factor(Country)) %>% dplyr::select(1,2,contains(scen[i])) %>% 
      rename_with(~all_of(c("Country","T0","T1","T2","T3"))) %>% 
      ggparcoord(columns = c(3:ncol(.)), groupColumn = 1,scale="globalminmax",showPoints = TRUE, boxplot=F,title=t_title[i], alphaLines = 1) + 
      ylim(c(1,0)) + ylab("Range Loss") +
      scale_y_continuous(limits=limites) +
      theme(
        legend.position="bottom",
        plot.title = element_text(size=13),
        panel.background=element_blank(),
        panel.border=element_rect(colour="black",fill=NA)
      ) +
      xlab("")
    
  }
}
wrap_plots(plotas) +
  plot_layout(ncol=3,byrow=T,guides = "collect",tag_level = 'new') + 
  plot_annotation(title = paste("Median range loss by Country"),
                  caption='Ureta et al. 2020') & theme(legend.position = 'bottom')
dev.off()




#### STANDARDIZE RICHNESS VALUES BY COUNTRY ####
#### The procedure considers the maximum and minimum richness values observed within a country #### 
#### across all scenarios and times (including present).                                       ####

list_arch <- list.files(pattern="SumTib",recursive = F)
list_tibs <- lapply(list_arch, function(x) mget(load(x))) %>% lapply(function(x)x[[1]]) %>% bind_rows() 
list_tibs %>% ungroup() %>% group_by(Country,Scenario,Time,CellID) -> sum_tib
sum_tib %>% dplyr::select(Longitud,Latitude,SR) %>% 
  summarise(SR=sum(SR),Longitud=first(Longitud),Latitude=first(Latitude)) %>% 
  ungroup() %>% group_by(Country) %>% mutate(SR_std = (SR-min(SR))/(max(SR)-min(SR))) -> df_tib
