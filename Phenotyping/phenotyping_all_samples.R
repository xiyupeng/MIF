library(readxl)
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(ggsci)

## read sample list-------------------------------------------------------------------

#sample_info<- read_excel("data/20230306_Melanoma_Slide_tissue_info_updated.xlsx")
sample_info<- read_excel("data/20230310_QC_cell_phenotyping.xlsx")

## list of all good samples 
selected_tissue<-sample_info %>% filter(Set %in% c(1,2,3,4,5,6,7,8)) %>% 
  filter(is.na(Status) | Status == "OK") %>% filter(is.na(Status_II)| Status_II == "OK" | regate_status == "OK" )
tissueID<-selected_tissue$TissueID

## add back samples that need re-threshold SOX10
# selected_tissue<-sample_info %>% filter(Set %in% c(1,2,3,4,5,6,8)) 
# tissueID<-selected_tissue$TissueID[grep("SOX10",selected_tissue$Status)]

## read Seurat object for each selected sample
## Note:
## two samples have two rds files
## The phenotyping will be done separately
#objectpath<-"~/flow_cytometry/multiplexed_image/melanoma_data/Seurat_object/no_tumor/"
objectpath<-"~/flow_cytometry/multiplexed_image/melanoma_data/seurat_object/"
list.files(objectpath,pattern = ".rds")
objectname<-paste0(objectpath,tissueID,".rds")


## read the celltype signature matrix--------------------------------------------------

marker_list<-c("panCK.SOX10","CD20","CD68","CD3","CD8","FOXP3","Ki67","TCF1.7","LAG3","PD.1","PD.L1","TOX")
celltype_signature <- read_excel("data/20230823 Cell Types v6 reformat.xlsx")


## identify major cell type with cell lineage marker "panCK.SOX10","CD20","CD68","CD3","CD8","FOXP3"
majortype_signature<-celltype_signature %>% 
                     dplyr::select(panCK.SOX10,CD20,CD68,CD3,CD8,FOXP3) %>%
                     unique()
### change it when adding new celltypes
#rownames(majortype_signature)<-c("1:tumor/epi","2:Macrophage","3:Bcell","4:CD4T","5:CD8T","6:Treg","7:other")
rownames(majortype_signature)<-c("3:Bcell","4:CD4T","5:CD8T","2:Macrophage","7:other","6:Treg","1:tumor/epi")
## the final matrix we used to phenotype major celltypes
majortype_signature

## identify functional subset with marker  "Ki67","TCF1.7","LAG3","PD.1","PD.L1","TOX"
celltype_signature$majortype<-as.factor(celltype_signature$Celltype)
#celltype_signature$majortype<-recode_factor(celltype_signature$majortype,
#                                           'Tumor v. Epithelium'= "1:tumor/epi", 'Macrophage' = "2:Macrophage",
#                                           'B cell'= "3:Bcell", 'CD4+ T cell'="4:CD4T",'CD8+ T cell'="5:CD8T",
#                                           'Treg' = "6:Treg", 'Other' = "7:other")
## new version
celltype_signature$majortype<-recode_factor(celltype_signature$majortype,
                                           'TumorEpithelial'= "1:tumor/epi", 'Macrophage' = "2:Macrophage",
                                           'B cell'= "3:Bcell", 'CD4 T cell'="4:CD4T",'CD8 T cell'="5:CD8T",
                                           'Treg' = "6:Treg", 'Other' = "7:other")

subtype_signature<-celltype_signature %>% 
                  dplyr::select(majortype,subtype, Ki67, PD.1, LAG3, TOX, TCF1.7, PD.L1)
subtype_signature$majortype<-as.character(subtype_signature$majortype)                  

## phenotyping each sample---------------------------------------------------------------

source("scripts/Cell_phenotyping.R")

cell_composition<-NULL
cell_composition_tumor<-NULL
cell_composition_stroma<-NULL
output_figure<-1
output_rds<-0
#select_dist<-1028 ## 100 microns // 286 pixel for 100 microns when scale = 0.35 // 1028 pixel for 360 microns

for(i in 1:length(objectname)){
  
  tissue<-readRDS(file = objectname[i])
  cat(tissueID[i],"\n")
  
  ## if no cell phenotyping result
  if(! "majortype" %in% colnames(tissue@meta.data)){
    ## phenotyping Major Celltypes
    result<-phenotyping_MajorCelltype(majortype_signature,tissue)
    tissue$majortype<-as.factor(result$celltype)
    Idents(tissue)<-"majortype"
  }
  
  ## Select cells within certain distances to tumor
  #tissue<-subset(tissue,subset = dist2tumor <= select_dist)
  
  
  if(output_figure){  ## output figures
    
    ###Check the quality of major cell type with visualization--------------------------------
    marker_list<-c("panCK.SOX10","CD20","CD3","CD8","CD68","FOXP3","Ki67","TCF1.7","LAG3","PD.1","PD.L1","TOX")
    
    # ## heatmap
    # DoHeatmap(subset(tissue, downsample = 1000),features = marker_list,slot = "scale.data",disp.max = 1, disp.min = 0)
    # ggsave(paste0("figures/heatmap/","heatmap_",tissueID[i],".png"),width = 7,height = 7,units = "in")
    # 
    # ## ridge plot
    # plist<-list()
    # for(m in marker_list){
    #   plist[[m]]<-RidgePlot(object = tissue, features = m,group.by = paste0("Classify_",m),slot = "scale.data")
    # }
    # wrap_plots(plist, ncol = 3)
    # ggsave(paste0("figures/ridge/","ridge_",tissueID[i],".png"),width = 12,height = 12,units = "in")
    # 
    
    ## tissue plot
    tissue@meta.data %>%
      ggplot(aes(X,Y,color = majortype,shape = region))+geom_point(size = 0.5)+
      scale_shape_manual(name='Region',breaks=c('glass', 'Stroma', 'Tumor'),values=c('glass'=1, 'Stroma'=3, 'Tumor'=16))+
      scale_color_manual(name='Celltype',breaks=c("1:tumor/epi","2:Macrophage","3:Bcell",
                                                "4:CD4T","5:CD8T","6:Treg","7:other"),
                         values=c("1:tumor/epi" = "red" ,"2:Macrophage"  = "#C49A00" ,"3:Bcell" = "#53B400",
                                                "4:CD4T" =  "#00C094" ,"5:CD8T" = "#00B6EB" ,"6:Treg" = "#A58AFF" ,
                                  "7:other" = "gray"))
      
    ggsave(paste0("figures/tissue/","tissue",tissueID[i],".png"),width = 36,height = 36,units = "in")
  }
  
  if(output_rds){
    ## keep all cells
    saveRDS(tissue, file = paste0("~/flow_cytometry/multiplexed_image/melanoma_data/Seurat_object/seurat_objects_wphenotypes/",tissueID[i],".rds"))
  }
  
  ## phenotyping all listed functional subsets----------------------------------------------
  ## The whole slide (all regions)
  counts<-Calc_FuncSet_Composition(tissue, subtype_signature)
  cell_composition<-cbind(cell_composition,counts)
  
  ## tumor region
  counts<-Calc_FuncSet_Composition(tissue, subtype_signature,"Tumor")
  cell_composition_tumor<-cbind(cell_composition_tumor,counts)
  
  ## stroma region
  counts<-Calc_FuncSet_Composition(tissue, subtype_signature,"Stroma")
  cell_composition_stroma<-cbind(cell_composition_stroma,counts)
  
  rm(tissue)
}

## all region
colnames(cell_composition)<-tissueID
rownames(cell_composition)<-celltype_signature$subtype
celltype_signature_composition<-as.data.frame(cbind(celltype_signature,cell_composition))
write.csv(celltype_signature_composition,file = "data/cell_phenotyping_result.csv",
          quote = FALSE,row.names = FALSE)

## tumor region
colnames(cell_composition_tumor)<-tissueID
rownames(cell_composition_tumor)<-celltype_signature$subtype
celltype_signature_composition<-as.data.frame(cbind(celltype_signature,cell_composition_tumor))
write.csv(celltype_signature_composition,file = "data/cell_phenotyping_result_tumor_region.csv",
          quote = FALSE,row.names = FALSE)

## stroma region
colnames(cell_composition_stroma)<-tissueID
rownames(cell_composition_stroma)<-celltype_signature$subtype
celltype_signature_composition<-as.data.frame(cbind(celltype_signature,cell_composition_stroma))
write.csv(celltype_signature_composition,file = "data/cell_phenotyping_result_stroma_region.csv",
          quote = FALSE,row.names = FALSE)
