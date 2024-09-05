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

length(tissueID)


## Only focus on sample we selected
mif_cohort_2023_08_21 <- read_excel("data/mif_cohort_2023-08-21.xlsx")
tissueID<-mif_cohort_2023_08_21$`Tissue ID`
tissueID<-sub("_"," ",tissueID)


## read Seurat object for each selected sample
## The phenotyping will be done separately
objectpath<-"~/flow_cytometry/multiplexed_image/melanoma_data/Seurat_object/no_tumor/"
#objectpath<-"~/flow_cytometry/multiplexed_image/melanoma_data/corrected_seurat_object/"
list.files(objectpath,pattern = ".rds")
objectname<-paste0(objectpath,tissueID,".rds")


## read the celltype signature matrix--------------------------------------------------

marker_list<-c("panCK.SOX10","CD20","CD68","CD3","CD8","FOXP3","Ki67","TCF1.7","LAG3","PD.1","PD.L1","TOX")
celltype_signature <- read_excel("data/20230823 Cell Types v6 reformat.xlsx")


## identify major cell type with cell lineage marker "panCK.SOX10","CD20","CD68","CD3","CD8","FOXP3"
majortype_signature<-celltype_signature %>% 
  dplyr::select(panCK.SOX10,CD20,CD68,CD3,CD8,FOXP3) %>%
  unique()
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

## phenotyping each sample-------------------------------------------------------------------

source("scripts/Cell_phenotyping.R")

feature_table_all_images<-NULL
#select_dist<-200 ## 100 microns
print_figure<-1
nbins<-36

###i<-52 A good example 

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
  
  if(sum(tissue$majortype == "3:Bcell")){
    
    ## Select cells within certain distances to tumor
    tissue<-subset(tissue,subset = (dist2Bcell <= 720 & majortype %in% c("2:Macrophage", "4:CD4T", "5:CD8T", "6:Treg")))
    #table(tissue$majortype)
    tissue_save<-tissue
    
    if(print_figure){
      
      ### all immune cells 
      p0<-gghistogram(tissue@meta.data, x = "dist2Bcell",color = "majortype",fill = "majortype",bins = nbins)
      
      ## Only Look at T cells
      if(sum(tissue$majortype %in% c("4:CD4T", "5:CD8T","6:Treg"))){
        
        tissue<-subset(tissue,subset = (majortype %in% c("4:CD4T", "5:CD8T","6:Treg")))
        
        p1<-gghistogram(tissue@meta.data, x = "dist2Bcell",color = "majortype",fill = "majortype",bins = nbins)
        
        tissue$Classify_TCF1.7<-as.factor(tissue$Classify_TCF1.7)
        p2<-gghistogram(tissue@meta.data, x = "dist2Bcell",color = "Classify_TCF1.7",fill = "Classify_TCF1.7",bins = nbins)
        
        tissue$Classify_LAG3<-as.factor(tissue$Classify_LAG3)
        p3<-gghistogram(tissue@meta.data, x = "dist2Bcell",color = "Classify_LAG3",fill = "Classify_LAG3", bins = nbins)
        
        tissue$Classify_Ki67<-as.factor(tissue$Classify_Ki67)
        p4<-gghistogram(tissue@meta.data, x = "dist2Bcell",color = "Classify_Ki67",fill = "Classify_Ki67", bins = nbins)
        
        tissue$Classify_PD.1<-as.factor(tissue$Classify_PD.1)
        p5<-gghistogram(tissue@meta.data, x = "dist2Bcell",color = "Classify_PD.1",fill = "Classify_PD.1", bins = nbins)
        
        tissue$Classify_PD.L1<-as.factor(tissue$Classify_PD.L1)
        p6<-gghistogram(tissue@meta.data, x = "dist2Bcell",color = "Classify_PD.L1",fill = "Classify_PD.L1", bins = nbins)
        
        tissue$Classify_TOX<-as.factor(tissue$Classify_TOX)
        p7<-gghistogram(tissue@meta.data, x = "dist2Bcell",color = "Classify_TOX",fill = "Classify_TOX", bins = nbins)
        
        wrap_plots(p0,p1,p2,p3,p4,p5,p6,p7,nrow = 4)
        ggsave(paste0("figures/dist2Bcells/","tissue",tissueID[i],".png"),width = 12,height = 12,units = "in")
      }else{
        wrap_plots(p0)
        ggsave(paste0("figures/dist2Bcells/","tissue",tissueID[i],".png"),width = 6,height = 3,units = "in")
      }
    
    }
    
    ## cut into bins
    tissue<-tissue_save
    tissue$dist2Bcell_cut<-cut(tissue$dist2Bcell, nbins)
    ## as.matrix(table(tissue$regions,tissue$majortype))
    ## as.matrix(table(tissue$regions,tissue$Classify_TCF1.7))
    tissue$tissueID<-tissueID[i]
    
    feature_table<-tissue@meta.data %>% select("tissueID","majortype","dist2Bcell","dist2Bcell_cut",
                                "Classify_Ki67","Classify_TCF1.7","Classify_LAG3",
                                "Classify_PD.1","Classify_PD.L1","Classify_TOX")
    ##dim(feature_table)
    feature_table_all_images<-rbind(feature_table_all_images,feature_table)
  }
  
}

## Output data
saveRDS(feature_table_all_images,file = "figures/dist2Bcell_all_images.rds")



