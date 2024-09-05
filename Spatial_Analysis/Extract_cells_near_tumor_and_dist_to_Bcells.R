# Calculate distance between each cell to tumor cells (SOX10+ in tumor region)
library(readxl)
library(Seurat)
library(tidyverse)
library(RANN)
source("scripts/Cell_phenotyping.R")
#source("scripts/spatial_func.R")

GetCoords <- function(tissue) {
  coords <- as.data.frame(cbind(tissue$X,tissue$Y))
  colnames(coords) <- c("X", "Y")
  return(coords)
}

## read the celltype signature matrix---------------------------------------------------------

marker_list<-c("panCK.SOX10","CD20","CD68","CD3","CD8","FOXP3","Ki67","TCF1.7","LAG3","PD.1","PD.L1","TOX")
celltype_signature <- read_excel("data/20230508 Cell Types v5_reformat.xlsx")


## identify major cell type with cell lineage marker "panCK.SOX10","CD20","CD68","CD3","CD8","FOXP3"
##majortype_signature<-as.matrix(celltype_signature[c(1,4,7,10,17,24,25),3:8])
majortype_signature<-celltype_signature %>% 
  dplyr::select(panCK.SOX10,CD20,CD68,CD3,CD8,FOXP3) %>%
  unique %>% as.matrix()
rownames(majortype_signature)<-c("1:tumor/epi","2:Macrophage","3:Bcell","4:CD4T","5:CD8T","6:Treg","7:other")

## identify functional subset with marker  "Ki67","TCF1.7","LAG3","PD.1","PD.L1","TOX"
celltype_signature$majortype<-as.factor(celltype_signature$Celltype)
celltype_signature$majortype<-recode_factor(celltype_signature$majortype,
                                            'Tumor v. Epithelium'= "1:tumor/epi", 'Macrophage' = "2:Macrophage",
                                            'B cell'= "3:Bcell", 'CD4+ T cell'="4:CD4T",'CD8+ T cell'="5:CD8T",
                                            'Treg' = "6:Treg", 'Other' = "7:other")
subtype_signature<-celltype_signature %>% 
  dplyr::select(majortype,subtype, Ki67, PD.1, LAG3, TOX, TCF1.7, PD.L1)
subtype_signature$majortype<-as.character(subtype_signature$majortype)      


####-------------------------------------------------------------------------------------------
### read sample list 
sample_info<- read_excel("data/20230310_QC_cell_phenotyping.xlsx")

selected_tissue<-sample_info %>% filter(Set %in% c(1,2,3,4,5,6,7,8)) %>% 
  filter(is.na(Status) | Status == "OK") %>% filter(is.na(Status_II)| Status_II == "OK" | regate_status == "OK" )

nrow(selected_tissue)
tissueID<-selected_tissue$TissueID

#objectpath<-"~/flow_cytometry/multiplexed_image/melanoma_data/Seurat_object/"
objectpath<-"Seurat_object/"
list.files(objectpath,pattern = ".rds")
objectname<-paste0(objectpath,tissueID,".rds")

majortype<-subtype_signature$majortype
signature<-subtype_signature %>% select(-majortype, -subtype)
markers<-colnames(signature)

##run parallelt
library(foreach)
library(doParallel)
registerDoParallel(15)

output_rds<-1

results<-foreach(i = 1:length(tissueID)) %dopar% {
  
  
  ## read the data
  cat(tissueID[i],"\n")
  tissue<-readRDS(file = objectname[i])
  
  ## if no cell phenotyping result
  if(! "majortype" %in% colnames(tissue@meta.data)){
    result<-phenotyping_MajorCelltype(majortype_signature,tissue)
    tissue$majortype<-as.factor(result$celltype)
    Idents(tissue)<-"majortype"
  }
  
  ###Check the quality of major cell type with visualization--------------------------------
  ##marker_list<-c("panCK.SOX10","CD20","CD3","CD8","CD68","FOXP3","Ki67","TCF1.7","LAG3","PD.1","PD.L1","TOX")
  
  ## split into two sets set1/2 (tumor in tumor region or not)
  tissue$is_tumor_cell<-ifelse(tissue$majortype == "1:tumor/epi" &  tissue$region =="Tumor", 1, 0)
  tissue_tumor<-subset(tissue,subset = is_tumor_cell == 1)
  tissue_other<-subset(tissue,subset = is_tumor_cell == 0)
  
  ## calculate the min distance of each cell in set 2 to tumor cells in set 1
  k_neighbor_centers <- nn2(GetCoords(tissue_tumor), GetCoords(tissue_other), k=1, treetype = "kd",
                            searchtype = "priority")
  tissue_other$dist2tumor<-k_neighbor_centers$nn.dists[,1]
  
  ## calculate the min distance of each cell in set 2 to B cells (Inf if no B cell)
  if(sum(tissue_other$majortype == "3:Bcell")){
    tissue_Bcells<-subset(tissue_other,subset = majortype == "3:Bcell")
    k_neighbor_centers <- nn2(GetCoords(tissue_Bcells), GetCoords(tissue_other), k=1, treetype = "kd",
                            searchtype = "priority")
    tissue_other$dist2Bcell<-k_neighbor_centers$nn.dists[,1]
  }else{
    tissue_other$dist2Bcell<-Inf
  }
  
  if(output_rds){
    saveRDS(tissue_other, file = paste0(objectpath,"no_tumor/",tissueID[i],".rds"))
  }
  
 
  
}
