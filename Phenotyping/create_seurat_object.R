library(tidyverse)
library(Seurat)
library(SeuratObject)

source("scripts/normalization.R")

## 1,2,3,4,6,8 staining set
filepath<- "~/flow_cytometry/multiplexed_image/melanoma_data/"
#filepath<-"W:/Staining Set 4 corrections/Halo archive 2023-02-24 10-32-23 - v3.4.2986/ObjectData/"
filenames<-list.files(filepath,pattern = "Count_object_results.csv")

## 5,7 staining set
filepath<-"C:/Users/pengx1/Box/Xiyu-Jim-Fiona Shared Melanoma Files/JWS HALO CSVs (Set 5 + 7)/7th Staining/"
filenames<-list.files(filepath,pattern = ".csv")

objectpath<-"~/flow_cytometry/multiplexed_image/melanoma_data/Seurat_object/"

head(filenames)
length(filenames)

## provide the list of markers
marker_list<-c("panCK.SOX10","CD20","CD3","CD8","CD68","FOXP3","Ki67","TCF1.7","LAG3","PD.1","PD.L1","TOX")

for(i in 1:length(filenames)){
  
  ## for set 5, 7
  sample<-gsub(".csv","",filenames[i]) #or
  ## for set 1,2,3,4,6,8
  ##sample<-gsub("_R1.*","",filenames[i]) 
  
  IMC_data<-read.csv(paste0(filepath,filenames[i]))
  
  IMC_data<-IMC_data %>%
    mutate(X = (XMin + XMax)/2, Y = (YMin + YMax)/2)
  
  ## transform the intensity values with asinh() 
  mat.trans<-IMC_data %>% dplyr::select(ends_with("Cell.Intensity")) %>% 
    select(paste0(marker_list,".Cell.Intensity")) %>%
    asinh(.)
  colnames(mat.trans)<-gsub(".Cell.Intensity","",colnames(mat.trans))
  
  ## obtain the pos/neg classifier
  mat.classify<-IMC_data %>% dplyr::select(ends_with("Positive.Classification"))
  colnames(mat.classify)<-gsub(".Positive.Classification","",colnames(mat.classify))
  
  ## check column for each marker.
  for( m in 1:length(marker_list)){
    marker<-marker_list[m]
    if(!marker %in% colnames(mat.classify)){
      mat.classify$marker<-0  ## cells are all negative if the column does not exist
      colnames(mat.classify)[ncol(mat.classify)]<-marker
      }
  }
  
  ## scale the intensity values 
  mat.trans %>% dplyr::select(all_of(marker_list)) ->mat.trans
  mat.classify %>% dplyr::select(all_of(marker_list))->mat.classify
  mat.norm<-vrescale(raw = mat.trans,label = mat.classify)

  ## set cell ids for each matrix
  mat.norm<-as.matrix(mat.norm)
  mat.trans<-as.matrix(mat.trans)
  rownames(mat.trans)<-rownames(IMC_data)
  rownames(mat.norm)<-rownames(IMC_data)
  
  ## create seurat object
  ## store the transformed intensities
  MIF <- CreateSeuratObject(counts = t(mat.trans), assay = "protein")
  
  ## Add scaled data in the slot "scale.data"
  MIF[["protein"]]<-SetAssayData(MIF[["protein"]], slot = "scale.data", new.data = t(mat.norm))
  
  ## Add meta data 
  ## Add pos/neg classification of each marker in meta data
  colnames(mat.classify)<-paste0("Classify_",colnames(mat.classify))
  MIF<-AddMetaData(MIF, metadata =  as.data.frame(mat.classify))
  MIF[['X']]<-IMC_data$X
  MIF[['Y']]<-IMC_data$Y
  MIF[['region']]<-IMC_data$Classifier.Label
  MIF[["sample"]]<-sample
  
  ## save the seurat object for each sample
  objectname<-paste0(objectpath,sample,".rds")
  saveRDS(MIF, file = objectname)
  
  cat(i,"\n")
}

