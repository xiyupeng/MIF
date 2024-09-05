## phenotyping_MajorCelltype
#
# This function phenotyping major celltypes based on transformed scaled intensity values
# majortype_signature is a celltype X marker matrix
# panCK.SOX10 CD20 CD3 CD8 CD68 FOXP3
# 1:tumor/epi            1    0   0   0    0     0
# 2:Macrophage           0    0   0   0    1     0
# 3:Bcell                0    1   0   0    0     0
# 4:CD4T                 0    0   1   0    0     0
# 5:CD8T                 0    0   1   1    0     0
# 6:Treg                 0    0   1   0    0     1
# 7:other                0    0   0   0    0     0
## In Seurat object, 
## the scale.data slot contains marker transformed intensity values scaled to [0,1]
phenotyping_MajorCelltype<-function(majortype_signature, SeuratObject){
  
  cells<-GetAssayData(SeuratObject,slot = "scale.data")
  majortype<-rownames(majortype_signature)
  marker<-colnames(majortype_signature)
  majortype_signature<-as.matrix(majortype_signature)
  
  ## CHECK 
  if(length(majortype) == 0 | length(marker) == 0)
      stop("Need input a valid majortype_signature matrix")
  
  for(m in marker){
    if(!(m %in% rownames(cells)))
       stop(paste0("marker ",m," was not found"))
  }
  
  if(any(cells>1) | any(cells < 0))
    stop("The scaled intensity values should be in [0,1].")
  
  cells_sub<-cells[marker,]
  assign_prob<-majortype_signature %*% log(cells_sub) + (1-majortype_signature) %*% log(1-cells_sub)
  assign_celltype<-majortype[apply(assign_prob,2,which.max)]
  return(list(celltype = assign_celltype, score = assign_prob))
}


## phenotyping_FunctionalSet
#
# count number of cells with positive markers in a specific cell type
# marker<-"Ki67"
# celltype<-"1:tumor/epi"
# phenotyping_FunctionalSet(tissue, marker, celltype)
phenotyping_FunctionalSet<-function(Seurat_object,marker,celltype, tissue_region = "All"){
  ## CHECK
  for(m in marker){
    if(!(m %in% rownames(Seurat_object)))
      stop(paste("marker",m,"was not found in the list of markers !\n"))
  }
  
    
  if(length(celltype) == 0){
    
    if(length(marker) == 0)
      stop("no input cell type nor marker!\n")
    
    mat<-Seurat_object@meta.data %>% 
      select(paste0("Classify_",marker))
    
  }else{
    
    if(length(celltype)>1) {
      stop("Please enter a single cell type!\n")
    }
    
    alltypes<-levels(as.factor(Seurat_object$majortype))
    
    if(!(celltype %in% alltypes))
        return(0)
  
    if(tissue_region == "All"){
      
      if(length(marker) == 0)
        return(sum(Seurat_object$majortype == celltype))
      
      mat<-Seurat_object@meta.data %>% 
        filter(majortype == celltype) %>% 
        select(paste0("Classify_",marker))
      
    }else{
      if(!(tissue_region %in% c("Tumor","Stroma","Glass")))
        stop("Region should be either Tumor, Stroma, Glass, or All!")
      
      mat<-Seurat_object@meta.data %>% 
        filter(majortype == celltype) %>% 
        filter(region == tissue_region) 
      
      if(length(marker) == 0)
        return(nrow(mat))
      
      mat<-mat %>% select(paste0("Classify_",marker))
    }
  }
  sum(apply(mat,1,sum)==length(marker))
}


Calc_FuncSet_Composition<-function(Seurat_object, subtype_signature,tissue_region = "All"){
  majortype<-subtype_signature$majortype
  signature<-subtype_signature %>% select(-majortype, -subtype)
  cellcounts<-c()
  markers<-colnames(signature)
  
  for(k in 1:length(majortype)){
    count<-phenotyping_FunctionalSet(Seurat_object,
                                     markers[which(!is.na(signature[k,]))],
                                     majortype[k],tissue_region)
    cellcounts<-c(cellcounts,count)
  }
  
  return(cellcounts)
}


get_func_subset_ID<-function(Seurat_object,mtype, markers, tissue_region = "All"){
  ## CHECK
  for(m in markers){
    if(!(m %in% rownames(Seurat_object)))
      stop(paste("marker",m,"was not found in the list of markers !\n"))
  }
  
  if(length(markers) == 0 & length(mtype) == 0)
    stop("no input cell type nor marker!\n")
  
  if(!(tissue_region %in% c("Tumor","Stroma","Glass","All")))
    stop("Region should be either Tumor, Stroma, Glass, or All!")
  
  if(length(markers) == 0){
    sets<-Seurat_object@meta.data %>% 
      mutate(subset = if_else((majortype == mtype), 1, 0) ) %>%
      filter(subset == 1)
    
    if(tissue_region == "All"){
      ids<-rownames(sets)
    }else{
      ids<-sets %>% filter(region == tissue_region) %>% rownames()
    }
    
    return(ids)
  }
  
  marker_colnames<-paste0("Classify_",markers)
  
  if(length(mtype) == 0){
    sets<-Seurat_object@meta.data %>% 
      mutate(ind = rowSums(across(marker_colnames))) %>%  ## check if it is 1 for selected markers
      mutate(subset = if_else(ind == length(markers), 1, 0) ) %>% ## select subset when all markers is 1
      filter(subset == 1) 
    
    if(tissue_region == "All"){
      ids<-rownames(sets)
    }else{
      ids<-sets %>% filter(region == tissue_region) %>% rownames()
    }
      
    return(ids)
  }
  
  sets<-Seurat_object@meta.data %>% 
    mutate(ind = rowSums(across(marker_colnames))) %>%
    mutate(subset = if_else((majortype == mtype & ind == length(markers)), 1, 0) ) %>%
    filter(subset == 1)
  
  if(tissue_region == "All"){
    ids<-rownames(sets)
  }else{
    ids<-sets %>% filter(region == tissue_region) %>% rownames()
  }
  
  return(ids) 
}
