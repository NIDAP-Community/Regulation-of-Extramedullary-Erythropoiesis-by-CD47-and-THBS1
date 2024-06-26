# Annotate Cell Types with SingleR [scRNA-seq][CCBR] (ba316b6a-e424-49ab-a8af-32af8c85529f): v125
RBCProg_ReclusteredCellTypes_SO <- function(RBCProg_Reclustered_SO,CellDex_snapshotDate_2021_10_19, RBCProg_Reclustered_MetadataTable) {

## --------- ##
## Libraries ##
## --------- ##

    library(nidapFunctions)
    library("SCWorkflow")
    library(grid)
    library(gridExtra)
    library(gridBase)
    library(cowplot)
    

## -------------------------------- ##
## User-Defined Template Parameters ##
## -------------------------------- ##

# Primary Inputs:
CombinedSO <-RBCProg_Reclustered_SO
Metadata <- RBCProg_Reclustered_MetadataTable
RDS <- CellDex_snapshotDate_2021_10_19

# Basic Parameters:
species <- "Mouse"

# Visualization Parameters:
Reduction <- "umap"
Legend_Dot_Size <- 2
imageType = "png"

# Advanced Parameters:
doFineTuning <- FALSE
useSpark <-FALSE
Number_of_cells_per_partition <- 400
useClusters <- FALSE
ClusterColumn <- c()[1]

## -------------------------------- ##
## Errors                           ##
## -------------------------------- ##

## -------------------------------- ##
## Functions                        ##
## -------------------------------- ##

annotateCellTypes <- function(object,
                              species = "Mouse",
                              reduction.type = "umap",
                              legend.dot.size = 2,
                              do.finetuning = FALSE,
                              local.celldex = NULL,
                              use.clusters = NULL) {
  ## -------------------------------- ##
  ## Functions                        ##
  ## -------------------------------- ##
  
  library(Seurat)
  library(cowplot)
  library(ggplot2)
  library(RColorBrewer)
  library(SingleR)
  library(celldex)
  
  
  .annotations <- function(so) {
    so.counts = GetAssayData(object = so)[, colnames(x = so)]
    if (species == "Human") {
      
      #HPCA block
      if (!is.null(local.celldex)) {
        HPCA <- local.celldex[[1]]
      } else {
        HPCA <- celldex::HumanPrimaryCellAtlasData(ensembl = FALSE)
      }
      
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = clusterList,
        fine.tune = do.finetuning,
        ref = HPCA,
        labels = HPCA$label.main
      )
      if(is.null(use.clusters))
      {
        so[["HPCA_main"]] <- 
          singler$labels[match(rownames(so[[]]), rownames(singler))]
      } else {
        so[["HPCA_main"]] <- 
          singler$labels[match(so[[]][[use.clusters]], rownames(singler))]  
      }
      
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = clusterList,
        fine.tune = do.finetuning,
        ref = HPCA,
        labels = HPCA$label.fine
      )
      if(is.null(use.clusters))
      {
        so[["HPCA"]] <-
          singler$labels[match(rownames(so[[]]), rownames(singler))]
      } else {
        so[["HPCA"]] <- 
          singler$labels[match(so[[]][[use.clusters]], rownames(singler))]  
      }
      
      #BP_encode block
      if (!is.null(local.celldex)) {
        BP <- local.celldex[[2]]
      } else {
        BP <- celldex::BlueprintEncodeData(ensembl = FALSE)
      }
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = clusterList,
        fine.tune = do.finetuning,
        ref = BP,
        labels = BP$label.main
      )
      if(is.null(use.clusters))
      {
        so[["BP_encode_main"]] <-
          singler$labels[match(rownames(so[[]]), rownames(singler))]
      } else {
        so[["BP_encode_main"]] <- 
          singler$labels[match(so[[]][[use.clusters]], rownames(singler))]  
      }
      
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = clusterList,
        fine.tune = do.finetuning,
        ref = BP,
        labels = BP$label.fine
      )
      if(is.null(use.clusters))
      {
        so[["BP_encode"]] <-
          singler$labels[match(rownames(so[[]]), rownames(singler))]
      } else {
        so[["BP_encode"]] <- 
          singler$labels[match(so[[]][[use.clusters]], rownames(singler))]  
      }      
    }
    if (species == "Mouse") {
      
      #mouseRNAseq block
      if (!is.null(local.celldex)) {
        mousernaseq <- local.celldex[[3]]
      } else {
        mousernaseq <- celldex::MouseRNAseqData(ensembl = FALSE)
      }
      
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = clusterList,
        fine.tune = do.finetuning,
        ref = mousernaseq,
        labels = mousernaseq$label.main
      )
      if(is.null(use.clusters))
      {
        so[["mouseRNAseq_main"]] <-
          singler$labels[match(rownames(so[[]]), rownames(singler))]
      } else {
        so[["mouseRNAseq_main"]] <- 
          singler$labels[match(so[[]][[use.clusters]], rownames(singler))]  
      }   
      
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = clusterList,
        fine.tune = do.finetuning,
        ref = mousernaseq,
        labels = mousernaseq$label.fine
      )
      if(is.null(use.clusters))
      {
        so[["mouseRNAseq"]] <-
          singler$labels[match(rownames(so[[]]), rownames(singler))]
      } else {
        so[["mouseRNAseq"]] <- 
          singler$labels[match(so[[]][[use.clusters]], rownames(singler))]  
      }       
      
      
      #ImmGen block
      if (!is.null(local.celldex)) {
        immgen <- local.celldex[[4]]
      } else {
        immgen <- celldex::ImmGenData(ensembl = FALSE)
      }
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = clusterList,
        fine.tune = do.finetuning,
        ref = immgen,
        labels = immgen$label.main
      )
      if(is.null(use.clusters))
      {
        so[["immgen_main"]] <- 
          singler$labels[match(rownames(so[[]]), rownames(singler))]
      } else {
        so[["immgen_main"]] <- 
          singler$labels[match(so[[]][[use.clusters]], rownames(singler))]  
      }     
      
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = clusterList,
        fine.tune = do.finetuning,
        ref = immgen,
        labels = immgen$label.fine
      )
      if(is.null(use.clusters))
      {
        so[["immgen"]] <- 
          singler$labels[match(rownames(so[[]]), rownames(singler))]
      } else {
        so[["immgen"]] <- 
          singler$labels[match(so[[]][[use.clusters]], rownames(singler))]  
      }     
      
    }
    return(so)
  }
  
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  # Getting Cluster List from MetaData of Seurat Object:
  if(!is.null(use.clusters))
  {
    clusterList <- object@meta.data[[use.clusters]]
  } else {
    clusterList <- NULL
  }
  
  # Running Annotation:
  object <- .annotations(object)
  print("done")
  
  # Assigning Colors (choice of datasets depends on species):
  if (species == "Human") {
    numColors = max(length(unique(object@meta.data$BP_encode_main)), length(unique(object@meta.data$HPCA_main)))
  } else {
    numColors = max(length(unique(object@meta.data$mouseRNAseq_main)), length(unique(object@meta.data$immgen_main)))
  }
  colpaired = colorRampPalette(brewer.pal(12, "Paired"))
  cols = c(
    "#e6194B",
    "#3cb44b",
    "#4363d8",
    "#f58231",
    "#911eb4",
    "#42d4f4",
    "#f032e6",
    "#bfef45",
    "#fabebe",
    "#469990",
    "#e6beff",
    "#9A6324",
    "#800000",
    "#aaffc3",
    "#808000",
    "#000075",
    colpaired(numColors)
  )
  
  # Creating plots (choice of datasets depends on species):
  if (species == "Human") {
    p1 = DimPlot(object, reduction = reduction.type, group.by = "HPCA_main") +
      scale_color_manual(values = cols) +
      theme(legend.position = "top") +
      guides(override.aes = list(size = legend.dot.size),
             colour = guide_legend(ncol = 4)) + ggtitle("HPCA Main Cell Type Annotations")
    p2 = DimPlot(object, reduction = reduction.type, group.by="BP_encode_main") +
      scale_color_manual(values = cols) +
      theme(legend.position = "top") +
      guides(override.aes = list(size = legend.dot.size),
             colour = guide_legend(ncol = 4)) + ggtitle("BP Encode Main Cell Type Annotations")
  } else {
    p1 = DimPlot(object, reduction = reduction.type, group.by = "immgen_main") +
      scale_color_manual(values = cols) +
      theme(legend.position = "top") +
      guides(override.aes = list(size = legend.dot.size),
             colour = guide_legend(ncol = 4)) + ggtitle("Immgen Main Cell Type Annotations")
    p2 = DimPlot(object, reduction=reduction.type, group.by="mouseRNAseq_main") +
      scale_color_manual(values = cols) +
      theme(legend.position = "top") +
      guides(override.aes = list(size = legend.dot.size),
             colour = guide_legend(ncol = 4)) + ggtitle("Mouse RNAseq Main Cell Type Annotations")
  }
  
  slot(object, "commands") <- list()
  #  cat("\nSingleR Object Checksum:\n")
  #  print(digest::digest(object))
  
  # Returning Seurat Object and 2 plots:
  return(list(
    "object" = object,
    "p1" = p1,
    "p2" = p2
  ))
}


## --------------- ##
## Main Code Block ##
## --------------- ##

# Loading Seurat object
    cat("1. Reading Seurat Object from dataset: seurat_object.rds\n\n")
# auto removed:     fs <- CombinedSO$fileSystem()
#path <- "./rds_output/var_NA.rds"
    so <- CombinedSO

# Loading CellDex annotation
#cell.dex <- list(HPCA, BP, mousernaseq, immgen)
    cat("2. Reading pre-saved CellDex Data from dataset: CellDexDatabase.rds\n\n")
# auto removed:     fs2 <- RDS$fileSystem()
#    path2 <- fs2$get_path("CellDexDatabase.rds", 'r')
#    CellDexData <- readRDS(path2)
    CellDexData <- CellDex_snapshotDate_2021_10_19

if (useClusters)
  {
      anno <- annotateCellTypes(object = so,
                            species = species,
                            reduction.type = Reduction,
                            legend.dot.size = Legend_Dot_Size,
                            do.finetuning = doFineTuning,
                            local.celldex = CellDexData,
                            use.clusters = ClusterColumn)
  } else {
      anno <- annotateCellTypes(object = so,
                            species = species,
                            reduction.type = Reduction,
                            legend.dot.size = Legend_Dot_Size,
                            do.finetuning = doFineTuning,
                            local.celldex = CellDexData)

  }
    gc()

## Print figures
    print(anno$p1)
    print(anno$p2)

## Save the annotated Seurat object
# auto removed:     output <- new.output()
# auto removed:     output_fs <- output$fileSystem()
return(anno$object)

# auto removed:     return(NULL)
}

print("template_function_RBCProg_ReclusteredCellTypes_SO.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_RBCProg_Reclustered_SO<-readRDS(paste0(rds_output,"/var_RBCProg_Reclustered_SO.rds"))
Input_is_Seurat_count <- 1
#for(item in var_RBCProg_Reclustered_SO){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_RBCProg_Reclustered_SO<-as.data.frame(var_RBCProg_Reclustered_SO)}else{var_RBCProg_Reclustered_SO <- var_RBCProg_Reclustered_SO}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_CellDex_snapshotDate_2021_10_19 <- readRDS("./nidap_downloads/CellDexDatabase.rds")
Input_is_Seurat_count <- 1
#for(item in var_CellDex_snapshotDate_2021_10_19){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_CellDex_snapshotDate_2021_10_19<-as.data.frame(var_CellDex_snapshotDate_2021_10_19)}else{var_CellDex_snapshotDate_2021_10_19 <- var_CellDex_snapshotDate_2021_10_19}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_RBCProg_Reclustered_MetadataTable<-NULL #(paste0(rds_output,"/var_RBCProg_Reclustered_MetadataTable.rds"))
Input_is_Seurat_count <- 0
#for(item in var_RBCProg_Reclustered_MetadataTable){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_RBCProg_Reclustered_MetadataTable<-as.data.frame(var_RBCProg_Reclustered_MetadataTable)}else{var_RBCProg_Reclustered_MetadataTable <- var_RBCProg_Reclustered_MetadataTable}
invisible(graphics.off())
var_RBCProg_ReclusteredCellTypes_SO<-RBCProg_ReclusteredCellTypes_SO(var_RBCProg_Reclustered_SO,var_CellDex_snapshotDate_2021_10_19,var_RBCProg_Reclustered_MetadataTable)
invisible(graphics.off())
saveRDS(var_RBCProg_ReclusteredCellTypes_SO, paste0(rds_output,"/var_RBCProg_ReclusteredCellTypes_SO.rds"))
