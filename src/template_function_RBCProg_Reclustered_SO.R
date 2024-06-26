# Recluster Seurat Object [scRNA-seq][CCBR] (576fe688-c445-48c6-a40a-033da171c149): v19
RBCProg_Reclustered_SO <- function(ccbr1072_Fig6a, RBCProg_Filtered_MetadataTable) {
    #image: png

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    library("SCWorkflow")
    
    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##
    
    #Basic Parameters:
    seurat_object = ccbr1072_Fig6a
    meta <- RBCProg_Filtered_MetadataTable
    reduction <- "umap"

    #Old Clustering Column Parameters:
    old_columns_to_save <- c("SCT_snn_res_0_2","SCT_snn_res_0_4","SCT_snn_res_0_6","SCT_snn_res_0_8","SCT_snn_res_1","SCT_snn_res_1_2")
    prepend_text = "old"

    #Reclustering Parameters:
    number_of_pcs = 30
    cluster_resolution_low_range <- 0.2
    cluster_resolution_high_range <- 1.2
    cluster_resolution_range_bins <- 0.2
    
    #Advanced Parameters:
    seurat_object_filename <- "seurat_object.rds"

    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##
    
    ## --------- ##
    ## Functions ##
    ## --------- ##

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

    ## Input SO.
    cat("Reading Seurat Object from dataset: seurat_object.rds\n\n")
    # path <- nidapGetPath(seurat_object,seurat_object_filename)
    SO <- seurat_object

    ## Recluster the SO.
    reclustered_SO <- reclusterSeuratObject(object = SO,
                                            prepend.txt = prepend_text,
                                            old.columns.to.save = old_columns_to_save,
                                            number.of.pcs = number_of_pcs,
                                            cluster.resolution.low.range = cluster_resolution_low_range,
                                            cluster.resolution.high.range = cluster_resolution_high_range,
                                            cluster.resolution.range.bins = cluster_resolution_range_bins,
                                            reduction.type = reduction
                                            )

    ## Print reclustered dimensionality reduction plot.
    print(reclustered_SO$plot)

    ## Output reclustered SO.
# auto removed:     output <- new.output()
# auto removed:     output_fs <- output$fileSystem()
return(reclustered_SO$object)
    return(output_fs)
}

print("template_function_RBCProg_Reclustered_SO.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ccbr1072_Fig6a<-readRDS(paste0(rds_output,"/var_ccbr1072_Fig6a.rds"))
Input_is_Seurat_count <- 1
#for(item in var_ccbr1072_Fig6a){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_ccbr1072_Fig6a<-as.data.frame(var_ccbr1072_Fig6a)}else{var_ccbr1072_Fig6a <- var_ccbr1072_Fig6a}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_RBCProg_Filtered_MetadataTable<-NULL #(paste0(rds_output,"/var_RBCProg_Filtered_MetadataTable.rds"))
Input_is_Seurat_count <- 0
#for(item in var_RBCProg_Filtered_MetadataTable){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_RBCProg_Filtered_MetadataTable<-as.data.frame(var_RBCProg_Filtered_MetadataTable)}else{var_RBCProg_Filtered_MetadataTable <- var_RBCProg_Filtered_MetadataTable}
invisible(graphics.off())
var_RBCProg_Reclustered_SO<-RBCProg_Reclustered_SO(var_ccbr1072_Fig6a,var_RBCProg_Filtered_MetadataTable)
invisible(graphics.off())
saveRDS(var_RBCProg_Reclustered_SO, paste0(rds_output,"/var_RBCProg_Reclustered_SO.rds"))
