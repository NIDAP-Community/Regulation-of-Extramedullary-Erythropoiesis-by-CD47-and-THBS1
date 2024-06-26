# DEG with Find Markers [scRNA-seq][CCBR] (54b6dd44-e233-4fb8-9bae-6ab0cb46e399): v93
R04C12_Filtered_DEG <- function(R04C12_Filtered_SO,R04C12_Filtered_SampleNames,R04C12_Filtered_MetadataTable) {

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    library("SCWorkflow")

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    seurat_object <- R04C12_Filtered_SO
    metadata_table <- R04C12_Filtered_MetadataTable
    samples <- 'c("CD47KO_1","CD47KO_2","CD47KO_3","TSP1KO_1","TSP1KO_2","TSP1KO_3","WT_1","WT_2","WT_3")'
    parameter_to_test <- "Treatment_Group"
    contrasts <- c("CD47KO-WT","TSP1KO-WT")
    test_to_use <- "MAST"
    log_fc_threshold <- 0
    use_spark <- FALSE
    assay_to_use <- "SCT"
    use_log_2 <- TRUE

    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##
    
    
    ## --------- ##
    ## Functions ##
    ## --------- ##

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

   ## Load SO 
        cat("Reading Seurat Object from dataset: seurat_object.rds\n\n")
# auto removed:         fs <- seurat_object$fileSystem()
path <- "./rds_output/var_R04C12_Filtered_SO,R04C12_Filtered_SampleNames,R04C12_Filtered_MetadataTable.rds"
        SO <- seurat_object
        print(SO)

results <- degGeneExpressionMarkers(object = SO,
                                     samples = samples,
                                     contrasts = contrasts,
                                     parameter.to.test = parameter_to_test,
                                     test.to.use = test_to_use,
                                     log.fc.threshold = log_fc_threshold,
                                     use.spark = use_spark,
                                     assay.to.use = assay_to_use)

    return(results$df)

}

print("template_function_R04C12_Filtered_DEG.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_R04C12_Filtered_SO<-readRDS(paste0(rds_output,"/var_R04C12_Filtered_SO.rds"))
Input_is_Seurat_count <- 1
#for(item in var_R04C12_Filtered_SO){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_R04C12_Filtered_SO<-as.data.frame(var_R04C12_Filtered_SO)}else{var_R04C12_Filtered_SO <- var_R04C12_Filtered_SO}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_R04C12_Filtered_SampleNames<-NULL #(paste0(rds_output,"/var_R04C12_Filtered_SampleNames.rds"))
Input_is_Seurat_count <- 0
#for(item in var_R04C12_Filtered_SampleNames){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_R04C12_Filtered_SampleNames<-as.data.frame(var_R04C12_Filtered_SampleNames)}else{var_R04C12_Filtered_SampleNames <- var_R04C12_Filtered_SampleNames}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_R04C12_Filtered_MetadataTable<-NULL #(paste0(rds_output,"/var_R04C12_Filtered_MetadataTable.rds"))
Input_is_Seurat_count <- 0
#for(item in var_R04C12_Filtered_MetadataTable){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_R04C12_Filtered_MetadataTable<-as.data.frame(var_R04C12_Filtered_MetadataTable)}else{var_R04C12_Filtered_MetadataTable <- var_R04C12_Filtered_MetadataTable}
invisible(graphics.off())
var_R04C12_Filtered_DEG<-R04C12_Filtered_DEG(var_R04C12_Filtered_SO,var_R04C12_Filtered_SampleNames,var_R04C12_Filtered_MetadataTable)
invisible(graphics.off())
saveRDS(var_R04C12_Filtered_DEG, paste0(rds_output,"/var_R04C12_Filtered_DEG.rds"))
