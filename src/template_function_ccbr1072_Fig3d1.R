# Color by Metadata [scRNA-Seq][CCBR] (0fac593c-01e9-4c02-96d0-cc4876c4dfef): v155
ccbr1072_Fig3d1 <- function(R04C12andC14_Filtered_SO,R04C12andC14_Filtered_SampleNames,R04C12andC14_Filtered_MetadataTable) {

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    library("SCWorkflow")
    library(grid)
    library(gridExtra)

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##
    
    #Basic Parameters:
    seurat_object <- R04C12andC14_Filtered_SO
    Sample_Table <- R04C12andC14_Filtered_SampleNames ### Double-check!!! Not in use??
    Metadata_Table <- R04C12andC14_Filtered_MetadataTable
    Samples_to_Include <- 'c("CD47KO_1","CD47KO_2","CD47KO_3","TSP1KO_1","TSP1KO_2","TSP1KO_3","WT_1","WT_2","WT_3")'
    Metadata_to_Plot <- 'c("SCT_snn_res_0_4")'

    #Advanced Parameters:
    Save_the_Entire_Dataset <- FALSE
    Use_CITE_seq <- FALSE

    #Visualization Parameters:
    Reduction_Type <- "tsne"
    Number_of_Columns_for_Final_Image <- 0
    Show_Labels <- FALSE
    Legend_Text_Size <- 1
    Legend_Position <- "top"
    Columns_to_Summarize <- c()
    Summarization_Cut_Off <- 5
    DPI <- 300
    #image: png
    imageType="png"
    Dot_Size <- 1

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
# auto removed:     fs <- seurat_object$fileSystem()
path <- "./rds_output/var_R04C12andC14_Filtered_SO,R04C12andC14_Filtered_SampleNames,R04C12andC14_Filtered_MetadataTable.rds"
    SO <- seurat_object
    print(SO)

results <- plotMetadata(object = SO,
  samples.to.include = Samples_to_Include,
  metadata.to.plot = Metadata_to_Plot,
  columns.to.summarize = Columns_to_Summarize,
  summarization.cut.off = Summarization_Cut_Off,
  reduction.type = Reduction_Type,
  use.cite.seq = Use_CITE_seq,
  show.labels = Show_Labels,
  legend.text.size = Legend_Text_Size,
  legend.position = Legend_Position,
  dot.size = Dot_Size
  ) 

## Print Graphic output

matched_commas <- gregexpr(",", Metadata_to_Plot, fixed = TRUE)
n_commas <- length(matched_commas[[1]])
    if (Number_of_Columns_for_Final_Image == 0) {
        n = ceiling((n_commas+1)^0.5)
    } else {
        n = Number_of_Columns_for_Final_Image
    }
 do.call("grid.arrange", c(results$plot, ncol=n))

## Save dataset if requested

if (Save_the_Entire_Dataset){
# auto removed:     output <- new.output()
# auto removed:     output_fs <- output$fileSystem()
return(results$object)
# auto removed:     return(NULL)
    } else {
    return(results$object@meta.data)
    }
}

print("template_function_ccbr1072_Fig3d1.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_R04C12andC14_Filtered_SO<-readRDS(paste0(rds_output,"/var_R04C12andC14_Filtered_SO.rds"))
Input_is_Seurat_count <- 1
#for(item in var_R04C12andC14_Filtered_SO){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_R04C12andC14_Filtered_SO<-as.data.frame(var_R04C12andC14_Filtered_SO)}else{var_R04C12andC14_Filtered_SO <- var_R04C12andC14_Filtered_SO}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_R04C12andC14_Filtered_SampleNames<-NULL #(paste0(rds_output,"/var_R04C12andC14_Filtered_SampleNames.rds"))
Input_is_Seurat_count <- 0
#for(item in var_R04C12andC14_Filtered_SampleNames){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_R04C12andC14_Filtered_SampleNames<-as.data.frame(var_R04C12andC14_Filtered_SampleNames)}else{var_R04C12andC14_Filtered_SampleNames <- var_R04C12andC14_Filtered_SampleNames}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_R04C12andC14_Filtered_MetadataTable<-NULL #(paste0(rds_output,"/var_R04C12andC14_Filtered_MetadataTable.rds"))
Input_is_Seurat_count <- 0
#for(item in var_R04C12andC14_Filtered_MetadataTable){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_R04C12andC14_Filtered_MetadataTable<-as.data.frame(var_R04C12andC14_Filtered_MetadataTable)}else{var_R04C12andC14_Filtered_MetadataTable <- var_R04C12andC14_Filtered_MetadataTable}
invisible(graphics.off())
var_ccbr1072_Fig3d1<-ccbr1072_Fig3d1(var_R04C12andC14_Filtered_SO,var_R04C12andC14_Filtered_SampleNames,var_R04C12andC14_Filtered_MetadataTable)
invisible(graphics.off())
saveRDS(var_ccbr1072_Fig3d1, paste0(rds_output,"/var_ccbr1072_Fig3d1.rds"))
