# Color by Metadata [scRNA-Seq][CCBR] (0fac593c-01e9-4c02-96d0-cc4876c4dfef): v155
ccbr1072_Fig3b <- function(AllCells_CellTypes_SO,AllCells_CellTypes_SampleNames,AllCells_CellTypes_MetadataTable) {

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    #nidapLoadPackages("SCWorkflow")
    library(SCWorkflow)
    library(grid)
    library(gridExtra)

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##
    
    #Basic Parameters:
    Seurat_Object <- AllCells_CellTypes_SO
    Sample_Table <- AllCells_CellTypes_SampleNames
    Metadata_Table <- AllCells_CellTypes_MetadataTable
    Samples_to_Include <- 'c("CD47KO_1","CD47KO_2","CD47KO_3","TSP1KO_1","TSP1KO_2","TSP1KO_3","WT_1","WT_2","WT_3")'
    Metadata_to_Plot <- 'c("mouseRNAseq_main","immgen_main")'

    #Advanced Parameters:
    Save_the_Entire_Dataset <- FALSE
    Use_CITE_seq <- FALSE

    #Visualization Parameters:
    Reduction_Type <- "tsne"
    Number_of_Columns_for_Final_Image <- 0
    Show_Labels <- FALSE
    Legend_Text_Size <- 1
    Legend_Position <- "right"
    Columns_to_Summarize <- c()
    Summarization_Cut_Off <- 5
    DPI <- 300
    #image: png
    imageType="png"
    Dot_Size <- 0.01

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
# auto removed:     fs <- Seurat_Object$fileSystem()
#path <- "./rds_output/var_AllCells_CellTypes_SO,AllCells_CellTypes_SampleNames,AllCells_CellTypes_MetadataTable.rds"
 #   SO <- readRDS(path)
    SO <- Seurat_Object
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

print("template_function_ccbr1072_Fig3b.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_AllCells_CellTypes_SO<-readRDS(paste0(rds_output,"/var_AllCells_CellTypes_SO.rds"))
Input_is_Seurat_count <- 1
#for(item in var_AllCells_CellTypes_SO){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_AllCells_CellTypes_SO<-as.data.frame(var_AllCells_CellTypes_SO)}else{var_AllCells_CellTypes_SO <- var_AllCells_CellTypes_SO}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
#var_AllCells_CellTypes_SampleNames<-readRDS(paste0(rds_output,"/var_AllCells_CellTypes_SampleNames.rds"))
var_AllCells_CellTypes_SampleNames<-NULL
Input_is_Seurat_count <- 0
#for(item in var_AllCells_CellTypes_SampleNames){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_AllCells_CellTypes_SampleNames<-as.data.frame(var_AllCells_CellTypes_SampleNames)}else{var_AllCells_CellTypes_SampleNames <- var_AllCells_CellTypes_SampleNames}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
#var_AllCells_CellTypes_MetadataTable<-readRDS(paste0(rds_output,"/var_AllCells_CellTypes_MetadataTable.rds"))
var_AllCells_CellTypes_MetadataTable<-NULL
Input_is_Seurat_count <- 0
#for(item in var_AllCells_CellTypes_MetadataTable){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_AllCells_CellTypes_MetadataTable<-as.data.frame(var_AllCells_CellTypes_MetadataTable)}else{var_AllCells_CellTypes_MetadataTable <- var_AllCells_CellTypes_MetadataTable}
invisible(graphics.off())
var_ccbr1072_Fig3b<-ccbr1072_Fig3b(var_AllCells_CellTypes_SO,var_AllCells_CellTypes_SampleNames,var_AllCells_CellTypes_MetadataTable)
invisible(graphics.off())
saveRDS(var_ccbr1072_Fig3b, paste0(rds_output,"/var_ccbr1072_Fig3b.rds"))
