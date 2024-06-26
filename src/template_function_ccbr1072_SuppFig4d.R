# Color by Gene [scRNA-Seq][CCBR] (82c9ec27-5cd8-4fd2-87a9-2a0f83e2ac3a): v103
ccbr1072_SuppFig4d <- function(AllCells_CellTypes_SO,AllCells_CellTypes_SampleNames) {
    
    
    ## --------- ##
    ## Libraries ##
    ## --------- ##
    
    library(nidapFunctions)
    #nidapLoadPackages("SCWorkflow")
    library("SCWorkflow")
    library(grid)
    library(gridExtra)
    library(scales)

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    # Primary Inputs:
    seurat_object <- AllCells_CellTypes_SO

    # Basic Parameters:
    samples_to_include <- 'c("CD47KO_1","CD47KO_2","CD47KO_3","TSP1KO_1","TSP1KO_2","TSP1KO_3","WT_1","WT_2","WT_3")'
    gene <- c("Ranbp1","Ranbp2","Spta1","Sptb","Tal1","Tfrc")

    # Visualization Parameters:
    reduction_type <- "tsne"
    number_of_rows <- 2
    color <- "red"
    point_size <- 1
    point_shape <- 16
    point_transparency <- 0.5
    image_type <- "png"

    # Advanced Parameters:
    save_the_entire_dataset <- FALSE
    use_cite_seq_data <- FALSE
    assay <- "SCT"

    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##
    
    
    ## --------- ##
    ## Functions ##
    ## --------- ##

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

    # Loading Seurat Object
    cat("Reading Seurat Object from dataset: seurat_object.rds\n\n")
# auto removed:     fs <- seurat_object$fileSystem()
#path <- "./rds_output/var_AllCells_CellTypes_SO,AllCells_CellTypes_SampleNames.rds"
    #SO <- readRDS(path)
    SO <- seurat_object
    print(SO)

ColorByGene.result <- colorByGene(object = SO,
                                    samples.to.include = samples_to_include,
                                    gene = gene,
                                    reduction.type = reduction_type,
                                    number.of.rows = number_of_rows,
                                    return.seurat.object = save_the_entire_dataset,
                                    color = color,
                                    point.size = point_size,
                                    point.shape = point_shape,
                                    point.transparency = point_transparency,
                                    use.cite.seq.data = use_cite_seq_data,
                                    assay = assay)

# Preparing Graphic Output
    if (number_of_rows == 0) {
        n = ceiling(length(ColorByGene.result$plot)^0.5)
    } else {
        n = number_of_rows
    }

do.call("grid.arrange", c(ColorByGene.result$plot, nrow=n))

# Saving dataset (if requested)
if(save_the_entire_dataset){
    return(ColorByGene.result$object)
  }
  else{
    gene = as.data.frame(gene) 
    return(gene)
  }
}

print("template_function_ccbr1072_SuppFig4d.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_AllCells_CellTypes_SO<-readRDS(paste0(rds_output,"/var_AllCells_CellTypes_SO.rds"))
Input_is_Seurat_count <- 1
##for(item in var_AllCells_CellTypes_SO){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_AllCells_CellTypes_SO<-as.data.frame(var_AllCells_CellTypes_SO)}else{var_AllCells_CellTypes_SO <- var_AllCells_CellTypes_SO}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
#var_AllCells_CellTypes_SampleNames<-readRDS(paste0(rds_output,"/var_AllCells_CellTypes_SampleNames.rds"))
var_AllCells_CellTypes_SampleNames <- NULL
Input_is_Seurat_count <- 0
##for(item in var_AllCells_CellTypes_SampleNames){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_AllCells_CellTypes_SampleNames<-as.data.frame(var_AllCells_CellTypes_SampleNames)}else{var_AllCells_CellTypes_SampleNames <- var_AllCells_CellTypes_SampleNames}
invisible(graphics.off())
var_ccbr1072_SuppFig4d<-ccbr1072_SuppFig4d(var_AllCells_CellTypes_SO,var_AllCells_CellTypes_SampleNames)
invisible(graphics.off())
saveRDS(var_ccbr1072_SuppFig4d, paste0(rds_output,"/var_ccbr1072_SuppFig4d.rds"))
