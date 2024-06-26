# Color by Gene [scRNA-Seq][CCBR] (82c9ec27-5cd8-4fd2-87a9-2a0f83e2ac3a): v103
ccbr1072_Fig7b <- function(RBCProg_ReclusteredCellTypes_SO,RBCProg_ReclusteredCellTypes_SampleNames) {
    
    
    ## --------- ##
    ## Libraries ##
    ## --------- ##
    
    library(nidapFunctions)
    library("SCWorkflow")
    library(grid)
    library(gridExtra)
    library(scales)

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    # Primary Inputs:
    seurat_object <- RBCProg_ReclusteredCellTypes_SO

    # Basic Parameters:
    samples_to_include <- 'c("CD47KO_1","CD47KO_2","CD47KO_3","TSP1KO_1","TSP1KO_2","TSP1KO_3","WT_1","WT_2","WT_3")'
    gene <- c("Klf1","Gata1","Ermap","Aqp1","Tfrc","Tmem56")

    # Visualization Parameters:
    reduction_type <- "umap"
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
path <- "./rds_output/var_RBCProg_ReclusteredCellTypes_SO,RBCProg_ReclusteredCellTypes_SampleNames.rds"
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

print("template_function_ccbr1072_Fig7b.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_RBCProg_ReclusteredCellTypes_SO<-readRDS(paste0(rds_output,"/var_RBCProg_ReclusteredCellTypes_SO.rds"))
Input_is_Seurat_count <- 1
#for(item in var_RBCProg_ReclusteredCellTypes_SO){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_RBCProg_ReclusteredCellTypes_SO<-as.data.frame(var_RBCProg_ReclusteredCellTypes_SO)}else{var_RBCProg_ReclusteredCellTypes_SO <- var_RBCProg_ReclusteredCellTypes_SO}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_RBCProg_ReclusteredCellTypes_SampleNames<-NULL #(paste0(rds_output,"/var_RBCProg_ReclusteredCellTypes_SampleNames.rds"))
Input_is_Seurat_count <- 0
#for(item in var_RBCProg_ReclusteredCellTypes_SampleNames){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_RBCProg_ReclusteredCellTypes_SampleNames<-as.data.frame(var_RBCProg_ReclusteredCellTypes_SampleNames)}else{var_RBCProg_ReclusteredCellTypes_SampleNames <- var_RBCProg_ReclusteredCellTypes_SampleNames}
invisible(graphics.off())
var_ccbr1072_Fig7b<-ccbr1072_Fig7b(var_RBCProg_ReclusteredCellTypes_SO,var_RBCProg_ReclusteredCellTypes_SampleNames)
invisible(graphics.off())
saveRDS(var_ccbr1072_Fig7b, paste0(rds_output,"/var_ccbr1072_Fig7b.rds"))
