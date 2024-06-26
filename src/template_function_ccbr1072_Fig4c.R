# Color by Gene [scRNA-Seq][CCBR] (82c9ec27-5cd8-4fd2-87a9-2a0f83e2ac3a): v103
ccbr1072_Fig4c <- function(R04C12andC14_Filtered_SO,R04C12andC14_Filtered_SampleNames) {
    
    
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
    seurat_object <- R04C12andC14_Filtered_SO

    # Basic Parameters:
    samples_to_include <- 'c("CD47KO_1","CD47KO_2","CD47KO_3","TSP1KO_1","TSP1KO_2","TSP1KO_3","WT_1","WT_2","WT_3")'
    gene <- c("Trim10","Gypa","Epb42","Spta1","Sptb","Xpo1")

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
path <- "./rds_output/var_R04C12andC14_Filtered_SO,R04C12andC14_Filtered_SampleNames.rds"
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

print("template_function_ccbr1072_Fig4c.R #########################################################################")
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
invisible(graphics.off())
var_ccbr1072_Fig4c<-ccbr1072_Fig4c(var_R04C12andC14_Filtered_SO,var_R04C12andC14_Filtered_SampleNames)
invisible(graphics.off())
saveRDS(var_ccbr1072_Fig4c, paste0(rds_output,"/var_ccbr1072_Fig4c.rds"))
