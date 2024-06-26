# Color by Gene [scRNA-Seq][CCBR] (82c9ec27-5cd8-4fd2-87a9-2a0f83e2ac3a): v103
ccbr1072_SuppFig7a <- function(AllCells_CellTypes_SO,AllCells_CellTypes_SampleNames) {
    
    
    ## --------- ##
    ## Libraries ##
    ## --------- ##
    
    library(nidapFunctions)
    #nidapLoadPackages("SCWorkflow")
    library("SCWorkflow")
    library(grid)
    library(gridExtra)
    library(scales)
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    # Primary Inputs:
    seurat_object <- AllCells_CellTypes_SO

    # Basic Parameters:
    samples_to_include <- 'c("CD47KO_1","CD47KO_2","CD47KO_3","TSP1KO_1","TSP1KO_2","TSP1KO_3","WT_1","WT_2","WT_3")'
    gene <- c("Gypa","Ermap","Klf1","Aqp1","Gata1")

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

colorByGene_mod <- function (object, samples.to.include, gene, reduction.type = "umap", 
    number.of.rows = 0, return.seurat.object = FALSE, color = "red", 
    point.size = 1, point.shape = 16, point.transparency = 0.5, 
    use.cite.seq.data = FALSE, assay = "SCT") 
{
    print(object)
    samples = eval(parse(text = gsub("\\[\\]", "c()", samples.to.include)))
    if (length(samples) == 0) {
        samples = unique(object@meta.data$orig.ident)
    }
    colnames(object@meta.data) <- gsub("orig_ident", "orig.ident", 
        colnames(object@meta.data))
    if ("active.ident" %in% slotNames(object)) {
        sample.name = as.factor(object@meta.data$orig.ident)
        names(sample.name) = names(object@active.ident)
        object@active.ident <- as.factor(vector())
        object@active.ident <- sample.name
        object.sub = subset(object, ident = samples)
    }
    else {
        sample.name = as.factor(object@meta.data$orig.ident)
        names(sample.name) = names(object@active.ident)
        object@active.ident <- as.factor(vector())
        object@active.ident <- sample.name
        object.sub = subset(object, ident = samples)
    }
    no.gene = gene[!gene %in% rownames(object.sub[[assay]]@scale.data)]
    if (!is.null(no.gene)) {
        print("Gene(s) missing from dataset:")
        print(no.gene)
    }
    gene = gene[gene %in% rownames(object.sub[[assay]]@scale.data)]
    if (length(gene) > 0) {
        .plotGene <- function(gene) {
            gene.mat = object.sub[[assay]]@scale.data[gene, ]
            gene.quant = quantile(gene.mat[gene.mat > 1], probs = c(0.1, 
                0.5, 0.9))
            gene.mat[gene.mat > gene.quant[3]] = gene.quant[3]
            gene.mat[gene.mat < gene.quant[1]] = 0
            if (!(use.cite.seq.data)) {
                if (reduction.type == "tsne") {
                  p1 <- DimPlot(object.sub, reduction = "tsne", 
                    group.by = "ident")
                  clus.mat = data.frame(umap1 = p1$data$tSNE_1, 
                    umap2 = p1$data$tSNE_2, gene = gene.mat)
                }
                else if (reduction.type == "umap") {
                  p1 <- DimPlot(object.sub, reduction = "umap", 
                    group.by = "ident")
                  clus.mat = data.frame(umap1 = p1$data$UMAP_1, 
                    umap2 = p1$data$UMAP_2, gene = gene.mat)
                }
                else {
                  p1 <- DimPlot(object.sub, reduction = "pca", 
                    group.by = "ident")
                  clus.mat = data.frame(umap1 = p1$data$PC_1, 
                    umap2 = p1$data$PC_2, gene = gene.mat)
                }
            }
            else {
                if (reduction.type == "tsne") {
                  p1 <- DimPlot(object.sub, reduction = "protein_tsne", 
                    group.by = "ident")
                  clus.mat = data.frame(umap1 = p1$data$protein_tsne_1, 
                    umap2 = p1$data$protein_tsne_2, gene = gene.mat)
                }
                else if (reduction.type == "umap") {
                  p1 <- DimPlot(object.sub, reduction = "protein_umap", 
                    group.by = "ident")
                  clus.mat = data.frame(umap1 = p1$data$protein_umap_1, 
                    umap2 = p1$data$protein_umap_2, gene = gene.mat)
                }
                else {
                  p1 <- DimPlot(object.sub, reduction = "protein_pca", 
                    group.by = "ident")
                  clus.mat = data.frame(umap1 = p1$data$protein_pca_1, 
                    umap2 = p1$data$protein_pca_2, gene = gene.mat)
                }
            }
            reduction.type.x <- paste0(reduction.type, "-1")
            reduction.type.y <- paste0(reduction.type, "-2")
            clus.mat <- clus.mat %>% dplyr::arrange(gene)
            g <- ggplot(clus.mat, aes(x = umap1, y = umap2)) + 
                theme_bw() + theme(legend.title = element_blank(), text = element_text(size = 18)) + ggtitle(gene) + geom_point(aes(colour = gene), 
                alpha = point.transparency, shape = point.shape, 
                size = point.size) + theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), panel.background = element_blank(), 
                legend.text = element_text(size = rel(0.5))) + 
                scale_color_gradient(limits = c(0, gene.quant[3]), 
                  low = "lightgrey", high = color) + xlab(reduction.type.x) + 
                ylab(reduction.type.y)
            return(g)
        }
        if (number.of.rows == 0) {
            n = ceiling(length(gene)^0.5)
        }
        else {
            n = number.of.rows
        }
        grob <- lapply(seq_along(gene), function(x) .plotGene(gene[x]))
        if (return.seurat.object) {
            result.list <- list(object = object, plot = grob)
            return(result.list)
        }
        else {
            gene = as.data.frame(gene)
            result.list <- list(object = gene, plot = grob)
            return(result.list)
        }
    }
    else {
        print("No genes found in dataset")
# auto removed:         return(NULL)
    }
}

ColorByGene.result <- colorByGene_mod(object = SO,
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

print("template_function_ccbr1072_SuppFig7a.R #########################################################################")
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
var_ccbr1072_SuppFig7a<-ccbr1072_SuppFig7a(var_AllCells_CellTypes_SO,var_AllCells_CellTypes_SampleNames)
invisible(graphics.off())
saveRDS(var_ccbr1072_SuppFig7a, paste0(rds_output,"/var_ccbr1072_SuppFig7a.rds"))
