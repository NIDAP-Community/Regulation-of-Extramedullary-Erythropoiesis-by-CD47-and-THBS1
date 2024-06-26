# Color by Gene Lists [scRNA-seq][CCBR] (d71ed4e6-a25d-4f66-a186-27c00a50a703): v104
AllCells_ModScoreRBCProgs_CBGL <- function(AllCells_CellTypes_SO,AllCells_CellTypes_SampleNames, ccbr1072_RBC_Progenitors_Markers) {
    
    #image: png
    
    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    #nidapLoadPackages("SCWorkflow")
    library(SCWorkflow)

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    #Input Parameters:
    #seurat_object <- AllCells_CellTypes_SO
    so <- AllCells_CellTypes_SO
    marker_list = ccbr1072_RBC_Progenitors_Markers

    #Basic Parameters:
    samples_to_include = eval(parse(text=gsub('\\[\\]','c()','c("CD47KO_1","CD47KO_2","CD47KO_3","TSP1KO_1","TSP1KO_2","TSP1KO_3","WT_1","WT_2","WT_3")')))
    sample_to_display <- c("CD47KO_1","CD47KO_2","CD47KO_3","TSP1KO_1","TSP1KO_2","TSP1KO_3","WT_1","WT_2","WT_3")
    cells_of_interest <- c("RBC_Progenitor")
    protein_presence <- FALSE
    save_the_entire_dataset <- FALSE

    #Plot Parameters:
    assay_to_plot <- "SCT"
    reduction_type = "tsne"
    point_transparency <- 0.5
    point_shape <- 16
    cite_seq <- FALSE

    #Advanced Parameters:
    seurat_object_filename <- "seurat_object.rds"

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

    #path <- nidapGetPath(seurat_object,seurat_object_filename)
    #so <- readRDS(path)
    
    colorByMarkerTable(object = so,
                        samples.subset = samples_to_include,
                        samples.to.display = sample_to_display,
                        marker.table = marker_list,
                        cells.of.interest = cells_of_interest,
                        protein.presence = protein_presence,
                        assay = assay_to_plot,
                        reduction.type = reduction_type,
                        point.transparency = point_transparency,
                        point.shape = point_shape,
                        cite.seq = FALSE)

# auto removed: return(NULL)
}

#################################################
## Global imports and functions included below ##
#################################################

colorByMarkerTable <- function (object, samples.subset, samples.to.display, marker.table, 
    cells.of.interest, protein.presence = FALSE, assay = "SCT", 
    reduction.type = "umap", point.transparency = 0.5, point.shape = 16, 
    cite.seq = FALSE){ 

        library(ggplot2)
        library(Seurat)
        
        .plotMarkers <- function(markers) {
            if (is.na(markers) == TRUE) {
                g <- ggplot() + theme_void()
                return(g)
            } else {
                markers.mat = object.sub[[assay]]@scale.data[markers, 
                    ]
                markers.quant = quantile(markers.mat[markers.mat > 
                    1], probs = c(0.1, 0.5, 0.9))
                markers.mat[markers.mat > markers.quant[3]] = markers.quant[3]
                markers.mat[markers.mat < markers.quant[1]] = 0
                if (!(cite.seq)) {
                    if (reduction.type == "tsne") {
                    p1 <- DimPlot(object.sub, reduction = "tsne", 
                        group.by = "ident")
                    clusmat = data.frame(umap1 = p1$data$tSNE_1, 
                        umap2 = p1$data$tSNE_2, markers = markers.mat, 
                        ident = as.factor(p1$data$ident))
                    } else if (reduction.type == "umap") {
                    p1 <- DimPlot(object.sub, reduction = "umap", 
                        group.by = "ident")
                    clusmat = data.frame(umap1 = p1$data$UMAP_1, 
                        umap2 = p1$data$UMAP_2, markers = markers.mat, 
                        ident = as.factor(p1$data$ident))
                    } else {
                    p1 <- DimPlot(object.sub, reduction = "pca", 
                        group.by = "ident")
                    clusmat = data.frame(umap1 = p1$data$PC_1, 
                        umap2 = p1$data$PC_2, markers = markers.mat, 
                        ident = as.factor(p1$data$ident))
                    }
                } else {
                    if (reduction.type == "tsne") {
                    p1 <- DimPlot(object.sub, reduction = "protein_tsne", 
                        group.by = "ident")
                    clusmat = data.frame(umap1 = p1$data$protein_tsne_1, 
                        umap2 = p1$data$protein_tsne_2, markers = markers.mat, 
                        ident = as.factor(p1$data$ident))
                    } else if (reduction.type == "umap") {
                    p1 <- DimPlot(object.sub, reduction = "protein_umap", 
                        group.by = "ident")
                    clusmat = data.frame(umap1 = p1$data$protein_umap_1, 
                        umap2 = p1$data$protein_umap_2, markers = markers.mat, 
                        ident = as.factor(p1$data$ident))
                    } else {
                    p1 <- DimPlot(object.sub, reduction = "protein_pca", 
                        group.by = "ident")
                    clusmat = data.frame(umap1 = p1$data$protein_pca_1, 
                        umap2 = p1$data$protein_pca_2, markers = markers.mat, 
                        ident = as.factor(p1$data$ident))
                    }
                }
                samples.caption <- paste(samples.to.display, sep = "", 
                    collapse = "\n")
                final_caption <- paste("Samples Displayed: ", samples.caption, 
                    sep = "", collapse = "\n")
                clusmat <- dplyr::mutate(clusmat, sample.markers = clusmat$markers * 
                    grepl(paste(samples.to.display, collapse = "|"), 
                    clusmat$ident))
                clusmat <- dplyr::arrange(clusmat,sample.markers)
                if (grepl("_neg", markers) == TRUE) {
                    clusmat <- dplyr::arrange(clusmat, desc(sample.markers))
                    g <- ggplot(clusmat, aes(x = umap1, y = umap2, 
                    group = ident)) + theme_bw() + theme(legend.title = element_blank()) + 
                    ggtitle(markers) + geom_point(aes(color = sample.markers, 
                    shape = ident), alpha = point.transparency, 
                    shape = point.shape, size = 1) + theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), panel.background = element_blank(), 
                    legend.text = element_text(size = rel(0.5))) + 
                    scale_color_gradient(limits = c(0, markers.quant[3]), 
                        low = "lightgrey", high = "red") + xlab("umap-1") + 
                    ylab("umap-2")
                    return(g)
                } else {
                    clusmat <- dplyr::arrange(clusmat, sample.markers)
                    g <- ggplot(clusmat, aes(x = umap1, y = umap2, 
                    group = ident)) + theme_bw() + theme(legend.title = element_blank()) + 
                    ggtitle(markers) + geom_point(aes(color = sample.markers, 
                    shape = ident), alpha = point.transparency, 
                    shape = point.shape, size = 1) + theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), panel.background = element_blank(), 
                    legend.text = element_text(size = rel(0.5))) + 
                    scale_color_gradient(limits = c(0, markers.quant[3]), 
                        low = "lightgrey", high = "red") + xlab("umap-1") + 
                    ylab("umap-2")
                    return(g)
                }
            }
        }
        if (length(samples.subset) == 0) {
            samples.subset = unique(object@meta.data$sample.name)
        }
        if ("active.ident" %in% slotNames(object)) {
            sample.name = as.factor(object@meta.data$orig.ident)
            names(sample.name) = names(object@active.ident)
            object@active.ident <- as.factor(vector())
            object@active.ident <- sample.name
            object.sub = subset(object, ident = samples.subset)
        } else {
            sample.name = as.factor(object@meta.data$orig.ident)
            names(sample.name) = names(object@active.ident)
            object@active.ident <- as.factor(vector())
            object@active.ident <- sample.name
            object.sub = subset(object, ident = samples.subset)
        }
        marker.table <- marker.table[cells.of.interest]
        present.marker.ls <- list()
        for (celltype in colnames(marker.table)) {
            print(names(marker.table[celltype]))
            present = lapply(marker.table[[celltype]], function(x) x %in% 
                rownames(object.sub$SCT@scale.data))
            absent.genes = unlist(marker.table[[celltype]])[present == 
                FALSE]
            present.genes = unlist(marker.table[[celltype]])[present == 
                TRUE]
            print(paste0("Genes not present: ", paste0(absent.genes, 
                collapse = ",")))
            print(paste0("Genes present: ", paste0(present.genes, 
                collapse = ",")))
            if (length(present.genes) == 0) {
                print(paste0(names(marker.table[celltype]), " genes were not found in object and will not be analyzed"))
            } else {
                present.marker.ls[[celltype]] <- present.genes
            }
        }
        
        padded.ls <- lapply(present.marker.ls, `length<-`, max(lengths(present.marker.ls)))
        markers.from.list <- do.call(cbind, padded.ls)
        markers.present = unlist(markers.from.list)

        if (!length(markers.present) > 0) {
            print("No markers found in dataset")
# auto removed:             return(NULL)
        }
        
        for (cell in colnames(markers.from.list)) {
            title <- cell
            markers.to.analyze <- as.character(markers.from.list[, 
                cell])
            grob <- lapply(markers.to.analyze, function(x) .plotMarkers(x))
            plot(gridExtra::arrangeGrob(grobs = grob, newpage = F, as.table = F, top = ggpubr::text_grob(title, 
                    size = 15, face = "bold")))
        }
       
# auto removed:         return(NULL)
    }

# Functions defined here will be available to call in
# the code for any table.

#install_bioconductor_package <- function(pkg) {
#}

print("template_function_AllCells_ModScoreRBCProgs_CBGL.R #########################################################################")
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
#if(Input_is_Seurat_count == 0 ){
#var_AllCells_CellTypes_SampleNames<-as.data.frame(var_AllCells_CellTypes_SampleNames)}else{var_AllCells_CellTypes_SampleNames <- var_AllCells_CellTypes_SampleNames}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ccbr1072_RBC_Progenitors_Markers<-readRDS(paste0(rds_output,"/var_ccbr1072_RBC_Progenitors_Markers.rds"))
Input_is_Seurat_count <- 0
##for(item in var_ccbr1072_RBC_Progenitors_Markers){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_ccbr1072_RBC_Progenitors_Markers<-as.data.frame(var_ccbr1072_RBC_Progenitors_Markers)}else{var_ccbr1072_RBC_Progenitors_Markers <- var_ccbr1072_RBC_Progenitors_Markers}
invisible(graphics.off())
var_AllCells_ModScoreRBCProgs_CBGL<-AllCells_ModScoreRBCProgs_CBGL(var_AllCells_CellTypes_SO,var_AllCells_CellTypes_SampleNames,var_ccbr1072_RBC_Progenitors_Markers)
invisible(graphics.off())
saveRDS(var_AllCells_ModScoreRBCProgs_CBGL, paste0(rds_output,"/var_AllCells_ModScoreRBCProgs_CBGL.rds"))
