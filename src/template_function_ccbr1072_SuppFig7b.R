# Module Score Cell Classification [scRNA-seq][CCBR] (10cf059e-0bd0-4a1c-9a3a-65b8688dab23): v191
ccbr1072_SuppFig7b <- function(AllCells_CellTypes_SO,AllCells_CellTypes_SampleNames,ccbr1072_RBC_Progenitors_Markers, ccbr1072_MultiLevelClassesForModScores) {

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    #nidapLoadPackages("SCWorkflow")
    library("SCWorkflow")

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    #Primary Inputs:
    seurat_object = AllCells_CellTypes_SO
    marker_table <- ccbr1072_RBC_Progenitors_Markers
    levels_data_frame <- ccbr1072_MultiLevelClassesForModScores

    #Basic Parameters:
    sample_names = eval(parse(text=gsub('\\[\\]','c()','c("CD47KO_1","CD47KO_2","CD47KO_3","TSP1KO_1","TSP1KO_2","TSP1KO_3","WT_1","WT_2","WT_3")')))
    sample_to_display <- c("CD47KO_1","CD47KO_2","CD47KO_3","TSP1KO_1","TSP1KO_2","TSP1KO_3","WT_1","WT_2","WT_3")
    celltypes_to_analyze <- c("RBC_Progenitor")
    manual_threshold <- as.numeric(c("0.06"))
    general_class <- c("RBC_Progenitor")
    multi_level_class <- FALSE

    #Plot Parameters:
    reduction = "tsne"
    nbins <- 24
    gradient_density_font_size <- 6
    violinplot_font_size <- 6
    step_size <- 0.1
    
    #Advanced Parameters:
    seurat_object_filename <- "seurat_object.rds"

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

    #path <- nidapGetPath(seurat_object,seurat_object_filename)
    #so <- readRDS(path)
    so <- seurat_object

    manual_threshold <- as.numeric(manual_threshold)

    MS_res <- modScore(object = so,
             samples.subset = sample_names,
             sample.to.display = sample_to_display,
             marker.table = marker_table,
#             cite.seq = proteins_presence,
             celltypes = celltypes_to_analyze,
             threshold = manual_threshold,
             general.class = general_class,
             multi.lvl = multi_level_class,
             lvl.df = levels_data_frame,
             reduction = reduction,
             nbins = nbins,
             gradient.ft.size = gradient_density_font_size,
             violin.ft.size = violinplot_font_size,
             step.size = step_size)

# auto removed:     output <- new.output()
# auto removed:     output_fs <- output$fileSystem()
return(MS_res)

# auto removed:     return(NULL)
}
#################################################
## Global imports and functions included below ##
#################################################
#install_bioconductor_package <- function(pkg) {
#  }
#install_bioconductor_package("GenomeInfoDbData_1.2.1_r351")
#suppressMessages(library(GenomeInfoDbData))
# Functions defined here will be available to call in
# the code for any table.

modScore <- function(object, samples.subset, sample.to.display, marker.table, 
    cite.seq = FALSE, celltypes, threshold = c(0), general.class, 
    multi.lvl = FALSE, lvl.df, reduction = "tsne", nbins = 10, 
    gradient.ft.size = 10, violin.ft.size = 10, step.size = 0.1) 
{
    library(Seurat)
    library(gridExtra)
    library(grid)
    library(dplyr)
    library(ggplot2)

    .modScoreCall <- function(ms.meta, threshold, reject) {
        thres.ls <- list()
        for (i in 1:ncol(ms.meta)) {
            thres.ls[[i]] <- rep(threshold[i], nrow(ms.meta))
        }
        thres.df <- data.frame(matrix(unlist(thres.ls), nrow = nrow(ms.meta)))

        # For negative markers, reverse selection to cells left of threshold
        if (grepl("_neg",colnames(ms.meta))){ 
            thres.filter <- ms.meta < thres.df
        } else {
            thres.filter <- ms.meta > thres.df
        }
        
        ms.meta.filt <- ms.meta * thres.filter
        max.col.vec <- max.col(ms.meta.filt)
        zero.filt <- as.integer(!apply(ms.meta.filt, 1, function(find_zero_rows) all(find_zero_rows == 
            0)))
        final.filt <- (max.col.vec * zero.filt) + 1
        append.name <- c(reject, names(ms.meta))
        dupl.data <- ms.meta
        dupl.data[, "Likely_CellType"] <- append.name[final.filt]
        return(dupl.data)
    }

    marker.tab <- unlist(marker.table)
    if (!"Barcode" %in% colnames(object@meta.data)) {
        object@meta.data$Barcode <- rownames(object@meta.data)
    }
    if (length(samples.subset) == 0) {
        samples.subset = unique(object@meta.data$sample.name)
    }
    colnames(object@meta.data) <- gsub("orig_ident", "orig.ident", 
        colnames(object@meta.data))
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
    rm(object)
    marker = select(marker.table, celltypes)
    marker.list = as.list(marker)
    if (sum(unlist(marker.list) %in% rownames(object.sub@assays$SCT@data)) == 
        0) {
        stop("No genes from list was found in data")
    }
    if (length(threshold) != length(celltypes)) {
        if (sum(threshold) == 0) {
            threshold <- rep(0, length(celltypes))
            print("Manual threshold set to zero - outputing preliminary data")
        } else {
            stop("Threshold length does not match # celltypes to analyze")
        }
    }
    names(threshold) <- celltypes
    figures <- list()
    exclude_cells <- c()
    h = 0
    j = 1
    for (h in seq_along(marker.list)) {
        print(names(marker.list[h]))
        present = lapply(marker.list[[h]], function(x) x %in% 
            rownames(object.sub@assays$SCT@data))
        absentgenes = unlist(marker.list[[h]])[present == FALSE]
        absentgenes = absentgenes[is.na(absentgenes) == F]
        presentgenes = unlist(marker.list[[h]])[present == TRUE]
        presentgenes = presentgenes[is.na(presentgenes) == F]
        print(paste0("Genes not present: ", paste0(absentgenes, 
            collapse = ",")))
        print(paste0("Genes present: ", paste0(presentgenes, 
            collapse = ",")))
        if (length(presentgenes) == 0) {
            print(paste0(names(marker.list[h]), " genes were not found in object and will not be analyzed"))
            exclude_cells[j] <- h
            j = j + 1
        }
    }
    if (length(exclude_cells) > 0) {
        marker.list <- marker.list[-exclude_cells]
    } else {
        marker.list <- marker.list
    }
    for (i in seq_along(marker.list)) {
        object.sub = AddModuleScore(object.sub, marker.list[i], 
            name = names(marker.list[i]), nbin = nbins, assay = "SCT")
        m = paste0(names(marker.list[i]), "1")
        object.sub@meta.data[[m]] <- scales::rescale(object.sub@meta.data[[m]], 
            to = c(0, 1))
        clusid = object.sub@meta.data[[m]]
        d <- density(clusid)
        if (reduction == "tsne") {
            p1 <- DimPlot(object.sub, reduction = "tsne", group.by = "ident")
        } else if (reduction == "umap") {
            p1 <- DimPlot(object.sub, reduction = "umap", group.by = "ident")
        } else {
            p1 <- DimPlot(object.sub, reduction = "pca", group.by = "ident")
        }
        if (reduction == "tsne") {
            clusmat = data.frame(ident = p1$data$ident, umap1 = p1$data$tSNE_1, 
                umap2 = p1$data$tSNE_2, clusid = as.numeric(object.sub@meta.data[[m]]))
        } else if (reduction == "umap") {
            clusmat = data.frame(ident = p1$data$ident, umap1 = p1$data$UMAP_1, 
                umap2 = p1$data$UMAP_2, clusid = as.numeric(object.sub@meta.data[[m]]))
        } else {
            clusmat = data.frame(ident = p1$data$ident, umap1 = p1$data$PC_1, 
                umap2 = p1$data$PC_2, clusid = as.numeric(object.sub@meta.data[[m]]))
        }
        
        clusmat <- mutate(clusmat, sample_clusid = clusmat$clusid * 
            grepl(paste(sample.to.display, collapse = "|"), clusmat$ident))
        umap.pos <- clusmat %>% group_by(clusid) %>% dplyr::summarise(umap1.mean = mean(umap1), 
            umap2.mean = mean(umap2))
        title = as.character(m)
        clusmat <- clusmat %>% dplyr::arrange(clusid)
        clusid.df <- data.frame(id = object.sub@meta.data$orig.ident, 
            ModuleScore = object.sub@meta.data[[m]])
        g <- ggplot(clusmat, aes(x = umap1, y = umap2)) + theme_bw() + 
            theme(legend.title = element_blank(), text = element_text(size = 24)) + geom_point(aes(colour = sample_clusid), 
            alpha = 0.5, shape = 20, size = 1) + scale_color_gradientn(colours = c("blue4", 
            "lightgrey", "red"), values = scales::rescale(c(0, 
            threshold[i]/2, threshold[i], (threshold[i] + 1)/2, 
            1), limits = c(0, 1))) + guides(colour = guide_legend(override.aes = list(size = 5, 
            alpha = 1))) + theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), panel.background = element_blank()) + 
            xlab("tsne-1") + ylab("tsne-2")
        g1 <- RidgePlot(object.sub, features = m, group.by = "orig.ident") + 
            theme(legend.position = "none", title = element_blank(), 
                text = element_text(size = 48)) + 
            geom_vline(xintercept = threshold[i], linetype = "dashed", 
                color = "red3") + scale_x_continuous(breaks = seq(0, 
            1, step.size))
        g2 <- ggplot(clusid.df, aes(x = id, y = ModuleScore)) + 
            geom_violin(aes(fill = id)) + theme_classic() + theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), legend.title = element_blank(), 
            panel.background = element_blank(), axis.text.x = element_blank(), 
            legend.text = element_text(size = rel(0.6)), legend.position = "top", 
            text = element_text(size = 20)) + 
            guides(colour = guide_legend(override.aes = list(size = 5, 
                alpha = 1))) + geom_hline(yintercept = threshold[i], 
            linetype = "dashed", color = "red3") + scale_y_continuous(breaks = seq(0, 
            1, step.size))
        g3 <- ggplot(data.frame(x = d$x, y = d$y), aes(x, y)) + 
            xlab("ModuleScore") + ylab("Density") + geom_line() + 
            geom_segment(aes(xend = d$x, yend = 0, colour = x)) + 
            scale_y_log10() + scale_color_gradientn(colours = c("blue4", 
            "lightgrey", "red"), values = scales::rescale(c(0, 
            threshold[i]/2, threshold[i], (threshold[i] + 1)/2, 
            1), limits = c(0, 1))) + geom_vline(xintercept = threshold[i], 
            linetype = "dashed", color = "red3") + geom_vline(xintercept = threshold[i], 
            linetype = "dashed", color = "red3") + scale_x_continuous(breaks = seq(0, 
            1, step.size)) + theme(legend.title = element_blank(), 
            text = element_text(size = 20))
        figures[[i]] = arrangeGrob(g, g1, g2, g3, ncol = 2, top = textGrob(names(marker.list[i]), 
            gp = gpar(fontsize = 14, fontface = "bold")))
    }
    colnames(object.sub@meta.data)[colnames(object.sub@meta.data) %in% 
        paste0(names(marker.list), 1)] <- names(marker.list)
    general.class <- general.class[general.class %in% colnames(object.sub@meta.data)]
    trunc.meta.gen <- object.sub@meta.data[general.class]
    gen.thrs.vec <- threshold[general.class]
    call.res <- .modScoreCall(trunc.meta.gen, gen.thrs.vec, reject = "unknown")
    call.res$Barcode <- rownames(call.res)
    if (multi.lvl) {
        for (k in 1:ncol(lvl.df)) {
            sub.class.call <- list()
            store.sub.class <- lvl.df[[k]][!is.na(lvl.df[[k]])]
            parent.class <- unique(gsub("(.*)-(.*)", "\\1", store.sub.class))
            for (parent in parent.class) {
                sub.class <- store.sub.class[grepl(parent, store.sub.class)]
                children_class <- gsub("(.*)-(.*)", "\\2", sub.class)
                parents <- call.res$Barcode[call.res$Likely_CellType == 
                  parent]
                trunc.meta.parent <- object.sub@meta.data[parents, 
                  ] %>% select(children_class)

                  # Stop hierarchical classification in case no parent cell can be called
                  if (nrow(trunc.meta.parent) == 0){
                      stop(paste0("No ",parent," can be called in ","level ",k-1," classification, try setting more lenient thresholds"))}

                for (child in children_class) {
                  plot.title <- paste("Density plot for", child, 
                    "Module Scores within", parent, "population", 
                    sep = " ")
                  figures[[length(figures) + 1]] <- ggplot(trunc.meta.parent, aes_string(x = child)) + 
                    geom_density() + ggtitle(plot.title) + geom_vline(xintercept = threshold[child], 
                    linetype = "dashed", color = "red3")
                }
                trunc.meta.no.parent <- call.res[!call.res$Likely_CellType == 
                  parent, ]
                non.parent <- rownames(trunc.meta.no.parent)
                child.thres.vec <- threshold[children_class]

                
                sub.class.call[[match(parent, parent.class)]] <- .modScoreCall(trunc.meta.parent, 
                  child.thres.vec, reject = parent) %>% select(Likely_CellType)
            }
            sub.class.call <- do.call(rbind, sub.class.call)
            sub.class.call$Barcode <- rownames(sub.class.call)
            call.res$temp.call <- sub.class.call$Likely_CellType[match(call.res$Barcode, 
                sub.class.call$Barcode)]
            call.res <- call.res %>% mutate(Likely_CellType = case_when(is.na(temp.call) ~ 
                Likely_CellType, TRUE ~ temp.call))
            call.res$temp.call <- NULL
        }
    }
    object.sub@meta.data$Likely_CellType <- call.res$Likely_CellType[match(object.sub@meta.data$Barcode, 
        call.res$Barcode)]
    
    lapply(figures, plot)

    return(object.sub)
}

print("template_function_ccbr1072_SuppFig7b.R #########################################################################")
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
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ccbr1072_RBC_Progenitors_Markers<-readRDS(paste0(rds_output,"/var_ccbr1072_RBC_Progenitors_Markers.rds"))
Input_is_Seurat_count <- 0
##for(item in var_ccbr1072_RBC_Progenitors_Markers){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_ccbr1072_RBC_Progenitors_Markers<-as.data.frame(var_ccbr1072_RBC_Progenitors_Markers)}else{var_ccbr1072_RBC_Progenitors_Markers <- var_ccbr1072_RBC_Progenitors_Markers}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ccbr1072_MultiLevelClassesForModScores<-readRDS(paste0(rds_output,"/var_ccbr1072_MultiLevelClassesForModScores.rds"))
Input_is_Seurat_count <- 0
#for(item in var_ccbr1072_MultiLevelClassesForModScores){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_ccbr1072_MultiLevelClassesForModScores<-as.data.frame(var_ccbr1072_MultiLevelClassesForModScores)}else{var_ccbr1072_MultiLevelClassesForModScores <- var_ccbr1072_MultiLevelClassesForModScores}
invisible(graphics.off())
var_ccbr1072_SuppFig7b<-ccbr1072_SuppFig7b(var_AllCells_CellTypes_SO,var_AllCells_CellTypes_SampleNames,var_ccbr1072_RBC_Progenitors_Markers,var_ccbr1072_MultiLevelClassesForModScores)
invisible(graphics.off())
saveRDS(var_ccbr1072_SuppFig7b, paste0(rds_output,"/var_ccbr1072_SuppFig7b.rds"))
