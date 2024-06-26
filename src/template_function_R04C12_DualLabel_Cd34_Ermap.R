# Dual Labeling [scRNA-seq][CCBR] (4452781a-7015-4dad-923c-61e2a183855f): v112
R04C12_DualLabel_Cd34_Ermap <- function(R04C12_Filtered_SO,R04C12_Filtered_SampleNames) {  # produces 2 color channels and the overlay

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
    
    #Primary inputs:
    seurat_object <- R04C12_Filtered_SO

    #Basic Parameters:
    samples <- c("CD47KO_1","CD47KO_2","CD47KO_3","TSP1KO_1","TSP1KO_2","TSP1KO_3","WT_1","WT_2","WT_3")
    marker1 <- "Cd34"
    marker_1_threshold <- 0.02
    marker_1_type <- "SCT"
    marker2 <- "Ermap"
    marker_2_threshold <- 0.02
    marker_2_type <- "SCT"
    
    #Filter Parameters:
    filter_data <- TRUE
    parameter_name <- "Cd34_Ermap"
    M1_filter_direction <- "greater than"
    M2_filter_direction <- "greater than"
    apply_filter_1 <- TRUE
    apply_filter_2 <- TRUE
    filter_condition <- TRUE
    
    #Visualization Parameters:
    data_reduction <- "tsne"
    add_marker_thresholds <- TRUE
    density_heatmap <- TRUE
    point_size <- 0.5
    point_shape <- 16
    point_transparency <- 0.5

    #Advanced Parameters:
    trim_marker_1 <- FALSE
    trim_marker_2 <- FALSE
    pre_scale_trim <- 0.99    
    display_unscaled_values <- FALSE
    seurat_object_filename <- "seurat_object.rds"

    ## -------------------------------- ##
    ## Functions                        ##
    ## -------------------------------- ##
    
    # path <- nidapGetPath(seurat_object, seurat_object_filename)
    so <- seurat_object
    
    #In case of NIDAPism:
    colnames(so@meta.data)[colnames(so@meta.data) == "orig_ident"] <- "orig.ident"

    so.result <- dualLabeling(object = so,
                           samples = samples,
                           marker.1 = marker1,
                           marker.2 = marker2,
                           marker.1.type = marker_1_type,
                           marker.2.type = marker_2_type,
                           data.reduction = data_reduction,
                           point.size = point_size,
                           point.shape = point_shape,
                           point.transparency = point_transparency,
                           add.marker.thresholds = add_marker_thresholds,
                           marker.1.threshold = marker_1_threshold,
                           marker.2.threshold = marker_2_threshold,
                           filter.data = filter_data,
                           M1.filter.direction = M1_filter_direction,
                           M2.filter.direction = M2_filter_direction,
                           apply.filter.1 = apply_filter_1,
                           apply.filter.2 = apply_filter_2,
                           filter.condition = filter_condition,
                           parameter.name = parameter_name,
                           trim.marker.1 = trim_marker_1,
                           trim.marker.2 = trim_marker_2,
                           pre.scale.trim = pre_scale_trim,
                           density.heatmap = density_heatmap,
                           display.unscaled.values = display_unscaled_values) 

    #First plot showing tSNE 
    grid.draw(so.result$plot)

    #Second plot showing numbers of cells annotated for filtering
    g <- so.result$plot2
    g$width <- g$width * 2
    g$height <- g$height * 2
    grid.newpage()
    grid.draw(g)

# auto removed:     output <- new.output()
# auto removed:     output_fs <- output$fileSystem()
return(so.result$object)

    return(so.result$object)
}

#################################################
## Global imports and functions included below ##
#################################################

dualLabeling <- function (object, samples, marker.1, marker.2, marker.1.type = "SCT", 
    marker.2.type = "SCT", data.reduction = "umap", point.size = 0.5, 
    point.shape = 16, point.transparency = 0.5, add.marker.thresholds = FALSE, 
    marker.1.threshold = 0.5, marker.2.threshold = 0.5, filter.data = FALSE, 
    M1.filter.direction = "greater than", M2.filter.direction = "greater than", 
    apply.filter.1 = TRUE, apply.filter.2 = TRUE, filter.condition = TRUE, 
    parameter.name = "Marker", trim.marker.1 = TRUE, trim.marker.2 = TRUE, 
    pre.scale.trim = 0.99, density.heatmap = FALSE, display.unscaled.values = FALSE) 
{
    library(grid)
    library(gridExtra)
    library(scales)
    library(Seurat)
    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(ggExtra)

    if (!(marker.1 %in% rownames(object))) {
        stop(sprintf("%s is not found in dataset", marker.1))
    }
    if (!(marker.2 %in% rownames(object))) {
        stop(sprintf("%s is not found in dataset", marker.2))
    }
    if (!(marker.1.type %in% names(object@assays))) {
        stop(sprintf("%s slot is not found in dataset", marker.1.type))
    }
    if (!(marker.2.type %in% names(object@assays))) {
        stop(sprintf("%s slot is not found in dataset", marker.2.type))
    }
    .ggOverlay <- function(so.sub, df, marker.1, marker.2) {
        df <- df %>% arrange(mark1.scale)
        xmin <- min(df$dr1) - 0.1 * min(df$dr1)
        xmax <- max(df$dr1) + 0.1 * min(df$dr1)
        gg.z1 <- ggplot(df, aes(dr1, dr2)) + geom_point(color = rgb(red = df$mark1.scale, 
            green = 0, blue = 0), shape = point.shape, size = point.size, 
            alpha = point.transparency) + theme_classic() + xlab(paste0(data.reduction, 
            "-1")) + ylab(paste0(data.reduction, "-2")) + ggtitle(marker.1) + 
            coord_fixed()
        df <- df %>% arrange(mark2.scale)
        gg.z2 <- ggplot(df, aes(dr1, dr2)) + geom_point(color = rgb(red = 0, 
            green = df$mark2.scale, blue = 0), shape = point.shape, 
            size = point.size, alpha = point.transparency) + 
            theme_classic() + xlab(paste0(data.reduction, "-1")) + 
            ylab(paste0(data.reduction, "-2")) + ggtitle(marker.2) + 
            coord_fixed()
        df <- df %>% mutate(avg = mark2.scale + mark1.scale) %>% 
            arrange(avg)
        gg <- ggplot(df, aes(dr1, dr2)) + geom_point(color = rgb(red = df$mark1.scale, 
            green = df$mark2.scale, blue = 0), shape = point.shape, 
            size = point.size, alpha = point.transparency) + 
            theme_classic() + xlab(paste0(data.reduction, "-1")) + 
            ylab(paste0(data.reduction, "-2")) + ggtitle("Combined") + 
            coord_fixed()
        return(list(gg.z1, gg.z2, gg))
    }
    .ggOverlay2 <- function(so.sub, df, marker.1, marker.2) {
        df <- df %>% arrange(mark1.scale)
        if (display.unscaled.values == TRUE) {
            label1.min <- paste("unscaled min:", round(min(mark1), 
                digits = 2))
            label1.max <- paste("unscaled max:", round(max(mark1), 
                digits = 2))
            label1 <- paste(as.character(marker.1), label1.min, 
                label1.max, sep = "\n")
            label2.min <- paste("unscaled min:", round(min(mark2), 
                digits = 2))
            label2.max <- paste("unscaled max:", round(max(mark2), 
                digits = 2))
            label2 <- paste(as.character(marker.2), label2.min, 
                label2.max, sep = "\n")
        }
        else {
            label1 <- as.character(marker.1)
            label2 <- as.character(marker.2)
        }
        gg.z1 <- ggplot(df, aes(mark1.scale, mark2.scale)) + 
            geom_point(color = rgb(red = df$mark1.scale, green = 0, 
                blue = 0), shape = 20, size = point.size) + theme_classic() + 
            xlab(label1) + ylab(label2) + coord_fixed()
        df <- df %>% arrange(mark2.scale)
        gg.z2 <- ggplot(df, aes(mark1.scale, mark2.scale)) + 
            geom_point(color = rgb(red = 0, green = df$mark2.scale, 
                blue = 0), shape = 20, size = point.size) + theme_classic() + 
            xlab(label1) + ylab(label2) + coord_fixed()
        df <- df %>% mutate(avg = mark2.scale + mark1.scale) %>% 
            arrange(avg)
        gg <- ggplot(df, aes(mark1.scale, mark2.scale)) + geom_point(color = rgb(red = df$mark1.scale, 
            green = df$mark2.scale, blue = 0), shape = 20, size = point.size) + 
            theme_classic() + xlab(label1) + ylab(label2) + coord_fixed()
        if (add.marker.thresholds == TRUE) {
            gg.z1 <- gg.z1 + geom_vline(xintercept = t1, linetype = "dashed") + 
                geom_hline(yintercept = t2, linetype = "dashed")
            gg.z2 <- gg.z2 + geom_vline(xintercept = t1, linetype = "dashed") + 
                geom_hline(yintercept = t2, linetype = "dashed")
            gg <- gg + geom_vline(xintercept = t1, linetype = "dashed") + 
                geom_hline(yintercept = t2, linetype = "dashed")
        }
        return(list(gg.z1, gg.z2, gg))
    }
    if ("active.ident" %in% slotNames(object)) {
        sample.name <- as.factor(object@meta.data$orig.ident)
        names(sample.name) <- names(object@active.ident)
        object@active.ident <- as.factor(vector())
        object@active.ident <- sample.name
        so.sub <- subset(object, ident = samples)
    }
    else {
        sample.name <- as.factor(object@meta.data$orig.ident)
        names(sample.name) <- names(object@active.ident)
        object@active.ident <- as.factor(vector())
        object@active.ident <- sample.name
        so.sub <- subset(object, ident = samples)
    }
    t1 <- marker.1.threshold
    t2 <- marker.2.threshold
    mark1 <- so.sub@assays[[marker.1.type]]@scale.data[marker.1, 
        ]
    if (trim.marker.1 == TRUE) {
        q1 <- quantile(mark1, pre.scale.trim)
        q0 <- quantile(mark1, 1 - pre.scale.trim)
        mark1[mark1 < q0] <- q0
        mark1[mark1 > q1] <- q1
    }
    mark1.scale <- rescale(mark1, to = c(0, 1))
    mark2 <- so.sub@assays[[marker.2.type]]@scale.data[marker.2, 
        ]
    if (trim.marker.2 == TRUE) {
        q1 <- quantile(mark2, pre.scale.trim)
        q0 <- quantile(mark2, 1 - pre.scale.trim)
        mark2[mark2 < q0] <- q0
        mark2[mark2 > q1] <- q1
    }
    mark2.scale <- rescale(mark2, to = c(0, 1))
    df <- data.frame(cbind(dr1 = so.sub@reductions[[data.reduction]]@cell.embeddings[, 
        1], dr2 = so.sub@reductions[[data.reduction]]@cell.embeddings[, 
        2], mark1.scale, mark2.scale))
    gg.list <- .ggOverlay(so.sub, df, marker.1, marker.2)
    gg.list2 <- .ggOverlay2(so.sub, df, marker.1, marker.2)
    if (density.heatmap == TRUE) {
        x = df$mark1.scale
        y = df$mark2.scale
        df_heatmap <- data.frame(x = x, y = y, d <- densCols(x, 
            y, nbin = 1000, bandwidth = 1, colramp <- colorRampPalette(rev(rainbow(10, 
                end = 4/6)))))
        p <- ggplot(df_heatmap) + geom_point(aes(x, y, col = d), 
            size = 1) + scale_color_identity() + xlab(marker.1) + 
            ylab(marker.2) + theme_bw() + geom_vline(xintercept = t1, 
            linetype = "dashed") + geom_hline(yintercept = t2, 
            linetype = "dashed") 
            
        p2 <- ggMarginal(p, df_heatmap, x = marker.1, y = marker.2, type = "density")

        grob <- arrangeGrob(gg.list[[1]], gg.list[[2]], gg.list[[3]], 
            gg.list2[[1]], gg.list2[[2]], gg.list2[[3]], p2, ncol = 3)
    } else {
        grob <- arrangeGrob(gg.list[[1]], gg.list[[2]], gg.list[[3]], 
            gg.list2[[1]], gg.list2[[2]], gg.list2[[3]], ncol = 3)
    } 
    
    if (filter.data == TRUE && (apply.filter.1 == TRUE | apply.filter.2 == 
        TRUE)) {
        df <- df %>% mutate(sample = so.sub@meta.data$orig.ident) %>% 
            mutate(cellbarcode = rownames(so.sub@meta.data))
        if (M1.filter.direction == "greater than") {
            ind1 <- df$mark1.scale > t1
        } else {
            ind1 <- df$mark1.scale < t1
        }
        cat("\n")
        print("Marker 1 filter:")
        print(sum(ind1))
        if (M2.filter.direction == "greater than") {
            ind2 <- df$mark2.scale > t2
        } else {
            ind2 <- df$mark2.scale < t2
        }
        cat("\n")
        print("Marker 2 filter:")
        print(sum(ind2))
        if (apply.filter.1 == TRUE) {
            if (apply.filter.2 == TRUE) {
                if (filter.condition == TRUE) {
                  df <- df[c(ind1 & ind2), ]
                } else {
                  df <- df[c(ind1 | ind2), ]
                }
            } else {
                df <- df[ind1, ]
            }
        } else {
            if (apply.filter.2) {
                df <- df[ind2, ]
            }
        }
        colnames(df)[3:4] <- c(marker.1, marker.2)
        so.sub.df <- so.sub@meta.data %>% mutate(x = case_when(rownames(so.sub@meta.data) %in% 
            df$cellbarcode ~ TRUE, TRUE ~ FALSE))
        colnames(so.sub.df) <- sub("x", parameter.name, colnames(so.sub.df))
        data.filt <- as.data.frame.matrix(table(so.sub.df[[parameter.name]], 
            so.sub.df$orig.ident))
        data.filt$Total <- rowSums(data.filt)
        data.filt <- data.filt %>% rownames_to_column("Passed Filter")
        if (filter.condition == TRUE) {
            cond = "and"
        } else {
            cond = "or"
        }
        if (apply.filter.1 == TRUE) {
            if (apply.filter.2 == TRUE) {
                titlename <- paste("Number of cells that pass filters:\n", 
                  marker.1, M1.filter.direction, marker.1.threshold, 
                  cond, marker.2, M2.filter.direction, marker.2.threshold)
            } else {
                titlename <- paste("Number of cells that pass filter:\n", 
                  marker.1, M1.filter.direction, marker.1.threshold)
            }
        } else {
            titlename <- paste("Number of cells that pass filter:\n", 
                marker.2, M2.filter.direction, marker.2.threshold)
        }
        title <- textGrob(titlename, y = 1, vjust = 1, gp = gpar(fontsize = 15))
        grid.table <- tableGrob(data.filt, rows = NULL)
        g <- arrangeGrob(grid.table, top = title)
        g$heights[[2]] <- unit(0.5, "npc") - max(g$grobs[[1]]$heights)
        rownames(so.sub.df) <- rownames(so.sub@meta.data)
        so.sub@meta.data <- so.sub.df
    } else {
        g <- textGrob("No filtering thresholds applied")
    }
    result.list <- list(object = so.sub, plot = grob, plot2 = g)
    return(result.list)
}

# Functions defined here will be available to call in
# the code for any table.

print("template_function_R04C12_DualLabel_Cd34_Ermap.R #########################################################################")
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
invisible(graphics.off())
var_R04C12_DualLabel_Cd34_Ermap<-R04C12_DualLabel_Cd34_Ermap(var_R04C12_Filtered_SO,var_R04C12_Filtered_SampleNames)
invisible(graphics.off())
saveRDS(var_R04C12_DualLabel_Cd34_Ermap, paste0(rds_output,"/var_R04C12_DualLabel_Cd34_Ermap.rds"))
