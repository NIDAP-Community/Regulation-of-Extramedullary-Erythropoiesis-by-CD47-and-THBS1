# Violin Plot from Seurat Object [scRNASeq][CCBR] (10e8db05-d402-4310-bb56-7e00ee139141): v91
R04C14_Filtered_ViolinPlot_Fig4_Xpo1 <- function(R04C14_Filtered_SO,R04C14_Filtered_MetadataTable) {

    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    
    library("SCWorkflow")

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##
    
    #Primary Inputs:
    seurat_object <- R04C14_Filtered_SO

    #Basic Parameters:
    assay <- 'SCT'
    slot <- 'scale.data'
    group <- 'Treatment_Group'
    facet_on <- TRUE
    facet_var <- 'Treatment_Group'
    group_subset <- c()
    genes <- c("Xpo1")

    ##Filter and Normalization:
    #filter_outliers <- FALSE
    #scale_data <- TRUE
    #log_scale_data <- FALSE
    #outlier_low <- 0.1
    #outlier_high <- 0.9

    # Visualziation:
    jitter_points <- TRUE
    jitter_dot_size <- 0.5

    #Advanced Parameters:
    seurat_object_filename <- "seurat_object.rds"

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

    # path <- nidapGetPath(seurat_object,seurat_object_filename)
    so <- seurat_object

    colnames(so@meta.data) <- gsub("\\.","_",colnames(so@meta.data))

    # Subset data (optional)
    if(!is.null(group_subset)){
        library(Seurat)
        
        Idents(so) <- so[[group]]
        so <- subset(so, idents = group_subset)
        }
    
    genes <- genes[genes %in% rownames(Seurat::GetAssayData(so, assay = assay, slot = slot))]

    for (indv_gene in genes){
        violins_res <- violinPlot(object = so, assay, slot, genes = indv_gene, group, facet_data = facet_on, facet_by = facet_var, jitter_points, jitter_dot_size)
        
        plot(violins_res$fig)
    }

return(violins_res$stat)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

violinPlot <- function (object, assay, slot, genes, group, facet_data = FALSE, facet_by = "", jitter_points, jitter_dot_size) 
{
    library(Seurat)
    library(ggplot2)
    library(gridExtra)
    library(tidyr)
    library(dplyr)
    library(broom)

    if (!assay %in% Assays(object)) {
        stop("expression data type was not found in Seurat object")
    } else if (!slot %in% slotNames(object[[assay]])) {
        stop("slot not found in Seurat[[assay]]")
    } else if (all(!genes %in% rownames(object[[assay]]))) {
        stop("no genes were found in Seurat object")
    } else if (!group %in% colnames(object@meta.data)) {
        stop("grouping parameter was not found in Seurat object")
    } else if (!is.null(facet_by)) {
        if (!facet_by %in% colnames(object@meta.data)) {
            stop("facet parameter was not found in Seurat object")
        }
    }

    # Scale to non-negative for visualization
    gene_mtx <- as.matrix(GetAssayData(object, assay = assay, slot = slot))

    ## For CCBR 1072 - comment out for downstream visualization ##
    #gene_mtx <- scales::rescale(gene_mtx, to = c(0,1))
    ## End ##

    print(paste0(genes[!genes %in% rownames(gene_mtx)], 
        " not found and will not be displayed"))

    genes.present <- genes[genes %in% rownames(gene_mtx)]

    meta_sub <- object@meta.data[,c(group,facet_by)]

    for (col in genes.present) {
            meta_sub[[col]] <- gene_mtx[col,]
            }

    data_df <- meta_sub %>% pivot_longer(genes.present, names_to = "Gene", values_to = "Expression")

    
    # set Gene as factor in data_df, so faceted plots will not be alphabetical 
    data_df$Gene <- factor(data_df$Gene, levels = genes.present)

    unique_facets <- unique(object@meta.data[,facet_by])
    available_linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")

        # If you have more unique values than available linetypes, this will recycle them
    linetype_mapping <- rep(available_linetypes, length.out = length(unique_facets))

    available_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")

# Map the colors to the unique values
# If there are more unique sets than available colors, this will recycle the colors
    color_mapping <- setNames(rep(available_colors, length.out = length(unique_facets)), unique_facets)

## For CCBR 1072 - add Matt's legacy code ##
    # Set up the common elements of the plot
  g <- ggplot(data_df, aes(x = .data[[group]], y = Expression, fill = .data[[facet_by]])) + 
    geom_violin(scale = "width", position = position_dodge(width = 0.9), trim = F) +
    #geom_jitter(height = 0, width = 0.05, size=0.1) +
    scale_fill_brewer(palette = "Set1") + 
    # scale_linetype_manual(values = linetype_mapping) + 
    facet_wrap(~ Gene, scales = "free_y", ncol = 3, strip.position = "left") + 
    theme_classic() + 
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 14, color = "black", face = "bold"),
          strip.text.y = element_text(size = 14, color = "black", face = "bold"),
          strip.placement = "outside")
  
  g <- g + ylim(0,max(data_df$Expression) + 20)
  
  g <- g + geom_jitter(height = 0)
  
  g <- g + geom_boxplot(width=0.1, fill="white")

  ## End of modification

    # Function to calculate p-values for a single gene within a cell type
calculate_p_values <- function(data, data_group, data_gene) {
  # Subset data for the specific cell type and gene
  data_sub <- data[data[,group] == data_group & data[,"Gene"] == data_gene,]
  
  # Perform ANOVA and Tukey HSD
  fit <- aov(as.formula(paste("Expression ~", facet_by)), data = data_sub)
  tukey_result <- TukeyHSD(fit)
  
  # Tidy up the results and add metadata
  tidy_tukey_result <- tidy(tukey_result)
  tidy_tukey_result$gene <- data_gene
  tidy_tukey_result$group <- data_group
  
  return(tidy_tukey_result)
}

# List unique cell types
unique_groups <- unique(data_df[[group]])

# Filter out 
facet_df <- table(data_df[[group]], data_df[[facet_by]])

# Find rows with more than one non-zero column
count_non_zero <- function(row) {
  sum(row != 0)
}
non_zero_counts <- apply(facet_df, 1, count_non_zero)

# Use rownames whose values are in more than 1 column 
unique_groups <- names(non_zero_counts)[non_zero_counts > 1]

# Calculate p-values for each cell type and gene
p_values_list <- list()
for (indv_group in unique_groups) {
  p_values_list[[indv_group]] <- do.call(rbind, lapply(genes.present, function(gene) calculate_p_values(data_df, data_group = indv_group, data_gene = gene)))
}

# Combine the results into a single data frame
p_values_df <- do.call(rbind, p_values_list)

    final_res <- list(fig = g, stat = p_values_df)

    return(final_res)
}

#install_bioconductor_package <- function(pkg) {
#}

print("template_function_R04C14_Filtered_ViolinPlot_Fig4_Xpo1.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_R04C14_Filtered_SO<-readRDS(paste0(rds_output,"/var_R04C14_Filtered_SO.rds"))
Input_is_Seurat_count <- 1
#for(item in var_R04C14_Filtered_SO){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_R04C14_Filtered_SO<-as.data.frame(var_R04C14_Filtered_SO)}else{var_R04C14_Filtered_SO <- var_R04C14_Filtered_SO}
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_R04C14_Filtered_MetadataTable<-NULL #(paste0(rds_output,"/var_R04C14_Filtered_MetadataTable.rds"))
Input_is_Seurat_count <- 0
#for(item in var_R04C14_Filtered_MetadataTable){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_R04C14_Filtered_MetadataTable<-as.data.frame(var_R04C14_Filtered_MetadataTable)}else{var_R04C14_Filtered_MetadataTable <- var_R04C14_Filtered_MetadataTable}
invisible(graphics.off())
var_R04C14_Filtered_ViolinPlot_Fig4_Xpo1<-R04C14_Filtered_ViolinPlot_Fig4_Xpo1(var_R04C14_Filtered_SO,var_R04C14_Filtered_MetadataTable)
invisible(graphics.off())
saveRDS(var_R04C14_Filtered_ViolinPlot_Fig4_Xpo1, paste0(rds_output,"/var_R04C14_Filtered_ViolinPlot_Fig4_Xpo1.rds"))
