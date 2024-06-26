# Combine & Normalize [scRNA-seq][CCBR] (7f919d9a-be13-43dd-967f-cbbc8c1cf813): v227
QC_CombNorm_SO <- function(ccbr1072_SuppFig3) {

## --------- ##
## Libraries ##
## --------- ##

    library(nidapFunctions)
#    nidapLoadPackages(c("SCWorkflow","Seurat","reshape2","magrittr","dplyr","ggplot2","ggpubr","gridExtra","RColorBrewer"))
    library("SCWorkflow")
        library("Seurat")
        library("reshape2")
        library("magrittr")
        library("dplyr")
        library("ggplot2")
        library("ggpubr")
        library("gridExtra")
        library("RColorBrewer")
        
    
## -------------------------------- ##
## User-Defined Template Parameters ##
## -------------------------------- ##              
    
    # Primary Inputs
    input = ccbr1072_SuppFig3

    # Nomralization Variables
    SCT.level= "Merged"
    npcs = 30
    vars.to.regress = c()
    
    # FindVariableFeatures
    nfeatures = 2000
    low.cut = 0.1
    high.cut = 8
    low.cut.disp = 1
    high.cut.disp = 100000
    selection.method = "vst"

    # Dim Reduction
    draw.umap = TRUE
    draw.tsne = TRUE
    only.var.genes = FALSE 

    seed.for.pca = 42
    seed.for.tsne = 1
    seed.for.umap = 42

    # Clustering Varables
    clust.res.low=0.2
    clust.res.high = 1.2
    clust.res.bin = 0.2
    
    # Select PCs
    methods.pca = c("none")
    var.threshold = 0.1
    pca.reg.plot= FALSE
    jackstraw = FALSE
    jackstraw.dims = 5
        
    # Advanced
    exclude.sample = c()
    cell.count.limit = 50000
    project.name = "scRNAProject"
    cell.hashing.data = FALSE
    

## -------------------------------- ##
## Load Seurat Object.              ##
## -------------------------------- ##    

  cat("1. Reading Seurat Object from dataset: seurat_object.rds\n\n")
# auto removed:   fs <- input$fileSystem()
#path <- "./rds_output/var_ccbr1072_SuppFig3.rds"
#  so <- readRDS(path)
so <- input

## -------------------------------- ##
## Run Fucnction.                   ##
## -------------------------------- ##  
    
    out=combineNormalize(object =so,
                            
                            # Nomralization variables
                            SCT.level= SCT.level,
                            npcs = npcs,
                            vars.to.regress = vars.to.regress,
                            
                            # FindVariableFeatures
                            nfeatures = nfeatures,
                            low.cut = low.cut,
                            high.cut = high.cut,
                            low.cut.disp = low.cut.disp,
                            high.cut.disp = high.cut.disp,
                            selection.method = selection.method,

                            # Dim Reduction
                            draw.umap = draw.umap,
                            draw.tsne = draw.tsne,
                            only.var.genes = only.var.genes,

                            seed.for.pca = seed.for.pca,
                            seed.for.tsne = seed.for.tsne,
                            seed.for.umap = seed.for.umap,

                            # Clustering Varables
                            clust.res.low = clust.res.low,
                            clust.res.high = clust.res.high,
                            clust.res.bin = clust.res.bin,
                            
                            # Select PCs
                            methods.pca = methods.pca,
                            var.threshold = var.threshold,
                            pca.reg.plot = FALSE,
                            jackstraw = jackstraw,
                            jackstraw.dims = jackstraw.dims,

                            # Advanced
                            exclude.sample = exclude.sample,
                            cell.count.limit=cell.count.limit,
                            project.name = project.name,
                            cell.hashing.data = cell.hashing.data
    )

## -------------------------------- ##
## Plot Figures.                    ##
## -------------------------------- ##  

    plot(out$plots$TSNE)
    plot(out$plots$UMAP)

    if ('none'%in%methods.pca==F){
    plot(out$plots$`Elbow Plot`)
    }
    if (is.null(vars.to.regress)==F&pca.reg.plot==T){
        for(x in 1:length(vars.to.regress)){
            plot(out$plots$`Regression Plots`[[x]])
        }
    }
    if(jackstraw==T){
        plot(out$plots$JackStraw)
    }

## -------------------------------- ##
## Save Seurat object               ##
## -------------------------------- ##  

  gc()
# auto removed:   output <- new.output()
# auto removed:   output_fs <- output$fileSystem()
return(out$object)

# auto removed:     return(NULL)
}

#################################################
## Global imports and functions included below ##
#################################################

combineNormalize <- function(object,
                             
                             # Nomralization variables
                             npcs = 30,
                             SCT.level="Merged",
                             vars.to.regress = NULL,
                             
                             # FindVariableFeatures
                             nfeatures = 2000,
                             low.cut = 0.1,
                             high.cut = 8,
                             low.cut.disp = 1,
                             high.cut.disp = 100000,
                             selection.method = 'vst',
                             
                             # Dim Reduction
                             only.var.genes = FALSE,
                             draw.umap = TRUE,
                             draw.tsne = TRUE,
                             
                             seed.for.pca = 42,
                             seed.for.tsne = 1,
                             seed.for.umap = 42,
                             
                             # Clustering Varables
                             clust.res.low = 0.2,
                             clust.res.high = 1.2,
                             clust.res.bin = 0.2,
                             
                             # Select PCs
                             methods.pca = 'none',
                             var.threshold = 0.1,
                             pca.reg.plot = FALSE,
                             jackstraw = FALSE,
                             jackstraw.dims=5,
                             
                             exclude.sample = "",
                             cell.count.limit= 35000,
                             project.name = 'scRNAProject',
                             cell.hashing.data = FALSE
                             
                             
){
  
  
  
  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##
  
  
  ## --------- ##
  ## Functions ####
  ## --------- ##
  
  .plotPCA <- function(so,m,sample){
    p1 <- DimPlot(so, reduction = "pca")
    
    ## Data for PCA plot
    clusmat <- data.frame(umap1=p1$data$PC_1,
                          umap2=p1$data$PC_2, 
                          clusid=so@meta.data[[m]])
    
    ## Calculate Percent Varriation for each PC
    sumpcsd <- sum(so@reductions$pca@stdev)
    
    pcvar <- (so@reductions$pca@stdev/sumpcsd)*100
    pcvar <- formatC(pcvar,format = "g",digits=3)
    
    pcvar1 <- pcvar[1] 
    pcvar2 <- pcvar[2]
    
    
    .runCateg <- function(mat,sample){
      
      ## PCA plot for to display catigorical data
      colors=c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4","#42d4f4",
               "#f032e6","#bfef45","#fabebe","#469990","#e6beff","#9A6324",
               "#800000","#aaffc3","#808000","#000075","#a9a9a9")
      g <- ggplot(mat, aes(x=umap1, y=umap2)) +
        geom_point(aes(colour=clusid),size=1) +
        scale_color_manual(values=colors) +
        xlab(paste0("PC-1 ",pcvar[1],"%")) + 
        ylab(paste0("PC-2 ",pcvar[2],"%"))+
        theme_bw() +
        theme(legend.title=element_blank(),
              aspect.ratio = 1,
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              legend.text=element_text(size=rel(1.5)),
              axis.title=element_text(size=16),
              plot.title = element_text(size=16,face='bold')) +
        guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
        ggtitle(paste0(sample))#,"_",m))  
      return(g)
    }
    
    ## PCA plot for to display continuous data
    .runCont <- function(mat,midpt,maxpt,sample){
      g <- ggplot(mat, aes(x=umap1, y=umap2)) +
        theme_bw() +
        theme(legend.title=element_blank(),aspect.ratio = 1) +
        geom_point(aes(colour=clusid),size=1) +
        scale_colour_gradient2(low = "blue",
                               mid="lightgrey",high = "red",
                               limits=c(0, maxpt),
                               midpoint = midpt) +
        xlab(paste0("PC-1 ",pcvar[1],"%")) + 
        ylab(paste0("PC-2 ",pcvar[2],"%")) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              legend.text=element_text(size=rel(1.5)),
              axis.title=element_text(size=16),
              plot.title = element_text(size=16,face='bold')) +
        guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
        ggtitle(paste0(sample))#,"_",m))  
      return(g)
      
      
    }
    ##Run PCA plots with selected metadata
    cls=class(clusmat$clusid)
    if(cls == "factor"|cls == "character"){
      g = .runCateg(clusmat,sample)
      return(g)
    }else{
      clusmat %>% arrange(clusid) -> clusmat
      if(m=="percent.mt"){
        mid=5
        max = 10
        clusmat$clusid[clusmat$clusid > max] <- max
      }else{
        mid = quantile(clusmat$clusid)[3]
        max = quantile(clusmat$clusid,probs=0.95)
        clusmat$clusid[clusmat$clusid > max] <- max
      }
      g = .runCont(clusmat,mid,max,sample)
      return(g)
    }
    
  }
  
  
  .plotElbow <- function(so,sample){
    
    ##Select plot To determine number of PCs
    if("Elbow" %in% methods.pca){
      #Find Elbow:
      #NC Add comments for context of specific actions
      sumpcsd = sum(so@reductions$pca@stdev)
      pcvar = (so@reductions$pca@stdev/sumpcsd)*100
      cumu <- cumsum(pcvar)
      co1 <- which(cumu > 80 & pcvar < 5)[1]
      co2 <- sort(which((pcvar[1:length(pcvar) - 1] - pcvar[2:length(pcvar)]) > 
                          var.threshold), 
                  decreasing = T)[1] + 1
      pcs <- min(co1,co2)
      lab <- paste0("Elbow = ", pcs)
      xpos <- pcs + 4
    }
    
    if("Marchenko-Pastur" %in% methods.pca){
      #Using code from URD (https://rdrr.io/github/farrellja/URD/src/R/pca.R)
      #NC Add comments for context of specific actions
      pcaMarchenkoPastur <- function(M, N, pca.sdev, factor=1, do.print=T) {
        pca.eigenvalue <- (pca.sdev)^2
        marchenko.pastur.max <- (1+sqrt(M/N))^2
        pca.sig <- pca.eigenvalue > (marchenko.pastur.max * factor)
        if (do.print) {
          print(paste("Marchenko-Pastur eigenvalue null upper bound:", 
                      marchenko.pastur.max))
          if (factor != 1) {
            print(paste(length(which(pca.sig)), 
                        "PCs have eigenvalues larger than", 
                        factor, "times null upper bound."))
          } else {
            print(paste(length(which(pca.eigenvalue > marchenko.pastur.max)), 
                        "PCs have larger eigenvalues."))
          }}
        pca.sig.length = length(pca.sig[pca.sig==TRUE])
        return(pca.sig.length)
      }
      
      ## Determine Dimentions of Expression data for MarchenkoPastur
      M <- dim(so$RNA@data)[1]
      N <- dim(so$RNA@data)[2]
      pca.sdev <- so@reductions$pca@stdev
      pca.sig.num <- pcaMarchenkoPastur(M=M,N=N,pca.sdev = pca.sdev)
      lab2 = paste0("MP = ", pca.sig.num)
      xpos2 = pca.sig.num+4
    }
    ep <- ElbowPlot(so,ndims=(npcs+10)) + theme_bw() + 
      ggtitle(paste0(sample)) +
      theme(plot.title = element_text(size=16,face='bold'))
    
    if(exists("lab")){
      ep <- ep + 
        geom_vline(xintercept = pcs, color="red") +
        annotate("text",  x=xpos, y = 4, label = lab, color="red",size=4) 
    }
    if(exists("lab2")){
      ep <- ep + 
        geom_vline(xintercept = pca.sig.num, color="blue") +
        annotate("text",  x=xpos2, y = 6, label = lab2, color="blue",size=4)
    }
    return(ep)
  }
  
  ## --------------- ##
  ## Main Code Block ####
  ## --------------- ##
  
  
  #in case you want to redo this on a merged SO
  if (class(object) =="Seurat") {
    x =list()
    x[[1]] <- object
    object <- x
  }
  
  # If exclude option is TRUE, filter out undesirable sample
  if (length(c(exclude.sample)) == 0){
    object <- object
  } else {
    # object <- object[-c(exclude.sample)]
    object <- object[!names(object)%in%exclude.sample]
  }
  
  
  ### Auto detect number of cells and turn on Conserve memory ####
  
  ## Calculate total number of cells in input SO.
  cell.count <- sum(unlist((lapply(object, function(x) dim(x)[2]))))
  
  ## Setting a limit for cell numbers
  # if (dim(object)[2] < cell.count.limit) {
  if (cell.count < cell.count.limit) {
    too.many.cells <- FALSE
  } else {
    too.many.cells <- TRUE
  }
  
  ## 
  if (too.many.cells || only.var.genes) {
    conserve.memory <- TRUE
  } else {
    conserve.memory <- FALSE
  }
  
  
  ### Normalize Data ####
  
  if (SCT.level=="Merged") {
    #### Merge and SCTransform ####
    ### Merge Samples ====
    dat = vector()
    if (length(object) > 1) {
      for(i in 2:length(object)){dat=c(dat,object[[i]]) }
      
      object.merge <- merge(object[[1]], 
                            y = dat, 
                            add.cell.ids = names(object), 
                            project = project.name, 
                            merge.data = TRUE)
      allgenes <- rownames(object.merge)
      
    } else {
      object.merge <- object[[1]]
      allgenes <- rownames(object.merge)
    }
    
    #### SCTransform ####
    object.merge <-SCTransform(object.merge,
                               do.correct.umi = TRUE,
                               vars.to.regress=vars.to.regress, 
                               conserve.memory = conserve.memory, 
                               return.only.var.genes = only.var.genes)
    
    if(is.null(vars.to.regress)==F){
      object.merge.nr <-SCTransform(object.merge,
                                    do.correct.umi = TRUE,
                                    conserve.memory = conserve.memory, 
                                    return.only.var.genes = only.var.genes)
    }
    #### rescaling to lowest median SCT score accross samples ====
    object.merge=PrepSCTFindMarkers(object.merge)
    
    
    object.merge<-FindVariableFeatures(
      object = object.merge, 
      nfeatures = nfeatures, 
      mean.cutoff = c(low.cut, high.cut), 
      dispersion.cutoff=c(low.cut.disp,high.cut.disp), 
      selection.method = selection.method, 
      verbose = FALSE)
    
  } else if (SCT.level=="Sample"){
    #### SCTransform and merge ####
    #### SCTransform ####
    object.merge <- lapply(object,function(x){
      SCTransform(x,
                  do.correct.umi = TRUE,
                  vars.to.regress=vars.to.regress, 
                  conserve.memory = conserve.memory, 
                  return.only.var.genes = only.var.genes)
    })
    
    
    #### Integration features to set as Varable Features ####
    integ.features <- SelectIntegrationFeatures(
      object.list = object.merge, 
      nfeatures = nfeatures, 
      mean.cutoff = c(low.cut, high.cut), 
      dispersion.cutoff=c(low.cut.disp,high.cut.disp),
      normalization.method="SCT",
      verbose = FALSE)
    
    #### Merge Samples ====
    dat = vector()
    if (length(object.merge) > 1) {
      
      for(i in 2:length(object.merge)){dat=c(dat,object.merge[[i]]) }
      object.merge <- merge(object.merge[[1]], 
                            y = dat, 
                            add.cell.ids = names(object.merge), 
                            project = project.name, merge.data = TRUE)
      allgenes <- rownames(object.merge)
      
    } else {
      object.merge <- object.merge[[1]]
      allgenes <- rownames(object.merge)
    }
    
    #### rescaling to lowest median SCT score accross samples ####
    object.merge=PrepSCTFindMarkers(object.merge)
    
    
    #### Set Variable Features ####
    VariableFeatures(object.merge) = integ.features
    
    
    
    ### non-Regression Test
    if(is.null(vars.to.regress)==F){
      #### Merge
      object.merge.nr <- lapply(object,function(x){
        SCTransform(x,
                    do.correct.umi = TRUE,
                    conserve.memory = conserve.memory, 
                    return.only.var.genes = only.var.genes)
      })
      
      #### Integration Features
      integ.features.nr <- SelectIntegrationFeatures(
        object.list = object.merge.nr, 
        nfeatures = nfeatures, 
        mean.cutoff = c(low.cut, high.cut), 
        dispersion.cutoff=c(low.cut.disp,high.cut.disp),
        normalization.method="SCT",
        verbose = FALSE)
      
      #### Merge Samples 
      dat = vector()
      if (length(object) > 1) {
        for(i in 2:length(object.merge.nr)){dat=c(dat,object.merge.nr[[i]]) }
        object.merge.nr <- merge(object.merge.nr[[1]], 
                                 y = dat, 
                                 add.cell.ids = names(object.merge.nr), 
                                 project = project.name, merge.data = TRUE)
        allgenes.nr <- rownames(object.merge.nr)
      } else {
        object.merge.nr <- object.merge.nr[[1]]
        allgenes.nr <- rownames(object.merge.nr)
      }
      
      #### rescaling to lowest median SCT score accross samples ####
      object.merge.nr=PrepSCTFindMarkers(object.merge.nr)
      
      
      #### Set Variable Features ====
      VariableFeatures(object.merge.nr) = integ.features.nr
      
    }
    
  } else {stop("SCT method should be either Merged or Sample")}
  
  ### QC samples ####
  grobsList = list()  
  
  ## Internal Error Check
  ## Turn off pca.reg.plot if no regression variables selected.
  if(is.null(vars.to.regress)&pca.reg.plot==T){
    pca.reg.plot==F
    print('pca.reg.plot set to false because no 
          Regression Variables were selected')
  }
  
  if ('none'%in%methods.pca==F|pca.reg.plot==T|jackstraw==T) {
    
    #### PCA on individual samples ####
    object.merge.split=SplitObject(object.merge, split.by = "orig.ident")
    n=names(object.merge.split)
    object.merge.split=lapply(n,function(x){
      RunPCA(object = object.merge.split[[x]],
             npcs = (npcs+10), verbose = FALSE,
             seed.use = seed.for.pca)})
    names(object.merge.split)=n
    
  } else {"No PCA plots created"}
  
  
  #### Create PCA regression plots  ####
  if (pca.reg.plot==T) {
    
    object.merge.nr.split=SplitObject(object.merge.nr, split.by = "orig.ident")
    nr=names(object.merge.nr.split)
    object.merge.nr.split=lapply(nr,function(x){
      RunPCA(object = object.merge.nr.split[[x]],
             npcs = (npcs+10), verbose = FALSE,
             seed.use = seed.for.pca)})
    names(object.merge.nr.split)=nr
    
    
    k <- 1
    pca.grob=list()
    for (i in 1:length(vars.to.regress)){
      v=vars.to.regress[i]
      print(v)
      r <- lapply(names(object.merge.split), 
                  function(x){.plotPCA(object.merge.split[[x]],v,x)})
      nr <- lapply(names(object.merge.nr.split), 
                   function(x){.plotPCA(object.merge.nr.split[[x]],v,x)})
      
      grob = ggarrange(plotlist=r,ncol=1,
                       common.legend = F,
                       legend = 'right')%>%
        annotate_figure(top=text_grob(paste0(v,' Regression'), 
                                      face = "bold", size = 20),
                        fig.lab.size = 20,
                        fig.lab.face = 'bold',
                        fig.lab.pos='top')    
      grob.nr = ggarrange(plotlist=nr
                          ,ncol=1,
                          common.legend = F,
                          legend = 'right')%>%
        annotate_figure(top=text_grob('No Regression', 
                                      face = "bold", size = 20),
                        fig.lab.size = 20,
                        fig.lab.face ='bold',
                        fig.lab.pos='top')
      
      pca.grob=ggarrange(grob,grob.nr)
      
      grobsList[['Regression Plots']][[v]] =pca.grob
      
    }
    
  } else {
    print("PCA regression plots not created")
  }
  
  #### create Elbow plot ####
  if ('none'%in%methods.pca==F) {
    
    elbow.grob=lapply(names(object.merge.split),function(x){
      gg=
        .plotElbow(object.merge.split[[x]],sample=x)
      return(gg+ylab("")+xlab(""))
    } 
    )  
    elbow.comb=ggarrange(plotlist=elbow.grob,ncol =1) %>%
      annotate_figure(left=text_grob("Standard deviation",size=16,rot=90),
                      bottom=text_grob("PC",size=16))
    
    grobsList[["Elbow Plot"]]=elbow.comb
    grobsList[["Elbow Plot Individual"]]=elbow.grob
    
  }  
  
  #### Create Jackstraw plot  ####
  ## Jackstraw does not work with SCT data so create ScaleData
  if (jackstraw==T) {
    # jackstraw.dims = (npcs+10)
    js <-lapply(names(object.merge.split), function(x){ 
      jso=object.merge.split[[x]]
      DefaultAssay(jso)="RNA"
      jso=ScaleData(jso)%>%
        FindVariableFeatures(nfeatures = nfeatures, 
                             mean.cutoff = c(low.cut, high.cut), 
                             dispersion.cutoff=c(low.cut.disp,high.cut.disp), 
                             selection.method = selection.method, 
                             verbose = F)%>%
        RunPCA(npcs = (npcs+10), 
               verbose = FALSE,
               seed.use = seed.for.pca)%>%
        JackStraw(reduction = "pca",
                  dims = jackstraw.dims,
                  num.replicate = 100,
                  prop.freq = 0.01,
                  verbose = F,
                  maxit = 1000)%>%suppressWarnings()
    })
    names(js)=names(object.merge.split)
    js <- lapply(names(js), function(x) ScoreJackStraw(js[[x]],
                                                       dims = 1:jackstraw.dims))
    names(js )=names(object.merge.split)
    
    grob4 <- lapply(names(js), function(x) {JackStrawPlot(js[[x]],
                                                          dims = 1:jackstraw.dims)+
        ggtitle(x)+
        ylab("")+xlab("")+
        theme(plot.title = element_text(size=16,face='bold'),
              axis.title=element_text(size=16))
    })
    names(grob4 )=names(object.merge.split)
    
    js.comb=ggarrange(plotlist=grob4,ncol =1) %>%
      annotate_figure(left=text_grob("Theoretical",size=16,rot=90),
                      bottom=text_grob("Empirical",size=16))
    
    
    grobsList[['JackStraw']] <- js.comb
  }
  
  #### Detect Citeseq ####
  
  #initialize Citeseq functionality as false, 
  #later the template will check for a Protein assay and run if it finds it
  
  do.cite.seq <- FALSE
  
  if ("Protein" %in% names(object.merge@assays)){
    do.cite.seq <-TRUE
  }
  
  ### HTO ####
  if(cell.hashing.data){
    object.merge <- ScaleData(object.merge, assay = "HTO")
  }
  
  ### Dimension reduction ####
  
  object.merge <- RunPCA(object = object.merge, 
                         npcs = npcs, verbose = FALSE,
                         seed.use = seed.for.pca)
  if(is.null(vars.to.regress)==F){
    object.merge.nr <- RunPCA(object = object.merge.nr, 
                              npcs = npcs, verbose = FALSE,
                              seed.use = seed.for.pca)
  }
  
  object.merge <- RunUMAP(object = object.merge, 
                          reduction = "pca", 
                          dims = 1:npcs, 
                          seed.use=seed.for.umap)
  object.merge <- RunTSNE(object = object.merge, 
                          reduction = "pca", 
                          dim.embed = 2, 
                          dims = 1:npcs, 
                          seed.use = seed.for.tsne)
  object.merge <- FindNeighbors(object.merge, dims = 1:npcs)
  
  
  ### Citeseq ####
  
  #check for CITE-seq data and if so, run reductions
  if(do.cite.seq) {
    object.merge <- ScaleData(object.merge, assay = "Protein")
    
    print("finding protein variable features...")
    VariableFeatures(object.merge,assay="Protein") <- 
      rownames(object.merge$Protein)
    
    #Support for partial
    if(all(sapply(seq.along(object),
                  function(i) "Protein" %in% names(object[[i]]@assays)))){
      print("running protein pca...")
      object.merge = ScaleData(object.merge, 
                               assay = 'Protein',
                               verbose = FALSE) 
      object.merge <- RunPCA(object = object.merge, 
                             assay="Protein",
                             npcs = npcs,
                             reduction.name="protein.pca",
                             seed.use = seed.for.pca,
                             verbose = FALSE)
      object.merge <- RunUMAP(object = object.merge, 
                              assay="Protein", 
                              features=rownames(object.merge$Protein), 
                              reduction.name="protein.umap",
                              seed.use=seed.for.umap)
      object.merge <- RunTSNE(object = object.merge, 
                              assay="Protein", 
                              features=rownames(object.merge$Protein),
                              seed.use = seed.for.tsne,
                              reduction.name="protein.tsne",
                              check.duplicates=F)
      object.merge <- FindNeighbors(object.merge, assay="Protein",
                                    graph.name="Protein.snn",
                                    features=rownames(object.merge$Protein))
    }else{
      do.cite.seq <- FALSE #set to false so we don't cluster protein
    }
    
  } else {
    do.cite.seq <- FALSE
  }
  
  
  ## Cluster ####
  
  for (i in seq(clust.res.low,clust.res.high,clust.res.bin)){
    object.merge <- FindClusters(object.merge, resolution = i, algorithm = 1)
    if(do.cite.seq==TRUE){
      object.merge <- FindClusters(object.merge, 
                                   graph.name="Protein_snn",
                                   resolution = i, 
                                   algorithm = 1)
    }
  }
  print("Clustering successful!")
  
  
  ## create plots ####
  
  n <- 60
  qual.col.pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  qual.col.pals = qual.col.pals[c(7,6,2,1,8,3,4,5),]
  cols = unlist(mapply(brewer.pal, 
                       qual.col.pals$maxcolors, 
                       rownames(qual.col.pals)))
  
  
  
  if(draw.tsne){
    p1 <- DimPlot(object.merge, 
                  reduction = "tsne",group.by = "orig.ident", 
                  repel = TRUE,pt.size=.75) + 
      theme_classic() + 
      scale_color_manual(values=cols) + 
      theme(legend.position="right", 
            legend.text=element_text(size=rel(1.5)),
            aspect.ratio = 1,
            plot.title = element_text(size=16,face='bold',hjust = 0.5),
            axis.title=element_text(size=16)) +
      guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +         
      ggtitle("RNA TSNE")
    
    grobsList[['TSNE']] <- p1
    print("Added RNA TSNE")
    print(length(grobsList))
    
  }
  if(draw.umap){
    p2 <- DimPlot(object.merge, 
                  reduction = "umap", 
                  group.by = "orig.ident", 
                  repel = TRUE, pt.size=0.75) + 
      theme_classic() + 
      scale_color_manual(values=cols) + 
      theme(legend.position="right", 
            legend.text=element_text(size=rel(1.5)),
            aspect.ratio = 1,
            plot.title = element_text(size=16,face='bold',hjust = 0.5),
            axis.title=element_text(size=16)) + 
      ggtitle("RNA UMAP")
    
    grobsList[['UMAP']] <- p2
    print("Added RNA UMAP")
    print(length(grobsList))
    
  }
  
  ### CITEseq Figures
  if(draw.tsne & do.cite.seq){ 
    p3 <- DimPlot(object.merge, 
                  reduction = "protein_tsne", 
                  group.by = "orig.ident", 
                  repel = TRUE, pt.size=0.75) + 
      theme_classic() + 
      scale_color_manual(values=cols) + 
      theme(legend.position="right", 
            legend.text=element_text(size=rel(1.5)),
            aspect.ratio = 1,
            plot.title = element_text(size=16,face='bold',hjust = 0.5),
            axis.title=element_text(size=16)) + 
      ggtitle("Antibody TSNE")
    
    grobsList[['CITEseq TSNE']] <- p3
    print("Added Antibody TSNE")
    print(length(grobsList))
    
  }
  if(draw.umap & do.cite.seq==TRUE){ 
    p4 <- DimPlot(object.merge, 
                  reduction = "protein_umap", 
                  group.by = "orig.ident", 
                  repel = TRUE, pt.size=0.75) + 
      theme_classic() + 
      scale_color_manual(values=cols) + 
      theme(legend.position="right", 
            legend.text=element_text(size=rel(1.5)),
            aspect.ratio = 1,
            plot.title = element_text(size=16,face='bold',hjust = 0.5),
            axis.title=element_text(size=16)) + 
      ggtitle("Antibody UMAP")
    grobsList[['CITEseq UMAP']] <- p4
    print("Added Antibody UMAP")
    print(length(grobsList))
    
  }
  
  return(list(object=object.merge,
              plots=grobsList))
}

print("template_function_QC_CombNorm_SO.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
currentdir <- getwd()
rds_output <- paste0(currentdir,'/rds_output')
var_ccbr1072_SuppFig3<-readRDS(paste0(rds_output,"/var_ccbr1072_SuppFig3.rds"))
Input_is_Seurat_count <- 1
#for(item in var_ccbr1072_SuppFig3){ if (class(item)=="Seurat"){Input_is_Seurat_count = Input_is_Seurat_count + 1}}
if(Input_is_Seurat_count == 0 ){
var_ccbr1072_SuppFig3<-as.data.frame(var_ccbr1072_SuppFig3)}else{var_ccbr1072_SuppFig3 <- var_ccbr1072_SuppFig3}
invisible(graphics.off())
var_QC_CombNorm_SO<-QC_CombNorm_SO(var_ccbr1072_SuppFig3)
invisible(graphics.off())
saveRDS(var_QC_CombNorm_SO, paste0(rds_output,"/var_QC_CombNorm_SO.rds"))
