
# RUN Seurat Heatmap with a list of specific genes

run_Seurat_heatmap <- function(TPM_file="TPM.txt",
                               group_file="group.txt",
                               group_col="GroupID",
                               list.gene.file='list_gene.txt',
                               group.order=NULL,
                               base_name='',
                               ifpng=T,
                               mt="^MT-", 
                               scale=F,
                               scale_factor=1e2,
                               add_missing_gene=FALSE,
                               col_blue_white_red =TRUE){
  
  # load library
  library(Seurat)
  
  pbmc.data <- read.table( TPM_file , sep= "\t" , row.names = 1 , check.names = F, header = T )
  sample <- read.table(group_file, sep = "\t" , header = T, row.names = 1 )
  
  pbmc.data.sub <- pbmc.data[,colnames(pbmc.data) %in% rownames(sample)]
  colnames(pbmc.data.sub) <- paste(sample[colnames(pbmc.data.sub),group_col],colnames(pbmc.data.sub),sep="_")
  
  pbmc.data.sub[1:5,1:5]
  pbmc.data=pbmc.data.sub
  
  if(!is.null(list.gene.file)){
    # load list.gene
    list.gene <- read.table(list.gene.file, sep="\t")
    list.gene <- list.gene[,1]
    head(list.gene)
    if(!is.factor(list.gene)){list.gene <- as.factor(list.gene)} # list.gene has to be a factor
  } else {
    list.gene <- rownames(pbmc.data)
  }
  
  # be sure that list.gene into rownames(tpm)
  print(paste(length(list.gene),'genes'))
  s <- sum(!list.gene %in% rownames(pbmc.data))
  if(s>0){
    if(s==1){
      print(paste(s,'gene not found :',paste(list.gene[!list.gene %in% rownames(pbmc.data)], collapse=', ')))
    } else {
      print(paste(s,'genes not found :',paste(list.gene[!list.gene %in% rownames(pbmc.data)], collapse=', ')))
    }
    # if missing genes, add to heatmap or not
    if(add_missing_gene==TRUE){
      n_missing <- length(list.gene[!list.gene %in% rownames(pbmc.data)])
      pbmc.data[c((dim(pbmc.data)[1]+1):(dim(pbmc.data)[1]+n_missing)),] <- 0
      rownames(pbmc.data) <- c(rownames(pbmc.data)[1:(dim(pbmc.data)[1]-n_missing)],as.character(list.gene[!list.gene %in% rownames(pbmc.data)]))
    }
  }
  
  pbmc <- CreateSeuratObject(pbmc.data, min.cells = 0, min.features = 0, project = base_name)
  
  #mito.genes <- grep(pattern = mt, x = rownames(x = pbmc@data), value = TRUE)
  #percent.mito <- colSums(pbmc@raw.data[mito.genes, ]) / colSums(pbmc@raw.data)
  #pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
  #pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = mt)
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = scale_factor)
  #pbmc@meta.data

 pbmc <- ScaleData(pbmc, features = list.gene)
  
  # be sure that length of group.order is the same that number of clustering 
  if(!is.null(group.order) & 
     ((length(group.order) != length(levels(pbmc))) | 
      (sum(group.order %in% levels(pbmc))!=length(levels(pbmc))))){
    print('Number or label of clusters in group.order is incorrect')
  } else {
    
    levels(pbmc) <- group.order
    
    # Heatmap
    if(col_blue_white_red==TRUE){
      pdf(paste(base_name,"_cluster_hm_list_genes_ordered_BWR.pdf",sep=""))
      print(DoHeatmap(object = pbmc, features = list.gene,
                      angle = 90,size = 3,draw.lines = T,lines.width = 2,raster = F)
            + NoLegend()
            + scale_fill_gradientn(colors = c("blue", "white", "red"))
            +theme(axis.text.y = element_text(size = 7)))
      dev.off()
    }else{
      pdf(paste(base_name,"_cluster_hm_list_genes_ordered.pdf",sep=""))
      print(DoHeatmap(object = pbmc, features = list.gene,
                      angle = 90,size = 3,draw.lines = T,lines.width = 2,raster = F)
            + NoLegend()
            + scale_fill_gradientn(colors = c("#FF00FF", "#000000", "#FFFF00"),na.value = "white")
            +theme(axis.text.y = element_text(size = 7)))
      dev.off()
    }
    
    if(ifpng==TRUE){
      if(col_blue_white_red==TRUE){
        png(paste(base_name,"_cluster_hm_BWR.png",sep=""),width = 2*480, height = 2*480, res = 2*72,type = "cairo")
        print(DoHeatmap(object = pbmc, features = list.gene,
                        angle = 90,size = 3,draw.lines = T,lines.width = 2,raster = F)
              + NoLegend()
              + scale_fill_gradientn(colors = c("blue", "white", "red"))
              +theme(axis.text.y = element_text(size = 7)))
        dev.off()
      }else{
        png(paste(base_name,"_cluster_hm.png",sep=""),width = 2*480, height = 2*480, res = 2*72,type = "cairo")
        print(DoHeatmap(object = pbmc, features = list.gene,
                        angle = 90,size = 3,draw.lines = T,lines.width = 2,raster = F)
              + NoLegend()
              + scale_fill_gradientn(colors = c("#FF00FF", "#000000", "#FFFF00"),na.value = "white")
              +theme(axis.text.y = element_text(size = 7)))
        dev.off()
      }
    }
  }
  
}
