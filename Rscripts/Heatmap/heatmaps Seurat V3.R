##### MINE 

seurat_analysis_deg <- function(TPM_file="ZhangColon_Scenic_MacroMono_TPM.txt",
                                sample_file="Scenic_ZhangColon_MacroMono_Dimred.txt",
                                group_col="Phenograph_md0.001_k150",
                                base_name="ZhangColon_MacroMono_Regulon_RNA_FC0.05", 
                                mt="^MT-", 
                                scale=F, 
                                scale_factor=1e2,
                                group1=NULL,
                                group2=NULL,
                                test.use="bimod",
                                col.low = "#FF00FF",
                                col.mid = "#000000", col.high = "#FFFF00",
                                logfc.threshold = 0.05,
                                plotDEGhm=T,topn=NULL,
                                out_folder="Output_RNA/"){
  
  library(Seurat)
  library(dplyr)
  library(Matrix)
  
  pbmc.data <- read.table(TPM_file,sep="\t",header=T,row.names=1,check.names = F)

  sample <- read.table(sample_file,sep="\t",header=T,row.names=1)

  pbmc.data <- pbmc.data[,colnames(pbmc.data) %in% rownames(sample)]

  colnames(pbmc.data) <- paste(sample[colnames(pbmc.data),group_col],colnames(pbmc.data),sep="_")

  pbmc.data <- read.table(TPM_file,sep="\t",header=T,row.names=1,check.names=F)

  if(!is.null(sample_file)){
    sample <- read.table(sample_file,sep="\t",header=T,row.names=1)
  }

  levels(sample$CellType)
  head(sample)
  #
  common_cell <- intersect(row.names(sample),colnames(pbmc.data))
  pbmc.data <- pbmc.data[,common_cell]
  if(ncol(sample)==1) sample <- cbind(sample,sample)
  sample <- sample[common_cell,]

  if((!is.null(group1)) & (!is.null(group2))){
    group12 <- row.names(sample)[as.character(sample[,group_col])==as.character(group1) | as.character(sample[,group_col])==as.character(group2)]
    pbmc.data <- pbmc.data[,group12]
  }

  if(!is.null(sample_file)){
    colnames(pbmc.data) <- paste(sample[colnames(pbmc.data),group_col],colnames(pbmc.data),sep="_")
  }

  pbmc <- CreateSeuratObject(counts = pbmc.data, project = base_name, min.cells = 3, min.features = 200)
  
  #pbmc <- readRDS("Srt_obJ_forDEG.RDS")
  #saveRDS(object = pbmc,file = "Srt_obJ_forDEG.RDS")
  
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = mt)
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = scale_factor)
  
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  
  if(is.null(group1)&is.null(group2)){
    pbmc.markers <- FindAllMarkers(object = pbmc, test.use=test.use, logfc.threshold = logfc.threshold)
  }else{
    pbmc.markers <- FindMarkers(object =pbmc, ident.1 = group1, ident.2 = group2,test.use=test.use,logfc.threshold = logfc.threshold)
  }
  
  if(sum(pbmc.markers$p_val_adj<0.05)==0){
    print("No DEGs with adjusted p values < 0.05")
    return(0)
  }
  ####################
  
  rotate = FALSE
  
  if(is.null(group1) & is.null(group2)){
    
    if(length(levels(pbmc.markers$cluster))>2){
      deg <- pbmc.markers[pbmc.markers$p_val_adj<0.05 & pbmc.markers$avg_logFC>0,]
      rotate=TRUE
    } else {
      deg <- pbmc.markers[pbmc.markers$p_val_adj<0.05,]
    }
    deg <- deg[order(deg$cluster,deg$avg_logFC,decreasing = T),]
    genelist <- unique(deg$gene)
    
    # keep topn 
    if(!is.null(topn)){
      
      degtop <- deg[order(deg$avg_logFC,decreasing = T),]
      if(length(grep(mt,rownames(degtop)))>0){
        degtop <- degtop[!grepl(mt,rownames(degtop)),]
      }
      if(length(grep("^Gm[0-9]",rownames(degtop)))>0){
        degtop <- degtop[!grepl("^Gm[0-9]",rownames(degtop)),]
      }
      if(length(grep("^RP[0-9]",rownames(degtop)))>0){
        degtop <- degtop[!grepl("^RP[0-9]",rownames(degtop)),]
      }
      degtop <- degtop[!duplicated(degtop$gene),]
      topnn <- topn
      if(topn>dim(degtop)[1]){topnn <- dim(degtop)[1]}
      degtop <- degtop[1:topnn,]
      
      degtop <- degtop[order(degtop$cluster,degtop$avg_logFC,decreasing = T),]
      top50deg <- unique(degtop$gene)
    }
    
  } else {
    
    deg <- pbmc.markers[pbmc.markers$p_val_adj<0.05,]
    deg <- deg[order(deg$avg_logFC,decreasing = T),]
    genelist <- rownames(deg)
    
    if(!is.null(topn)){
      degtop <- deg[order(abs(deg$avg_logFC),decreasing = T),]
      if(length(grep(mt,rownames(degtop)))>0){
        degtop <- degtop[!grepl(mt,rownames(degtop)),]
      }
      if(length(grep("^Gm[0-9]",rownames(degtop)))>0){
        degtop <- degtop[!grepl("^Gm[0-9]",rownames(degtop)),]
      }
      if(length(grep("^RP[0-9]",rownames(degtop)))>0){
        degtop <- degtop[!grepl("^RP[0-9]",rownames(degtop)),]
      }
      topnn <- topn
      if(topn>dim(degtop)[1]){topnn <- dim(degtop)[1]}
      degtop <- degtop[1:topnn,]
      
      degtop <- degtop[order(degtop$avg_logFC,decreasing = T),]
      top50deg <- rownames(degtop)
      
    }
  }
  
  # library(xlsx)
  # write.xlsx(pbmc.markers,paste(base_name,"_seurat_bimod.xls",sep=""))
  write.table(pbmc.markers,paste(out_folder,"/",base_name,"_seurat_bimod.txt",sep=""),sep="\t",col.names=NA, quote=F)
  # write.xlsx(deg,paste(base_name,"_seurat_bimod_DEG.xls",sep=""))
  write.table(deg,paste(out_folder,"/",base_name,"_seurat_bimod_DEG.txt",sep=""),sep="\t",col.names=NA, quote=F)
  
  # setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
  if(plotDEGhm){
    
    # adjust cex.row
    cex.row <- ifelse(length(genelist)>=80,5,
                      ifelse(length(genelist)>=70,6,
                             ifelse(length(genelist)>=60,7,
                                    ifelse(length(genelist)>=50,8,
                                           ifelse(length(genelist)>=40,9,
                                                  ifelse(length(genelist)>=30,10,
                                                         ifelse(length(genelist)>=20,11,12)))))))
    pdf(paste(out_folder,"/",base_name,"_DEG_hm.pdf",sep=""))
    print(DoHeatmap(object = pbmc, features = genelist,
                    angle = 45,size = 3))
    dev.off()
    
    png(paste(out_folder,"/",base_name,"_DEG_hm.png",sep=""), width = 2*480, height = 2*480, res = 2*72)
    print(DoHeatmap(object = pbmc, features = genelist,
                    angle = 45,size = 3))
    dev.off()
    
    if(!is.null(topn)){
      
      # adjust cex.row
      cex.row <- ifelse(length(top50deg)>=80,5,
                        ifelse(length(top50deg)>=70,6,
                               ifelse(length(top50deg)>=60,7,
                                      ifelse(length(top50deg)>=50,8,
                                             ifelse(length(top50deg)>=40,9,
                                                    ifelse(length(top50deg)>=30,10,
                                                           ifelse(length(top50deg)>=20,11,12)))))))
      
      pdf(paste(base_name,out_folder,"/","_DEG_top",topn,"_hm.pdf",sep=""))
      print(DoHeatmap(object = pbmc, features = top50deg,
                      angle = 45,size = 3))
      dev.off()
      
      png(paste(base_name,out_folder,"/","_DEG_top",topn,"_hm.png",sep=""), width = 2*480, height = 2*480, res = 2*72)
      print(DoHeatmap(object = pbmc, features = top50deg,
                      angle = 45,size = 3))
      dev.off()
    }
  }
  
  if(!is.null(group1)&!is.null(group2)){
    # return up and down regulated 
    return(paste0(dim(deg)[1],' DEGs (adjusted pvalue < 0.05, logFC=',logfc.threshold,') between ',group1,' and ',group2))
  } 
  
}

seurat_analysis_deg()
