library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
read.table('FILE.txt',sep = "/t",header = T,row.names = 1)


#Tonsil
Tonsil1 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Tonsil2 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)


#Tissue Tonsil
Tonsil.list = list()
Tonsil.list[[1]] <- CreateSeuratObject(Cillo)
Tonsil.list[[2]] <- CreateSeuratObject(Tang)


#Colon
Colon1 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Colon2 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Colon3 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Colon4 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)



Colon.list = list()
Colon.list[[1]] <- CreateSeuratObject(Colon1)
Colon.list[[2]] <- CreateSeuratObject(Colon2)
Colon.list[[3]] <- CreateSeuratObject(Colon3)
Colon.list[[4]] <- CreateSeuratObject(Colon4)



#Lung

Lung1 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Lung2 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Lung3 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Lung4 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Lung5 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Lung6 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Lung7 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)

Lung.list = list()

Lung.list[[1]] <- CreateSeuratObject(Lung1)
Lung.list[[2]] <- CreateSeuratObject(Lung2)
Lung.list[[3]] <- CreateSeuratObject(Lung3)
Lung.list[[4]] <- CreateSeuratObject(Lung4)
Lung.list[[5]] <- CreateSeuratObject(Lung5)
Lung.list[[6]] <- CreateSeuratObject(Lung6)
Lung.list[[7]] <- CreateSeuratObject(Lung7)

#Liver
Liver1 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Liver2 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Liver3 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Liver4 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Liver5 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)

Liver.list = list()

Liver.list[[1]] <- CreateSeuratObject(Liver1)
Liver.list[[2]] <- CreateSeuratObject(Liver2)
Liver.list[[3]] <- CreateSeuratObject(Liver3)
Liver.list[[4]] <- CreateSeuratObject(Liver4)
Liver.list[[5]] <- CreateSeuratObject(Liver5)



#Pancreas
Pancreas = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)

Srt_Peng = CreateSeuratObject(Pancreas)


#Spleen
Spleen = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)

Srt_Rudensky = CreateSeuratObject(Spleen)

#Skin
Skin1 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Skin2 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Skin3 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Skin4 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)

Skin.list <- list()
Skin.list[[1]] <- CreateSeuratObject(Skin1)
Skin.list[[2]] <- CreateSeuratObject(Skin2)
Skin.list[[3]] <- CreateSeuratObject(Skin3)
Skin.list[[4]] <- CreateSeuratObject(Skin4)




#Gasteric Cancer

Gasteric = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)

Srt_Gasteric = CreateSeuratObject(Gasteric)

#Breast
Breast1 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Breast2 = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
                       
list.breast <- list()
list.breast[[1]] <- CreateSeuratObject(Breast1)
list.breast[[2]] <- CreateSeuratObject(Breast2)

Global.list = c(Tonsil.list,Colon.list,Lung.list,Liver.list,Srt_Peng,Srt_Rudensky,Srt_Gasteric,list.breast,Skin.list)

#naming files 
names(Global.list) = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27") 

# normalizing data sets
for (i in 1:length(Global.list)) {
  Global.list[[i]] <- NormalizeData(Global.list[[i]], verbose = FALSE)
  Global.list[[i]] <- FindVariableFeatures(Global.list[[i]], selection.method = "vst", 
                                             nfeatures = 3000, verbose = T)
}

# running anchors between data sets
reference.list <- Global.list[c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","Sato","Emma","19","20","21","22","23","24","25","26","27")] #
file.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30 ,k.filter = 30)
file.integrated <- IntegrateData(anchorset = file.anchors, dims = 1:30)
saveRDS(file.integrated, file = "Step3.RDS", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

#running PCA
DefaultAssay(file.integrated) <- "integrated"
file.integrated <- ScaleData(file.integrated, verbose = T)
file.integrated <- RunPCA(file.integrated, npcs = 50, verbose = T)

#jackstraw
file.integrated = JackStraw(file.integrated, dims = max(1:50))
file.integrated = ScoreJackStraw(file.integrated, dims = 1:50)
dims.use = file.integrated@reductions$pca@jackstraw@overall.p.values
dims.use = dims.use[dims.use[, 2] < 0.05, 1] #taking significant PCA and assigning dims for UMAP

#running tsne and Umap 
file.integrated <- RunTSNE(file.integrated, reduction = "pca", dims = dims.use,
                      do.fast = T, k.seed = 10, check_duplicates = FALSE,
                      perplexity = 30)
file.integrated <- RunUMAP(file.integrated, dims = dims.use)

saveRDS(file.integrated, file = "Integration.RDS", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
#print 
pca_res = file.integrated@reductions$pca@cell.embeddings

tsne_res = file.integrated@reductions$tsne@cell.embeddings

umap_res = file.integrated@reductions$umap@cell.embeddings


pca_res = cbind(rownames(pca_res), pca_res)
write.table(pca_res,file = "_pca_PC50.txt", sep = "\t", row.names = F, col.names = TRUE)
tsne_res = cbind(rownames(tsne_res), tsne_res)
write.table(tsne_res,file = "tsne_PC50.txt", sep = "\t", row.names = F, col.names = TRUE)
umap_res = cbind(rownames(umap_res), umap_res)
write.table(umap_res,file = "umap_PC50.txt", sep = "\t", row.names = F, col.names = TRUE)

sig_pc_score = file.integrated@reductions$pca@jackstraw$overall.p.values
select_pc = sig_pc_score[sig_pc_score[, 2] < 0.05, 1]
sig_pc_score = sig_pc_score[select_pc, ]
fast_save_table(sig_pc_score, "", "-50PC_SigPC.txt")
