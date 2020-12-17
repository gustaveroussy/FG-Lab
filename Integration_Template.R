library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
read.table('FILE.txt',sep = "/t",header = T,row.names = 1)


#Tonsil
Cillo = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Tang = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)


#Tissue Tonsil
Tonsil.list = list()
Tonsil.list[[1]] <- CreateSeuratObject(Cillo)
Tonsil.list[[2]] <- CreateSeuratObject(Tang)


#Colon
Smillie = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Zhang10x = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
James = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Sam = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)



Colon.list = list()
Colon.list[[1]] <- CreateSeuratObject(Smillie)
Colon.list[[2]] <- CreateSeuratObject(Zhang10x)
Colon.list[[3]] <- CreateSeuratObject(James)
Colon.list[[4]] <- CreateSeuratObject(Sam)



#Lung

Zilionis = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Reyfman = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Maier = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Lambrechts = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
CAD = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Kid = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Kiml = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)

Lung.list = list()

Lung.list[[1]] <- CreateSeuratObject(Zilionis)
Lung.list[[2]] <- CreateSeuratObject(Reyfman)
Lung.list[[3]] <- CreateSeuratObject(Maier)
Lung.list[[4]] <- CreateSeuratObject(Lambrechts)
Lung.list[[5]] <- CreateSeuratObject(CAD)
Lung.list[[6]] <- CreateSeuratObject(Kid)
Lung.list[[7]] <- CreateSeuratObject(Kiml)

#Liver
Sharma = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
ZhangLiver = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Grun = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Rama = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Macp = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)

Liver.list = list()

Liver.list[[1]] <- CreateSeuratObject(Sharma)
Liver.list[[2]] <- CreateSeuratObject(ZhangLiver)
Liver.list[[3]] <- CreateSeuratObject(Grun)
Liver.list[[4]] <- CreateSeuratObject(Rama)
Liver.list[[5]] <- CreateSeuratObject(Macp)



#Pancreas
Peng = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)

Srt_Peng = CreateSeuratObject(Peng)


#Spleen
Rudensky = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)

Srt_Rudensky = CreateSeuratObject(Rudensky)

#Skin
Cheng = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Kim = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Sato = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Asc = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)

Skin.list <- list()
Skin.list[[1]] <- CreateSeuratObject(Cheng)
Skin.list[[2]] <- CreateSeuratObject(Kim)
Skin.list[[3]] <- CreateSeuratObject(Sato)
Skin.list[[4]] <- CreateSeuratObject(Asc)




#Gasteric Cancer

Gasteric = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)

Srt_Gasteric = CreateSeuratObject(Gasteric)

#Breast
swar = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
Azizi = read.table('FILE.txt',sep = "/t",header = T,row.names = 1)
                       
list.breast <- list()
list.breast[[1]] <- CreateSeuratObject(swar)
list.breast[[2]] <- CreateSeuratObject(Azizi)

Global.list = c(Tonsil.list,Colon.list,Lung.list,Liver.list,Srt_Peng,Srt_Rudensky,Srt_Gasteric,list.breast,Skin.list)

#naming files 
names(Global.list) = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","Sato","Emma","19","20","21","22","23","24","25","26","27") 

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
save
#running tsne and Umap 
file.integrated <- RunTSNE(file.integrated, reduction = "pca", dims = dims.use,
                      do.fast = T, k.seed = 10, check_duplicates = FALSE,
                      perplexity = 30)
file.integrated <- RunUMAP(file.integrated, dims = dims.use)

saveRDS(file.integrated, file = "20201103_DC3_Integration_V2.RDS", ascii = FALSE, version = NULL,
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
