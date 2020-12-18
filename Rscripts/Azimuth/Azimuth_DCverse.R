

remotes::install_github("satijalab/seurat", ref = "release/4.0.0")
remotes::install_github("jlmelville/uwot")
remotes::install_github("mojaveazure/seurat-disk")

#Test making Azimuth

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(Matrix)


RenameIdents(reference)

reference = readRDS('../20201120_Reference_Azimuth.RDS')
anchors = readRDS('anchors.RDS')
pbmc3k2 = readRDS('SuccessAzimuth.RDS')
reference = UpdateSeuratObject(reference)
rm(pmbc.updated)
reference[["umap"]] <- CreateDimReducObject(embeddings = Umap, key = "UMAP_", assay = DefaultAssay(reference))

reference <- AddMetaData(reference,Meta,col.name = NULL)




#PC selection
dims.use = file.integrated@reductions$pca@jackstraw@overall.p.values
dims.use = dims.use[dims.use[, 2] < 10e-100, 1]

#Try Run Umap
reference <- RunUMAP(file.integrated, dims.use,min.dist = 0.05,return.model = T,reduction.key = "UMAP_")

#Phenograph
reference <- FindNeighbors(reference,dims = dims.use)


reference_TPM <- fast_read_table('IDO_TPM_ALL_Macroverse_TOp164Genes.txt')
Meta <- read.table('20201013_DCverse_Meta_top25PCs_Full_Azimuth.txt',sep = '\t',header = T, row.names = 1)
PCA <- fast_read_table('20200624_Global_Macro_Mono_pca_PC100.txt')
Umap <- Meta[,c(0,14:15)]
Umap <- as.data.frame(Umap)
reference = CreateSeuratObject(reference_TPM,meta.data = Meta,)
reference = file.integrated
reference@meta.data <- Meta

reference@reductions$umap@key <- "UMAP_"

#PLotting

DimPlot(object = reference, reduction = "umap", group.by = 'DC_Phenograph_k100', label = T, label.size = 3, repel = TRUE) + NoLegend()
DimPlot(object = reference, reduction = "umap", group.by = "Phenograph", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
DimPlot(object = reference, reduction = "umap", split.by = "Phenograph", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
DimPlot(object = reference, reduction = "umap", group.by = "Phenograph", label = TRUE, label.size = 3, repel = TRUE,
        cols = c('white', 'white', 'white', 'darkgreen', 'white', 'Purple', 'white', 'white', 'white', 'white', 'white',
                 'white', 'white', 'white', 'white', 'white', 'white')) + NoLegend()
DimPlot(object = reference, reduction = "umap", group.by = "Phenograph", label = TRUE, label.size = 3, repel = TRUE,
        cols = c('LightGrey', 'LightGrey', 'LightGrey','darkgreen', 'LightGrey', 'Purple', 'LightGrey',  'LightGrey', 'LightGrey', 'LightGrey', 'LightGrey',
                 'LightGrey', 'LightGrey', 'LightGrey', 'LightGrey', 'LightGrey', 'LightGrey')) + NoLegend()

#Change UMAP to one we use

Umap = as.matrix(Umap)
reference@reductions$umap@cell.embeddings = Umap
colnames(x = reference[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:2)
reference@reductions$umap@key <- "UMAP_"



#Load Data
pbmc3k2 = readRDS('COVID_BLood_Myeloid.RDS')
pbmc3k2 = read.table('metaD_TPM_monocytesonly/TPM_monocytesonly.txt',sep = '\t',header = T, row.names = 1)
pbmc3ksmeta = read.table('metaD_TPM_monocytesonly/metaD_monocytesonly.txt',sep = '\t',header = T,row.names = 1)
pbmc3k2 = CreateSeuratObject(pbmc3k2,meta.data = pbmc3ksmeta)
pbmc3k2 <- NormalizeData(pbmc3k2, verbose = FALSE)
#pbmc3k2 = readRDS('20201120_Ankur_Azimuth.RDS')

anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc3k2,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:100,
  npcs = 100
)
?FindTransferAnchors

pbmc3k2 <- MapQuery(
  anchorset = anchors,
  query = pbmc3k2,
  reference = reference,
  refdata = list(
    meta.data= "Phenograph"
  ),
  reference.reduction = "pca",
  reduction.model = "umap"
)

saveRDS(anchors, "anchors.RDS")
?MapQuery

p1 = DimPlot(pbmc3k2, reduction = "ref.umap", group.by = "predicted.meta.data",split.by = 'Status', label = TRUE, label.size = 3, repel = TRUE)
p1

reference$id <- 'reference'
pbmc3k2$id <- 'query'
refquery <- merge(reference, pbmc3k2)
refquery[["pca"]] <- merge(reference[["pca"]], pbmc3k2[["ref.pca"]])
refquery <- RunUMAP(refquery, reduction = 'pca', dims = 1:50)
DimPlot(refquery, group.by = 'id')

test = pbmc3k

UpdateSeuratObject(reference)


#Mapping

pbmc3k2 <- MapQuery(
  anchorset = anchors,
  query = pbmc3k2,
  reference = reference,
  refdata = list(
    celltype.l1 = "Phenograph",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "pca",
  reduction.model = "umap"
)


TPMCOVID = as.data.frame(pbmc3k2@assays$RNA@data)

write.table(TPMCOVID,sep = '\t',col.names = NA,row.names = T,"TPM_COVID_BLOOD_Myeloid.txt")

saveRDS(reference, file = "20201217_Reference_DCverse_Azimuth.RDS")
saveRDS(pbmc3k2, file = "20201202_COVID_BLOOD_Azimuth.RDS")
#UMAP from Azimuth

umap_res = pbmc3k2@reductions$ref.umap@cell.embeddings
umap_res = cbind(rownames(umap_res), umap_res)
write.table(umap_res,file = "20201204_Azimuth_COVID_Blood_UMAP_new.txt", sep = "\t", row.names = F, col.names = TRUE)

#UMAP from Reference
umap_ref = reference@reductions$umap@cell.embeddings
umap_ref = cbind(rownames(umap_ref), umap_ref)
write.table(umap_ref,file = "20201204_Azimuth_Ref_UMAP_MD0.01_0.05 PC.txt", sep = "\t", row.names = F, col.names = TRUE)

FeaturePlot(pbmc3k2, features = c('C1QA'),split.by = 'disease', cols = c('blue','green','red'))

Pheno = pbmc3k2@meta.data
write.table(Pheno, file = "20201204_COVID_BLood_Metadata.txt",sep = '\t', row.names = F,col.names = T)
