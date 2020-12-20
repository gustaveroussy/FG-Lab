
remotes::install_github("satijalab/seurat", ref = "release/4.0.0")
remotes::install_github("jlmelville/uwot")
remotes::install_github("mojaveazure/seurat-disk")

#Download Azimuth packages

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(Matrix)

#Load Reference

reference = readRDS('Reference_Azimuth.RDS')

#PLotting

DimPlot(object = reference, reduction = "umap", group.by = 'Phenograph', label = T, label.size = 3, repel = TRUE) + NoLegend()

#Load Query Data
Query = read.table('metaD_TPM_monocytesonly/TPM_monocytesonly.txt',sep = '\t',header = T, row.names = 1)
Querymeta = read.table('metaD_TPM_monocytesonly/metaD_monocytesonly.txt',sep = '\t',header = T,row.names = 1)
Query = CreateSeuratObject(Query,meta.data = Querymeta)
Query <- NormalizeData(pbmc3k2, verbose = FALSE)

# Generate anchor point for query dataset, use for dims and npcs the amount of PC that you have used in your reference
anchors <- FindTransferAnchors(
  reference = reference,
  query = Query,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:100,
  npcs = 100
)

Query <- MapQuery(
  anchorset = anchors,
  query = Query,
  reference = reference,
  refdata = list(
    meta.data= "Phenograph"
  ),
  reference.reduction = "pca",
  reduction.model = "umap"
)

#Plot predicted UMAP of the query dataset

p1 = DimPlot(pbmc3k2, reduction = "ref.umap", group.by = "predicted.meta.data", label = TRUE, label.size = 3, repel = TRUE)
p1

#Merge datasets
reference$id <- 'reference'
pbmc3k2$id <- 'query'
refquery <- merge(reference, pbmc3k2)
refquery[["pca"]] <- merge(reference[["pca"]], pbmc3k2[["ref.pca"]])
refquery <- RunUMAP(refquery, reduction = 'pca', dims = 1:50)
DimPlot(refquery, group.by = 'id')

#Save RDS
saveRDS(reference, file = "Reference_Azimuth.RDS")
saveRDS(pbmc3k2, file = "Query_Azimuth.RDS")

#UMAP from Azimuth
umap_res = Query@reductions$ref.umap@cell.embeddings
umap_res = cbind(rownames(umap_res), umap_res)
write.table(umap_res,file = "Azimuth_Query_UMAP.txt", sep = "\t", row.names = F, col.names = TRUE)

#UMAP from Reference
umap_ref = reference@reductions$umap@cell.embeddings
umap_ref = cbind(rownames(umap_ref), umap_ref)
write.table(umap_ref,file = "Azimuth_Ref_UMAP.txt", sep = "\t", row.names = F, col.names = TRUE)

#Phenograph from Azimuth
Pheno = pbmc3k2@meta.data
write.table(Pheno, file = "Query_Metadata.txt",sep = '\t', row.names = F,col.names = T)

