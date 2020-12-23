
library(Seurat)


seurat3_integration = function (..., sample_list = NULL, method = c("standard", "sctransform"), 
                                base_names, integration_base_name, nfeatures = 56000, dr_dim.use = 1:50, 
                                use_sig_pc_only = T, selection.method = "vst", normalization.method = "SCT", 
                                max.features = 200, k.score = 30, k.filter = 50, integration_dims = 1:30,
                                base_name = "Project") 
{
  # browser()
  if (packageVersion("Seurat") < "3") {
    stop(paste0("Require Seurat Version 3, now Seurat Version is: ", 
                packageVersion("Seurat")))
  }
  pars = c(list(...), sample_list)
  if (length(pars) != length(base_names)) {
    stop("Length of sample and base names doesn't match.")
  }
  for (i in 1:length(pars)) {
    colnames(pars[[i]]) = paste0(base_names[i], "_", colnames(pars[[i]]))
  }
  item_so = list()
  method = method[1]
  for (i in seq_along(pars)) {
    if (class(pars[[i]]) != "Seurat") {
      item_so[[i]] = CreateSeuratObject(pars[[i]], project = base_names[i])
      if (tolower(method) == "sctransform") {
        print("sctransform")
        item_so[[i]] <- SCTransform(item_so[[i]], verbose = T, 
                                    variable.features.n = nfeatures)
      }
      else if (tolower(method) == "standard") {
        print("standard")
        item_so[[i]] <- NormalizeData(item_so[[i]], verbose = T)
        item_so[[i]] <- FindVariableFeatures(item_so[[i]], 
                                             selection.method = selection.method, nfeatures = nfeatures, 
                                             verbose = T)
      }
    }
    else {
      item_so[[i]] = pars[[i]]
    }
  }
  print("FindIntegrationAnchors")
  if (tolower(method) == "sctransform") {
    features <- SelectIntegrationFeatures(object.list = item_so, 
                                          nfeatures = nfeatures)
    item_so <- PrepSCTIntegration(object.list = item_so, 
                                  anchor.features = features, verbose = T)
    anchors <- FindIntegrationAnchors(object.list = item_so, 
                                      normalization.method = normalization.method, anchor.features = features, 
                                      verbose = T, max.features = max.features, k.score = k.score, 
                                      k.filter = k.filter, dims = integration_dims)
  } else if (tolower(method) == "standard") {
    anchors <- FindIntegrationAnchors(object.list = item_so, 
                                      anchor.features = nfeatures, max.features = max.features, 
                                      k.score = k.score, k.filter = k.filter, dims = integration_dims)
  }
  print("IntegrateData")
  integrated <- IntegrateData(anchorset = anchors, normalization.method = normalization.method, 
                              verbose = FALSE)
  print("RunPCA")
  integrated <- RunPCA(object = integrated, npcs = max(dr_dim.use), 
                       verbose = T)
  dims.use = dr_dim.use
  print(ezp("use_sig_pc_only: ", use_sig_pc_only))
  if (use_sig_pc_only) {
    integrated = JackStraw(integrated, dims = max(dr_dim.use))
    integrated = ScoreJackStraw(integrated, dims = dr_dim.use)
    dims.use = integrated@reductions$pca@jackstraw@overall.p.values
    dims.use = dims.use[dims.use[, 2] < 0.05, 1]
  }
  print("RunTSNE")
  integrated <- RunTSNE(integrated, reduction = "pca", dims = dims.use, 
                        do.fast = T, k.seed = 10, check_duplicates = FALSE, perplexity = 30)
  print("RunUMAP")
  integrated <- RunUMAP(integrated, dims = dims.use)
  f_pca = integrated@reductions$pca@cell.embeddings
  f_tsne = integrated@reductions$tsne@cell.embeddings
  f_umap = integrated@reductions$umap@cell.embeddings
  # browser()
  res = c(integrated, item_so)
  saveRDS(res, ezp(base_name, "_integrated_so.RDS"))
  write.table(as.matrix(res[[1]]@assays$integrated@data), base_name, "_integrated_expression.txt",sep = '\t',row.names = 1,col.names = NA)
  write.table(as.matrix(res[[1]]@meta.data), base_name, "_integrated_metadata.txt",sep = '\t',row.names = 1,col.names = NA)
  
  write.table(res[[1]]@reductions$pca@cell.embeddings, base_name, '_PCA.txt',sep = '\t',row.names = 1,col.names = NA)
  write.table(res[[1]]@reductions$tsne@cell.embeddings, base_name, '_tSNE.txt',sep = '\t',row.names = 1,col.names = NA)
  write.table(res[[1]]@reductions$umap@cell.embeddings, base_name, '_UMAP.txt',sep = '\t',row.names = 1,col.names = NA)
  
  return(res[[1]])
}

integrate_data = function(..., base_name){
  # browser()
  print("integrating")
  so = seurat3_integration(sample_list = list(...), base_names = paste0("Dt", 1:length(list(...)))
                           , integration_base_name = paste0("Dt", 1:length(list(...)))
                           , k.filter = 30, base_name = base_name)
  # saveRDS(so, paste0(base_name, "_so.RDS"))
  # so = so[[3]]
  # sig_pc = so@reductions$pca@jackstraw$empirical.p.values
  # fast_save_table(sig_pc, base_name, '_sigPCA.txt')
  # fast_save_table(so@reductions$pca@cell.embeddings, base_name, '_PCA.txt')
  # fast_save_table(so@reductions$tsne@cell.embeddings, base_name, '_PCA.txt')
  # fast_save_table(so@reductions$umap@cell.embeddings, base_name, '_PCA.txt')
  # saveRDS(so@assays$integrated@data, paste0(base_name, "_data.RDS"))
}

res = integrate_data(Acsites, CAD1 , CAD2 , Villani, Azizi, Swarrick, CAD , Zhang10x , Zhangsmrt2 ,Smillie  , Stomach,
                     Arazi ,ZhangSmrt2Liver , Rama, Sharma , Grun, Kimlung , Maier , Lambrechts
                     , Zilionis , Reyfman , Pancreas , Cheng , Kim , Sato,Emma,Spleen,Cillo,TangHuau,Samsung, base_name = "DCverse_Selected_2020-06-18")

TPM = fast_read_table("TPM_Batch1_PSGONLY.txt")
TPZ = fast_read_table("TPM_Batch2_PSGONLY.txt")
#TPF = fast_read_table("TPM_BlisterOnly.txt")

test = fast_read_table("test2_integrated_expression.txt")

#file name change
old_name = colnames(tpm)
View(old_name)
new_name = ezp(get_string_part(old_name, "", -1), "", gsub("^dt(.*)", "\\1", get_string_part(old_name, "_", 1), ignore.case = T))
View(new_name)
colnames(tpm) = new_name

