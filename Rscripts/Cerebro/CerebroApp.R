

library(Seurat)

Covid = CreateSeuratObject(tissue_TPM,meta.data = selected_tpm,min.cells = 3)


BiocManager::install('romanhaa/cerebroApp')

if ( 'BiocManager' %in% installed.packages() == FALSE ) install.packages('BiocManager')

library(BiocManager)

if ( 'cerebroApp' %in% installed.packages() == FALSE ) install('cerebroApp')
if ( 'SingleR' %in% installed.packages() == FALSE ) install('SingleR')
if ( 'monocle' %in% installed.packages() == FALSE ) install('monocle')

library(cerebroApp)

exportFromSeurat(
  reference,
  assay = "integrated",
  slot = "data",
  file = "Macroverse_Try.crb",
  experiment_name = "Macroverse",
  organism = 'hg',
  groups = c('Tissue','Phenograph','Global Healthy-1_Cancer-2_Other-3','Condition per tissue'),
  cell_cycle = NULL,
  nUMI = "Patient No",
  nGene = "Patient No",
  add_all_meta_data = TRUE,
  use_delayed_array = FALSE,
  verbose = FALSE
)

cerebroApp::launchCerebroV1.3(maxFileSize = 20000)

