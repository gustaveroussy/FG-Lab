source('Seurat_Heatmap_V3.R')
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ezTools)
              
run_Seurat_heatmap( TPM_file="../TPM.txt",                             # TPM_file.txt
                    group_file="../Group.txt",                          # Sample_file.txt (cam be based on cluster file)
                    group_col="CellType",  # Name of the concerned Sample_file.txt column
                    list.gene.file='GeneList.txt',  # Can be DEGs file (list.gene.file=NULL, keep all the genes in tpm.txt)
                    group.order=c('Macro','cMo','CD16posMono','cDC1','cDC2','MigDC'),                      # Change the group order (group.order=NULL, keep the default sorted order)
                    base_name='Heatmap',                                    # Name given to all output files
                    ifpng=T,
                    mt="^MT-", 
                    scale=F,
                    scale_factor=1e2,
                    add_missing_gene=FALSE,
                    col_blue_white_red =TRUE)
