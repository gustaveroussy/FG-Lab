library(nichenetr)
library(Seurat)
library(tidyverse)
library(ezTools)

########## (1) Choose receiver population (phenograph cluster), 
########## (2) Load RDS object for the region of interest (Tumor Core, periphery, healthy tissue)
########## and (3) select genes of interest (enriched in the phenograph of interest)

# Example here: population of interest is phenograph cluster #2, Region of interest: Tumor Core

BaseName <- "Pheno2Receiver_Core_Sharma_"

seuratObj <- readRDS("../../Sharma_TumorCore_HCC_withPhenoAnnot.RDS")
seuratObj@meta.data$CT_mye_and_Tcells_and_pheno %>% table() # Different populations of the dataset

# Decide on the receiver and sender population
receiver = "Pheno.2"

sender_celltypes = c("Tcells_B3GNT7","Tcells_CCR7Tcells","Tcells_CD4","Tcells_CRTAMTcells",
                     "Tcells_Cytotoxic","Tcells_GEMTcells","Tcells_HBA2Tcells","Tcells_Tregs")

# Decide on the genes of interest
Pheno_of_interset <- 2

DEG_list_0.25 <- fast_read_table("../../nonnorm_DEG_allcells_pheno_logFC0.25.txt")

DEG_geneset_oi <- DEG_list_0.25[DEG_list_0.25$cluster == Pheno_of_interset, 7] # we want DEG of pheno of interset, and select the column 7 which corresponds to the genes. 

DER_list <- read.table("../../DER_Scenic.txt", sep='\t', header = T,check.names = F, stringsAsFactors = F)
head(DER_list)
DER_geneset_oi <- DER_list[DER_list$cluster == Pheno_of_interset, 2] # we want DEG of pheno of interset, and the column 2 which corresponds to the genes. 

geneset_oi <- c(DEG_geneset_oi,DER_geneset_oi)

############################################################################################################################################
############################################################################################################################################

##### Load Ligand-target prior model, ligand-receptor network, and weighted integrated networks

ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

########################################## Start NicheNet ##########################################

Idents(seuratObj) <- "CT_mye_and_Tcells_and_pheno"

#### STEP 1 : Define the receiver and the sender populations

## receiver: expressed genes

expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.15)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender : expressed genes

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.15) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()



#### STEP 3 : Define a set of potential ligands: expressed by the sender niche, 
# and that bind a (putative) receptor expressed by the received/target population

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)

expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

#### STEP 4 : Perform NicheNet ligand activity analysis: rank the potentail ligands based on their target
# genes in the gene set of interset (as compared to the background set of genes)

ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes,
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
fast_save_table(ligand_activities,"",paste0(BaseName,"Ligands_activities_pearson_correlation.txt"))
ligand_activities

## Here the top 20 ligands are kept
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

# To see which cells of the sender cells are expressing the top 20 ligands
senderSrtObj <- subset(seuratObj, idents = sender_celltypes)

pdf(paste0(BaseName,"best_upstream_ligands.pdf"))
DotPlot(senderSrtObj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
dev.off()

#### STEP 5 : infer the receptors and top-predicted target genes of ligands that are top-ranked in the ligand analysis

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))

pdf(paste0(BaseName,"p_ligand_target_network.pdf"))
p_ligand_target_network
dev.off()

## Receptors of top-ranked ligands, but after conisdering only bona fide 
## ligand-receptor interactions dicuments in literature and publicly available databases

lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")

pdf(paste0(BaseName,"p_ligand_receptor_network.pdf"))
p_ligand_receptor_network
dev.off()

## Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions 
##documented in literature and publicly available databases

lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")

pdf(paste0(BaseName,"p_ligand_receptor_network_STRICT.pdf"))
p_ligand_receptor_network_strict
dev.off()


