library(Matrix)

mtx <- readMM(file = "GSM3972009_69_matrix.mtx")
mtx <- as.data.frame(mtx)

col <- read.table("GSM3972009_69_barcodes.tsv.gz")
row <- read.table("GSM3972009_69_genes.tsv.gz")

colnames(mtx) <- col$V1
rownames(mtx) <- row$V1
mtx$gene_name <- row$V2

# remove duplicated genes
library(gdata)
dup <- unique(mtx$gene_name[duplicated(mtx$gene_name)])
tab_dup <- lapply(dup,function(l){
  tab <- mtx[mtx$gene_name %in% l,-dim(mtx)[2]]
  tab <- tab[which.max(rowSums(tab)),]
  tab$gene_name <- as.character(l)
  return(tab)
})
tab_dup <- do.call(rbind,tab_dup)
tab_non_dup <- mtx[!mtx$gene_name %in% dup,]
tab <- rbind(tab_non_dup,tab_dup)

rownames(tab) <- tab$gene_name
tab <- tab[,-dim(tab)[2]]

write.table(tab,"GSM69_TPM.txt", col.names=NA, row.names=T, quote=F, sep='\t')