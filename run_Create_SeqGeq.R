source('Create_SeqGeq_file.R')

# Convert mtx to txt
#TPM_Aymeric <- readMM("GSE131535_matrix.mtx")
#col <- read.table("GSE131535_barcodes.tsv")
#row <- read.table("GSE131535_genes.tsv")
#TPM_Aymeric <- as.matrix(TPM_Aymeric)
#colnames(TPM_Aymeric) <- col$V1
#rownames(TPM_Aymeric) <- row$V2
#write.table(TPM_Aymeric,"GSE131535_TPM.txt",quote=F, sep='\t',row.names = T, col.names=NA)

# Read the TPM or rc table
TPM_Aymeric <- read.table("tpmlungALL.txt",sep="\t",header=T,row.names=1,check.names = F)


# Read the file containing the Dim Red and other information. The first column contains the cell IDs.
Dim_Red_all <- read.table("transposed_infomat_healty lung.txt"
                          ,sep="\t",header=T,row.names=1,check.names = F
                          , stringsAsFactors = F)

SG_Header <- read.table("Header_SeqGeq.txt",sep="\t",header=F,row.names=1,check.names = F)

# the name of the output SeqGeq file 
basename = "2020_18_03_SG_Healtht lung Braga CAD.txt"

CreateSeqGeq(TPM=TPM_Aymeric, Dim_Red=Dim_Red_all, Header = SG_Header, basename=basename)
