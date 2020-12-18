CreateSeqGeq <- function(TPM, Dim_Red, Header, basename){
  
  # Take the cells that are in common
  common_cell <- intersect(row.names(Dim_Red),colnames(TPM))
  TPM <- TPM[,colnames(TPM) %in% common_cell]
  Dim_Red <- Dim_Red[rownames(Dim_Red) %in% common_cell,]
  
  # Transpose the GLOBAL file before merging 
  TR_Global <- t(Dim_Red)
  
  # Merge with the TPM table
  
  Aymeric_merge <- rbind(TPM[, sort(colnames(TPM))], TR_Global[, sort(colnames(TPM))])
  
  # Merge with the TPM table
 
  Aymeric_merge[dim(Aymeric_merge)[1]+1,] <- colnames(Aymeric_merge)
  rownames(Aymeric_merge)[dim(Aymeric_merge)[1]] <- 'CellID'
  Aymeric_merge <- Aymeric_merge[c(nrow(Aymeric_merge), 1:(nrow(Aymeric_merge) - 1)),] # this puts the last col with the gene names at the beginning
  
  # Merge with Header
  
  nb_of_col_to_add <- ncol(Aymeric_merge) - ncol(Header)
  m <- matrix(NA, ncol = nb_of_col_to_add, nrow = nrow(Header))
  m<- data.frame(m)
  Header_w_NA <- cbind(Header,m)
  colnames(Header_w_NA) <- colnames(Aymeric_merge)
  
  SG <- rbind(Header_w_NA,Aymeric_merge)
  
  
  write.table(SG, basename, na ="", sep="\t",row.names = T,col.names =F )

} 
