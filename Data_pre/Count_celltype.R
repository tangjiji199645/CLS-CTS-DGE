library(data.table)
library(Seurat)
library(BisqueRNA)

#setwd("/home/jyang51/YangLabData/stang/Layer/sample_data_code/")
###load cell type annotation with layer######
cell_annotation<-fread("cell_annotation_10sample.txt",header=TRUE)

Individual_id<-unique(cell_annotation$individualID)

###get expression matrix and cell type annotation for each cell type#####
expression_matrix_Ast<-c()
expression_matrix_Ex<-c()
expression_matrix_In<-c()
expression_matrix_Micro<-c()
expression_matrix_Oli<-c()
expression_matrix_Opc<-c()

cell_annotation_Ast<-c()
cell_annotation_Ex<-c()
cell_annotation_In<-c()
cell_annotation_Micro<-c()
cell_annotation_Oli<-c()
cell_annotation_Opc<-c()

for (i in 1:length(Individual_id)){
  load(paste0(Individual_id[i],"_expression_matrix.rdata"))
  
  cell_annotation_sample<-cell_annotation[cell_annotation$individualID==Individual_id[i],]
  
  print(i)
  #Cell annotation by cell type
  cell_annotation_Ast_1<-cell_annotation_sample[cell_annotation_sample$cell.type=="Astrocyte",]
  cell_annotation_Ex_1<-cell_annotation_sample[cell_annotation_sample$cell.type=="Excitatory Neurons",]
  cell_annotation_In_1<-cell_annotation_sample[cell_annotation_sample$cell.type=="Inhibitory Neurons",]
  cell_annotation_Micro_1<-cell_annotation_sample[cell_annotation_sample$cell.type=="Microglia",]
  cell_annotation_Oli_1<-cell_annotation_sample[cell_annotation_sample$cell.type=="Oligodendrocytes",]
  cell_annotation_Opc_1<-cell_annotation_sample[cell_annotation_sample$cell.type=="OPCs",]
  
  
  
  cell_annotation_Ast<-rbind(cell_annotation_Ast,cell_annotation_Ast_1)
  cell_annotation_Ex<-rbind(cell_annotation_Ex,cell_annotation_Ex_1)
  cell_annotation_In<-rbind(cell_annotation_In,cell_annotation_In_1)
  cell_annotation_Micro<-rbind(cell_annotation_Micro,cell_annotation_Micro_1)
  cell_annotation_Oli<-rbind(cell_annotation_Oli,cell_annotation_Oli_1)
  cell_annotation_Opc<-rbind(cell_annotation_Opc,cell_annotation_Opc_1)
  
  
  
  #Gene expression by cell type
  expression_matrix_Ast_1<-expression_matrix[,cell_annotation_Ast_1$cellBarcode]
  expression_matrix_Ex_1<-expression_matrix[,cell_annotation_Ex_1$cellBarcode]
  expression_matrix_In_1<-expression_matrix[,cell_annotation_In_1$cellBarcode]
  expression_matrix_Micro_1<-expression_matrix[,cell_annotation_Micro_1$cellBarcode]
  expression_matrix_Oli_1<-expression_matrix[,cell_annotation_Oli_1$cellBarcode]
  expression_matrix_Opc_1<-expression_matrix[,cell_annotation_Opc_1$cellBarcode]
  
  
  expression_matrix_Ast<-cbind(expression_matrix_Ast,expression_matrix_Ast_1)
  expression_matrix_Ex<-cbind(expression_matrix_Ex,expression_matrix_Ex_1)
  expression_matrix_In<-cbind(expression_matrix_In,expression_matrix_In_1)
  expression_matrix_Micro<-cbind(expression_matrix_Micro,expression_matrix_Micro_1)
  expression_matrix_Oli<-cbind(expression_matrix_Oli,expression_matrix_Oli_1)
  expression_matrix_Opc<-cbind(expression_matrix_Opc,expression_matrix_Opc_1)
  
  
  
}


save(cell_annotation_Ast, file="cell_annotation_Ast.rdata")
save(cell_annotation_Ex, file="cell_annotation_Ex.rdata")
save(cell_annotation_In, file="cell_annotation_In.rdata")
save(cell_annotation_Micro, file="cell_annotation_Micro.rdata")
save(cell_annotation_Oli, file="cell_annotation_Oli.rdata")
save(cell_annotation_Opc, file="cell_annotation_Opc.rdata")



save(expression_matrix_Ast, file="expression_matrix_Ast.rdata")
save(expression_matrix_Ex, file="expression_matrix_Ex.rdata")
save(expression_matrix_In, file="expression_matrix_In.rdata")
save(expression_matrix_Micro, file="expression_matrix_Micro.rdata")
save(expression_matrix_Oli, file="expression_matrix_Oli.rdata")
save(expression_matrix_Opc, file="expression_matrix_Opc.rdata")


###count the read counts for each sample on cell type#####

count_cell<-array(0, dim = c(36601,10,6))
celltype<-c("Ast","Ex","In","Micro","Oli","Opc")


rownames(count_cell)<-rownames(expression_matrix_Ast)
colnames(count_cell)<-Individual_id
dimnames(count_cell)[[3]] =celltype


scRNA_count<-function(ScRNA,meta_data,individualID){
  cts_cell<-array(0, dim = c(36601,length(individualID)))
  for(j in 1:length(individualID)) {
    
    print(j)
    cts_cell[,j] =apply(as.matrix(ScRNA[, meta_data$individualID == individualID[j]]),1,sum)

  }
  
  return(cts_cell)
}

scRNA_count_Ast<-scRNA_count(expression_matrix_Ast,cell_annotation_Ast,Individual_id)
scRNA_count_Ex<-scRNA_count(expression_matrix_Ex,cell_annotation_Ex,Individual_id)
scRNA_count_In<-scRNA_count(expression_matrix_In,cell_annotation_In,Individual_id)
scRNA_count_Micro<-scRNA_count(expression_matrix_Micro,cell_annotation_Micro,Individual_id)
scRNA_count_Oli<-scRNA_count(expression_matrix_Oli,cell_annotation_Oli,Individual_id)
scRNA_count_Opc<-scRNA_count(expression_matrix_Opc,cell_annotation_Opc,Individual_id)

count_cell[,,1]<-scRNA_count_Ast
count_cell[,,2]<-scRNA_count_Ex
count_cell[,,3]<-scRNA_count_In
count_cell[,,4]<-scRNA_count_Micro
count_cell[,,5]<-scRNA_count_Oli
count_cell[,,6]<-scRNA_count_Opc

save(count_cell, file="count_celltype.rdata")

###normalization#####

cpm_scRNA<-function(sc){
  
  cpm_sc<-sc
  
  for (i in 1:dim(sc)[3]){
    cpm_sc[,,i]<-cpm(sc[,,i]+1,log=TRUE)
  }
  
  return(cpm_sc)
}


count_cell_cpm<-cpm_scRNA(count_cell)


save(count_cell_cpm, file="count_cell_cpm.rdata")




###count the read counts for each sample on cell type, layer#####

scRNA_count<-function(ScRNA,meta_data,individualID){
  cts_cell<-array(0, dim = c(nrow(ScRNA),length(individualID),7))
  for(j in 1:length(individualID)) {
    print(j)
    cts_cell[,j,1] =rowSums(as.matrix(ScRNA[, meta_data$individualID == individualID[j] &meta_data$layer ==1]))
    cts_cell[,j,2] =rowSums(as.matrix(ScRNA[, meta_data$individualID == individualID[j] &meta_data$layer ==2]))
    cts_cell[,j,3] =rowSums(as.matrix(ScRNA[, meta_data$individualID == individualID[j] &meta_data$layer ==3]))
    cts_cell[,j,4] =rowSums(as.matrix(ScRNA[, meta_data$individualID == individualID[j] &meta_data$layer ==4]))
    cts_cell[,j,5] =rowSums(as.matrix(ScRNA[, meta_data$individualID == individualID[j] &meta_data$layer ==5]))
    cts_cell[,j,6] =rowSums(as.matrix(ScRNA[, meta_data$individualID == individualID[j] &meta_data$layer ==6]))
    cts_cell[,j,7] =rowSums(as.matrix(ScRNA[, meta_data$individualID == individualID[j] &meta_data$layer ==7]))
    
    
    #cts[,j,k] = cpm(cts[,j,k])
    
  }
  rownames(cts_cell)<-rownames(ScRNA)
  colnames(cts_cell)<-individualID
  dimnames(cts_cell)[[3]]<-c("1","2","3","4","5","6","7")
  return(cts_cell)
}

Ast_layer_count<-scRNA_count(expression_matrix_Ast,cell_annotation_Ast,Individual_id)
Ex_layer_count<-scRNA_count(expression_matrix_Ex,cell_annotation_Ex,Individual_id)
In_layer_count<-scRNA_count(expression_matrix_In,cell_annotation_In,Individual_id)
Micro_layer_count<-scRNA_count(expression_matrix_Micro,cell_annotation_Micro,Individual_id)
Oli_layer_count<-scRNA_count(expression_matrix_Oli,cell_annotation_Oli,Individual_id)
Opc_layer_count<-scRNA_count(expression_matrix_Opc,cell_annotation_Opc,Individual_id)

#normalization
Ast_layer_count_cpm<-cpm_scRNA(Ast_layer_count)
Ex_layer_count_cpm<-cpm_scRNA(Ex_layer_count)
In_layer_count_cpm<-cpm_scRNA(In_layer_count)
Micro_layer_count_cpm<-cpm_scRNA(Micro_layer_count)
Oli_layer_count_cpm<-cpm_scRNA(Oli_layer_count)
Opc_layer_count_cpm<-cpm_scRNA(Opc_layer_count)

save(Ast_layer_count_cpm, file="Ast_layer_count_cpm.rdata")
save(Ex_layer_count_cpm, file="Ex_layer_count_cpm.rdata")
save(In_layer_count_cpm, file="In_layer_count_cpm.rdata")
save(Micro_layer_count_cpm, file="Micro_layer_count_cpm.rdata")
save(Oli_layer_count_cpm, file="Oli_layer_count_cpm.rdata")
save(Opc_layer_count_cpm, file="Opc_layer_count_cpm.rdata")





