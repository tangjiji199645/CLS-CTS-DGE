library(data.table)

#load phenotype
phenotype<-read.table("phenotype.txt",header=TRUE)
phenotype$pmi<-as.numeric(phenotype$pmi)
#load count for each celltype
load("count_cell_cpm.rdata")

#prepare the data for conducting gemma
phenotype$Study<-ifelse(phenotype$study=="ROS",1,0)

phenotype_cov=phenotype[,c("msex","age_death","pmi","Study")]

phenotype_bimbam=cbind(rep(1,nrow(phenotype_cov)),phenotype_cov)

write.table(phenotype_bimbam,"gemma/cov_phenotype.txt",row.names=F, col.names = F, quote = F, sep = "\t")

write.table(phenotype$amyloid_sqrt,"gemma/phenotype.txt",row.names=F, col.names = F, quote = F, sep = "\t")


#prepare the expression level for each cell type
cell_type<-dimnames(count_cell_cpm)[[3]]

bim_celltype<-function(i){
  celltype_bimbam <-data.frame(genename=rownames(count_cell_cpm),A="A",T="T",count_cell_cpm[,,i])
  write.table(celltype_bimbam,paste0("gemma/bim_",cell_type[i],".txt"),row.names=F, col.names = F, quote = F, sep = "\t")
  
}

bim_celltype(1)
bim_celltype(2)
bim_celltype(3)
bim_celltype(4)
bim_celltype(5)
bim_celltype(6)




#load count for each celltype and layer

load("Ast_layer_count_cpm.rdata")
load("Ex_layer_count_cpm.rdata")
load("In_layer_count_cpm.rdata")
load("Micro_layer_count_cpm.rdata")
load("Oli_layer_count_cpm.rdata")
load("Opc_layer_count_cpm.rdata")

#prepare the expression level for each cell type, layer
bim_celltype_layer<-function(data,celltype){
  layer_bimbam_1 <-data.frame(genename=rownames(data),A="A",T="T",data[,,1])
  write.table(layer_bimbam_1,paste0("gemma/bim_",celltype,"_layer1.txt"),row.names=F, col.names = F, quote = F, sep = "\t")
  
  layer_bimbam_2 <-data.frame(genename=rownames(data),A="A",T="T",data[,,2])
  write.table(layer_bimbam_2,paste0("gemma/bim_",celltype,"_layer2.txt"),row.names=F, col.names = F, quote = F, sep = "\t")
  
  layer_bimbam_3 <-data.frame(genename=rownames(data),A="A",T="T",data[,,3])
  write.table(layer_bimbam_3,paste0("gemma/bim_",celltype,"_layer3.txt"),row.names=F, col.names = F, quote = F, sep = "\t")
  
  layer_bimbam_4 <-data.frame(genename=rownames(data),A="A",T="T",data[,,4])
  write.table(layer_bimbam_4,paste0("gemma/bim_",celltype,"_layer4.txt"),row.names=F, col.names = F, quote = F, sep = "\t")
  
  layer_bimbam_5 <-data.frame(genename=rownames(data),A="A",T="T",data[,,5])
  write.table(layer_bimbam_5,paste0("gemma/bim_",celltype,"_layer5.txt"),row.names=F, col.names = F, quote = F, sep = "\t")
  
  layer_bimbam_6 <-data.frame(genename=rownames(data),A="A",T="T",data[,,6])
  write.table(layer_bimbam_6,paste0("gemma/bim_",celltype,"_layer6.txt"),row.names=F, col.names = F, quote = F, sep = "\t")
  
  layer_bimbam_7 <-data.frame(genename=rownames(data),A="A",T="T",data[,,7])
  write.table(layer_bimbam_7,paste0("gemma/bim_",celltype,"_layer7.txt"),row.names=F, col.names = F, quote = F, sep = "\t")
  
}



bim_celltype_layer(Ast_layer_count_cpm,"Ast")
bim_celltype_layer(Ex_layer_count_cpm,"Ex")
bim_celltype_layer(In_layer_count_cpm,"In")
bim_celltype_layer(Micro_layer_count_cpm,"Mic")
bim_celltype_layer(Oli_layer_count_cpm,"Oli")
bim_celltype_layer(Opc_layer_count_cpm,"Opc")




