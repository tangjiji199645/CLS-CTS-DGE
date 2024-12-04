library(MASS)
library(sfsmisc)
library(data.table)
library(R.utils)
library(edgeR)

#function to conduct DGE by linear regression


linear_reg<-function(cov,RNA,phenotype){
  cov_stat <- matrix(NA, nrow = nrow(RNA), ncol = 4)
  for(i in 1:nrow(RNA)){
    if(i %% 1000 == 0) {
      print(paste("Percentage ", 100 * i/nrow(RNA), "%") )
    }
    fit_data<-data.frame(y=phenotype, raw_RNA=unlist(RNA[i,]),cov)
    fit_mod = lm(y ~., data = fit_data)
    fit=summary(fit_mod)
    cov_stat[i,] = fit$coefficients[2,]
  }
  
  colnames(cov_stat)<-c("beta","se","stat","pval")
  
  return(cov_stat)
}





############################################################################
############################################################################

#load count for each cell type and layer

load("Ast_layer_count_cpm.rdata")
load("Ex_layer_count_cpm.rdata")
load("In_layer_count_cpm.rdata")
load("Micro_layer_count_cpm.rdata")
load("Oli_layer_count_cpm.rdata")
load("Opc_layer_count_cpm.rdata")


cov_amyloid<-read.table("phenotype.txt",header=TRUE)
cov_amyloid$pmi<-as.numeric(cov_amyloid$pmi)


reg_celltype_layer<-function(count,celltype,cov,phenotype){
  
  layer<-c("1","2","3","4","5","6","WM")
  output=c()
  for (i in 1:7){
    output_layer=linear_reg(cov[,1:4],count[,,i],phenotype)
    output_layer=as.data.frame(output_layer)
    output_layer$layer=layer[i]
    output_layer$gene_name=rownames(count)
    output_layer$cell_type=celltype
    output<-rbind(output,output_layer)
    
  }
  
  return(output)
}

output_Ast_layer<-reg_celltype_layer(Ast_layer_count_cpm,"Ast",cov_amyloid,cov_amyloid$amyloid_sqrt)

output_Ex_layer<-reg_celltype_layer(Ex_layer_count_cpm,"Ex",cov_amyloid,cov_amyloid$amyloid_sqrt)

output_In_layer<-reg_celltype_layer(In_layer_count_cpm,"In",cov_amyloid,cov_amyloid$amyloid_sqrt)

output_Micro_layer<-reg_celltype_layer(Micro_layer_count_cpm,"Micro",cov_amyloid,cov_amyloid$amyloid_sqrt)

output_Oli_layer<-reg_celltype_layer(Oli_layer_count_cpm,"Oli",cov_amyloid,cov_amyloid$amyloid_sqrt)

output_Opc_layer<-reg_celltype_layer(Opc_layer_count_cpm,"Opc",cov_amyloid,cov_amyloid$amyloid_sqrt)



write.table(output_Ast_layer,"output/output_Ast_layer.txt", row.names=F, col.names = T, quote = F, sep = "\t")

write.table(output_Ex_layer,"output/output_Ex_layer.txt", row.names=F, col.names = T, quote = F, sep = "\t")

write.table(output_In_layer,"output/output_In_layer.txt", row.names=F, col.names = T, quote = F, sep = "\t")

write.table(output_Micro_layer,"output/output_Micro_layer.txt", row.names=F, col.names = T, quote = F, sep = "\t")

write.table(output_Oli_layer,"output/output_Oli_layer.txt", row.names=F, col.names = T, quote = F, sep = "\t")

write.table(output_Opc_layer,"output/output_Opc_layer.txt", row.names=F, col.names = T, quote = F, sep = "\t")




