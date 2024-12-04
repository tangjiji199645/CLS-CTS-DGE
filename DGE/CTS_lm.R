library(MASS)
library(sfsmisc)
library(data.table)
library(R.utils)

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

#load phenotype
cov_amyloid<-read.table("phenotype.txt",header=TRUE)
cov_amyloid$pmi<-as.numeric(cov_amyloid$pmi)



#count for each cell

load("count_cell_cpm.rdata")




output_Ast=linear_reg(cov_amyloid[,1:4],count_cell_cpm[,,1],cov_amyloid$amyloid_sqrt)
output_Ex=linear_reg(cov_amyloid[,1:4],count_cell_cpm[,,2],cov_amyloid$amyloid_sqrt)
output_In=linear_reg(cov_amyloid[,1:4],count_cell_cpm[,,3],cov_amyloid$amyloid_sqrt)
output_Mic=linear_reg(cov_amyloid[,1:4],count_cell_cpm[,,4],cov_amyloid$amyloid_sqrt)
output_Oli=linear_reg(cov_amyloid[,1:4],count_cell_cpm[,,5],cov_amyloid$amyloid_sqrt)
output_Opc=linear_reg(cov_amyloid[,1:4],count_cell_cpm[,,6],cov_amyloid$amyloid_sqrt)




write.table(output_Ast,"output/output_Ast.txt", row.names=F, col.names = T, quote = F, sep = "\t")

write.table(output_Ex,"output/output_Ex.txt", row.names=F, col.names = T, quote = F, sep = "\t")

write.table(output_In,"output/output_In.txt", row.names=F, col.names = T, quote = F, sep = "\t")

write.table(output_Mic,"output/output_Mic.txt", row.names=F, col.names = T, quote = F, sep = "\t")

write.table(output_Oli,"output/output_Oli.txt", row.names=F, col.names = T, quote = F, sep = "\t")

write.table(output_Opc,"output/output_Opc.txt", row.names=F, col.names = T, quote = F, sep = "\t")







