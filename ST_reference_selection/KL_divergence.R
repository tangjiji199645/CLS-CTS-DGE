library(philentropy)
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

#CosMx as reference
CosMx_layer<-fread("/Users/tangji/Library/CloudStorage/OneDrive-EmoryUniversity/backup/research/Spatial/github/CLS_CTS_DGE/CosMx_layer.txt",header=TRUE)

table(CosMx_layer$layer)
table(CosMx_layer$cell_type)

CosMx_layer$layer<-factor(CosMx_layer$layer,levels = c("L1","L2","L3","L4","L5","L6","WM"))
CosMx_layer$cell_type<-factor(CosMx_layer$cell_type,levels=c("Ast","Ex","Inh","Micro","OPC","oligodendrocyte"))
#Get the cell type proportion for each layer for CosMx
layer_freq_CosMx<-CosMx_layer %>%
  group_by(layer,cell_type) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

layer_freq_CosMx<-data.frame(layer_freq_CosMx)

#Inferred layer by sample 151673
layer_pred<-fread("/Users/tangji/Library/CloudStorage/OneDrive-EmoryUniversity/backup/research/Spatial/ROSMAP/ROSMAP_400/layer_predict_400_151673_all.txt",header=TRUE)
  
layer_pred$cell.type<-factor(layer_pred$cell.type,levels=c("Astrocyte","Excitatory Neurons","Inhibitory Neurons","Microglia","OPCs","Oligodendrocytes"))
#Get the cell type proportion for each layer by sample 151673
layer_freq<-layer_pred %>%
    group_by(layer,cell.type,.drop = FALSE) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
layer_freq<-data.frame(layer_freq)
layer_freq$layer<-as.factor(layer_freq$layer)
  
#Get KL divergence between sample 151673 and CosMx from L1 to L6 and white matter
dis_73<-sapply(1:7,function(i){KL(rbind(layer_freq_CosMx$freq[((i-1)*6+1):(6*i)],layer_freq$freq[((i-1)*6+1):(6*i)]))})

