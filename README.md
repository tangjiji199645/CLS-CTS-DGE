# Cortical-layer specific (CLS) and Cell type specific (CTS) DGE analyses

Paper preprint: https://www.medrxiv.org/content/10.1101/2024.11.18.24317499v1
Under construction.

## Workflow for CLS-CTS DGE Analyses of AD-related phenotypes.
![alt text](https://www.medrxiv.org/content/medrxiv/early/2024/11/18/2024.11.18.24317499/F1.large.jpg?width=800&height=600&carousel=1)

A) CeLEry tool was first used to infer the spatial locations of cells from the ROS/MAP snRNA-Seq data of DLPFC. 

B) Pseudo bulk read counts were generated from the scRNA-Seq data for each cell type in each cortical layer, which were used to conduct CLS-CTS DGE analyses of three AD-related phenotypes, Beta-amyloid, tangle density, and cognitive decline. Top significant genes from the CLS-CTS DGE analyses were further used to conduct gene set enrichment analyses.

**For LMM-DGE, please see detail at https://github.com/tangjiji199645/LMM_DGE_Pipeline**

## Infer cortical layer of snRNA-Seq data using spatial transcriptome as reference
Please find code under folder Layer_prediction. 

LIBD data can be download from http://spatial.libd.org/spatialLIBD/

For CeLEry installment, please see detail at https://github.com/QihuangZhang/CeLEry/

Prepoessing for LIBD data: data_preproessing.ipynb
Training model: LIBD_train.ipynb
Predicting: LIBD_predict.ipynb

## Select the best spatial transcriptome reference by KL divergence
Please find code under folder ST_reference_selection.

Select the best spatial transcriptome reference by KL divergence: KL_divergence.R

## Prepare the data for conducting Cortical-layer specific (CLS) and Cell type specific (CTS) DGE analyses

Please find code under folder Data_pre.

Count the snRNA-Seq data into layer level: Count_celltype.R

Prepare the data for conducting LMM-DGE by Gemma: gemma_pre.R

## Conducting Cortical-layer specific (CLS) and Cell type specific (CTS) DGE analyses

Please find code under folder DGE.

Cell type specific (CTS) DGE analyses by standard linear model: CTS_lm.R

Cortical-layer specific (CLS) and Cell type specific (CTS) DGE analyses by standard linear model: CLS_CTS_lm.R

Cell type specific (CTS) DGE analyses by linear mixed model: gemma_celltype.sh

Cortical-layer specific (CLS) and Cell type specific (CTS) DGE analyses by standard linear model: gemma_celltype_layer.sh

For more detail for LMM-DGE, please see detail at https://github.com/tangjiji199645/LMM_DGE_Pipeline















