#!/bin/bash

#Get Kinship matrix
./gemma-0.98.5-linux-static-AMD64 -g gemma -g bim_Ast_layer1.txt.gz -p amyloid_phenotype.txt -c cov_amyloid.txt -gk 2 -notsnp -o cov_amyloid_Ast_layer1

./gemma-0.98.5-linux-static-AMD64 -g gemma -g bim_Ast_layer2.txt.gz -p amyloid_phenotype.txt -c cov_amyloid.txt -gk 2 -notsnp -o cov_amyloid_Ast_layer2

./gemma-0.98.5-linux-static-AMD64 -g gemma -g bim_Ast_layer3.txt.gz -p amyloid_phenotype.txt -c cov_amyloid.txt -gk 2 -notsnp -o cov_amyloid_Ast_layer3

./gemma-0.98.5-linux-static-AMD64 -g gemma -g bim_Ast_layer4.txt.gz -p amyloid_phenotype.txt -c cov_amyloid.txt -gk 2 -notsnp -o cov_amyloid_Ast_layer4

./gemma-0.98.5-linux-static-AMD64 -g gemma -g bim_Ast_layer5.txt.gz -p amyloid_phenotype.txt -c cov_amyloid.txt -gk 2 -notsnp -o cov_amyloid_Ast_layer5

./gemma-0.98.5-linux-static-AMD64 -g gemma -g bim_Ast_layer6.txt.gz -p amyloid_phenotype.txt -c cov_amyloid.txt -gk 2 -notsnp -o cov_amyloid_Ast_layer6

./gemma-0.98.5-linux-static-AMD64 -g gemma -g bim_Ast_layer7.txt.gz -p amyloid_phenotype.txt -c cov_amyloid.txt -gk 2 -notsnp -o cov_amyloid_Ast_layer7


#Conduct DGE for each layer
./gemma-0.98.5-linux-static-AMD64 -g gemma -g bim_Ast_layer1.txt.gz -p amyloid_phenotype.txt -k output/cov_amyloid_Ast_layer1.sXX.txt -c cov_amyloid.txt -lmm 4 -notsnp -o amyloid_Ast_layer1

./gemma-0.98.5-linux-static-AMD64 -g gemma -g bim_Ast_layer2.txt.gz -p amyloid_phenotype.txt -k output/cov_amyloid_Ast_layer2.sXX.txt -c cov_amyloid.txt -lmm 4 -notsnp -o amyloid_Ast_layer2

./gemma-0.98.5-linux-static-AMD64 -g gemma -g bim_Ast_layer3.txt.gz -p amyloid_phenotype.txt -k output/cov_amyloid_Ast_layer3.sXX.txt -c cov_amyloid.txt -lmm 4 -notsnp -o amyloid_Ast_layer3

./gemma-0.98.5-linux-static-AMD64 -g gemma -g bim_Ast_layer4.txt.gz -p amyloid_phenotype.txt -k output/cov_amyloid_Ast_layer4.sXX.txt -c cov_amyloid.txt -lmm 4 -notsnp -o amyloid_Ast_layer4

./gemma-0.98.5-linux-static-AMD64 -g gemma -g bim_Ast_layer5.txt.gz -p amyloid_phenotype.txt -k output/cov_amyloid_Ast_layer5.sXX.txt -c cov_amyloid.txt -lmm 4 -notsnp -o amyloid_Ast_layer5

./gemma-0.98.5-linux-static-AMD64 -g gemma -g bim_Ast_layer6.txt.gz -p amyloid_phenotype.txt -k output/cov_amyloid_Ast_layer6.sXX.txt -c cov_amyloid.txt -lmm 4 -notsnp -o amyloid_Ast_layer6

./gemma-0.98.5-linux-static-AMD64 -g gemma -g bim_Ast_layer7.txt.gz -p amyloid_phenotype.txt -k output/cov_amyloid_Ast_layer7.sXX.txt -c cov_amyloid.txt -lmm 4 -notsnp -o amyloid_Ast_layer7


done


exit


