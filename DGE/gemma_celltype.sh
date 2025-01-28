#!/bin/bash

#Get Kinship matrix
./gemma-0.98.5-linux-static-AMD64 -g gemma -g bim_Ast.txt.gz -p amyloid_phenotype.txt -c cov_amyloid.txt -gk 2 -notsnp -o cov_amyloid_Ast


#Conduct DGE for each layer
./gemma-0.98.5-linux-static-AMD64 -g gemma -g bim_Ast.txt.gz -p amyloid_phenotype.txt -k output/cov_amyloid_Ast.sXX.txt -c cov_amyloid.txt -lmm 4 -notsnp -o amyloid_Ast


done


exit


