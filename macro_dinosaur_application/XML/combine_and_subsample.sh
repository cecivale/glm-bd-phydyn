#!/bin/bash


# Combine and subsample chains
for analysis in GLM_sumprior_3p GLM_sumprior_12p GLM_sumprior_error_distr_3p no_GLM 
do 
    for output in log trees
    do
    logcombiner \
        -log  results/${analysis}/Dinosaurs_${analysis}_A.${output} -log results/${analysis}/Dinosaurs_${analysis}_B.${output}  \
        -log  results/${analysis}/Dinosaurs_${analysis}_C.${output} -log results/${analysis}/Dinosaurs_${analysis}_D.${output}  \
        -log  results/${analysis}/Dinosaurs_${analysis}_E.${output}  \
        -b 10 -resample 10000000  \
        -o results/Dinosaurs_${analysis}_combined.${output} 
    done
done

