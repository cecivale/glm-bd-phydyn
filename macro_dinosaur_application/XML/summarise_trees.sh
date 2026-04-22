#!/bin/bash

for analysis in GLM_sumprior_3p GLM_sumprior_12p GLM_sumprior_error_distr_3p no_GLM 
do
    for topology in CCD0 MCC
do
    treeannotator -topology CCD0 -height median -b 0 -lowMem \
        -file results/Dinosaurs_${analysis}_combined.trees results/Dinosaurs_${analysis}_combined.${topology}.tree
    done
done

