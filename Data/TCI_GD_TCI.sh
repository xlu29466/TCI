#!/bin/bash

# vairables of file pathnames
SGAMatrix=./SGAmatrix_brca.csv
DEGMatrix=./DEGmatrix_brca.csv
globalPrior=./strcpriors.poplv.brca.csv
SGAprior=./strcpriors.tumorspecific.brca.csv

dictPrefix=./globDriver.dict

# search population-wide drivers for each phenotype in DEGMatrix
../TCI_GenerateGlobalDriver/TCI_GD -p $globalPrior -f $SGAMatrix -d $DEGMatrix -o "$dictPrefix.csv";

# TCI analysis
mkdir -p ./Results

../TCI/TCI -p $SGAprior -f $SGAMatrix -d $DEGMatrix -g $dictPrefix".csv" -o ./Results



