#!/bin/bash

# File path
SGAori=./DataSource/snp5e4wBP.csv # SGA matrix
DEGori=./DataSource/BPcommon3.csv # DEG matrix
globalPrior=./DataSource/global_prior.csv # global prior
prior=./DataSource/individual_prior.csv # tumor specific prior

# TCI global drivers
dictPrefix=./DataSource/dict
mkdir -p ./DataSource
./TDIC_GD_exeOMP -p $globalPrior -f $SGAori -d $DEGori -o "$dictPrefix.csv";

# TCI tumor specific drivers
mkdir -p ./Results/TDI
./PanCanTDICexeOOMP -p $prior -f $SGAori -d $DEGori -g $dictPrefix".csv" -o ./Results/TDI

