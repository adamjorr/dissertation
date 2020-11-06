#!/usr/bin/env bash

BASEDIR=ajorr1@hines.scit.us:/home/ajorr1/variant-standards/CHM-eval/hg19/chr1/reads/2020-09-13/variants/
scp ${BASEDIR}/\*.snp_roc.tsv.gz .
#FPRS=$(seq 0 20 100)
#for i in $FPRS; do
    #for j in $FPRS; do
        #sleep 2
        #scp ${BASEDIR}/${i}_${j}/snp_roc.tsv.gz ./${i}_${j}.snp_roc.tsv.gz
    #done
#done
