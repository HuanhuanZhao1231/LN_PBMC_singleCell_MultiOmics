#!/bin/bash
#PBS -N SLE_swiss
#PBS -l nodes=node1:ppn=4
#PBS -l walltime=168:00:00
#PBS -e oe
#PBS -o SLE_swiss.log

cd /node1/liyunna/coloc_lxy/eGFR/
swiss --assoc inputforswissGWAS.txt --ld-clump --ld-clump-source 1000G_2014-11_ASN --clump-ld-thresh 0.8 --clump-p 5e-05 --out myresults