#!/bin/bash

# load the white_british and chinese list
white_british_list="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/pheno/covar/covars_white_british_qced.txt"
chinese_list="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/pheno/covar/covars_chinese_qced.txt"

for group in "white_british" "chinese"; do
    if [ "$group" == "white_british" ]; then
        list=$white_british_list
        mkdir -p /data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/geno/White_British/0_sample_qc
        for chr in {1..22}; do
            echo "Processing $group chromosome $chr"
            /data1/jiapl_group/lishuhua/software/general/plink \
                --bfile /data1/jiapl_group/lishuhua/UKB/Genotype/White_British/eur_chr${chr} \
                --keep $list \
                --geno 0.02 \
                --maf 0.01 \
                --hwe 1e-10 midp \
                --make-bed \
                --out /data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/geno/White_British/0_sample_qc/chr${chr}
        done
    else
        list=$chinese_list
        mkdir -p /data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/geno/Chinese/0_sample_qc
        for chr in {1..22}; do
            echo "Processing $group chromosome $chr"
            /data1/jiapl_group/lishuhua/software/general/plink \
                --bfile /data1/jiapl_group/lishuhua/UKB/Genotype/Chinese/chr${chr} \
                --keep $list \
                --geno 0.02 \
                --maf 0.05 \
                --hwe 1e-6 midp \
                --make-bed \
                --out /data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/geno/Chinese/0_sample_qc/chr${chr}
        done
    fi
done