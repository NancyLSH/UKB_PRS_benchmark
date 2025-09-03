#!/bin/bash
# This script runs GWAS for the UKB White British population using plink2.

PHENO_FILE=/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/pheno/trait/White_British/p50_int.txt
PHENO_NAME=p50_int
COVAR_FILE=/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/pheno/covar/covars_white_british_final.tsv

for CHR in {1..22}; do
    echo "Running GWAS for chromosome $CHR"
     /data1/jiapl_group/lishuhua/software/general/plink2 \
        --bfile /data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/geno/White_British/0_sample_qc/chr${CHR} \
        --pheno $PHENO_FILE \
        --pheno-name $PHENO_NAME \
        --covar $COVAR_FILE \
        --covar-name age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
        --glm hide-covar cols=+a1freq \
        --threads 24 \
        --no-input-missing-phenotype \
        --covar-variance-standardize \
        --out /data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/gwas/White_British/${PHENO_NAME}_chr${CHR}
done