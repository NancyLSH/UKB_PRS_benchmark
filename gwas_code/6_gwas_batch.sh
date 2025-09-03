#!/bin/bash

PHENO_FILE="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/pheno/trait/White_British"
# PHENO_LIST="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/pheno/trait/trait_list_white_british.txt"
PHENO_LIST="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/pheno/trait/trait_list_white_british_2.txt"
COVAR_FILE="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/pheno/covar/covars_white_british_final.tsv"
GENO_DIR="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/geno/White_British/0_sample_qc"
OUTDIR="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/gwas/White_British"
PLINK="/data1/jiapl_group/lishuhua/software/general/plink2"
MAX_JOBS=8
THREADS_PER_JOB=24

while read -r PHENO_NAME; do
    echo "Processing phenotype: ${PHENO_NAME}"
    JOBS=0
    for CHR in {1..22}; do
        echo "Submitting ${CHR} for ${PHENO_NAME}"
        ${PLINK} --bfile ${GENO_DIR}/chr${CHR} --pheno ${PHENO_FILE}/${PHENO_NAME}.txt --pheno-name ${PHENO_NAME} --covar ${COVAR_FILE} --covar-name age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 --glm hide-covar cols=+a1freq --threads ${THREADS_PER_JOB} --no-input-missing-phenotype --covar-variance-standardize --out ${OUTDIR}/${PHENO_NAME}_${CHR} &
        JOBS=$((JOBS+1))
        if [ ${JOBS} -ge ${MAX_JOBS} ]; then
            wait
            JOBS=0
        fi
    done
    wait
    echo "Finished ${PHENO_NAME} GWAS"
done < "${PHENO_LIST}"
