#!/bin/bash

PHENO_FILE_EUR="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/pheno/trait/White_British"
COVAR_FILE_EUR="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/pheno/covar/covars_white_british_final.tsv"
GENO_DIR_EUR="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/geno/White_British/0_sample_qc"
OUTDIR_EUR="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/gwas/White_British"
PHENO_FILE_EAS="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/pheno/trait/Chinese"
COVAR_FILE_EAS="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/pheno/covar/covars_chinese_final.tsv"
GENO_DIR_EAS="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/geno/Chinese/0_sample_qc"
OUTDIR_EAS="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/gwas/Chinese"
PLINK="/data1/jiapl_group/lishuhua/software/general/plink2"
MAX_JOBS=8
THREADS_PER_JOB=24

# while read -r PHENO_NAME; do
#     echo "Processing phenotype: ${PHENO_NAME}"
#     JOBS=0
#     for CHR in {1..22}; do
#         echo "Submitting ${CHR} for ${PHENO_NAME} (EUR)"
#         ${PLINK} --bfile ${GENO_DIR_EUR}/chr${CHR} --pheno ${PHENO_FILE_EUR}/${PHENO_NAME}.txt --pheno-name ${PHENO_NAME} --covar ${COVAR_FILE_EUR} --covar-name age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 --glm hide-covar no-firth cols=+a1freq --threads ${THREADS_PER_JOB} --no-input-missing-phenotype --covar-variance-standardize --out ${OUTDIR_EUR}/${PHENO_NAME}_chr${CHR} &
#         JOBS=$((JOBS+1))
#         if [ ${JOBS} -ge ${MAX_JOBS} ]; then
#             wait
#             JOBS=0
#         fi
#     done
#     wait
#     echo "Finished ${PHENO_NAME} GWAS (EUR)"
# done <<< $'p20116_int\np20117_int'

while read -r PHENO_NAME; do
    echo "Processing phenotype: ${PHENO_NAME}"
    JOBS=0
    for CHR in {1..22}; do
        echo "Submitting ${CHR} for ${PHENO_NAME} (EAS)"
        ${PLINK} --bfile ${GENO_DIR_EAS}/chr${CHR} --pheno ${PHENO_FILE_EAS}/${PHENO_NAME}.txt --pheno-name ${PHENO_NAME} --covar ${COVAR_FILE_EAS} --covar-name age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --glm hide-covar no-firth cols=+a1freq --threads ${THREADS_PER_JOB} --no-input-missing-phenotype --covar-variance-standardize --out ${OUTDIR_EAS}/${PHENO_NAME}_${CHR} &
        JOBS=$((JOBS+1))
        if [ ${JOBS} -ge ${MAX_JOBS} ]; then
            wait
            JOBS=0
        fi
    done
    wait
    echo "Finished ${PHENO_NAME} GWAS (EAS)"
done <<< $'p50_int\np20116_int\np20117_int'