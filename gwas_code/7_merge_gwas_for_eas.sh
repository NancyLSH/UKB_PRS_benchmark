#!/bin/bash

# ==============================================================================
#           A Robust Script to Merge Per-Chromosome GWAS Results
# ==============================================================================
#
# This script iterates through a list of phenotypes, checks for the existence
# of all 22 per-chromosome GWAS files, and merges them if all are present.
#

# --- 1. CONFIGURATION SECTION (请根据您的实际情况修改此部分) ---

# 定义您的性状列表 (用空格分隔)
# 例如: TRAIT_LIST=("BMI" "Height" "T2D_risk" "HDL")
TRAIT_LIST=("p48_int" "p50_int" "p102_int" "p4079_int" "p4080_int" "p21001_int" "p21002_int" "p20116_int" "p20117_int" "p30000_int" "p30010_int" "p30020_int" "p30080_int" "p30120_int" "p30130_int" "p30140_int" "p30150_int" "p30620_int" "p30650_int" "p30670_int" "p30690_int" "p30700_int" "p30730_int" "p30740_int" "p30760_int" "p30780_int" "p30870_int" "p30880_int")

# 定义您的GWAS结果文件存放的目录
GWAS_DIR="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/gwas/Chinese"

# 定义合并后文件的输出目录
OUTPUT_DIR="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/merged_gwas/Chinese"

# 定义您的文件名结构 (这是关键！)
# 脚本会用 性状名 和 染色体编号 来拼凑出完整的文件名
# 例如，对于BMI的1号染色体，文件名是: ${GWAS_DIR}/BMI_chr1.BMI.glm.linear
FILE_PREFIX="" # 如果您的文件名开头没有共同前缀，此项可留空，例如 "BMI_chr1..."
FILE_SUFFIX="_" # 性状名和染色体编号之间的部分
FILE_EXTENSION=".glm.linear" # 文件扩展名
FILE_EXTENSION2=".glm.logistic"

# --- 2. SCRIPT LOGIC (通常无需修改此部分) ---

# 确保输出目录存在
mkdir -p "$OUTPUT_DIR"

echo "=================================================="
echo "Starting GWAS file merge process."
echo "GWAS Directory: $GWAS_DIR"
echo "Output Directory: $OUTPUT_DIR"
echo "=================================================="

# 对性状列表中的每一个性状进行循环
for PHENO in "${TRAIT_LIST[@]}"
do
    echo # 输出一个空行，让格式更清晰
    echo "--- Processing phenotype: [ $PHENO ] ---"

    # 初始化一个标志，用来追踪文件是否齐全
    all_files_exist=true

    # **安全检查**: 检查1到22号染色体的文件是否都存在
    echo "Checking for presence of all 22 chromosome files..."
    for i in {1..22}
    do
        # 根据您定义的规则拼凑出预期的文件名
        if [ "$PHENO" == "p20116_int" ] || [ "$PHENO" == "p20117_int" ]; then
            FILE_SUFFIX="_chr"
            FINAL_FILE_EXTENSION="${FILE_EXTENSION2}"
        else
            FILE_SUFFIX="_chr"
            FINAL_FILE_EXTENSION="${FILE_EXTENSION}"
        fi
        EXPECTED_FILE="${GWAS_DIR}/${FILE_PREFIX}${PHENO}${FILE_SUFFIX}${i}.${PHENO}${FINAL_FILE_EXTENSION}"

        # 使用 [ ! -f "文件名" ] 来判断文件是否存在
        if [ ! -f "$EXPECTED_FILE" ]; then
            echo "  [ERROR] Missing file: $EXPECTED_FILE"
            all_files_exist=false
            break # 只要有一个文件缺失，就没必要再检查了，直接跳出这个内层循环
        fi
    done

    # **条件合并**: 只有在所有文件都存在的情况下才执行合并操作
    if [ "$all_files_exist" = true ]; then
        echo "  All 22 files found. Proceeding with merge."

        # 准备用于合并的文件名模式 (带通配符*)
        if [ "$PHENO" == "p20116_int" ] || [ "$PHENO" == "p20117_int" ]; then
            FILE_PATTERN="${GWAS_DIR}/${FILE_PREFIX}${PHENO}${FILE_SUFFIX}*.${PHENO}${FILE_EXTENSION2}"
        else
            FILE_PATTERN="${GWAS_DIR}/${FILE_PREFIX}${PHENO}${FILE_SUFFIX}*.${PHENO}${FILE_EXTENSION}"
        fi
        # 定义最终的输出文件名
        OUTPUT_FILE="${OUTPUT_DIR}/${PHENO}.merged.glm.linear"

        # 使用我们之前讨论过的强大awk单行命令进行合并
        # ls -v 能保证染色体按 1, 2, ..., 10, 11 的自然顺序排列
        awk 'FNR==1 && NR!=1 {next} {print}' $(ls -v ${FILE_PATTERN}) > "$OUTPUT_FILE"

        # 检查输出文件是否成功生成且不为空
        if [ -s "$OUTPUT_FILE" ]; then
            echo "  Successfully merged. Output file created:"
            echo "  -> $OUTPUT_FILE"
        else
            echo "  [ERROR] Merging failed for phenotype [ $PHENO ]. Output file is empty."
        fi

    else
        # 如果文件不齐全，则打印跳过信息
        echo "  Skipping merge for phenotype [ $PHENO ] due to missing files."
    fi
done

echo # 输出一个空行
echo "=================================================="
echo "Process finished."
echo "=================================================="