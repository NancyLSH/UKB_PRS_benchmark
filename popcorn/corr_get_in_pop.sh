#!/bin/bash

# 设置输出文件名
# OUTPUT_FILE="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/CAS/gwas/corr/genetic_correlation_summary.csv"
OUTPUT_FILE="/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/merged_gwas/White_British/corr/genetic_correlation_summary.csv"

# 为CSV文件写入表头
# 这个表头是根据log文件中的列名手动整理的
echo "source_file,p1,p2,rg,se,z,p,h2_obs,h2_obs_se,h2_int,h2_int_se,gcov_int,gcov_int_se" > "$OUTPUT_FILE"

# 遍历当前目录下所有以 .log 结尾的文件
# for LOG_FILE in /data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/CAS/gwas/corr/*.log; do
for LOG_FILE in /data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/merged_gwas/White_British/corr/*.log; do
  # 检查文件是否存在，以防没有匹配的log文件而出错
  if [ -f "$LOG_FILE" ]; then
    # 使用 awk 工具来精准提取数据
    # /Summary.../ { ... } 是一个模式匹配，当awk找到包含"Summary..."的行时，执行花括号内的命令
    # 第一个 getline 命令会读取并跳过结果的表头行
    # 第二个 getline 命令会读取我们真正需要的数据行
    # BEGIN { OFS="," } 设置输出字段的分隔符为逗号
    # 最后，打印来源文件名(FILENAME)和数据行的12个字段
    awk '
      BEGIN { OFS="," }
      /Summary of Genetic Correlation Results/ {
        getline; # 读取并跳过标题行
        getline; # 读取数据行
        # 打印来源文件名和该行的12个字段
        print FILENAME, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12;
      }
    ' "$LOG_FILE" >> "$OUTPUT_FILE"
  fi
done

echo "提取完成！"
echo "所有结果已汇总到文件: $OUTPUT_FILE"