import pandas as pd
import sys
import os

def process_single_file(file_path, output_file):
    if file_path.endswith('.glm.linear'):
        # 读取文件
        df = pd.read_csv(file_path, sep='\t', usecols=["ID", "ALT", "REF", "BETA", "SE"])
        # 重命名列
        df.columns = ['SNP', 'A1', 'A2', 'BETA', 'SE']
        # 保存文件
        df.to_csv(output_file, sep='\t', index=False, header=True)
        print(f"Processed and saved: {output_file}")
    elif file_path.endswith('.glm.logistic'):
        # 读取文件
        df = pd.read_csv(file_path, sep='\t', usecols=["ID", "ALT", "REF", "OR", "LOG(OR)_SE"])
        # 重命名列
        df.columns = ['SNP', 'A1', 'A2', 'OR', 'SE']
        # 保存文件
        df.to_csv(output_file, sep='\t', index=False, header=True)
        print(f"Processed and saved: {output_file}")
    else:
        print(f"Skipping unsupported file type: {file_path}")

def process_directory(dir_path, output_dir):
    # 读取dir_path文件夹下的所有文件
    for d in os.listdir(dir_path):
        file_path = os.path.join(dir_path, d)
        output_prefix = d.replace('.glm.linear', '').replace('.glm.logistic', '')
        output_file = os.path.join(output_dir, f"{output_prefix}.txt")
        process_single_file(file_path, output_file)

if __name__ == "__main__":
    gwas_dir = "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/CAS/gwas/gwas/"
    output_dir = "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/prscsx/train/EAS/"
    process_directory(gwas_dir, output_dir)