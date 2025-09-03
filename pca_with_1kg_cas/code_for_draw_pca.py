import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

def plot_pca():
    """
    Reads PLINK PCA output and visualizes population stratification.
    """
    # --- 文件路径配置 ---
    eigenvec_file = '/data1/jiapl_group/lishuhua/project/PRS_benchmark/code/real_data/pca_code/1kg_ref_pca.eigenvec'
    pca_proj_file = '/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/Merge_CAS_UKB_1kg/CAS_UKB_pcs_projected.sscore'
    panel_1kg_file = '/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/1kg/20130606_g1k.ped'
    cohort_a_list = '/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/Merge_CAS_UKB_1kg/CAS_samples.txt'
    cohort_b_list = '/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/Merge_CAS_UKB_1kg/UKB_EAS_samples.txt'
    output_png = '/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/Merge_CAS_UKB_1kg/pca_plot.png'

    print("--- Starting PCA Visualization with Python ---")

    try:
        # 1. 读取PLINK计算出的PCA结果
        print(f"Reading PCA results from {eigenvec_file}...")
        pca_df = pd.read_csv(eigenvec_file, sep="\t")
        # make sure the IID column is string type
        pca_df['IID'] = pca_df['IID'].astype(str)
        pca_projected = pd.read_csv(pca_proj_file, sep="\t", header=0)
        pca_projected['IID'] = pca_projected['IID'].astype(str)
        score_cols = [col for col in pca_projected.columns if '_AVG' in col]
        rename_dict = {old_col: f'PC{i+1}' for i, old_col in enumerate(score_cols)}
        pca_projected = pca_projected.rename(columns=rename_dict)[['#FID', 'IID'] + list(rename_dict.values())]
        pca_data_full = pd.concat([pca_df, pca_projected], ignore_index=True)
        
        # 2. 准备样本标签
        print("Preparing sample labels...")
        # 1KG 标签
        panel_1kg = pd.read_csv(panel_1kg_file, sep="\t")
        panel_1kg = panel_1kg[["Family ID", "Individual ID", "Population"]]
        panel_1kg.columns = ["FID", "IID", "Population"]
        panel_1kg["Superpopulation"]="Unknown_1kg"
        panel_1kg.loc[panel_1kg["Population"].isin(["JPT","CHB","CDX","CHS","JPT","KHV","CHD"]),"Superpopulation"]="EAS"
        panel_1kg.loc[panel_1kg["Population"].isin(["BEB","GIH","ITU","PJL","STU"]),"Superpopulation"]="SAS"
        panel_1kg.loc[panel_1kg["Population"].isin(["CLM","MXL","PEL","PUR"]),"Superpopulation"]="AMR"
        panel_1kg.loc[panel_1kg["Population"].isin(["CEU","FIN","GBR","IBS","TSI"]),"Superpopulation"]="EUR"
        panel_1kg.loc[panel_1kg["Population"].isin(["ACB","ASW","ESN","GWD","LWK","MSL","YRI"]),"Superpopulation"]="AFR"
        panel_1kg = panel_1kg[['IID', 'Superpopulation']]
        # 将Superpopulation列重命名为Cohort
        panel_1kg.rename(columns={'Superpopulation': 'Cohort'}, inplace=True)
        # 队列A 标签
        fam_a = pd.read_csv(cohort_a_list, sep=' ', header=None, usecols=[1], names=['IID'])
        fam_a['Cohort'] = 'CAS Cohort'
        # 队列B 标签
        fam_b = pd.read_csv(cohort_b_list, sep=' ', header=None, usecols=[1], names=['IID'])
        fam_b['IID'] = fam_b['IID'].astype(str)  # 确保IID是字符串类型
        fam_b['Cohort'] = 'UK Biobank (EAS)'
        
        # 合并所有标签信息
        labels_df = pd.concat([panel_1kg, fam_a, fam_b], ignore_index=True)
        
        # 3. 将PCA坐标与样本标签合并
        print("Merging PCA data with labels...")
        plot_data = pd.merge(pca_data_full, labels_df, on='IID', how='left')
        
        # 处理可能未匹配上的样本
        plot_data['Cohort'].fillna('Unknown', inplace=True)
        if 'Unknown' in plot_data['Cohort'].unique():
            print(f"Warning: Found {plot_data['Cohort'].value_counts()['Unknown']} samples with no cohort label.")
        if 'Unknown_1kg' in plot_data['Cohort'].unique():
            print(f"Warning: Found {plot_data['Cohort'].value_counts()['Unknown_1kg']} samples with no 1KG label.")

        # 4. 绘图
        print("Generating plot...")
        # 设置因子水平，让图例顺序更美观，并把您的队列点放在最上层
        cohort_order = ['CAS Cohort', 'UK Biobank (EAS)', 'EAS', 'EUR', 'SAS', 'AMR', 'AFR', 'Unknown']
        plot_data['Cohort'] = pd.Categorical(plot_data['Cohort'], categories=cohort_order, ordered=True)
        
        plt.figure(figsize=(12, 9))
        
        # 使用 seaborn 绘图，它会自动处理颜色和图例
        sns.scatterplot(
            data=plot_data,
            x='PC1',
            y='PC2',
            hue='Cohort',      # 根据Cohort列来分配颜色
            style='Cohort',    # 根据Cohort列来分配形状
            s=50,              # 点的大小
            alpha=0.7          # 透明度
        )
        
        plt.title('Population Stratification using PCA', fontsize=18)
        plt.xlabel('Principal Component 1', fontsize=12)
        plt.ylabel('Principal Component 2', fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.legend(title='Population / Cohort', bbox_to_anchor=(1.02, 1), loc='upper left')
        
        # 5. 保存图像
        plt.tight_layout() # 调整布局，防止图例被截断
        plt.savefig(output_png, dpi=300)
        
        print(f"\nSuccess! Plot saved to {output_png}")

    except FileNotFoundError as e:
        print(f"\nError: A required file was not found.")
        print(f"Details: {e}")
        print("Please make sure all input files are in the correct directory.")
        sys.exit(1)


if __name__ == '__main__':
    plot_pca()