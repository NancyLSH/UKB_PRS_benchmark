import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_heritability_comparison(df):
    """
    绘制三个人群遗传力的比较图（森林图样式）。

    参数:
    df (pd.DataFrame): 包含遗传力数据的DataFrame。
                       需要包含以下列: 'trait', 
                                     'h2_eur', 'se_eur',
                                     'h2_eas', 'se_eas',
                                     'h2_afr', 'se_afr'
    """
    # --- 数据准备和计算 ---
    # 计算95%置信区间
    df['eur_ci_low'] = df['h2_eur'] - 1.96 * df['se_eur']
    df['eur_ci_high'] = df['h2_eur'] + 1.96 * df['se_eur']
    df['eas_ci_low'] = df['h2_eas'] - 1.96 * df['se_eas']
    df['eas_ci_high'] = df['h2_eas'] + 1.96 * df['se_eas']
    df['afr_ci_low'] = df['h2_afr'] - 1.96 * df['se_afr']
    df['afr_ci_high'] = df['h2_afr'] + 1.96 * df['se_afr']

    # 关键步骤：根据三组遗传力的平均值对DataFrame进行排序
    df['h2_mean'] = df[['h2_eur', 'h2_eas', 'h2_afr']].mean(axis=1)
    df_sorted = df.sort_values('h2_mean', ascending=True).reset_index(drop=True)

    # --- 开始绘图 ---
    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(12, 10))

    # 为每个点定义一个小的垂直偏移量，以避免重叠
    offsets = {'eur': -0.15, 'eas': 0, 'afr': 0.15}
    colors = {'eur': 'royalblue', 'eas': 'darkorange', 'afr': 'seagreen'}
    labels = {'eur': 'European', 'eas': 'East Asian', 'afr': 'African'}

    # 绘制误差棒和散点
    for pop in ['eur', 'eas', 'afr']:
        y_pos = df_sorted.index + offsets[pop]
        ax.errorbar(df_sorted[f'h2_{pop}'], y_pos,
                    xerr=(df_sorted[f'h2_{pop}'] - df_sorted[f'{pop}_ci_low'], 
                          df_sorted[f'{pop}_ci_high'] - df_sorted[f'h2_{pop}']),
                    fmt='none', color=colors[pop], capsize=5, alpha=0.8,
                    linewidth=1.5)
        
        ax.scatter(df_sorted[f'h2_{pop}'], y_pos, 
                   color=colors[pop], s=150, zorder=2, 
                   label=labels[pop], edgecolors='white', linewidth=1)

    # --- 美化图表 ---
    ax.set_yticks(df_sorted.index)
    ax.set_yticklabels(df_sorted['trait'], fontsize=16)
    ax.set_xlabel('SNP-based Heritability ($h^2_{SNP}$)', fontsize=16)
    ax.set_ylabel('')
    ax.set_title('Heritability Comparison Across Three Ancestries', fontsize=18, pad=20)
    
    # 创建图例
    # 为了避免图例中出现重复项，我们只为第一个循环中的标签创建图例
    handles, labels_map = ax.get_legend_handles_labels()
    # 根据标签顺序创建一个有序的图例
    by_label = dict(sorted(zip(labels_map, handles), key=lambda x: list(labels.values()).index(x[0])))
    ax.legend(by_label.values(), by_label.keys(), title='Population', fontsize=12, title_fontsize=13, loc='lower right')
    
    ax.grid(True, which='major', axis='x', linestyle='--', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)

    plt.tight_layout()
    # 您可以取消下面一行的注释，将图片保存到本地
    # plt.savefig("heritability_comparison_3_populations.png")
    plt.show()

# --- 创建一个示例DataFrame来演示 ---
# 在实际使用中，您应该加载自己的数据
data = {
    'trait': ['Height', 'BMI', 'Systolic BP', 'Diastolic BP', 'Fasting Glucose', 'HDL Cholesterol', 'LDL Cholesterol', 'Triglycerides'],
    'h2_eur': [0.5, 0.3, 0.25, 0.22, 0.18, 0.35, 0.3, 0.15],
    'se_eur': [0.03, 0.02, 0.025, 0.02, 0.015, 0.03, 0.025, 0.01],
    'h2_eas': [0.48, 0.32, 0.23, 0.21, 0.2, 0.33, 0.28, 0.16],
    'se_eas': [0.035, 0.022, 0.028, 0.022, 0.018, 0.032, 0.028, 0.012],
    'h2_afr': [0.52, 0.28, 0.28, 0.25, 0.16, 0.38, 0.32, 0.18],
    'se_afr': [0.04, 0.025, 0.03, 0.025, 0.02, 0.035, 0.03, 0.015]
}
df_sample = pd.DataFrame(data)

# 调用函数进行绘图
plot_heritability_comparison(df_sample)