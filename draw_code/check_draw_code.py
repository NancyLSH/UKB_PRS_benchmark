import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.lines import Line2D

# --- 1. 准备数据 ---
# 在您的实际应用中，您会从一个文件中加载计算好的结果
# 例如: df = pd.read_csv('my_prs_results.csv')
#
# 为了演示，我们这里创建一个模拟的DataFrame。
# 您的真实数据应该整理成下面这种“长格式”。

print("--- 用于绘图的模拟数据 ---")
print(merge_res.head())


# --- 2. 开始绘图 ---

# 设置绘图风格和字体，确保中文能够正确显示
sns.set_style("whitegrid")
# plt.rcParams['font.sans-serif'] = ['SimHei'] # 'SimHei' 是常用的黑体
plt.rcParams['axes.unicode_minus'] = False # 解决负号显示问题

# 创建一个图形和坐标轴，设置合适的尺寸
fig, ax = plt.subplots(figsize=(14, 8))

# 使用seaborn绘制分组柱状图
# x: x轴类别 (表型)
# y: y轴数值 (增量R²)
# hue: 分组的类别 (方法)
# data: 输入的DataFrame
barplot = sns.barplot(
    x='pheno',
    y='r2_full_median',
    hue='method',
    data=merge_res,
    palette='muted',  # 'viridis', 'muted', 'colorblind' 都是不错的色板
    ax=ax,
    # errorbar=None  # 不直接传递误差棒
)

# --- 3. 手动添加自定义的误差棒 ---
# 这种方法能精确控制非对称的误差棒，是更严谨的做法
# 获取条形的位置
x_coords = [p.get_x() + 0.5 * p.get_width() for p in ax.patches]
y_coords = [p.get_height() for p in ax.patches]

# 重新排序数据以匹配条形图的顺序
plot_order = merge_res.groupby(['pheno', 'method']).groups.keys()

errors = {}
for trait, method in plot_order:
    row = merge_res[(merge_res['pheno'] == trait) & (merge_res['method'] == method)]
    if row.shape[0] == 1:
    # if not row.empty:
        # print(f"Processing trait: {trait}, method: {method}")
        median = row['r2_full_median'].iloc[0]
        ci_lower = row['r2_full_CI_lower'].iloc[0]
        ci_upper = row['r2_full_CI_upper'].iloc[0]
        # print(f"Median: {median}, CI Lower: {ci_lower}, CI Upper: {ci_upper}")
        errors[(trait, method)] = [median - ci_lower, ci_upper - median]

# 从ax.get_xticklabels()获取正确的trait顺序
trait_order = [tick.get_text() for tick in ax.get_xticklabels()]
software_order = [text.get_text() for text in ax.legend_.get_texts()]
print(f"Trait order: {trait_order}")
print(f"Software order: {software_order}")

err_values = []
for trait in trait_order:
    for software in software_order:
        err = errors.get((trait, software), [0, 0])
        err_values.append(err)
err_values = np.array(err_values).T

ax.errorbar(
    x=x_coords,
    y=y_coords,
    yerr=err_values,
    fmt='none',  # 不显示数据点
    capsize=3,  # 误差线顶端横线长度
    ecolor='black',
    elinewidth=1,  # 误差线宽度
    linestyle='none',  # 不显示误差线的线型
    alpha=0.7,  # 误差线透明度
)


# --- 4. 美化与调整 ---

# 设置x轴的标签旋转角度，避免重叠
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

# 对于不同的category，对x轴标签进行分组显示
# 这里假设您的DataFrame中有一个'category'列来区分不同
# 根据分组标签改变Y轴标签的颜色
for label in ax.get_xticklabels():
    trait_name = label.get_text()
    # 根据类别设置标签颜色
    category = merge_res.loc[merge_res['pheno'] == trait_name, 'category'].values[0]
    label.set_color('#073B4C' if category == 'Blood biochemistry' else '#118AB2' if category == 'Body size measures' else '#E07A5F' if category == "Blood count" else "#EF476F" if category == "Blood pressure" else "#06D6A0")

legend_elements = [
    Line2D([0], [0], color='#073B4C', lw=4, label='Blood biochemistry'),
    Line2D([0], [0], color='#118AB2', lw=4, label='Body size measures'),
    Line2D([0], [0], color='#E07A5F', lw=4, label='Blood count'),
    Line2D([0], [0], color='#EF476F', lw=4, label='Blood pressure'),
]
ax.add_artist(ax.legend(handles=legend_elements, loc='center left', fontsize=10, title='Trait Categories', title_fontsize='11', bbox_to_anchor=(0, 0.70)))
# 设置Y轴的格式，转为百分比显示更直观
# ax.set_yticklabels([f'{y*100:.1f}%' for y in ax.get_yticks()])


ax.set_title(r'Full Model $R^2$ Comparison Across Phenotypes', fontsize=18, pad=20)
ax.set_xlabel('Phenotype', fontsize=14)
ax.set_ylabel(r'$R^2$', fontsize=14)

# # 调整图例
# # 再加一张图例，显示不同的PRS算法
# legend_elements = [
#     Line2D([0], [0], color='blue', lw=4, label='PRS-CSx'),
#     Line2D([0], [0], color='orange', lw=4, label='XPASS'),
# ]
# ax.legend(handles=legend_elements, title='PRS Algorithm', fontsize=11, title_fontsize='13')
ax.legend(title='PRS Algorithm', title_fontsize='11', fontsize='10', loc='upper left')

# 调整坐标轴刻度标签的字体大小
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=11)

# 确保布局紧凑
plt.tight_layout()

# 保存图形（可选）
# plt.savefig('incremental_r2_comparison.png', dpi=300)

# 显示图形
plt.show()