# (load_ldsc_rg_data 和 load_popcorn_pge_data 函数保持不变)

def create_comprehensive_heatmap(df_pop1, df_pop2, pge_data):
    """
    创建并绘制最终的合并热图。
    上三角: 人群1内部相关性 (df_pop1)
    下三角: 人群2内部相关性 (df_pop2)
    对角线: 跨人群相关性 (pge_data)
    NA值: 灰色斜线
    """
    # 1. 获取所有性状的唯一且排序的列表
    all_traits = sorted(list(
        set(df_pop1['p1']) | set(df_pop1['p2']) |
        set(df_pop2['p1']) | set(df_pop2['p2']) |
        set(pge_data.keys())
    ))
    
    # 2. 创建一个主数据矩阵
    matrix = pd.DataFrame(index=all_traits, columns=all_traits, dtype=np.float64)
    
    # 填充上三角 (人群1)
    for _, row in df_pop1.iterrows():
        matrix.loc[row['p1'], row['p2']] = row['rg']
    
    # 填充下三角 (人群2)
    for _, row in df_pop2.iterrows():
        matrix.loc[row['p2'], row['p1']] = row['rg']
        
    # 填充对角线 (PGE)
    for trait, pge_value in pge_data.items():
        if trait in all_traits:
            matrix.loc[trait, trait] = pge_value

    # 3. 创建掩码
    # 我们创建一个完整的矩阵，然后用颜色区分不同部分
    # 这里我们不再需要复杂的掩码，因为我们将根据值来决定颜色
    
    # --- 开始绘图 ---
    fig, ax = plt.subplots(figsize=(24, 20))

    # 绘制基础热图，NA值将是透明的
    sns.heatmap(matrix, cmap="vlag", annot=True, fmt=".2f", 
                ax=ax, cbar=True, annot_kws={"size": 9},
                linewidths=0.5, linecolor='white')

    # --- 覆盖NA单元格为灰色斜线 ---
    from matplotlib.patches import Rectangle
    for (i, j), val in np.ndenumerate(matrix):
        if pd.isna(val):
            # i是行索引, j是列索引
            # Rectangle的坐标是(x, y)，对应(j, i)
            ax.add_patch(Rectangle((j, i), 1, 1, fill=True, facecolor='lightgrey', hatch='//', edgecolor='white', lw=0.5))

    # --- 美化图表 ---
    ax.set_title('Comprehensive Genetic Correlation Map', fontsize=22, fontweight='bold', pad=20)
    ax.set_xlabel('Traits', fontsize=16)
    ax.set_ylabel('Traits', fontsize=16)
    plt.xticks(rotation=90, fontsize=12)
    plt.yticks(rotation=0, fontsize=12)
    
    # 反转Y轴，使原点在左下角
    ax.invert_yaxis()

    # 手动创建图例
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], color='#5A8EAF', lw=4, label='Population 1 (e.g., CAS) - rg'),
                       Line2D([0], [0], color='#E67F83', lw=4, label='Population 2 (e.g., EUR) - rg'),
                       Patch(facecolor='lightgrey', hatch='//', label='NA Value')]
    # 这里我们简化图例，因为对角线的值还是在同一个颜色条上
    ax.legend(handles=legend_elements, loc='upper right', fontsize=14)

    fig.tight_layout(rect=[0, 0, 0.9, 1])

    # 保存图表
    output_image_path = "comprehensive_genetic_heatmap.png"
    plt.savefig(output_image_path, dpi=300, bbox_inches='tight')
    print(u"\n综合热图已成功保存到: {}".format(output_image_path))
    
    plt.show()

# (The rest of the script remains the same, but you might need to adjust the main block
# if you used the previous version's separate matrices logic)
if __name__ == '__main__':
    # --- 请在这里配置您的文件名 ---
    pop1_rg_file = 'genetic_correlation_summary_cas.csv'
    pop2_rg_file = 'genetic_correlation_summary_eur.csv' # <--- 确保您有这个文件
    popcorn_pge_file = 'popcorn_results.tsv'

    data_pop1 = load_ldsc_rg_data(pop1_rg_file, 'Population 1')
    data_pop2 = load_ldsc_rg_data(pop2_rg_file, 'Population 2')
    pge_values = load_popcorn_pge_data(popcorn_pge_file)

    if data_pop1 is not None and data_pop2 is not None and pge_values is not None:
        print(u"\n所有数据加载成功，正在生成热图...")
        create_comprehensive_heatmap(data_pop1, data_pop2, pge_values)
    else:
        print(u"\n因文件缺失或读取失败，无法生成图表。请检查文件名和文件内容。")