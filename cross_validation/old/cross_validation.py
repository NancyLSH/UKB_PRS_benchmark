# =============================================================================
# 最终版 V2：更灵活的交叉验证ID列表生成脚本
#
# 版本特性：
# - 新增 `is_continuous` 参数，可灵活选择分层抽样或随机抽样。
# - 完美适配连续性状（如身高）和分类性状（如疾病状态）。
# =============================================================================

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, KFold # Import KFold
import os
import io

def create_plink_split_files(
    fam_filepath, 
    pheno_filepath, 
    output_prefix: str,
    is_continuous: bool = False, # <<< 新增参数
    n_splits: int = 10, 
    random_state: int = 42
):
    """
    读取.fam和表型文件，根据表型是分类还是连续，选择合适的抽样方法生成ID列表。

    Args:
        ...
        is_continuous (bool, optional): 
            标记当前表型是否为连续性状。
            - False (默认): 使用 StratifiedKFold (分层抽样)。
            - True: 使用 KFold (标准随机抽样)。
        ...
    """
    print(f"\n{'='*25}\n处理表型文件: {os.path.basename(str(pheno_filepath))}\n{'='*25}")
    
    # ... (数据加载和合并部分与之前完全相同)
    master_samples = pd.read_csv(fam_filepath, sep='\s+', header=None, usecols=[0, 1], names=['FID', 'IID'])
    pheno_data = pd.read_csv(pheno_filepath, sep='\s+', header=None, names=['FID', 'IID', 'Pheno'])
    merged_data = pd.merge(master_samples, pheno_data, on=['FID', 'IID'], how='left')
    merged_data.dropna(subset=['Pheno'], inplace=True)
    
    y_pheno = merged_data['Pheno'].values
    sample_ids_for_split = merged_data[['FID', 'IID']]
    dummy_X = np.zeros((len(y_pheno), 1))

    # --- 步骤 4: 根据表型类型选择抽样器 ---
    if is_continuous:
        # 如果是连续性状，使用标准 KFold
        print(f"检测到连续性状，使用标准KFold进行随机抽样...")
        splitter = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)
        # KFold的split方法不需要y作为参数
        all_fold_indices = [test_idx for _, test_idx in splitter.split(dummy_X)]
    else:
        # 如果是分类性状，使用 StratifiedKFold
        print(f"检测到分类性状，使用StratifiedKFold进行分层抽样...")
        # 检查y是否为整数类型，以避免潜在错误
        if pd.api.types.is_float_dtype(y_pheno):
             y_pheno = y_pheno.astype(int)
        splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
        all_fold_indices = [test_idx for _, test_idx in splitter.split(dummy_X, y_pheno)]
    
    # ... (文件生成和保存部分与之前完全相同)
    output_dir = os.path.dirname(output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i in range(n_splits):
        test_fold_num = i
        tune_fold_num = (i + 1) % n_splits
        train_fold_nums = [j for j in range(n_splits) if j != test_fold_num and j != tune_fold_num]

        test_indices = all_fold_indices[test_fold_num]
        tune_indices = all_fold_indices[tune_fold_num]
        train_indices = np.concatenate([all_fold_indices[k] for k in train_fold_nums])

        sample_ids_for_split.iloc[train_indices].to_csv(f"{output_prefix}_{i+1}_train_ids.txt", sep='\t', index=False, header=False)
        sample_ids_for_split.iloc[tune_indices].to_csv(f"{output_prefix}_{i+1}_tune_ids.txt", sep='\t', index=False, header=False)
        sample_ids_for_split.iloc[test_indices].to_csv(f"{output_prefix}_{i+1}_test_ids.txt", sep='\t', index=False, header=False)
    
    print(f"成功为 '{os.path.basename(str(pheno_filepath))}' 生成了所有 {n_splits} 折的ID文件。")

# =============================================================================
# 主程序：演示如何为不同类型的表型调用
# =============================================================================
if __name__ == "__main__":
    
    # --- 模拟文件环境 ---
    fam_content = """
    FAM0 SAMPLE_0 ...
    FAM1 SAMPLE_1 ...
    FAM1 SAMPLE_2 ...
    """
    mock_fam_file = io.StringIO(fam_content)

    # 模拟 AD.pheno (二元分类)
    pheno_AD_content = "FAM0 SAMPLE_0 0\nFAM1 SAMPLE_1 1\nFAM1 SAMPLE_2 1"
    mock_pheno_AD_file = io.StringIO(pheno_AD_content)
    
    # 模拟 Height.pheno (连续)
    pheno_Height_content = "FAM0 SAMPLE_0 176.5\nFAM1 SAMPLE_1 165.2\nFAM1 SAMPLE_2 180.1"
    mock_pheno_Height_file = io.StringIO(pheno_Height_content)

    # --- 为分类性状 (AD) 调用 ---
    create_plink_split_files(
        fam_filepath=mock_fam_file,
        pheno_filepath=mock_pheno_AD_file,
        output_prefix="results/eur_AD_cv",
        is_continuous=False # 明确告知这是分类性状 (默认值)
    )
    mock_fam_file.seek(0) # 重置指针

    # --- 为连续性状 (Height) 调用 ---
    create_plink_split_files(
        fam_filepath=mock_fam_file,
        pheno_filepath=mock_pheno_Height_file,
        output_prefix="results/eur_Height_cv",
        is_continuous=True # <<< 关键：明确告知这是连续性状
    )