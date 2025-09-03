import pandas as pd
import numpy as np
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score, average_precision_score, mean_squared_error
from sklearn.calibration import calibration_curve
from scipy.stats import pearsonr
import os
import warnings

# 忽略statsmodels在某些拟合中可能产生的警告
from statsmodels.tools.sm_exceptions import ConvergenceWarning
warnings.simplefilter('ignore', ConvergenceWarning)

# --- 1. 指标计算函数 (从原脚本逻辑封装而来) ---

def calculate_continuous_metrics(df, base_covars, full_covars):
    """为给定的数据集(df)计算连续性状的所有性能指标。"""
    # 增量R²
    model_base = sm.OLS(df["trait"], sm.add_constant(df[base_covars])).fit()
    model_full = sm.OLS(df["trait"], sm.add_constant(df[full_covars])).fit()
    r2_incremental = model_full.rsquared - model_base.rsquared

    # 皮尔逊相关系数 (SCORE vs. 表型残差)
    pheno_residuals = model_base.resid
    corr, _ = pearsonr(df["SCORE"], pheno_residuals)

    # RMSE
    prediction_full = model_full.predict(sm.add_constant(df[full_covars]))
    rmse = np.sqrt(mean_squared_error(df["trait"], prediction_full))

    # 分位数均值
    df['quantile'] = pd.qcut(df['SCORE'], 5, labels=False, duplicates='drop')
    quantile_means = df.groupby('quantile')['trait'].mean()
    
    return {
        "r2_incremental": r2_incremental,
        "r2_full": model_full.rsquared,
        "rmse": rmse,
        "pearson_r": corr,
        "top_quintile_mean": quantile_means.iloc[-1],
        "bottom_quintile_mean": quantile_means.iloc[0]
    }

def calculate_binary_metrics(df, base_covars, full_covars):
    """为给定的数据集(df)计算二元性状的所有性能指标。"""
    # AUC 和 PR-AUC
    logit_model = sm.Logit(df["trait"], sm.add_constant(df[full_covars])).fit(disp=0)
    pred_prob = logit_model.predict(sm.add_constant(df[full_covars]))
    auc = roc_auc_score(df["trait"], pred_prob)
    pr_auc = average_precision_score(df["trait"], pred_prob)

    # 每1-SD的OR
    df["prs_scaled"] = (df["SCORE"] - df["SCORE"].mean()) / df["SCORE"].std()
    logit_model_scaled = sm.Logit(df["trait"], sm.add_constant(df[base_covars + ["prs_scaled"]])).fit(disp=0)
    or_per_sd = np.exp(logit_model_scaled.params["prs_scaled"])

    # 分位数OR
    df['prs_quintile'] = pd.qcut(df['SCORE'], 5, labels=False, duplicates='drop')
    reference_quintile = 2
    or_quintiles = {}
    for q in range(5):
        if q == reference_quintile:
            or_quintiles[f'OR_Quintile_{q+1}'] = 1.0
            continue
        temp_df = df[df['prs_quintile'].isin([q, reference_quintile])].copy()
        temp_df['is_current_quintile'] = (temp_df['prs_quintile'] == q).astype(int)
        X_quintile = sm.add_constant(temp_df[['is_current_quintile'] + base_covars])
        model_q = sm.Logit(temp_df["trait"], X_quintile).fit(disp=0)
        or_quintiles[f'OR_Quintile_{q+1}'] = np.exp(model_q.params['is_current_quintile'])
        
    results = {
        "auc": auc,
        "pr_auc": pr_auc,
        "or_per_sd": or_per_sd,
    }
    results.update(or_quintiles) # 将分位数OR合并到结果字典
    return results

# --- 2. Bootstrap核心分析函数 ---

def bootstrap_analysis(df, n_bootstrap, analysis_func, base_covars, full_covars):
    """
    对给定的数据集执行Bootstrap分析。
    
    参数:
    - df: 完整的数据集 (DataFrame)。
    - n_bootstrap: Bootstrap重复次数。
    - analysis_func: 用于在每个样本上计算指标的函数 (例如, calculate_continuous_metrics)。
    - base_covars, full_covars: 协变量列表。

    返回:
    - 一个包含点估计和95%置信区间的字典。
    """
    bootstrap_results = []
    for i in range(n_bootstrap):
        # 创建自助样本 (有放回抽样)
        sample_df = df.sample(n=len(df), replace=True)
        
        try:
            # 在自助样本上计算指标
            metrics = analysis_func(sample_df, base_covars, full_covars)
            bootstrap_results.append(metrics)
        except Exception as e:
            # 在某些自助样本中，由于数据特殊性(如二元性状只有一类)，模型可能无法拟合
            # 此时可以跳过这次失败的抽样
            print(f"Bootstrap iteration {i+1} failed with error: {e}. Skipping.")
            continue
    
    # 将结果列表转换为DataFrame
    results_df = pd.DataFrame(bootstrap_results)
    
    # 计算点估计 (中位数) 和 95% CI
    final_report = {}
    for col in results_df.columns:
        point_estimate = results_df[col].median()
        ci_lower = results_df[col].quantile(0.025)
        ci_upper = results_df[col].quantile(0.975)
        final_report[f'{col}_median'] = point_estimate
        final_report[f'{col}_CI_lower'] = ci_lower
        final_report[f'{col}_CI_upper'] = ci_upper
        
    return final_report

# --- 3. 主执行流程 ---

def main():
    # --- 参数设置 ---
    # !! 注意: 请根据您的实际路径修改下面的变量 !!
    cleaned_prs_path = "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/tlprs/res/test_prs_cleaned"
    covar_path = "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/pheno/covar/covars_chinese_final.tsv"
    output_dir = "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/tlprs/res/full_res"
    pheno_dir = "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/pheno/trait/Chinese"
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    n_bootstrap = 1000 # 推荐值, 可设为100进行快速测试
    
    # --- 初始化 ---
    final_results_continuous = []
    final_results_binary = []
    covar_cols = ["FID", "IID", "age", "sex"] + [f"PC{i}" for i in range(1, 11)]
    base_covars = ["age", "sex"] + [f"PC{i}" for i in range(1, 11)]
    full_covars = base_covars + ["SCORE"]

    # --- 数据加载与处理 ---
    covars = pd.read_csv(covar_path, sep='\t', usecols=covar_cols)
    
    for prs_file in os.listdir(cleaned_prs_path):
        trait_id_from_prs = prs_file.split("_")[1]
        trait = trait_id_from_prs + '_int' # 构造匹配的表型文件名
        
        # 寻找对应的表型文件
        pheno_file_found = None
        for pheno_file in os.listdir(pheno_dir):
            if pheno_file.startswith(trait) and pheno_file.endswith(".txt"):
                pheno_file_found = pheno_file
                break
        
        if not pheno_file_found:
            print(f"Warning: No phenotype file found for PRS trait {trait_id_from_prs}. Skipping.")
            continue

        print(f"\nProcessing Trait: {trait_id_from_prs}")
        pheno = pd.read_csv(os.path.join(pheno_dir, pheno_file_found), sep='\t', header=None)
        pheno.columns = ["FID", "IID", "trait"]
        prs = pd.read_csv(os.path.join(cleaned_prs_path, prs_file), sep='\t')
        
        # 合并数据
        # make sure FID and IID are in the same format
        pheno["FID"] = pheno["FID"].astype(str)
        pheno["IID"] = pheno["IID"].astype(str)
        prs["FID"] = prs["FID"].astype(str)
        prs["IID"] = prs["IID"].astype(str)
        covars["FID"] = covars["FID"].astype(str)
        covars["IID"] = covars["IID"].astype(str)
        merged_data = pd.merge(pheno, prs, on=["FID", "IID"], how="inner")
        merged_data = pd.merge(merged_data, covars, on=["FID", "IID"], how="inner")

        # --- NEW: Add a defensive data cleaning and type conversion step ---
        # 1. Define all columns that MUST be numeric for the analysis
        numeric_cols = ["trait", "SCORE", "age", "sex"] + [f"PC{i}" for i in range(1, 11)]

        # 2. Loop through the columns and force them to be numeric.
        # The 'coerce' option is key: it will turn any problematic non-numeric
        # value (e.g., a string 'NA') into a proper NaN.
        for col in numeric_cols:
            if col in merged_data.columns:
                merged_data[col] = pd.to_numeric(merged_data[col], errors='coerce')

        # 3. Now, drop any rows that contain NaN values (either original or newly created).
        original_rows = len(merged_data)
        merged_data.dropna(inplace=True)
        new_rows = len(merged_data)

        # 4. Optional: Print a warning if rows were dropped, so you are aware of data quality issues.
        if original_rows > new_rows:
            print(f"--> Warning: Dropped {original_rows - new_rows} rows due to non-numeric data or NaNs.")

        print(f"Data merged and cleaned for trait {trait_id_from_prs}. Total samples: {len(merged_data)}")

        # --- 执行分析 ---
        if trait == "p20116_int" or trait == "p20117_int":
            # 二元性状分析
            # 确保二元性状是0/1编码
            unique_vals = sorted(merged_data["trait"].unique())
            if set(unique_vals).issubset({0, 1}):
                pass # 已经是0/1
            elif len(unique_vals) == 2:
                print(f"Converting binary trait from {unique_vals} to 0/1.")
                merged_data["trait"] = (merged_data["trait"] == unique_vals[1]).astype(int)
            else:
                print(f"Error: Binary trait column for {trait_id_from_prs} contains unexpected values: {unique_vals}. Skipping.")
                continue

            analysis_report = bootstrap_analysis(merged_data, n_bootstrap, calculate_binary_metrics, base_covars, full_covars)
            analysis_report['trait'] = trait_id_from_prs
            final_results_binary.append(analysis_report)
        else:
            # 连续性状分析
            analysis_report = bootstrap_analysis(merged_data, n_bootstrap, calculate_continuous_metrics, base_covars, full_covars)
            analysis_report['trait'] = trait_id_from_prs
            final_results_continuous.append(analysis_report)

    # --- 4. 保存最终结果 ---
    if final_results_continuous:
        continuous_df = pd.DataFrame(final_results_continuous)
        continuous_df.to_csv(os.path.join(output_dir, "prs_continuous_metrics_with_ci.csv"), index=False)
        print("\nContinuous trait results saved to prs_continuous_metrics_with_ci.csv")
        print(continuous_df)

    if final_results_binary:
        binary_df = pd.DataFrame(final_results_binary)
        binary_df.to_csv(os.path.join(output_dir, "prs_binary_metrics_with_ci.csv"), index=False)
        print("\nBinary trait results saved to prs_binary_metrics_with_ci.csv")
        print(binary_df)

if __name__ == '__main__':
    main()