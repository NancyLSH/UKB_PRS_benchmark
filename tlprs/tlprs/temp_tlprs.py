import pandas as pd
import numpy as np
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score, average_precision_score, mean_squared_error
from scipy.stats import pearsonr
import os
import warnings

# Ignore warnings that may arise from certain fits in statsmodels
from statsmodels.tools.sm_exceptions import ConvergenceWarning
warnings.simplefilter('ignore', ConvergenceWarning)

# --- 1. Metric Calculation Functions ---

def calculate_continuous_metrics(df, base_covars, full_covars):
    """Calculates all performance metrics for a continuous trait for a given dataset (df)."""
    # Incremental R²
    model_base = sm.OLS(df["trait"], sm.add_constant(df[base_covars])).fit()
    model_full = sm.OLS(df["trait"], sm.add_constant(df[full_covars])).fit()
    r2_incremental = model_full.rsquared - model_base.rsquared

    # Pearson correlation coefficient (SCORE vs. phenotype residuals)
    pheno_residuals = model_base.resid
    corr, _ = pearsonr(df["SCORE"], pheno_residuals)

    # RMSE
    prediction_full = model_full.predict(sm.add_constant(df[full_covars]))
    rmse = np.sqrt(mean_squared_error(df["trait"], prediction_full))

    # Quantile means
    df['quantile'] = pd.qcut(df['SCORE'], 5, labels=False, duplicates='drop')
    quantile_means = df.groupby('quantile')['trait'].mean()
    
    return {
        "r2_incremental": r2_incremental,
        "r2_full": model_full.rsquared,
        "rmse": rmse,
        "pearson_r": corr,
        "top_quintile_mean": quantile_means.iloc[-1] if not quantile_means.empty else np.nan,
        "bottom_quintile_mean": quantile_means.iloc[0] if not quantile_means.empty else np.nan
    }

def calculate_binary_metrics(df, base_covars, full_covars):
    """Calculates all performance metrics for a binary trait for a given dataset (df)."""
    # AUC and PR-AUC
    logit_model = sm.Logit(df["trait"], sm.add_constant(df[full_covars])).fit(disp=0)
    pred_prob = logit_model.predict(sm.add_constant(df[full_covars]))
    auc = roc_auc_score(df["trait"], pred_prob)
    pr_auc = average_precision_score(df["trait"], pred_prob)

    # OR per 1-SD
    df["prs_scaled"] = (df["SCORE"] - df["SCORE"].mean()) / df["SCORE"].std()
    logit_model_scaled = sm.Logit(df["trait"], sm.add_constant(df[base_covars + ["prs_scaled"]])).fit(disp=0)
    or_per_sd = np.exp(logit_model_scaled.params["prs_scaled"])

    # Quantile OR
    df['prs_quintile'] = pd.qcut(df['SCORE'], 5, labels=False, duplicates='drop')
    reference_quintile = 2 # Middle quintile
    or_quintiles = {}
    for q in range(5):
        if q == reference_quintile:
            or_quintiles[f'OR_Quintile_{q+1}'] = 1.0
            continue
        
        # Check if both current and reference quintiles exist in the data
        if not df['prs_quintile'].isin([q, reference_quintile]).all():
             or_quintiles[f'OR_Quintile_{q+1}'] = np.nan
             continue
             
        temp_df = df[df['prs_quintile'].isin([q, reference_quintile])].copy()
        
        # Check for sufficient data in both groups for stable model fitting
        if temp_df['trait'].nunique() < 2 or temp_df['prs_quintile'].nunique() < 2:
            or_quintiles[f'OR_Quintile_{q+1}'] = np.nan
            continue
            
        temp_df['is_current_quintile'] = (temp_df['prs_quintile'] == q).astype(int)
        X_quintile = sm.add_constant(temp_df[['is_current_quintile'] + base_covars])
        try:
            model_q = sm.Logit(temp_df["trait"], X_quintile).fit(disp=0)
            or_quintiles[f'OR_Quintile_{q+1}'] = np.exp(model_q.params['is_current_quintile'])
        except Exception:
            or_quintiles[f'OR_Quintile_{q+1}'] = np.nan
            
    results = {
        "auc": auc,
        "pr_auc": pr_auc,
        "or_per_sd": or_per_sd,
    }
    results.update(or_quintiles) # Merge quantile ORs into the results dictionary
    return results

# --- 2. Main Execution Flow ---

def main():
    # --- Parameter Settings ---
    # !! Note: Please modify the variables below according to your actual paths !!
    cleaned_prs_path = "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/tlprs/res/test_prs_cleaned"
    covar_path = "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/pheno/covar/covars_chinese_final.tsv"
    output_dir = "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/tlprs/res/full_res"
    pheno_dir = "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/pheno/trait/Chinese"
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # --- Initialization ---
    final_results_continuous = []
    final_results_binary = []
    covar_cols = ["FID", "IID", "age", "sex"] + [f"PC{i}" for i in range(1, 11)]
    base_covars = ["age", "sex"] + [f"PC{i}" for i in range(1, 11)]
    full_covars = base_covars + ["SCORE"]

    # --- Data Loading and Processing ---
    covars = pd.read_csv(covar_path, sep='\t', usecols=covar_cols)
    
    # Iterate through existing PRS files
    for prs_file in os.listdir(cleaned_prs_path):
        # Check if it is a file and has the expected format, not a directory
        if not os.path.isfile(os.path.join(cleaned_prs_path, prs_file)):
            continue

        trait_id_from_prs = prs_file.split("_")[1]
        trait_pheno_name = trait_id_from_prs + '_int' # Construct matching phenotype filename
        
        # Find the corresponding phenotype file
        pheno_file_found = None
        for pheno_file in os.listdir(pheno_dir):
            if pheno_file.startswith(trait_pheno_name) and pheno_file.endswith(".txt"):
                pheno_file_found = pheno_file
                break
        
        if not pheno_file_found:
            print(f"Warning: No phenotype file found for PRS trait {trait_id_from_prs}. Skipping.")
            continue

        print(f"\nProcessing Trait: {trait_id_from_prs}")
        pheno_path = os.path.join(pheno_dir, pheno_file_found)
        current_prs_path = os.path.join(cleaned_prs_path, prs_file)

        pheno = pd.read_csv(pheno_path, sep='\t', header=None)
        pheno.columns = ["FID", "IID", "trait"]
        prs = pd.read_csv(current_prs_path, sep='\t')
        
        # Merge data
        # Ensure FID and IID are in the same format before merging
        for df in [pheno, prs, covars]:
            df["FID"] = df["FID"].astype(str)
            df["IID"] = df["IID"].astype(str)
            
        merged_data = pd.merge(pheno, prs, on=["FID", "IID"], how="inner")
        merged_data = pd.merge(merged_data, covars, on=["FID", "IID"], how="inner")

        # Defensive data cleaning and type conversion
        numeric_cols = ["trait", "SCORE", "age", "sex"] + [f"PC{i}" for i in range(1, 11)]
        for col in numeric_cols:
            if col in merged_data.columns:
                merged_data[col] = pd.to_numeric(merged_data[col], errors='coerce')

        original_rows = len(merged_data)
        merged_data.dropna(subset=numeric_cols, inplace=True)
        new_rows = len(merged_data)

        if original_rows > new_rows:
            print(f"--> Warning: Dropped {original_rows - new_rows} rows due to non-numeric data or NaNs.")

        if new_rows == 0:
            print(f"--> Error: No valid samples remained for trait {trait_id_from_prs} after cleaning. Skipping.")
            continue
            
        print(f"Data merged and cleaned for trait {trait_id_from_prs}. Total samples: {len(merged_data)}")

        # --- Perform Analysis ---
        if trait_pheno_name in ["p20116_int", "p20117_int"]:
            # Binary trait analysis
            # Ensure binary trait is 0/1 coded
            unique_vals = sorted(merged_data["trait"].unique())
            if not set(unique_vals).issubset({0, 1}):
                if len(unique_vals) == 2:
                    print(f"Converting binary trait from {unique_vals} to 0/1.")
                    merged_data["trait"] = (merged_data["trait"] == unique_vals[1]).astype(int)
                else:
                    print(f"Error: Binary trait column for {trait_id_from_prs} contains unexpected values: {unique_vals}. Skipping.")
                    continue
            
            # Direct calculation of metrics
            analysis_report = calculate_binary_metrics(merged_data, base_covars, full_covars)
            analysis_report['trait'] = trait_id_from_prs
            final_results_binary.append(analysis_report)
        else:
            # Continuous trait analysis
            # Direct calculation of metrics
            analysis_report = calculate_continuous_metrics(merged_data, base_covars, full_covars)
            analysis_report['trait'] = trait_id_from_prs
            final_results_continuous.append(analysis_report)

    # --- 3. Save Final Results ---
    if final_results_continuous:
        continuous_df = pd.DataFrame(final_results_continuous)
        # Reorder columns to have 'trait' first
        cols = ['trait'] + [col for col in continuous_df.columns if col != 'trait']
        continuous_df = continuous_df[cols]
        continuous_df.to_csv(os.path.join(output_dir, "prs_continuous_metrics.csv"), index=False)
        print("\nContinuous trait results saved to prs_continuous_metrics.csv")
        print(continuous_df)

    if final_results_binary:
        binary_df = pd.DataFrame(final_results_binary)
        # Reorder columns to have 'trait' first
        cols = ['trait'] + [col for col in binary_df.columns if col != 'trait']
        binary_df = binary_df[cols]
        binary_df.to_csv(os.path.join(output_dir, "prs_binary_metrics.csv"), index=False)
        print("\nBinary trait results saved to prs_binary_metrics.csv")
        print(binary_df)

if __name__ == '__main__':
    main()