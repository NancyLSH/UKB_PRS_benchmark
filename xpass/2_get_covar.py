import pandas as pd

def get_covar(file_path, output_path):
    """
    Read a file and extract the covariance matrix.
    
    Args:
        file_path (str): Path to the file containing covariance data.
        
    Returns:
        pd.DataFrame: DataFrame containing the covariance matrix.
    """
    covar = pd.read_csv(file_path, sep='\t', usecols=[f'PC{i}' for i in range(1, 10)], header=0)
    covar.to_csv(output_path, sep=' ', index=False, header=False)

if __name__ == "__main__":
    covar_path = "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/CAS/geno/CAS_final/CAS_pca_2.eigenvec"
    output_path = "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/CAS/pheno/covar_for_xpass.txt"
    get_covar(covar_path, output_path)
    print(f"Covariance matrix saved to {output_path}")