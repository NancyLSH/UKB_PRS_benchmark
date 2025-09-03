package <- "/data1/jiapl_group/lishuhua/software/PRS/PROSPER/PROSPER-main/"
base_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/prosper/"
eur_base_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/Cross_Validation/UKB_EUR/"

# trait_list <- c("waist", "height", "pulse", "dbp")
# trait_list <- c("sbp", "smoke", "drink", "bmi")
# trait_list <- c("wbc", "rbc", "hb", "plt")
# trait_list <- c("lymph", "mono", "neut", "eos")
# trait_list <- c("alt", "ast", "bun", "cholesterol")
# trait_list <- c("creatinine", "ggt", "glucose", "hdl")
trait_list <- c("ldl", "triglycerides", "ua")

for (trait in trait_list) {
  print(paste("Processing trait:", trait))
  for (i in 1:10){
    print(paste("  Processing fold:", i))
    output_dir <- paste0(base_dir, "res/", trait, "/group_", i, "/lassosum2/")
    if (!dir.exists(output_dir)){
      dir.create(output_dir, recursive = TRUE)
    }
    sst_eas <- paste0(base_dir, "train/EAS/", trait, "/group_", i, ".txt")
    if (!file.exists(sst_eas)){
      print(paste("  EAS summary statistic file does not exist for group", i, "of trait", trait))
      next
    }
    for (chrom in 1:22){
        print(paste("    Processing chromosome:", chrom))
        # summary statistic
        sst_eur <- paste0(base_dir, "train/EUR/", trait, "/group_", i, "_chr", chrom, ".txt")
        if (!file.exists(sst_eur)){
          print(paste("    EUR summary statistic file does not exist for group", i, "of trait", trait, "chromosome", chrom))
          next
        }
        # tuning and testing data
        tuning_eas <- paste0(base_dir, "tuning/EAS/", trait, "/group_", i)
        tuning_eur <- paste0(eur_base_dir, "tune/fold_", i, "/chr", chrom)
        tuning_eas_pheno <- paste0(base_dir, "tuning/EAS/", trait, "/group_", i, ".fam")
        tuning_eur_pheno <- paste0(base_dir, "tuning/EUR/", trait, "/group_", i, "_chr", chrom, ".fam")
        if (!file.exists(tuning_eas_pheno)){
          print(paste("    EAS tuning genotype file does not exist for group", i, "of trait", trait, "chromosome", chrom))
          next
        }
        if (!file.exists(tuning_eur_pheno)){
          print(paste("    EUR tuning genotype file does not exist for group", i, "of trait", trait, "chromosome", chrom))
          next
        }
        file_tuning_input <- paste0(tuning_eur, ",", tuning_eas)
        file_tuning_pheno <- paste0(tuning_eur_pheno, ",", tuning_eas_pheno)
        test_eas <- paste0(base_dir, "valid/EAS/", trait, "/group_", i)
        test_eur <- paste0(eur_base_dir, "test/fold_", i, "/chr", chrom)
        test_eas_pheno <- paste0(base_dir, "valid/EAS/", trait, "/group_", i, ".fam")
        test_eur_pheno <- paste0(base_dir, "valid/EUR/", trait, "/group_", i, "_chr", chrom, ".fam")
        if (!file.exists(test_eas_pheno)){
          print(paste("    EAS testing genotype file does not exist for group", i, "of trait", trait, "chromosome", chrom))
          next
        }
        if (!file.exists(test_eur_pheno)){
          print(paste("    EUR testing genotype file does not exist for group", i, "of trait", trait, "chromosome", chrom))
          next
        }
        file_test_pheno <- paste0(test_eur_pheno, ",", test_eas_pheno)
        file_test_input <- paste0(test_eur, ",", test_eas)
        output_file_prefix <- paste0(output_dir, "chr", chrom)
        if (file.exists(paste0(output_file_prefix, "/EUR/R2.txt")) && file.exists(paste0(output_file_prefix, "/EAS/R2.txt"))){
          print(paste("      Output file already exists for group", i, "of trait", trait, "chromosome", chrom, ", skipping..."))
          next
        }
        # Run LassoSum
        command <- paste0("Rscript ", package, "/scripts/lassosum2.R --PATH_package ", package, " --PATH_out ", output_file_prefix, " --PATH_plink /data1/jiapl_group/lishuhua/software/general/plink2 --FILE_sst ", sst_eur, ",", sst_eas, " --pop EUR,EAS --chrom ", chrom, " --bfile_tuning ", file_tuning_input, " --pheno_tuning ", file_tuning_pheno, " --bfile_testing ", file_test_input, " --pheno_testing ", file_test_pheno, " --testing TRUE --NCORES 5")
        print(paste("      Running command:", command))
        system(command)
    }
  }
}