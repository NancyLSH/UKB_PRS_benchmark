package <- "/data1/jiapl_group/lishuhua/software/PRS/PROSPER/PROSPER-main/"
base_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/prosper/"

# trait_list <- c("waist")
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
    output_dir <- paste0(base_dir, "res/", trait, "/group_", i, "/prosper/")
    if (!dir.exists(output_dir)){
      dir.create(output_dir, recursive = TRUE)
    }
    for (chrom in 1:22){
        print(paste("    Processing chromosome:", chrom))
        # tuning and testing data
        tuning_eas <- paste0(base_dir, "tuning/EAS/", trait, "/group_", i)
        tuning_eas_pheno <- paste0(base_dir, "tuning/EAS/", trait, "/group_", i, ".fam")
        if (!file.exists(tuning_eas_pheno)){
          print(paste("    EAS tuning genotype file does not exist for group", i, "of trait", trait, "chromosome", chrom))
          next
        }
        test_eas <- paste0(base_dir, "valid/EAS/", trait, "/group_", i)
        test_eas_pheno <- paste0(base_dir, "valid/EAS/", trait, "/group_", i, ".fam")
        if (!file.exists(test_eas_pheno)){
          print(paste("    EAS testing genotype file does not exist for group", i, "of trait", trait, "chromosome", chrom))
          next
        }

        # Run ensemble PROSPER
        output_prefix <- paste0(output_dir, "chr", chrom)
        if (!dir.exists(output_prefix)){
          dir.create(output_prefix, recursive = TRUE)
        }
        # check if output already exists
        # /data1/jiapl_group/lishuhua/project/PRS_benchmark/software/prosper/res/alt/group_1/prosper/chr1/after_ensemble_EAS/PROSPER_prs_file.txt
        if (file.exists(paste0(output_prefix, "/after_ensemble_EAS/PROSPER_prs_file.txt"))){
          print(paste("    Output already exists for group", i, "of trait", trait, "chromosome", chrom, ", skipping..."))
          next
        }
        command <- paste0("Rscript ", package, "/scripts/tuning_testing.R --PATH_plink /data1/jiapl_group/lishuhua/software/general/plink2 --PATH_out ", output_prefix, " --prefix EAS --bfile_tuning ", tuning_eas, " --pheno_tuning ", tuning_eas_pheno, " --bfile_testing ", test_eas, "  --pheno_testing ", test_eas_pheno, " --NCORES 5")
        print(paste("      Running command:", command))
        system(command)
    }
  }
}