package <- "/data1/jiapl_group/lishuhua/software/PRS/PROSPER/PROSPER-main/"
base_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/prosper/"

# trait_list <- c("waist")
# trait_list <- c("waist", "height", "pulse", "dbp")
# trait_list <- c("sbp", "smoke", "drink", "bmi")
# trait_list <- c("wbc", "rbc", "hb", "plt")
# trait_list <- c("lymph", "mono", "neut", "eos")
# trait_list <- c("alt", "ast", "bun", "cholesterol")
# trait_list <- c("creatinine", "ggt", "glucose", "hdl")
# trait_list <- c("ldl", "triglycerides", "ua")
trait_list <- c("smoke", "drink")

for (trait in trait_list) {
  print(paste("Processing trait:", trait))
  for (i in 1:10){
    print(paste("  Processing fold:", i))
    output_dir <- paste0(base_dir, "res/", trait, "/group_", i, "/prosper/")
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
        # load lassosum2 parameters
        lassosum_param_eas <- paste0(base_dir, "res/", trait, "/group_", i, "/lassosum2/chr", chrom , "/EAS/optimal_param.txt")
        lassosum_param_eur <- paste0(base_dir, "res/", trait, "/group_", i, "/lassosum2/chr", chrom , "/EUR/optimal_param.txt")
        if (!file.exists(lassosum_param_eas) || !file.exists(lassosum_param_eur)){
          print(paste("    Lassosum2 parameter file does not exist for group", i, "of trait", trait, "chromosome", chrom))
          next
        }
        
        # Run PROSPER
        output_prefix <- paste0(output_dir, "chr", chrom)
        if (!dir.exists(output_prefix)){
          dir.create(output_prefix, recursive = TRUE)
        }
        # check if output already exists
        # /data1/jiapl_group/lishuhua/project/PRS_benchmark/software/prosper/res/alt/group_1/prosper/chr1/before_ensemble/score_param.txt
        if (file.exists(paste0(output_prefix, "/before_ensemble/score_param.txt"))){
          print(paste("    Output already exists for group", i, "of trait", trait, "chromosome", chrom, ", skipping..."))
          next
        }

        command <- paste0("Rscript ", package, "/scripts/PROSPER.R --PATH_package ", package, " --PATH_out ", output_prefix, " --FILE_sst ", sst_eur, ",", sst_eas, " --pop EUR,EAS --chrom ", chrom, " --lassosum_param ", lassosum_param_eur, ",", lassosum_param_eas, " --NCORES 5")
        print(paste("      Running command:", command))
        system(command)
    }
  }
}