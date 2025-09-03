library(XPASS)
library(data.table)
library(RhpcBLASctl)
# blas_set_num_threads(30)

ref_EAS <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/CAS/geno/CAS_final/CAS_merged_qc_final"
ref_EUR <- "/data1/jiapl_group/lishuhua/software/PRS/XPASS/reference/1000G.EUR.QC.hm3.ind"
covar_EAS <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/CAS/pheno/covar_for_xpass.txt"
covar_EUR <- "/data1/jiapl_group/lishuhua/software/PRS/XPASS/reference/1000G.EUR.QC.hm3.ind.pc20.txt"
test <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/geno/Chinese/1_merged/merged"
sst_eas_path <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/xpass/train/EAS/"
sst_eur_path <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/xpass/train/EUR/"
output_path <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/xpass/res/data/"

trait_dict = list(
    p20116 = "smoke",
    p20117 = "drink",
    p30620 = "alt",
    p30650 = "ast",
    p30670 = "bun",
    p30690 = "cholesterol",
    p30700 = "creatinine",
    p30730 = "ggt",
    p30740 = "glucose",
    p30760 = "hdl",
    p30780 = "ldl",    
    p30870 = "triglycerides",
    p30880 = "ua"
)

# for each trait, run XPASS
for (trait in names(trait_dict)) {
    # set the output file name
    print(paste0("Running trait: ", trait_dict[[trait]]))
    print(paste0("Trait code: ", trait))
    if (trait == "p20116" || trait == "p20117") {
        # for smoking and drinking, use the 3M format
        sst_eas <- paste0(sst_eas_path, trait_dict[[trait]], "_raw.", trait_dict[[trait]], "_raw.txt")
        sst_eur <- paste0(sst_eur_path, trait, "_int.merged.txt")
    } else {
        # for other traits, use the standard format
        sst_eas <- paste0(sst_eas_path, trait_dict[[trait]], "_int.", trait_dict[[trait]], "_int.txt")
        sst_eur <- paste0(sst_eur_path, trait, "_int.merged.txt")
    }
    output_file_prefix <- paste0(output_path, trait_dict[[trait]], "_", trait)
    # run XPASS
    fit_test <- XPASS(file_z1 = sst_eas, 
                      file_z2 = sst_eur, 
                      file_ref1 = ref_EAS, 
                      file_ref2 = ref_EUR, 
                      file_cov1 = covar_EAS, 
                      file_cov2 = covar_EUR, 
                      file_predGeno = test, 
                      pop = "EAS", 
                      compPRS = TRUE, 
                      sd_method = "LD_block", 
                      compPosMean = TRUE, 
                      file_out = output_file_prefix)
    print(paste0("Finished running trait: ", trait_dict[[trait]]))
}