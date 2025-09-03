library(data.table)
library(lassosum)
library(TLPRS)
library(parallel)

ref_EAS_base <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/tlprs/train/EAS/train/"
gwas_EAS_base <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/tlprs/train/EAS/train/"
valid_EAS_base <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/tlprs/train/EAS/valid/"
beta_EUR_base <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/tlprs/train/EUR/lasso_beta/"
LDblocks <- "ASN.hg19"
output_base <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/tlprs/res/beta/"
Y_name <- "Pheno"
Covar_name <- c('age', 'sex', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')

trait_dict = list(
    # p48 = "waist"
    # p50 = "height",
    # p102 = "pulse",
    # p4079 = "dbp",
    # p4080 = "sbp",
    # p20116 = "smoke",
    # p20117 = "drink",
    # p21001 = "bmi",
    # p21002 = "weight",
    # p30000 = "wbc",
    # p30010 = "rbc",
    # p30020 = "hb",
    # p30080 = "plt",
    # p30120 = "lymph",
    # p30130 = "mono",
    # p30140 = "neut",
    # p30150 = "eos",
    # p30620 = "alt",
    # p30650 = "ast",
    # p30670 = "bun",
    # p30690 = "cholesterol",
    # p30700 = "creatinine",
    # p30730 = "ggt",
    p30740 = "glucose"
    # p30760 = "hdl",
    # p30780 = "ldl",    
    # p30870 = "triglycerides",
    # p30880 = "ua"
)

for (trait in names(trait_dict)) {
    cat("Processing trait:", trait_dict[[trait]], "\n")
    # Add your processing code here
    output_prefix <- paste0(output_base, trait_dict[[trait]])
    if (!dir.exists(output_prefix)) {
        dir.create(output_prefix, recursive = TRUE)
    }
    if (trait == "p20116" || trait == "p20117") {
        Ytype <- "B"
    } else {
        Ytype <- "C"
    }

    for (i in 9:10){
        print(paste0("Processing P: ", i))
        gwas_EAS <- paste0(gwas_EAS_base, trait_dict[[trait]], "/group_", i, "/train_gwas.txt")
        ref_EAS <- paste0(ref_EAS_base, trait_dict[[trait]], "/group_", i, "/train")
        valid_EAS <- paste0(valid_EAS_base, trait_dict[[trait]], "/group_", i, "/valid")
        valid_y_EAS <- paste0(valid_EAS_base, trait_dict[[trait]], "/group_", i, "/pheno_covar.txt")
        beta_EUR <- paste0(beta_EUR_base, trait_dict[[trait]], "/group_", i, "/lasso_beta.txt")
        output_prefix <- paste0(output_base, trait_dict[[trait]], "/group_", i, "/")
        if (!dir.exists(output_prefix)) {
            dir.create(output_prefix, recursive = TRUE)
        }
        output_path <- paste0(output_prefix, "res")
        # check if all files exist
        if (!file.exists(gwas_EAS)) {
            stop(paste0("GWAS file not found: ", gwas_EAS))
        }
        if (!file.exists(valid_y_EAS)) {
            stop(paste0("Validation phenotype file not found: ", valid_y_EAS))
        }
        if (!file.exists(beta_EUR)) {
            stop(paste0("Beta file not found: ", beta_EUR))
        }
        temp <- fread(beta_EUR, nrow=5)
        temp_gwas <- fread(gwas_EAS, nrow=5)
        print(temp)
        print(temp_gwas)
        # run TL_PRS
        cl <- makeCluster(14, type = "FORK")
        system.time({
            out.beta <- TL_PRS(valid_y_EAS, Covar_name, Y_name, Ytype = Ytype, ref_EAS, valid_EAS, beta_EUR, gwas_EAS, LDblock = LDblocks, output_path, cluster = cl)
        })
        stopCluster(cl)
        summary(out.beta)
        # write.table(out.beta, file = paste0(output_prefix, "/beta.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
}