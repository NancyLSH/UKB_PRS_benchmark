library(lassosum)
library(data.table)
library(parallel)

ref_EAS <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/tlprs/reference/EAS_1kg/1000G.EAS.QC.hm3.ind"
# ref_EAS <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/CAS/geno/CAS_final/CAS_merged_qc_final"
ref_EAS_base <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/tlprs/train/EAS/tune/"
gwas_EUR_base <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/merged_gwas/White_British/gwas/"
output_EUR_base <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/tlprs/beta/1_lassosum/"
LDblocks <- "ASN.hg19"

trait_dict = list(
    # p48 = "waist"
    # p50 = "height",
    # p102 = "pulse",
    # p4079 = "dbp",
    # p4080 = "sbp",
    # p20116 = "smoke",
    # p20117 = "drink",
    # p21001 = "bmi",
    # # p21002 = "weight",
    # p30000 = "wbc",
    # p30010 = "rbc",
    # p30020 = "hb",
    # p30080 = "plt",
    # p30120 = "lymph",
    # p30130 = "mono",
    p30140 = "neut",
    p30150 = "eos",
    # p30620 = "alt",
    # p30650 = "ast",
    # p30670 = "bun",
    # p30690 = "cholesterol",
    # p30700 = "creatinine",
    p30730 = "ggt",
    # p30740 = "glucose",
    # p30760 = "hdl",
    p30780 = "ldl",    
    p30870 = "triglycerides"
    # p30880 = "ua"
)

# 遍历每个 trait
for (trait in names(trait_dict)) {
    print(paste0("Processing trait: ", trait_dict[[trait]]))
    print(paste0("Trait code: ", trait))
    
    # 设置输出文件路径
    output_prefix <- paste0(output_EUR_base, trait_dict[[trait]])
    if (!dir.exists(output_prefix)) {
        dir.create(output_prefix, recursive = TRUE)
    }
    
    # 读取 GWAS summary statistics: p102_int.merged.glm.linear
    gwas_file_path <- paste0(gwas_EUR_base, trait, "_int.merged.glm.linear")
    gwas_file <- fread(gwas_file_path, header = TRUE)
    setnames(gwas_file, "#CHROM", "CHR")
    # change columns type of P
    gwas_file$P <- as.numeric(gwas_file$P)
    # remove P == 0 or P == NA
    gwas_file <- gwas_file[gwas_file$P != 0 & !is.na(gwas_file$P), ]
    # remove OBS_CT == NA
    gwas_file <- gwas_file[!is.na(gwas_file$OBS_CT), ]
    if (!file.exists(gwas_file_path)) {
        stop(paste0("GWAS file not found: ", gwas_file_path))
    }
    if (trait == 'p20116' || trait == 'p20117') {
        # remove OR == NA
        gwas_file <- gwas_file[!is.na(gwas_file$OR), ]
        cor <- p2cor(p = gwas_file$P, n = gwas_file$OBS_CT, sign = log(gwas_file$OR))
    }else{
        # remove BETA == NA
        gwas_file <- gwas_file[!is.na(gwas_file$BETA), ]
        cor <- p2cor(p = gwas_file$P, n = gwas_file$OBS_CT, sign = gwas_file$BETA)
    }
    if (any(is.na(cor))){
        # delete the NA value in cor
        cor_na <- which(is.na(cor))
        cor <- cor[!is.na(cor)]
        # delete the same index in ss
        gwas_file <- gwas_file[-cor_na, ]
    }
    for (i in 1:10) {
        print(paste0("Processing P", i))
        # 读取参考基因组数据
        test_EAS <- paste0(ref_EAS_base, trait_dict[[trait]], "/group_", i, "/tune")
        pheno_path <- paste0(ref_EAS_base, trait_dict[[trait]], "/group_", i, "/pheno.txt")
        covar_path <- paste0(ref_EAS_base, trait_dict[[trait]], "/group_", i, "/covar.txt")
        pheno_data <- fread(pheno_path, header = TRUE)
        pheno_data <- as.data.frame(pheno_data)
        covar_data <- fread(covar_path, header = TRUE)
        covar_data <- as.data.frame(covar_data)
        # make sure FID and IID are character, other columns are numeric
        pheno_data$FID <- as.character(pheno_data$FID)
        pheno_data$IID <- as.character(pheno_data$IID)
        pheno_data$Pheno <- as.numeric(pheno_data$Pheno)
        covar_data$FID <- as.character(covar_data$FID)
        covar_data$IID <- as.character(covar_data$IID)
        if (ncol(covar_data) > 2) {
            covar_data[, -c(1, 2)] <- lapply(covar_data[, -c(1, 2)], as.numeric)
        }
        print(head(pheno_data))
        print(head(covar_data))
        cl <- parallel::makePSOCKcluster(2)
        out <- lassosum.pipeline(
            cor = cor,
            LDblocks = LDblocks,
            snp = gwas_file$ID,
            chr = gwas_file$CHR,
            pos = gwas_file$POS,
            A1 = gwas_file$ALT,
            A2 = gwas_file$REF,
            ref.bfile = ref_EAS,
            test.bfile = test_EAS,
            cluster = cl
        )
        v <- validate(out, pheno = pheno_data, covar = covar_data)
        sumstats <- cbind(out$sumstats, v$best.beta)
        # 保存结果
        fwrite(sumstats, paste0(output_prefix, "/group_", i, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        stopCluster(cl)
        print(paste0("Finished processing P", i, " for trait: ", trait_dict[[trait]]))
    }
}