library(lassosum)
library(data.table)
library(parallel)

ref_EAS <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/tlprs/reference/EAS_1kg/1000G.EAS.QC.hm3.ind"
ref_EAS_base <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/ct/test/"
gwas_EAS_base <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/Cross_Validation/CAS/"
lasso_EAS_base <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/lassosum/test/EAS/"
output_EAS_base <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/lassosum/res/beta/EAS/"
LDblocks <- "ASN.hg19"

trait_dict = list(
    p48 = "waist",
    p50 = "height",
    p102 = "pulse",
    p4079 = "dbp",
    p4080 = "sbp",
    p20116 = "smoke",
    p20117 = "drink",
    p21001 = "bmi",
    p30000 = "wbc",
    p30010 = "rbc",
    p30020 = "hb",
    p30080 = "plt",
    p30120 = "lymph",
    p30130 = "mono",
    p30140 = "neut",
    p30150 = "eos",
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

for (trait in names(trait_dict)) {
    print(paste0("Processing trait: ", trait_dict[[trait]]))
    print(paste0("Trait code: ", trait))
    output_prefix <- paste0(output_EAS_base, trait_dict[[trait]])
    if (!dir.exists(output_prefix)) {
        dir.create(output_prefix, recursive = TRUE)
    }

    for (group in 1:10) {
        print(paste0("Processing group: ", group))
        gwas_file_path <- paste0(gwas_EAS_base, trait_dict[[trait]], "/group_", group, "/gwas/train.Pheno.glm.linear")
        gwas_file <- fread(gwas_file_path, header = TRUE)
        setnames(gwas_file, "#CHROM", "CHR")
        gwas_file$P <- as.numeric(gwas_file$P)
        gwas_file <- gwas_file[!is.na(gwas_file$P) & gwas_file$P > 0, ]
        gwas_file <- gwas_file[!is.na(gwas_file$OBS_CT),]
        gwas_file <- gwas_file[!is.na(gwas_file$BETA),]
        cor <- p2cor(p = gwas_file$P, n = gwas_file$OBS_CT, sign = gwas_file$BETA)
        if (any(is.na(cor))) {
            cor_na <- which(is.na(cor))
            cor <- cor[!is.na(cor)]
            gwas_file <- gwas_file[-cor_na, ]
        }
        # /data1/jiapl_group/lishuhua/project/PRS_benchmark/software/ct/test/alt/group_1/combine
        test_EAS <- paste0(ref_EAS_base, trait_dict[[trait]], "/group_", group, "/combine")
        pheno_path <- paste0(lasso_EAS_base, trait_dict[[trait]], "/group_", group, "/pheno.txt")
        covar_path <- paste0(lasso_EAS_base, trait_dict[[trait]], "/group_", group, "/covar.txt")
        pheno <- fread(pheno_path, header = TRUE)
        covar <- fread(covar_path, header = TRUE)
        pheno <- as.data.frame(pheno)
        covar <- as.data.frame(covar)
        pheno$FID <- as.character(pheno$FID)
        pheno$IID <- as.character(pheno$IID)
        pheno$Pheno <- as.numeric(pheno$Pheno)
        covar$FID <- as.character(covar$FID)
        covar$IID <- as.character(covar$IID)
        if (ncol(covar) > 2) {
            covar[, -c(1, 2)] <- lapply(covar[, -c(1, 2)], as.numeric)
        }
        print(head(pheno))
        print(head(covar))

        out <- lassosum.pipeline(
            cor = cor,
            LDblocks = LDblocks,
            snp = gwas_file$ID,
            chr = gwas_file$CHR,
            pos = gwas_file$POS,
            A1 = gwas_file$ALT,
            A2 = gwas_file$REF,
            ref.bfile = ref_EAS,
            test.bfile = test_EAS
        )

        v <- validate(out, pheno = pheno, covar = covar)
        sumstats <- cbind(out$sumstats, v$best.beta)
        write.table(sumstats, file = paste0(output_prefix, "/group_", group, "_beta.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
        print(paste0("Saved results to ", output_prefix, "/group_", group, "_beta.txt"))
    }
}