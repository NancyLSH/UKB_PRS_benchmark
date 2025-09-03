# 1. 加载所有必需的库
library(CTSLEB)
library(data.table)
library(dplyr)
library(caret)
library(SuperLearner)
library(ranger)
library(glmnet)
# 加载用于并行的库
library(foreach)
library(doParallel)

# 2. 设置并行处理环境
# 使用系统可用核心数减1，以保留一个核心用于系统操作
# 您也可以手动设置核心数，例如 num_cores <- 8
num_cores <- 4
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# print(paste("已注册并行后端，使用", num_cores, "个核心。"))
print(paste0("Have registered parallel backend with ", num_cores, " cores."))


# --- 全局变量和路径定义 (与原脚本相同) ---
trait_list <- c("waist", "height", "pulse", "dbp", "sbp", "smoke", "drink", "bmi", "wbc", "rbc", "hb", "plt", "lymph", "mono", "neut", "eos", "alt", "ast", "bun", "cholesterol", "creatinine", "ggt", "glucose", "hdl", "ldl", "triglycerides", "ua")

EUR_sst_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/ctsleb/train/EUR/"
EAS_sst_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/ctsleb/train/EAS/"
EAS_test_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/ctsleb/tuning/EAS/"
EUR_ref_plinkfile <- "/data1/jiapl_group/lishuhua/software/PRS/CT_SLEB/reference/EUR/merged/chr_all"
EAS_ref_plinkfile <- "/data1/jiapl_group/lishuhua/software/PRS/CT_SLEB/reference/EAS/merged/chr_all"
plink19_exec <- "/data1/jiapl_group/lishuhua/software/general/plink"
plink2_exec <- "/data1/jiapl_group/lishuhua/software/general/plink2"
output_base_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/ctsleb/res/ct_res/"
beta_output_base_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/ctsleb/res/beta/"
y_base_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/Cross_Validation/CAS/"

# 3. 使用 foreach 并行执行对每个表型的处理
# .packages 参数确保每个核心都加载了必要的库
# .errorhandling = 'pass' 确保一个表型的失败不会中断其他表型
results <- foreach(
    trait = trait_list,
    .packages = c("CTSLEB", "data.table", "dplyr", "SuperLearner", "glmnet", "caret", "ranger"),
    .errorhandling = 'pass'
) %dopar% {
    # check if the trait has prs_mat.txt, if not, then skip
    prs_mat_path <- paste0(output_base_dir, trait, "/prs_mat.txt")
    if (!file.exists(prs_mat_path)) {
        return(paste("PRS matrix does not exist for trait:", trait))
    }
    # check if eb_mat exists
    eb_mat_path <- paste0(output_base_dir, trait, "/prs_mat_eb.txt")
    if (file.exists(eb_mat_path)) {
        return(paste("EB matrix already exists for trait:", trait))
    }
    # Step1: read the summary statistic of EUR
    # 针对每个表型的内部循环 (groups 1-10) 保持不变
    EUR_sst_path <- paste0(EUR_sst_dir, trait, "_final.txt")
    EUR_sst <- fread(EUR_sst_path, header = TRUE)
    if (is.null(EUR_sst) || nrow(EUR_sst) == 0) {
        # 在 foreach 中，使用 return() 退出当前迭代，而不是 next
        return(paste("EUR summary statistic file is empty or does not exist for trait:", trait))
    }
    for (i in 1:10){
        EAS_sst_path <- paste0(EAS_sst_dir, trait, "/group_", i, "_final.txt")
        EAS_tune_path <- paste0(y_base_dir, trait, "/group_", i, "/pheno/tune_pheno.txt")
        EAS_test_path <- paste0(y_base_dir, trait, "/group_", i, "/pheno/test_pheno.txt")
        EAS_tune <- fread(EAS_tune_path, header = TRUE)
        EAS_test <- fread(EAS_test_path, header = TRUE)
        EAS_test_plinkfile <- paste0(EAS_test_dir, trait, "/group_", i, "_final")

        if (is.null(EAS_tune) || nrow(EAS_tune) == 0) {
            message(paste("EAS tune phenotype file is empty or does not exist for group", i, "of trait", trait))
            next
        }
        if (is.null(EAS_test) || nrow(EAS_test) == 0) {
            message(paste("EAS test phenotype file is empty or does not exist for group", i, "of trait", trait))
            next
        }
        if (!file.exists(EAS_sst_path)) {
            message(paste("EAS summary statistic file does not exist for group", i, "of trait", trait))
            next
        }
        EAS_sst <- fread(EAS_sst_path, header = TRUE)
        if (is.null(EAS_sst) || nrow(EAS_sst) == 0) {
            message(paste("EAS summary statistic file is empty or does not exist for group", i, "of trait", trait))
            next
        }
        
        print(paste("Processing trait:", trait, "for group:", i))
        print(head(EAS_sst))
        print(head(EUR_sst))

        output_dir <- paste0(output_base_dir, trait, "/group_", i, "/")
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        prs_mat <- read.table(file = paste0(output_dir, "prs_mat.txt"), header = TRUE, sep = "\t")
        # replace the colname #FID in prs_mat with FID
        colnames(prs_mat)[colnames(prs_mat) == "#FID"] <- "FID"
        prs_tune <- left_join(EAS_tune, prs_mat, by = c("FID", "IID"))
        
        n.total.prs <- length(pthres)^2*length(r2_vec)*length(wc_base_vec)
        prs_r2_vec_test <- rep(0, n.total.prs)
        
        for (p_ind in 1:n.total.prs){
            model <- lm(as.formula(paste("Pheno ~", paste0("PRS", p_ind))), data = prs_tune)
            prs_r2_vec_test[p_ind] <- summary(model)$r.squared
        }

        max_ind <- which.max(prs_r2_vec_test)
        print(paste("Max R2 for trait", trait, "in group", i, "is", max(prs_r2_vec_test), "for PRS", max_ind))

        best_snps <- colnames(prs_mat)[max_ind + 2]
        prs_mat_eb <- CalculateEBEffectSize(bfile = EAS_test_plinkfile, snp_ind = best_snps, plink_list = plink_list, out_prefix = paste0(output_dir, "prs_eb"), params_farm = PRS_farm)

        prs_mat_eb <- as.data.frame(prs_mat_eb)
        if (is.null(prs_mat_eb) || nrow(prs_mat_eb) == 0) {
            message(paste("PRS matrix with EB effect size is empty or does not exist for group", i, "of trait", trait))
            next
        }
        write.table(prs_mat_eb, file = paste0(output_dir, "prs_mat_eb.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }

    # 返回成功信息，用于最终的日志/结果汇总
    return(paste("Successfully processed trait prs mat:", trait))
}

# 4. 停止并行集群，释放资源
stopCluster(cl)

# (可选) 打印并行任务的结果摘要
print("Parallel processing completed for all traits.")
print(results)