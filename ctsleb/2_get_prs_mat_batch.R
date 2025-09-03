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
    
    # 整个原始 for 循环的主体部分都放在这里
    # 对于每个表型，这段代码将在一个单独的核心上运行

    # Step1: read the summary statistic of EUR
    EUR_sst_path <- paste0(EUR_sst_dir, trait, "_final.txt")
    EUR_sst <- fread(EUR_sst_path, header = TRUE)
    if (is.null(EUR_sst) || nrow(EUR_sst) == 0) {
        # 在 foreach 中，使用 return() 退出当前迭代，而不是 next
        return(paste("EUR summary statistic file is empty or does not exist for trait:", trait))
    }

    # 针对每个表型的内部循环 (groups 1-10) 保持不变
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

        PRS_farm <- SetParamsFarm(plink19_exec = plink19_exec, plink2_exec = plink2_exec, mem=60000, threads = 4)
        output_dir <- paste0(output_base_dir, trait, "/group_", i, "/")
        beta_output_dir <- paste0(beta_output_base_dir, trait, "/group_", i, "/")
        if (!dir.exists(beta_output_dir)) {
            dir.create(beta_output_dir, recursive = TRUE)
        }
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }

        prs_mat <- dimCT(results_dir = output_dir, sum_target = EAS_sst, sum_ref = EUR_sst, ref_plink = EUR_ref_plinkfile, target_plink = EAS_ref_plinkfile, test_target_plink= EAS_test_plinkfile, out_prefix = "res", params_farm = PRS_farm)
        
        if (is.null(prs_mat) || nrow(prs_mat) == 0) {
            message(paste("PRS matrix is empty or does not exist for group", i, "of trait", trait))
            next
        }

        prs_mat <- as.data.frame(prs_mat)
        write.table(prs_mat, file = paste0(output_dir, "prs_mat.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        colnames(prs_mat) <- c("FID", "IID", paste0("PRS", 1:(ncol(prs_mat) - 2)))
        prs_tune <- left_join(EAS_tune, prs_mat, by = c("FID", "IID"))
        
        n.total.prs <- length(pthres)^2*length(r2_vec)*length(wc_base_vec)
        prs_r2_vec_test <- rep(0, n.total.prs)
        
        for (p_ind in 1:n.total.prs){
            model <- lm(as.formula(paste("Pheno ~", paste0("PRS", p_ind))), data = prs_tune)
            prs_r2_vec_test[p_ind] <- summary(model)$r.squared
        }

        max_ind <- which.max(prs_r2_vec_test)
        print(paste("Max R2 for trait", trait, "in group", i, "is", max(prs_r2_vec_test), "for PRS", max_ind))
    }

    # 返回成功信息，用于最终的日志/结果汇总
    return(paste("Successfully processed trait prs mat:", trait))
}

# 4. 停止并行集群，释放资源
stopCluster(cl)

# (可选) 打印并行任务的结果摘要
print("Parallel processing completed for all traits.")
print(results)