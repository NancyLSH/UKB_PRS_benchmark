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
num_cores <- 6
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# print(paste("已注册并行后端，使用", num_cores, "个核心。"))
print(paste0("Have registered parallel backend with ", num_cores, " cores."))


# --- 全局变量和路径定义 (与原脚本相同) ---
# trait_list <- c("waist", "height", "pulse", "dbp", "sbp", "smoke", "drink", "bmi", "wbc", "rbc", "hb", "plt", "lymph", "mono", "neut", "eos", "alt", "ast", "bun", "cholesterol", "creatinine", "ggt", "glucose", "hdl", "ldl", "triglycerides", "ua")
trait_list <- c("waist")

EUR_sst_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/ctsleb/train/EUR/"
EAS_sst_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/ctsleb/train/EAS/"
EAS_test_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/ctsleb/tuning/EAS/"
EUR_ref_plinkfile <- "/data1/jiapl_group/lishuhua/software/PRS/CT_SLEB/reference/EUR/merged/chr_all"
EAS_ref_plinkfile <- "/data1/jiapl_group/lishuhua/software/PRS/CT_SLEB/reference/EAS/merged/chr_all"
plink19_exec <- "/data1/jiapl_group/lishuhua/software/general/plink"
plink2_exec <- "/data1/jiapl_group/lishuhua/software/general/plink2"
output_base_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/ctsleb/res/ct_res/"
beta_output_base_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/ctsleb/res/beta/"
r2_output_base_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/software/ctsleb/res/r2/"
y_base_dir <- "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/Cross_Validation/CAS/"

# 3. 使用 foreach 并行执行对每个表型的处理
# .packages 参数确保每个核心都加载了必要的库
# .errorhandling = 'pass' 确保一个表型的失败不会中断其他表型

for(trait in trait_list) {
    print(paste("Processing trait:", trait))   
    # 整个原始 for 循环的主体部分都放在这里
    # 对于每个表型，这段代码将在一个单独的核心上运行
    # create a ctsleb r2 result dataframe
    best_snp_result <- data.frame(trait = character(), group = integer(), best_snp = character(), best_r2 = numeric())
    ctsleb_r2_result <- data.frame(trait = character(), group = integer(), r2 = numeric())
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

        PRS_farm <- SetParamsFarm(plink19_exec = plink19_exec, plink2_exec = plink2_exec, mem=60000, threads = 4, r2_vec = c(0.01,0.05), wc_base_vec = c(50), pthres = c(1))
        output_dir <- paste0(output_base_dir, trait, "/group_", i, "/")
        beta_output_dir <- paste0(beta_output_base_dir, trait, "/group_", i, "/")
        if (!dir.exists(beta_output_dir)) {
            dir.create(beta_output_dir, recursive = TRUE)
        }
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }

        # Run dimCT and capture output for logging
        prs_mat <- dimCT(
            results_dir = output_dir,
            sum_target = EAS_sst,
            sum_ref = EUR_sst,
            ref_plink = EUR_ref_plinkfile,
            target_plink = EAS_ref_plinkfile,
            test_target_plink = EAS_test_plinkfile,
            out_prefix = "res",
            params_farm = PRS_farm
        )

        print(paste("dimCT completed for trait:", trait, "in group:", i))
        
        if (is.null(prs_mat) || nrow(prs_mat) == 1) {
            message(paste("PRS matrix is empty or does not exist for group", i, "of trait", trait))
            next
        }

        prs_mat <- as.data.frame(prs_mat)
        write.table(prs_mat, file = paste0(output_dir, "prs_mat.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        colnames(prs_mat)[colnames(prs_mat) == "#FID"] <- "FID"
        EAS_tune_ids <- EAS_tune[, c("FID", "IID")]
        # change FID and IID into character to avoid factor issues
        EAS_tune_ids$FID <- as.character(EAS_tune_ids$FID)
        EAS_tune_ids$IID <- as.character(EAS_tune_ids$IID)
        prs_tune <- left_join(EAS_tune_ids, prs_mat, by = c("FID", "IID"))
        print(head(prs_tune))
        # check if EAS_tune and prs_tune have the same number of rows
        if (nrow(EAS_tune) != nrow(prs_tune)) {
            message(paste("EAS tune and PRS tune have different number of rows for group", i, "of trait", trait))
            next
        }
        n.total.prs <- length(pthres)^2*length(r2_vec)*length(wc_base_vec)
        n_col_prs <- ncol(prs_tune) - 2  # -2 for FID and IID
        if (n.total.prs != n_col_prs) {
            message(paste("Number of PRS columns does not match expected count for group", i, "of trait", trait))
            next
        }
        prs_r2_vec_test <- rep(0, n.total.prs)

        print(paste("Calculating R2 for each PRS in group", i, "of trait", trait))
        
        for (p_ind in 1:n.total.prs){
            model <- lm(EAS_tune$Pheno ~ prs_tune[[p_ind + 2]])
            prs_r2_vec_test[p_ind] <- summary(model)$r.squared
        }

        max_ind <- which.max(prs_r2_vec_test)
        print(paste("Max R2 for trait", trait, "in group", i, "is", max(prs_r2_vec_test), "for PRS", max_ind))
        best_snp_result <- rbind(best_snp_result, data.frame(trait = trait, group = i, best_snp = colnames(prs_mat)[max_ind + 2], best_r2 = max(prs_r2_vec_test)))
        
        best_snps <- colnames(prs_mat)[max_ind + 2]
        prs_mat_eb <- CalculateEBEffectSize(
            bfile = EAS_test_plinkfile,
            snp_ind = best_snps,
            plink_list = plink_list,
            out_prefix = "prs_eb",
            results_dir = output_dir,
            params_farm = PRS_farm
        )

        print(paste("CalculateEBEffectSize completed for trait:", trait, "in group:", i))
        
        prs_mat_eb_2 <- as.data.frame(prs_mat_eb)
        if (is.null(prs_mat_eb_2) || nrow(prs_mat_eb_2) == 0) {
            message(paste("PRS matrix with EB effect size is empty or does not exist for group", i, "of trait", trait))
            next
        }
        colnames(prs_mat_eb_2)[colnames(prs_mat_eb_2) == "#FID"] <- "FID"
        prs_mat_eb_2$FID <- as.character(prs_mat_eb_2$FID)
        prs_mat_eb_2$IID <- as.character(prs_mat_eb_2$IID)
        write.table(prs_mat_eb_2, file = paste0(output_dir, "prs_mat_eb.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        EAS_test_ids <- EAS_test[, c("FID", "IID")]
        EAS_test_ids$FID <- as.character(EAS_test_ids$FID)
        EAS_test_ids$IID <- as.character(EAS_test_ids$IID)
        prs_tune <- left_join(EAS_tune_ids, prs_mat_eb_2, by = c("FID", "IID"))
        prs_validation <- left_join(EAS_test_ids, prs_mat_eb_2, by = c("FID", "IID"))
        prs_tune$FID <- as.character(prs_tune$FID)
        prs_tune$IID <- as.character(prs_tune$IID)
        prs_validation$FID <- as.character(prs_validation$FID)
        prs_validation$IID <- as.character(prs_validation$IID)

        cols_to_convert_tune <- setdiff(colnames(prs_tune), c("FID", "IID"))
        for (col in cols_to_convert_tune) {
            prs_tune[, (col) := as.numeric(as.character(prs_tune[[col]]))] # Use data.table syntax for efficiency
        }

        cols_to_convert_validation <- setdiff(colnames(prs_validation), c("FID", "IID"))
        for (col in cols_to_convert_validation) {
            prs_validation[, (col) := as.numeric(as.character(prs_validation[[col]]))] # Use data.table syntax for efficiency
        }

        Cleaned_Data <- PRS_Clean(Tune_PRS = prs_tune, Tune_Y = EAS_tune, Validation_PRS = prs_validation)
        prs_tune_sl <- Cleaned_Data$Cleaned_Tune_PRS
        prs_valid_sl <- Cleaned_Data$Cleaned_Validation_PRS
        
        if (is.null(prs_tune_sl) || nrow(prs_tune_sl) == 0) {
            message(paste("Cleaned tune PRS is empty or does not exist for group", i, "of trait", trait))
            next
        }
        if (is.null(prs_valid_sl) || nrow(prs_valid_sl) == 0) {
            message(paste("Cleaned validation PRS is empty or does not exist for group", i, "of trait", trait))
            next
        }

        SL.library <- c("SL.glmnet", "SL.ridge")
        sl <- SuperLearner(Y = EAS_tune$Pheno, X = prs_tune_sl[, -c(1, 2)], family = gaussian(), SL.library = SL.library)
        y_pred_valid <- predict(sl, prs_valid_sl[, -c(1, 2)], onlySL = TRUE)
        model <- lm(EAS_test$Pheno ~ y_pred_valid$pred)
        r2_ctsleb <- summary(model)$r.squared
        print(paste("R2 for trait", trait, "in group", i, "is", r2_ctsleb))

        ctsleb_r2_result <- rbind(ctsleb_r2_result, data.frame(trait = trait, group = i, r2 = r2_ctsleb))
        y_pred_tune <- predict(sl, prs_tune_sl[, -c(1, 2)], onlySL = TRUE)
        Predicted_Tune_Y <- y_pred_tune$pred
        Tune_PRS <- prs_tune_sl[, -c(1, 2)]

        # ... (Your previous code to predict y_pred_tune, etc.)

        # =======================================================================
        # --- FINAL FIX: Proactively handle NA values in P-value columns ---
        # =======================================================================
        message("--- FINAL FIX: Handling NA values in P-value columns before the final step ---")

        na_p_before <- sum(is.na(unique_infor_post$P))
        na_pref_before <- sum(is.na(unique_infor_post$P_ref))

        if (na_p_before > 0 || na_pref_before > 0) {
            message(paste("Found", na_p_before, "NAs in P column and", na_pref_before, "NAs in P_ref column. Replacing with 1."))
            # Replace any NA p-values with 1. This is a safe, non-significant value.
            unique_infor_post$P[is.na(unique_infor_post$P)] <- 1
            unique_infor_post$P_ref[is.na(unique_infor_post$P_ref)] <- 1
            message("... NA replacement complete.")
        } else {
            message("... No NA values found in P-value columns. Proceeding.")
        }

        # =======================================================================
        # Now, call the function with the cleaned data
        # =======================================================================

        Final_Betas <- NULL # Initialize to NULL

        # Use a final tryCatch as a safety net
        tryCatch({
            Final_Betas <- ExtractFinalBetas(
            Tune_PRS = prs_tune_sl[, -c(1, 2)],
            Predicted_Tune_Y = y_pred_tune$pred,
            prs_mat_eb = prs_mat_eb,
            unique_infor_post = unique_infor_post, # Now cleaned of NAs
            pthres = pthres # Use the full pthres vector here
        )
        }, error = function(e) {
            message("!!!!!!!! ERROR during ExtractFinalBetas call even after NA fix. !!!!!!!!")
            message("Original R Error Message:")
            message(e)
        })

        if (is.null(Final_Betas)) {
            message("ExtractFinalBetas failed to produce betas. Skipping file write for this group.")
            next
        }
    
        # If successful, write the betas to the file
        output_beta_path <- paste0(beta_output_dir, "final_betas.txt")
        write.table(Final_Betas, file = output_beta_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        # Final_Betas <- ExtractFinalBetas(Tune_PRS = prs_tune_sl[, -c(1,2)], Predicted_Tune_Y = y_pred_tune$pred, prs_mat_eb = prs_mat_eb, unique_infor_post = unique_infor_post, pthres = pthres)

        # output_beta_path <- paste0(beta_output_dir, "final_betas.txt")
        # write.table(Final_Betas, file = output_beta_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
    # 将每个 trait 的结果保存到一个文件中
    output_r2_path <- paste0(r2_output_base_dir, trait)
    if (!dir.exists(dirname(output_r2_path))) {
        dir.create(dirname(output_r2_path), recursive = TRUE)
    }
    write.table(ctsleb_r2_result, file = paste0(output_r2_path, "/insample_r2.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    # 将每个 trait 的 best_snp_result 保存到一个文件中
    output_best_snp_path <- paste0(output_base_dir, trait)
    if (!dir.exists(dirname(output_best_snp_path))) {
        dir.create(dirname(output_best_snp_path), recursive = TRUE)
    }
    write.table(best_snp_result, file = paste0(output_best_snp_path, "/best_snps.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    # 返回成功信息，用于最终的日志/结果汇总
    return(paste("Successfully processed trait:", trait))
}

# 4. 停止并行集群，释放资源
stopCluster(cl)

# (可选) 打印并行任务的结果摘要
print("Parallel processing completed for all traits.")