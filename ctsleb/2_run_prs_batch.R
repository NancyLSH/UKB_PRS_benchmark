library(CTSLEB)
library(data.table)
library(dplyr)
library(caret)
library(SuperLearner)
library(ranger)
library(glmnet)

# /data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/Cross_Validation/CAS/alt/group_1/pheno
# /data1/jiapl_group/lishuhua/project/PRS_benchmark/software/ctsleb/tuning/EAS/alt

# trait_list <- c("height", "pulse", "dbp", "sbp", "smoke", "drink", "bmi", "weight", "wbc", "rbc", "hb", "plt", "lymph", "mono", "neut", "eos", "alt", "ast", "bun", "cholesterol", "creatinine", "ggt", "glucose", "hdl", "ldl", "triglycerides", "ua")
# trait_list <- c("height", "pulse", "dbp")
# trait_list <- c("sbp", "smoke", "drink")
# trait_list <- c("bmi")
# trait_list <- c("plt")
# trait_list <- c("wbc", "lymph", "mono", "neut")
# trait_list <- c("alt", "ast")
# trait_list <- c("bun", "cholesterol", "creatinine")
# trait_list <- c("ggt", "glucose", "hdl")
# trait_list <- c("ldl", "triglycerides", "ua")
# trait_list <- c('waist')
trait_list <- c('drink')

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

for (trait in trait_list) {
    # Step1: read the summary statistic of EUR and EAS
    EUR_sst_path <- paste0(EUR_sst_dir, trait, "_final.txt")
    EUR_sst <- fread(EUR_sst_path, header = TRUE)
    if (is.null(EUR_sst) || nrow(EUR_sst) == 0) {
        message(paste("EUR summary statistic file is empty or does not exist for trait:", trait))
        next
    }
    for (i in 3:10){
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
        EAS_tune$FID <- as.character(EAS_tune$FID)
        EAS_tune$IID <- as.character(EAS_tune$IID)
        EAS_test$FID <- as.character(EAS_test$FID)
        EAS_test$IID <- as.character(EAS_test$IID)
        colnames(EAS_tune)[colnames(EAS_tune) == "FID"] <- "#FID"
        colnames(EAS_test)[colnames(EAS_test) == "FID"] <- "#FID"
        print(paste("Processing trait:", trait, "for group:", i))
        # print the first few rows of EAS_sst
        print(head(EAS_sst))
        print(head(EUR_sst))
        # PRS_farm <- SetParamsFarm(plink19_exec = plink19_exec, plink2_exec = plink2_exec, mem=60000, threads = 4, r2_vec = c(0.01,0.05), wc_base_vec = c(50), pthres = c(1))
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
        # prs_mat look like: ["FID","IID","PRS1", "PRS2", ...] but have no header
        if (is.null(prs_mat) || nrow(prs_mat) == 0) {
            message(paste("PRS matrix is empty or does not exist for group", i, "of trait", trait))
            next
        }
        prs_mat <- as.data.frame(prs_mat)
        prs_mat$`#FID` <- as.character(prs_mat$`#FID`)
        prs_mat$IID <- as.character(prs_mat$IID)
        write.table(prs_mat, file = paste0(output_dir, "prs_mat.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        EAS_tune_ids <- EAS_tune[, c("#FID", "IID")]
        prs_tune <- left_join(EAS_tune_ids, prs_mat, by = c("#FID", "IID"))
        n.total.prs <- length(pthres)^2*length(r2_vec)*length(wc_base_vec)
        prs_r2_vec_test <- rep(0, n.total.prs)
        for (p_ind in 1:n.total.prs){
            model <- lm(EAS_tune$Pheno ~ prs_tune[[p_ind + 2]])
            prs_r2_vec_test[p_ind] <- summary(model)$r.squared
        }

        max_ind <- which.max(prs_r2_vec_test)
        print(paste("Max R2 for trait", trait, "in group", i, "is", max(prs_r2_vec_test), "for PRS", max_ind))
        best_snps <- colnames(prs_mat)[max_ind + 2]
        prs_mat_eb <- CalculateEBEffectSize(bfile = EAS_test_plinkfile, snp_ind = best_snps, plink_list = plink_list, out_prefix = "prs_eb", results_dir = output_dir, params_farm = PRS_farm)
        EAS_test_ids <- EAS_test[, c("#FID", "IID")]
        prs_mat_eb <- as.data.frame(prs_mat_eb)
        prs_mat_eb$`#FID` <- as.character(prs_mat_eb$`#FID`)
        prs_mat_eb$IID <- as.character(prs_mat_eb$IID)
        if (is.null(prs_mat_eb) || nrow(prs_mat_eb) == 0) {
            message(paste("PRS matrix with effect size is empty or does not exist for group", i, "of trait", trait))
            next
        }
        write.table(prs_mat_eb, file = paste0(output_dir, "prs_mat_eb.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        prs_tune <- left_join(EAS_tune_ids, prs_mat_eb, by = c("#FID", "IID"))
        prs_validation <- left_join(EAS_test_ids, prs_mat_eb, by = c("#FID", "IID"))


        prs_tune_eb <- inner_join(EAS_tune[, c("#FID", "IID")], prs_mat_eb, by = c("#FID", "IID"))
        y_tune_eb <- inner_join(EAS_tune[, c("#FID", "IID", "Pheno")], prs_mat_eb, by = c("#FID", "IID"))
        prs_validation_eb <- inner_join(EAS_test[, c("#FID", "IID")], prs_mat_eb, by = c("#FID", "IID"))
        Tune_Y_cleaned <- y_tune_eb[, c("Pheno")]
        print(head(prs_tune_eb))
        print(head(prs_validation_eb))
        print(head(Tune_Y_cleaned))
        
        Cleaned_Data <- PRS_Clean(
            Tune_PRS = prs_tune_eb, 
            Tune_Y = Tune_Y_cleaned, 
            Validation_PRS = prs_validation_eb
        )
        if (is.null(Cleaned_Data)) {
            message(paste("Cleaned data is NULL for group", i, "of trait", trait))
            next
        }
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

        y_pred_tune <- predict(sl, prs_tune_sl[, -c(1, 2)], onlySL = TRUE)
        Predicted_Tune_Y <- y_pred_tune$pred
        Tune_PRS <- prs_tune_sl[, -c(1, 2)]

       

        Final_Betas <- ExtractFinalBetas(Tune_PRS = Tune_PRS, Predicted_Tune_Y = Predicted_Tune_Y, prs_mat_eb = prs_mat_eb, unique_infor_post = unique_infor_post, pthres = pthres)

        # Save the results
        output_beta_path <- paste0(beta_output_dir, "final_betas.txt")
        write.table(Final_Betas, file = output_beta_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
}