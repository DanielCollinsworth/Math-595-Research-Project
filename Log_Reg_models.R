# start Logging
sink("run_log.txt", split = TRUE)

# ===============================
# 1) Libraries
# ===============================
library(survey)
library(pROC)
library(PRROC)
library(dplyr)
library(caret)
library(splines)
library(ggeffects)


# ===============================
# 2) Load data
# ===============================
data.complete <- readRDS("data_complete_cleaned.rds")

cat("\n[INFO] Data loaded:\n")
cat(" - Rows:", nrow(data.complete), "  Cols:", ncol(data.complete), "\n")
cat(" - Cycles present:", paste(unique(as.character(data.complete$Cycle)), collapse=", "), "\n")


# ===============================
# 3) Split by cycle
# ===============================
train_df <- subset(data.complete, Cycle == "P")  # 2017–2020
test_df  <- subset(data.complete, Cycle == "L")  # 2021–2023

# Quartiles on train and cut both sets
ibi_quartiles <- quantile(train_df$IBI, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
cut_ibi <- function(x, breaks) cut(x, breaks = breaks, labels = c("Q1","Q2","Q3","Q4"), include.lowest = TRUE)

train_df$IBI_Category <- cut_ibi(train_df$IBI, ibi_quartiles)
test_df$IBI_Category  <- cut_ibi(test_df$IBI,  ibi_quartiles)
train_df$IBI_Category <- stats::relevel(train_df$IBI_Category, ref = "Q1")

# Target encodings
train_df$CVD <- factor(train_df$CVD, levels = c("No","Yes"))
test_df$CVD  <- factor(test_df$CVD,  levels = c("No","Yes"))
train_df$CVD_num <- as.integer(train_df$CVD == "Yes")
test_df$CVD_num  <- as.integer(test_df$CVD  == "Yes")


# ===============================
# 4) Build survey design on Training data
# ===============================
train_design <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR,
  nest = TRUE, survey.lonely.psu = "adjust",
  data = train_df
)


# ===============================
# 5) Train weighted logistic models
# ===============================
cat("\n[FIT] Model 1: IBI only\n")
model1 <- svyglm(CVD ~ IBI_Category, design = train_design, family = quasibinomial("logit"))
print(summary(model1))

cat("\n[FIT] Model 2: + Demographics\n")
model2 <- svyglm(CVD ~ IBI_Category + Age + Gender + Ethnicity + Education,
                 design = train_design, family = quasibinomial("logit"))
print(summary(model2))

cat("\n[FIT] Model 3: Fully adjusted\n")
model3 <- svyglm(CVD ~ IBI_Category + Age + Gender + Ethnicity + Education +
                   Smoking_status + Alcohol + BMI + Diabetes +
                   TOTAL_CHOLESTEROL + HDL_CHOLESTEROL + Hypertension,
                 design = train_design, family = quasibinomial("logit"))
print(summary(model3))


# ===============================
# 6) Predict on Test data (2021–2023)
# ===============================
test_df$pred1 <- predict(model1, newdata = test_df, type = "response")
test_df$pred2 <- predict(model2, newdata = test_df, type = "response")
test_df$pred3 <- predict(model3, newdata = test_df, type = "response")

# Some Example Predictions
cat("\n[INFO] Example predictions (first 5 rows):\n")
print(head(test_df[, c("CVD","pred1","pred2","pred3")], 5))

# Get Weighted Mean
cat("\n[TEST prevalence]\n")
prev <- with(test_df, stats::weighted.mean(CVD_num, WTMEC2YR, na.rm = TRUE))
cat(" - Weighted prevalence of CVD (Yes):", round(prev, 3), "\n")

# Weighted baseline accuracy, i.e. majority class
baseline_w_acc <- max(prev, 1 - prev)
cat(" - Weighted baseline accuracy (majority class):", round(baseline_w_acc, 3), "\n")



# Function to get Weighted metrics at a given threshold
  # Takes in predictions, target, weights, threshold
weighted_metrics <- function(pred, y, w, thr) {
  pred_cls <- ifelse(pred >= thr, 1L, 0L)
  keep <- is.finite(pred_cls) & is.finite(y) & is.finite(w)
  pred_cls <- pred_cls[keep]; y <- y[keep]; w <- as.numeric(w[keep])
  
  TPw <- sum(w * (pred_cls == 1L) * (y == 1L), na.rm = TRUE)
  FPw <- sum(w * (pred_cls == 1L) * (y == 0L), na.rm = TRUE)
  TNw <- sum(w * (pred_cls == 0L) * (y == 0L), na.rm = TRUE)
  FNw <- sum(w * (pred_cls == 0L) * (y == 1L), na.rm = TRUE)
  Wtot <- TPw + FPw + TNw + FNw
  
  acc  <- (TPw + TNw) / Wtot
  sens <- ifelse((TPw + FNw) > 0, TPw / (TPw + FNw), NA_real_)
  spec <- ifelse((TNw + FPw) > 0, TNw / (TNw + FPw), NA_real_)
  ppv  <- ifelse((TPw + FPw) > 0, TPw / (TPw + FPw), NA_real_)
  npv  <- ifelse((TNw + FNw) > 0, TNw / (TNw + FNw), NA_real_)
  
  lr_pos_raw <- ifelse((1 - spec) > 0, sens / (1 - spec), NA_real_)  
  lr_neg_raw <- ifelse(spec > 0, (1 - sens) / spec, NA_real_) 
  
  f2 <- ifelse(is.na(ppv) | is.na(sens) | (4*ppv + sens) == 0,
               NA_real_, (5 * ppv * sens) / (4 * ppv + sens))
  
  tibble::tibble(
    Threshold            = thr,
    Weighted_Accuracy    = round(acc, 3),
    Weighted_Recall      = round(sens, 3),
    Weighted_Specificity = round(spec, 3),
    Weighted_Precision   = round(ppv, 3),
    Weighted_NPV         = round(npv, 3),
    LR_Positive          = round(lr_pos_raw, 3),
    LR_Negative          = round(lr_neg_raw, 3),
    Weighted_F2          = round(f2, 3)
  )
}


metrics1_test <- with(test_df, weighted_metrics(pred1, CVD_num, WTMEC2YR, 0.5))
metrics2_test <- with(test_df, weighted_metrics(pred2, CVD_num, WTMEC2YR, 0.5))
metrics3_test <- with(test_df, weighted_metrics(pred3, CVD_num, WTMEC2YR, 0.5))

metrics_05 <- dplyr::bind_rows(metrics1_test, metrics2_test, metrics3_test)
cat("\n[WEIGHTED METRICS @ 0.5]\n"); print(metrics_05)

# ===============================
# 7) ROC/AUC/Metrics (Test data) – weighted
# ===============================
roc1_test <- pROC::roc(test_df$CVD, test_df$pred1, levels = c("No","Yes"),
                       weights = test_df$WTMEC2YR, quiet = TRUE)
roc2_test <- pROC::roc(test_df$CVD, test_df$pred2, levels = c("No","Yes"),
                       weights = test_df$WTMEC2YR, quiet = TRUE)
roc3_test <- pROC::roc(test_df$CVD, test_df$pred3, levels = c("No","Yes"),
                       weights = test_df$WTMEC2YR, quiet = TRUE)

auc1_test <- pROC::auc(roc1_test)
auc2_test <- pROC::auc(roc2_test)
auc3_test <- pROC::auc(roc3_test)

cat("\n[AUC on TEST]\n")
cat(" - Model 1 (IBI only):      ", round(auc1_test, 3), "\n")
cat(" - Model 2 (Demographics):  ", round(auc2_test, 3), "\n")
cat(" - Model 3 (Full):          ", round(auc3_test, 3), "\n")


# ===============================
# 7b) Youden thresholds on TEST (weighted)
# ===============================

youden_info <- function(roc_obj) {
  # Returns threshold, sensitivity, specificity, and J = sens + spec - 1
  cs <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  thr  <- as.numeric(cs[1, "threshold"])
  sens <- as.numeric(cs[1, "sensitivity"])
  spec <- as.numeric(cs[1, "specificity"])
  data.frame(
    Threshold   = thr,
    Sensitivity = sens,
    Specificity = spec,
    J           = sens + spec - 1
  )
}

youden1 <- youden_info(roc1_test)
youden2 <- youden_info(roc2_test)
youden3 <- youden_info(roc3_test)

youden_tbl <- dplyr::bind_rows(
  dplyr::mutate(youden1, Model = "Model 1", .before = 1),
  dplyr::mutate(youden2, Model = "Model 2", .before = 1),
  dplyr::mutate(youden3, Model = "Model 3", .before = 1)
) |>
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4)))

cat("\n[YOUDEN THRESHOLDS on TEST]\n"); print(youden_tbl)

# Evaluate weighted metrics at each model's Youden cutoff
metrics1_youden <- with(test_df, weighted_metrics(pred1, CVD_num, WTMEC2YR, thr = youden1$Threshold))
metrics2_youden <- with(test_df, weighted_metrics(pred2, CVD_num, WTMEC2YR, thr = youden2$Threshold))
metrics3_youden <- with(test_df, weighted_metrics(pred3, CVD_num, WTMEC2YR, thr = youden3$Threshold))

metrics_youden <- dplyr::bind_rows(
  dplyr::mutate(metrics1_youden, Model = "Model 1", .before = 1),
  dplyr::mutate(metrics2_youden, Model = "Model 2", .before = 1),
  dplyr::mutate(metrics3_youden, Model = "Model 3", .before = 1)
)

cat("\n[WEIGHTED METRICS @ YOUDEN THRESHOLD]\n"); print(metrics_youden)

# Save both tables
readr::write_csv(youden_tbl, "youden_thresholds_test.csv")
readr::write_csv(metrics_youden, "metrics_at_youden_threshold_test.csv")
cat("[FILES SAVED] youden_thresholds_test.csv, metrics_at_youden_threshold_test.csv\n")

# ===============================
# 7c) Model 3 cutoff calibration: Youden, Rule-in, Rule-out, Closest.topleft
# ===============================

# Helper - ROC coords
roc_all <- pROC::coords(
  roc3_test,
  x = "all",
  ret = c("threshold", "sensitivity", "specificity"),
  transpose = FALSE
)
roc_all <- dplyr::mutate(
  as.data.frame(roc_all),
  J = sensitivity + specificity - 1,
  dist_topleft = sqrt((1 - sensitivity)^2 + (1 - specificity)^2) # distance to (0,1)
)

# Get a row for each index 
# 1) Youden J (maximize J)
youden_row <- roc_all[which.max(roc_all$J), , drop = FALSE]

# 2) Rule-in, choose the max sensitivity with specificity >= 0.95
rulein_pool <- dplyr::filter(roc_all, specificity >= 0.95)
rulein_row <- if (nrow(rulein_pool)) {
  rulein_pool[which.max(rulein_pool$sensitivity), , drop = FALSE]
} else {
  # if none >= 0.95 take the row with maximum specificity
  roc_all[which.max(roc_all$specificity), , drop = FALSE]
}

# 3) Rule-out,  choose the with max specificity with sensitivity >= 0.95
ruleout_pool <- dplyr::filter(roc_all, sensitivity >= 0.95)
ruleout_row <- if (nrow(ruleout_pool)) {
  ruleout_pool[which.max(ruleout_pool$specificity), , drop = FALSE]
} else {
  # if none >= 0.95 take the row with maximum sensitivity
  roc_all[which.max(roc_all$sensitivity), , drop = FALSE]
}

# 4) Closest to (0,1)
closest_row <- roc_all[which.min(roc_all$dist_topleft), , drop = FALSE]

# table of all 4 thresholds
thr_tbl_m3 <- dplyr::bind_rows(
  dplyr::mutate(youden_row,  Index = "Youden J",        .before = 1),
  dplyr::mutate(rulein_row,  Index = "Rule-in (Spec≥0.95)", .before = 1),
  dplyr::mutate(ruleout_row, Index = "Rule-out (Sens≥0.95)", .before = 1),
  dplyr::mutate(closest_row, Index = "Closest to (0,1)", .before = 1)
) |>
  dplyr::transmute(
    Index,
    Threshold = as.numeric(threshold),
    Sensitivity = sensitivity,
    Specificity = specificity,
    J = J
  ) |>
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4)))

cat("\n[MODEL 3: Candidate thresholds]\n"); print(thr_tbl_m3)

# Evaluate weighted metrics at each threshold
m3_metrics_list <- lapply(thr_tbl_m3$Threshold, function(thr) {
  with(test_df, weighted_metrics(pred3, CVD_num, WTMEC2YR, thr = thr))
})
metrics_m3_cal <- dplyr::bind_rows(m3_metrics_list) |>
  dplyr::mutate(Index = thr_tbl_m3$Index, .before = 1)

cat("\n[MODEL 3: Weighted metrics at candidate thresholds]\n"); print(metrics_m3_cal)

# Save
readr::write_csv(thr_tbl_m3, "model3_thresholds_calibration.csv")
readr::write_csv(metrics_m3_cal, "model3_metrics_at_calibrated_thresholds.csv")
cat("[FILES SAVED] model3_thresholds_calibration.csv, model3_metrics_at_calibrated_thresholds.csv\n")


# Single view for thresholds + full metrics for Model 3
m3_summary <- thr_tbl_m3 |>
  dplyr::rename(Cutoff = Threshold) |>
  dplyr::left_join(
    dplyr::rename(metrics_m3_cal, Cutoff = Threshold),
    by = c("Index", "Cutoff")
  ) |>
  dplyr::select(
    Index, Cutoff,
    Weighted_Accuracy, Weighted_Recall, Weighted_Specificity,
    Weighted_Precision, Weighted_NPV, LR_Positive, LR_Negative, Weighted_F2
  ) |>
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))

cat("\n[MODEL 3: Calibrated thresholds + full metrics]\n"); print(m3_summary)
readr::write_csv(m3_summary, "model3_calibrated_thresholds_full_metrics.csv")
cat("[FILE SAVED] model3_calibrated_thresholds_full_metrics.csv\n")


# ===============================
# 8) Train-set ROC/AUC (unweighted)
# ===============================
train_df$pred1_tr <- predict(model1, type = "response")
train_df$pred2_tr <- predict(model2, type = "response")
train_df$pred3_tr <- predict(model3, type = "response")

roc1_train <- pROC::roc(train_df$CVD, train_df$pred1_tr, levels = c("No","Yes"), quiet = TRUE)
roc2_train <- pROC::roc(train_df$CVD, train_df$pred2_tr, levels = c("No","Yes"), quiet = TRUE)
roc3_train <- pROC::roc(train_df$CVD, train_df$pred3_tr, levels = c("No","Yes"), quiet = TRUE)

results_compare <- data.frame(
  Model = c("Model 1","Model 2","Model 3"),
  Train_AUC = c(auc(roc1_train), auc(roc2_train), auc(roc3_train)),
  Test_AUC  = c(auc1_test,       auc2_test,       auc3_test)
)

results_compare_print <- dplyr::mutate(results_compare,
                                       dplyr::across(where(is.numeric), ~ round(.x, 3)))
cat("\n[AUC TRAIN vs TEST]\n"); print(results_compare_print)
readr::write_csv(results_compare_print, "auc_train_vs_test.csv")
cat("[FILE SAVED] auc_train_vs_test.csv\n")

# ===============================
# 9) Plot Test ROC curves
# ===============================
# All 3 models 
plot(roc1_test, col = "red", main = "Test ROC (2021–2023) — models trained on 2017–2020",
     xlab = "1 - Specificity", ylab = "Sensitivity", lwd = 2)
lines(roc2_test, col = "blue", lwd = 2)
lines(roc3_test, col = "green", lwd = 2)
legend("bottomright",
       legend = c(paste("Model 1:", round(auc1_test, 3)),
                  paste("Model 2:", round(auc2_test, 3)),
                  paste("Model 3:", round(auc3_test, 3))),
       col = c("red","blue","green"), lty = 1, lwd = 2, cex = 0.85)

# Save PNG
png("ROC_Curves_TrainToTest.png", width = 800, height = 600)
plot(roc1_test, col = "red", main = "Test ROC (2021–2023) — models trained on 2017–2020",
     xlab = "1 - Specificity", ylab = "Sensitivity", lwd = 2)
lines(roc2_test, col = "blue", lwd = 2)
lines(roc3_test, col = "green", lwd = 2)
legend("bottomright",
       legend = c(paste("Model 1:", round(auc1_test, 3)),
                  paste("Model 2:", round(auc2_test, 3)),
                  paste("Model 3:", round(auc3_test, 3))),
       col = c("red","blue","green"), lty = 1, lwd = 2, cex = 0.85)
dev.off()
cat("[FILE SAVED] ROC_Curves_TrainToTest.png\n")



# Weighted precision–recall curve for Model 3
wpr_auc <- PRROC::pr.curve(
  scores.class0 = test_df$pred3[test_df$CVD_num == 0],
  scores.class1 = test_df$pred3[test_df$CVD_num == 1],
  weights.class0 = test_df$WTMEC2YR[test_df$CVD_num == 0],
  weights.class1 = test_df$WTMEC2YR[test_df$CVD_num == 1],
  curve = FALSE # took up too much memory
)$auc.integral
cat("Weighted PR-AUC (Model 3):", round(wpr_auc, 3), "\n")


# ===============================
# 10) Continuous IBI and model with splines
# ===============================

# Continuous IBI
model3_cont <- svyglm(CVD ~ IBI + Age + Gender + Ethnicity + Education +
                        Smoking_status + Alcohol + BMI + Diabetes +
                        TOTAL_CHOLESTEROL + HDL_CHOLESTEROL + Hypertension,
                      design = train_design, family = quasibinomial("logit"))

# Splines for IBI
model3_spline <- svyglm(CVD ~ ns(IBI, df = 3) + Age + Gender + Ethnicity + Education +
                          Smoking_status + Alcohol + BMI + Diabetes +
                          TOTAL_CHOLESTEROL + HDL_CHOLESTEROL + Hypertension,
                        design = train_design, family = quasibinomial("logit"))

# Compare TEST AUCs for these
test_df$pred3_cont   <- predict(model3_cont,   newdata=test_df, type="response")
test_df$pred3_spline <- predict(model3_spline, newdata=test_df, type="response")

roc3_cont   <- pROC::roc(test_df$CVD, test_df$pred3_cont,   levels=c("No","Yes"), weights=test_df$WTMEC2YR, quiet=TRUE)
roc3_spline <- pROC::roc(test_df$CVD, test_df$pred3_spline, levels=c("No","Yes"), weights=test_df$WTMEC2YR, quiet=TRUE)

cat("\n[Model 3 variants — AUC on TEST]\n")
cat("  IBI continuous: ", round(pROC::auc(roc3_cont),   3), "\n")
cat("  IBI spline:     ", round(pROC::auc(roc3_spline), 3), "\n")


# ===============================
# 10b) IBI as binary (top quartile vs others), fully adjusted model
# ===============================

# Create binary on train and test using existing quartiles
train_df$IBI_high <- as.integer(train_df$IBI >= ibi_quartiles[4])
test_df$IBI_high  <- as.integer(test_df$IBI  >= ibi_quartiles[4])

train_design <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR,
  nest = TRUE, survey.lonely.psu = "adjust",
  data = train_df
)

# Fit with IBI_high
model3_binary <- svyglm(
  CVD ~ IBI_high + Age + Gender + Ethnicity + Education +
    Smoking_status + Alcohol + BMI + Diabetes +
    TOTAL_CHOLESTEROL + HDL_CHOLESTEROL + Hypertension,
  design = train_design, family = quasibinomial("logit")
)
cat("\n[FIT] Model 3 (IBI binary: Q4 vs others)\n"); print(summary(model3_binary))


# Predict on TRAIN data
train_df$pred3_binary <- predict(model3_binary, type = "response")

roc3_binary_train <- pROC::roc(
  train_df$CVD, train_df$pred3_binary,
  levels = c("No","Yes"), quiet = TRUE
)
auc3_binary_train <- pROC::auc(roc3_binary_train)



# Predict & ROC/AUC on TEST
test_df$pred3_binary <- predict(model3_binary, newdata = test_df, type = "response")

roc3_binary <- pROC::roc(
  test_df$CVD, test_df$pred3_binary,
  levels = c("No","Yes"),
  weights = test_df$WTMEC2YR, quiet = TRUE
)
auc3_binary <- pROC::auc(roc3_binary)

cat("\n[Model 3 (binary IBI) — AUC on TEST]\n")
cat("  IBI binary (Q4 vs others): ", round(auc3_binary, 3), "\n")


# Compare continuous, spline, and binary models
results_compare2 <- results_compare |>
  dplyr::bind_rows(tibble::tibble(
    Model     = "Model 3 (IBI binary)",
    Train_AUC = as.numeric(auc3_binary_train),
    Test_AUC  = as.numeric(auc3_binary)
  ))


results_compare_print <- results_compare2 |>
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))

cat("\n[AUC TRAIN vs TEST] (updated with binary IBI)\n"); print(results_compare_print)
readr::write_csv(results_compare_print, "auc_train_vs_test.csv")

# ===============================
# 11)  Dose–response for IBI
# ===============================


lo <- as.numeric(quantile(train_df$IBI, 0.01, na.rm = TRUE))
hi <- as.numeric(quantile(train_df$IBI, 0.99, na.rm = TRUE))
step <- (hi - lo) / 200  # ~200 points

gg_ibi <- ggpredict(
  model3_cont,
  terms = sprintf("IBI [%.6f:%.6f by=%.6f]", lo, hi, step)
)

plot(gg_ibi) +
  theme_minimal() +
  labs(
    title = "Continuous IBI dose–response (trimmed to 1st–99th pct)",
    x = "IBI",
    y = "Predicted Probability of CVD"
  )


# ===============================
# 12) Save models & other objects
# ===============================
saveRDS(model1, "model1_ibi_only_trainP.rds")
saveRDS(model2, "model2_demographic_adjusted_trainP.rds")
saveRDS(model3, "model3_fully_adjusted_trainP.rds")
saveRDS(roc1_test, "roc_model1_testL.rds")
saveRDS(roc2_test, "roc_model2_testL.rds")
saveRDS(roc3_test, "roc_model3_testL.rds")
saveRDS(ibi_quartiles, "ibi_quartiles_trainP.rds")
saveRDS(model3_cont,   "model3_cont_trainP.rds")
saveRDS(model3_spline, "model3_spline_trainP.rds")
saveRDS(model3_binary, "model3_binary_trainP.rds")
saveRDS(roc3_cont,     "roc_model3_cont_testL.rds")
saveRDS(roc3_spline,   "roc_model3_spline_testL.rds")
saveRDS(roc3_binary,   "roc_model3_binary_testL.rds")
saveRDS(roc3_binary_train, "roc_model3_binary_trainP.rds")
readr::write_csv(
  dplyr::transmute(test_df, SEQN, CVD, CVD_num, pred1, pred2, pred3, WTMEC2YR),
  "test_predictions_L.csv"
)

cat("\n[FILES SAVED]\n",
    "- model1_ibi_only_trainP.rds\n",
    "- model2_demographic_adjusted_trainP.rds\n",
    "- model3_fully_adjusted_trainP.rds\n",
    "- roc_model1_testL.rds\n",
    "- roc_model2_testL.rds\n",
    "- roc_model3_testL.rds\n",
    "- ibi_quartiles_trainP.rds\n",
    "- test_predictions_L.csv\n",
    "- model3_cont_trainP.rds\n",
    "- model3_spline_trainP.rds\n",
    "- model3_binary_trainP.rds\n",
    "- roc_model3_cont_testL.rds\n",
    "- roc_model3_spline_testL.rds\n",
    "- roc_model3_binary_testL.rds\n",
    "- roc_model3_binary_trainP.rds"
    )



writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat("[FILE SAVED] sessionInfo.txt\n")

utils::savehistory("Rhistory_this_run.Rhistory")


# End logging 
sink()
