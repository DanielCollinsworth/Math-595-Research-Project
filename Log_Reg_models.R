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

# Example Predictions
cat("\n[INFO] Example predictions (first 5 rows):\n")
print(head(test_df[, c("CVD","pred1","pred2","pred3")], 5))

# Get Weighted Mean
cat("\n[TEST prevalence]\n")
prev <- with(test_df, stats::weighted.mean(CVD_num, WTMEC2YR, na.rm = TRUE))
cat(" - Weighted prevalence of CVD (Yes):", round(prev, 3), "\n")

# Weighted baseline accuracy -> majority class
baseline_w_acc <- max(prev, 1 - prev)
cat(" - Weighted baseline accuracy (majority class):", round(baseline_w_acc, 3), "\n")



# Function to get Weighted metrics at a given threshold
  # Takes in predictions, target, weights, threshold
weighted_metrics <- function(pred, y, w, thr) {
  pred_cls <- ifelse(pred >= thr, 1L, 0L) # as int
  keep <- is.finite(pred_cls) & is.finite(y) & is.finite(w)
  pred_cls <- pred_cls[keep]; y <- y[keep]; w <- as.numeric(w[keep])
  
  TPw <- sum(w * (pred_cls == 1L) * (y == 1L), na.rm = TRUE)
  FPw <- sum(w * (pred_cls == 1L) * (y == 0L), na.rm = TRUE)
  TNw <- sum(w * (pred_cls == 0L) * (y == 0L), na.rm = TRUE)
  FNw <- sum(w * (pred_cls == 0L) * (y == 1L), na.rm = TRUE)
  Wtot <- TPw + FPw + TNw + FNw
  
  acc  <- (TPw + TNw) / Wtot
  rec  <- ifelse((TPw + FNw) > 0, TPw / (TPw + FNw), NA_real_)  # sensitivity
  spec <- ifelse((TNw + FPw) > 0, TNw / (TNw + FPw), NA_real_)
  prec <- ifelse((TPw + FPw) > 0, TPw / (TPw + FPw), NA_real_)
  f2   <- ifelse(is.na(prec) | is.na(rec) | (4*prec + rec) == 0,
                 NA_real_, (5 * prec * rec) / (4 * prec + rec))
  
  tibble::tibble(
    Threshold = thr,
    Weighted_Accuracy    = round(acc, 3),
    Weighted_Recall      = round(rec, 3),
    Weighted_Specificity = round(spec, 3),
    Weighted_Precision   = round(prec, 3),
    Weighted_F2          = round(f2, 3)
  )
}


metrics1_test <- with(test_df, weighted_metrics(pred1, CVD_num, WTMEC2YR, 0.5))
metrics2_test <- with(test_df, weighted_metrics(pred2, CVD_num, WTMEC2YR, 0.5))
metrics3_test <- with(test_df, weighted_metrics(pred3, CVD_num, WTMEC2YR, 0.5))

metrics_5 <- dplyr::bind_rows(metrics1_test, metrics2_test, metrics3_test)
cat("\n[WEIGHTED METRICS @ 0.5]\n"); print(metrics_5)


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
# 8) Train-set ROC/AUC (unweighted)
# ===============================
train_df$pred1_tr <- predict(model1, type = "response")
train_df$pred2_tr <- predict(model2, type = "response")
train_df$pred3_tr <- predict(model3, type = "response")

roc1_train <- roc(train_df$CVD, train_df$pred1_tr, levels = c("No","Yes"))
roc2_train <- roc(train_df$CVD, train_df$pred2_tr, levels = c("No","Yes"))
roc3_train <- roc(train_df$CVD, train_df$pred3_tr, levels = c("No","Yes"))

results_compare <- data.frame(
  Model = c("Model 1","Model 2","Model 3"),
  Train_AUC = c(auc(roc1_train), auc(roc2_train), auc(roc3_train)),
  Test_AUC  = c(auc1_test,       auc2_test,       auc3_test)
)

results_compare_print <- dplyr::mutate(results_compare,
                                       dplyr::across(where(is.numeric), ~ round(.x, 3)))
cat("\n[AUC TRAIN vs TEST]\n"); print(results_compare_print)
write.csv(results_compare_print, "auc_train_vs_test.csv", row.names = FALSE)
cat("[FILE SAVED] auc_train_vs_test.csv\n")

# ===============================
# 9) Plot Test ROC curves
# ===============================
# All 3models 
plot(roc1_test, col = "red", main = "ROC: Train (2017–2020) → Test (2021–2023)",
     xlab = "1 - Specificity", ylab = "Sensitivity", lwd = 2)
lines(roc2_test, col = "blue", lwd = 2)
lines(roc3_test, col = "green", lwd = 2)
legend("bottomright",
       legend = c(paste("Model 1:", round(auc1_test, 3)),
                  paste("Model 2:", round(auc2_test, 3)),
                  paste("Model 3:", round(auc3_test, 3))),
       col = c("red","blue","green"), lty = 1, lwd = 2, cex = 0.85)
cat("\n[INFO] ROC plot displayed in Plots panel\n")

# Save PNG
png("ROC_Curves_TrainToTest.png", width = 800, height = 600)
plot(roc1_test, col = "red", main = "ROC: Train (2017–2020) → Test (2021–2023)",
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
# 10) Continuous IBI and model  with splines
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

# Compare TEST AUCs for these specs
test_df$pred3_cont   <- predict(model3_cont,   newdata=test_df, type="response")
test_df$pred3_spline <- predict(model3_spline, newdata=test_df, type="response")

roc3_cont   <- pROC::roc(test_df$CVD, test_df$pred3_cont,   levels=c("No","Yes"), weights=test_df$WTMEC2YR, quiet=TRUE)
roc3_spline <- pROC::roc(test_df$CVD, test_df$pred3_spline, levels=c("No","Yes"), weights=test_df$WTMEC2YR, quiet=TRUE)

cat("\n[Model 3 variants — AUC on TEST]\n")
cat("  IBI continuous: ", round(pROC::auc(roc3_cont),   3), "\n")
cat("  IBI spline:     ", round(pROC::auc(roc3_spline), 3), "\n")



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
    "- test_predictions_L.csv\n"
    )



writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat("[FILE SAVED] sessionInfo.txt\n")

utils::savehistory("Rhistory_this_run.Rhistory")


# End logging 
sink()
