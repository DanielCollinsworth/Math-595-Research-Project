# ===============================
# 1. Load Libraries
# ===============================
library(survey)
library(pROC)
library(dplyr)
library(ggplot2)

# ===============================
# 2. Load Data
# ===============================

data.complete <- readRDS("data_complete_cleaned.rds")
data.design <- readRDS("data_design_cleaned.rds")

# ===============================
# 3. Create IBI Categories (Q1, Q2, Q3, Q4)
# ===============================

# Calculate quartiles for IBI
ibi_quartiles <- quantile(data.complete$IBI, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
print("IBI Quartiles:")
print(ibi_quartiles)

# Create IBI categories
data.complete$IBI_Category <- cut(data.complete$IBI,
                                  breaks = ibi_quartiles,
                                  labels = c("Q1", "Q2", "Q3", "Q4"),
                                  include.lowest = TRUE)

# Set Q1 as reference
data.complete$IBI_Category <- relevel(data.complete$IBI_Category, ref = "Q1")

# Update survey design with new IBI category
data.design <- svydesign(
  id                = ~SDMVPSU,
  strata            = ~SDMVSTRA,
  weights           = ~WTMEC2YR,
  nest              = TRUE,
  survey.lonely.psu = "adjust",
  data              = data.complete
)

# ===============================
# 4. Weighted Logistic Regression Models
# ===============================

# -------------------------------
# Model 1: Non-adjusted model (IBI only)
# -------------------------------
model1 <- svyglm(CVD ~ IBI_Category,
                 design = data.design,
                 family = quasibinomial(link = "logit"))

print("=== MODEL 1: Non-adjusted (IBI only) ===")
summary(model1)

# ROC for Model 1
roc_model1 <- roc(response = data.complete$CVD,
                   predictor = fitted(model1))
auc_model1 <- auc(roc_model1)
print(paste("Model 1 AUC:", round(auc_model1, 3)))

# -------------------------------
# Model 2: Adjusted for demographic variables
# -------------------------------
model2 <- svyglm(CVD ~ IBI_Category + Age + Gender + Ethnicity + Education,
                 design = data.design,
                 family = quasibinomial(link = "logit"))

print("=== MODEL 2: Adjusted for demographic variables ===")
summary(model2)

# ROC for Model 2
roc_model2 <- roc(response = data.complete$CVD,
                   predictor = fitted(model2))
auc_model2 <- auc(roc_model2)
print(paste("Model 2 AUC:", round(auc_model2, 3)))

# -------------------------------
# Model 3: Fully adjusted model
# -------------------------------
model3 <- svyglm(CVD ~ IBI_Category + Age + Gender + Ethnicity + Education + 
                 Smoking_status + Alcohol + BMI + Diabetes + 
                 TOTAL_CHOLESTEROL + HDL_CHOLESTEROL + Hypertension,
                 design = data.design,
                 family = quasibinomial(link = "logit"))

print("=== MODEL 3: Fully adjusted model ===")
summary(model3)

# ROC for Model 3
roc_model3 <- roc(response = data.complete$CVD,
                   predictor = fitted(model3))
auc_model3 <- auc(roc_model3)
print(paste("Model 3 AUC:", round(auc_model3, 3)))

# ===============================
# 5. Create Combined ROC Plot
# ===============================

# Create a comprehensive ROC plot
png("ROC_Curves_Comparison.png", width = 800, height = 600)
plot(roc_model1, col = "red", main = "ROC Curves Comparison", 
     xlab = "1 - Specificity", ylab = "Sensitivity")
lines(roc_model2, col = "blue")
lines(roc_model3, col = "green")

# Add legend
legend("bottomright", 
       legend = c(paste("Model 1 (AUC =", round(auc_model1, 3), ")"),
                  paste("Model 2 (AUC =", round(auc_model2, 3), ")"),
                  paste("Model 3 (AUC =", round(auc_model3, 3), ")")),
       col = c("red", "blue", "green"),
       lty = 1, cex = 0.8)

dev.off()

# ===============================
# 6. Generate Results Summary
# ===============================

# Create results summary
results_summary <- data.frame(
  Model = c("Model 1", "Model 2", "Model 3"),
  Description = c("Non-adjusted (IBI only)", 
                  "Adjusted for demographic variables", 
                  "Fully adjusted model"),
  AUC = c(round(auc_model1, 3), round(auc_model2, 3), round(auc_model3, 3)),
  Covariates = c("IBI Category only",
                 "IBI + Age + Gender + Ethnicity + Education",
                 "IBI + Demographics + Behavioral + Clinical factors")
)

print("=== RESULTS SUMMARY ===")
print(results_summary)

# Save results to file
write.csv(results_summary, "logistic_regression_results.csv", row.names = FALSE)

# ===============================
# 7. Save Models
# ===============================
saveRDS(model1, "model1_ibi_only.rds")
saveRDS(model2, "model2_demographic_adjusted.rds")
saveRDS(model3, "model3_fully_adjusted.rds")

# Save ROC objects
saveRDS(roc_model1, "roc_model1.rds")
saveRDS(roc_model2, "roc_model2.rds")
saveRDS(roc_model3, "roc_model3.rds")

print("=== ANALYSIS COMPLETE ===")
print("Files saved:")
print("- ROC_Curves_Comparison.png")
print("- logistic_regression_results.csv")
print("- model1_ibi_only.rds")
print("- model2_demographic_adjusted.rds")
print("- model3_fully_adjusted.rds")
print("- roc_model1.rds")
print("- roc_model2.rds")
print("- roc_model3.rds")
  
