# ===============================
# 1. Load Libraries
# ===============================
library(survey)
library(pROC)
library(dplyr)
library(recipes)

# ===============================
# 2. Load Data
# ===============================

data.complete <- readRDS("data_complete_cleaned.rds")
data.design <- readRDS("data_design_cleaned.rds")


# ===============================
# 3. Make models adn plot ROC Curves
# ===============================

# -------------------------------
# 3a. Plain logistic regression
basic_glm <- glm(CVD ~ IBI,
               data = data.complete,
               family = binomial(link = "logit"))
summary(basic_glm)

  
  # make ROC for true labels vs predicted
  roc_obj <- roc(response = data.complete$CVD,
                 predictor = fitted(basic_glm))
  plot(roc_obj, col = "blue", main = "ROC Curve")
  auc(roc_obj)   # Area under the curve
  

  
# -------------------------------
# 3b. Logistic regression fitted to more fields
fit_glm <- glm(CVD ~ IBI + Age + Gender + BMI,
               data = data.complete, # use full data
               family = binomial(link = "logit"))
#summary(fit_glm)

  
  # make ROC for true labels vs predicted
  roc_fit_obj <- roc(response = data.complete$CVD,
                 predictor = fitted(fit_glm))
  plot(roc_fit_obj, col = "blue", main = "ROC Curve")
  auc(roc_fit_obj)


# -------------------------------
# 3c. Logistic regression fitted to survey

# Survey-weighted logistic regression
fit_svy <- svyglm(CVD ~ IBI + Age + Gender + BMI,
                  design = data.design, # use survey data
                  family = quasibinomial(link = "logit"))
#summary(fit_svy)

  # make ROC for true labels vs predicted
  roc_svy_obj <- roc(response = data.complete$CVD,
                 predictor = fitted(fit_svy))
  plot(roc_svy_obj, col = "blue", main = "ROC Curve")
  auc(roc_svy_obj)
  
  
# ===============================
# 4. Save models
# ===============================
saveRDS(basic_glm, "basic_glm.rds")
saveRDS(fit_glm, "fit_glm.rds")
saveRDS(fit_svy, "fit_svy.rds")
  
