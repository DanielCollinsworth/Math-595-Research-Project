
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
# str(data.complete)


cat("\n[INFO] Data loaded:\n")
cat(" - Rows:", nrow(data.complete), "  Cols:", ncol(data.complete), "\n")


# ===============================
# 3) Split by cycle & build IBI quartiles on TRAIN
# ===============================
#train_df <- subset(data.complete, Cycle == "P")  # 2017–2020
#test_df  <- subset(data.complete, Cycle == "L")  # 2021–2023

# look at everything
train_df <- data.complete


# Quartiles on  IBI
ibi_quartiles <- quantile(train_df$IBI, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)

# function to cut into quartiles
cut_ibi <- function(x, breaks) {
  cut(x, breaks = breaks, labels = c("Q1","Q2","Q3","Q4"),
      include.lowest = TRUE, right = TRUE, ordered_result = TRUE)
}

train_df$IBI_Category <- cut_ibi(train_df$IBI, ibi_quartiles)
test_df$IBI_Category  <- cut_ibi(test_df$IBI,  ibi_quartiles)

# factor and revel quartiles
train_df$IBI_Category <- factor(as.character(train_df$IBI_Category),
                                levels = c("Q1","Q2","Q3","Q4"), ordered = FALSE)
train_df$IBI_Category <- relevel(train_df$IBI_Category, ref = "Q1")


test_df$IBI_Category <- factor(as.character(test_df$IBI_Category),
                               levels = c("Q1","Q2","Q3","Q4"), ordered = FALSE)
test_df$IBI_Category <- relevel(test_df$IBI_Category, ref = "Q1")


# Quartiles as Numeric 1-4 
train_df$IBI_QuartileNum <- as.numeric(train_df$IBI_Category)
test_df$IBI_QuartileNum  <- as.numeric(test_df$IBI_Category)


# Target encodings
train_df$CVD <- factor(train_df$CVD, levels = c("No","Yes"))
test_df$CVD  <- factor(test_df$CVD,  levels = c("No","Yes"))
train_df$CVD_num <- as.integer(train_df$CVD == "Yes")
test_df$CVD_num  <- as.integer(test_df$CVD  == "Yes")

# size and cutoffs
cat("\n[SECTION 3 SUMMARY]\n")
cat(" - Train rows (Cycle P):", nrow(train_df), "\n")
cat(" - Test  rows (Cycle L):", nrow(test_df),  "\n")
cat(" - Quartile cutpoints (IBI, raw on TRAIN):\n"); print(ibi_quartiles)


# ===============================
# 4) Build survey design on Training data
# ===============================


train_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights  = ~WTMEC2YR,
  nest  = TRUE,
  survey.lonely.psu = "adjust",
  data = train_df
)


# ===============================
# 5) build Models and find OPs and CIs
# ===============================
cat("\n[FIT] Model 1: IBI only\n")
model1 <- svyglm(CVD ~ IBI_Category, design = train_design, family = binomial("logit"))
print(summary(model1))

cat("\n[FIT] Model 2: + Demographics\n")
model2 <- svyglm(CVD ~ IBI_Category + Age + Gender + Ethnicity + Education,
                 design = train_design, family = binomial("logit"))
print(summary(model2))

cat("\n[FIT] Model 3: Fully adjusted\n")
model3 <- svyglm(CVD ~ IBI_Category + Age + Gender + Ethnicity + Education +
                   Smoking_status + Alcohol + BMI + Diabetes +
                   TOTAL_CHOLESTEROL + HDL_CHOLESTEROL + Hypertension ,
                 design = train_design, family = binomial("logit"))
print(summary(model3))



cat("\n[FIT] Model 3: Fully adjusted\n")
model3_Q_Num <- svyglm(CVD ~ IBI_QuartileNum + Age + Gender + Ethnicity + Education +
                         Smoking_status + Alcohol + BMI + Diabetes +
                         TOTAL_CHOLESTEROL + HDL_CHOLESTEROL + Hypertension,
                       design = train_design, family = binomial("logit"))
print(summary(model3_Q_Num))






