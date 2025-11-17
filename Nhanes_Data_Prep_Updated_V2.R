setwd('C:/Users/danie/repos/595 Research/Code')
# ===============================
# 1. Load Required Libraries
# ===============================
if (!require("nhanesA")) install.packages("nhanesA")
library(nhanesA)
library(plyr)
library(dplyr)
library(tidyr)
library(knitr)
library(purrr)
library(sjPlot)
library(jstable)
library(broom)
library(arsenal)

# ===============================
# 2. Helper Function
# ===============================
get_nhanes_data <- function(nhanes_table, vars, translate = TRUE) {
  data <- nhanes(nhanes_table)[, vars, drop = FALSE]
  if (translate) {
    data <- nhanesTranslate(
      nh_table = nhanes_table,
      names(data),
      data = data
    )
  }
  data
}

# ===============================
# 3. Load Cycle L Data
# ===============================
cycle_L_demo <- get_nhanes_data(
  "DEMO_L",
  c("SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH3", "DMDEDUC2",
    "INDFMPIR", "WTINT2YR", "SDMVSTRA", "WTMEC2YR", "SDMVPSU")
) %>%
  rename(
    Gender        = RIAGENDR,
    Age           = RIDAGEYR,
    Ethnicity_raw = RIDRETH3,
    Education_raw = DMDEDUC2,
    Poverty_Ratio = INDFMPIR,
    SAMPLEWEIGHT  = WTINT2YR
  ) %>%
  mutate(
    # Ethnicity: American / Asian/Other / Hispanic
    Ethnicity = case_when(
      Ethnicity_raw %in% c("Non-Hispanic White") ~ "American",
      Ethnicity_raw %in% c("Mexican American", "Other Hispanic") ~ "Hispanic",
      Ethnicity_raw %in% c("Non-Hispanic Asian",
                           "Non-Hispanic Black",
                           "Other Race - Including Multi-Racial") ~ "Asian/Other",
      TRUE ~ NA_character_
    ),
    Ethnicity = factor(Ethnicity,
                       levels = c("American","Asian/Other","Hispanic")),
    # Education: Less than HS / HS+GED / Some college or higher
    Education = case_when(
      Education_raw %in% c("Less than 9th grade",
                           "9-11th grade (Includes 12th grade with no diploma)") ~
        "Less than High School",
      Education_raw %in% c("High school graduate/GED or equivalent") ~
        "HighSchoolGraduate/GED",
      Education_raw %in% c("Some college or AA degree",
                           "College graduate or above") ~
        "SomeCollegeorHigher",
      TRUE ~ NA_character_
    ),
    Education = factor(
      Education,
      levels = c("Less than High School",
                 "HighSchoolGraduate/GED",
                 "SomeCollegeorHigher")
    ),
    Gender = factor(Gender, levels = c("Male","Female"))
  )

cycle_L_BMI <- get_nhanes_data("BMX_L", c("SEQN", "BMXBMI", "BMXWAIST")) %>%
  rename(BMI = BMXBMI, Waist_Circ = BMXWAIST)

cycle_L_BP <- nhanes("BPXO_L")[, c(
  "SEQN",
  "BPXOSY1","BPXOSY2","BPXOSY3",
  "BPXODI1","BPXODI2","BPXODI3"
)] %>%
  mutate(
    SYSTOLIC_BP  = rowMeans(select(., BPXOSY1, BPXOSY2, BPXOSY3), na.rm = TRUE),
    DIASTOLIC_BP = rowMeans(select(., BPXODI1, BPXODI2, BPXODI3), na.rm = TRUE)
  ) %>%
  select(SEQN, SYSTOLIC_BP, DIASTOLIC_BP)

# Smoking: 3-level, binary derived later
cycle_L_Smoke <- get_nhanes_data("SMQ_L", c("SEQN", "SMQ020", "SMQ040")) %>%
  mutate(
    Smoking_status = case_when(
      SMQ020 != "Yes" ~ "Never Smoker",
      SMQ020 == "Yes" & SMQ040 %in% c("Every day", "Some days") ~ "Current Smoker",
      SMQ020 == "Yes" & SMQ040 == "Not at all" ~ "Former Smoker",
      TRUE ~ NA_character_
    )
  ) %>%
  select(SEQN, Smoking_status)

# Alcohol: detailed categories; binary later
cycle_L_Alch <- get_nhanes_data("ALQ_L", c("SEQN", "ALQ121")) %>%
  rename(Alcohol_freq = ALQ121) %>%
  filter(!Alcohol_freq %in% c("Refused","Don't know")) %>%
  mutate(
    Alcohol = case_when(
      Alcohol_freq %in% c("Every day", "Nearly every day",
                          "3 to 4 times a week", "2 times a week",
                          "Once a week") ~ "Frequent Drinker",
      Alcohol_freq %in% c("2 to 3 times a month", "Once a month",
                          "7 to 11 times in the last year",
                          "3 to 6 times in the last year",
                          "1 to 2 times in the last year") ~ "Occasional Drinker",
      Alcohol_freq %in% c("Never in the last year") ~ "Non-Drinker",
      TRUE ~ NA_character_
    )
  )

cycle_L_Diabetes <- get_nhanes_data("DIQ_L", c("SEQN", "DIQ010")) %>%
  rename(Diabetes = DIQ010) %>%
  filter(Diabetes != "Don't know")

cycle_L_CProt <- get_nhanes_data("HSCRP_L", c("SEQN", "LBXHSCRP")) %>%
  rename(C_REATIVE_PROTEIN = LBXHSCRP)

cycle_L_TChol <- get_nhanes_data("TCHOL_L", c("SEQN", "LBXTC")) %>%
  rename(TOTAL_CHOLESTEROL = LBXTC)

cycle_L_HDL <- get_nhanes_data("HDL_L", c("SEQN", "LBDHDD")) %>%
  rename(HDL_CHOLESTEROL = LBDHDD)

cycle_L_CVD <- get_nhanes_data(
  "MCQ_L",
  c("SEQN", "MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160F")
) %>%
  mutate(across(c(MCQ160B, MCQ160C, MCQ160D, MCQ160E, MCQ160F), as.character)) %>%
  mutate(across(c(MCQ160B, MCQ160C, MCQ160D, MCQ160E, MCQ160F),
                ~ ifelse(. == "Don't know", NA, .))) %>%
  mutate(
    CVD = ifelse(
      MCQ160B == "Yes" | MCQ160C == "Yes" |
        MCQ160D == "Yes" | MCQ160E == "Yes" |
        MCQ160F == "Yes", "Yes", "No"
    )
  ) %>%
  rename(
    HEART_FAILURE          = MCQ160B,
    CORONARY_HEART_DISEASE = MCQ160C,
    ANGINA                 = MCQ160D,
    HEART_ATTACK           = MCQ160E,
    STROKE                 = MCQ160F
  )

cycle_L_Depres <- get_nhanes_data("DPQ_L", c(
  "SEQN","DPQ010","DPQ020","DPQ030","DPQ040","DPQ050",
  "DPQ060","DPQ070","DPQ080","DPQ090"
)) %>%
  mutate(across(starts_with("DPQ0"), ~ ifelse(. %in% c(7,9), NA, .))) %>%
  mutate(
    PHQ9_Score = rowSums(select(., starts_with("DPQ0")), na.rm = TRUE),
    Depression = case_when(
      PHQ9_Score == 0 ~ NA_character_,
      PHQ9_Score >= 10 ~ "Yes",
      TRUE ~ "No"
    )
  ) %>%
  select(SEQN, Depression)

cycle_L_CBC <- get_nhanes_data("CBC_L",
                               c("SEQN","LBDLYMNO","LBDNENO","LBDMONO","LBXPLTSI")) %>%
  rename(
    Lymphocyte = LBDLYMNO,
    Neutrophil = LBDNENO,
    Monocyte   = LBDMONO,
    Platelet   = LBXPLTSI
  )

cycle_L_list <- list(
  cycle_L_demo, cycle_L_BMI, cycle_L_BP, cycle_L_Smoke, cycle_L_Alch,
  cycle_L_Diabetes, cycle_L_CProt, cycle_L_TChol, cycle_L_HDL,
  cycle_L_CVD, cycle_L_Depres, cycle_L_CBC
)

cycle_L_data <- reduce(cycle_L_list, full_join, by = "SEQN") %>%
  mutate(Cycle = "L")

# ===============================
# 4. Load Cycle P Data
# ===============================

cycle_P_demo <- get_nhanes_data(
  "P_DEMO",
  c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH3","DMDEDUC2",
    "INDFMPIR","WTINTPRP","WTMECPRP","SDMVPSU","SDMVSTRA")
) %>%
  rename(
    Gender        = RIAGENDR,
    Age           = RIDAGEYR,
    Ethnicity_raw = RIDRETH3,
    Education_raw = DMDEDUC2,
    Poverty_Ratio = INDFMPIR,
    SAMPLEWEIGHT  = WTINTPRP,
    WTMEC2YR      = WTMECPRP
  ) %>%
  mutate(
    Ethnicity = case_when(
      Ethnicity_raw %in% c("Non-Hispanic White") ~ "American",
      Ethnicity_raw %in% c("Mexican American", "Other Hispanic") ~ "Hispanic",
      Ethnicity_raw %in% c("Non-Hispanic Asian",
                           "Non-Hispanic Black",
                           "Other Race - Including Multi-Racial") ~ "Asian/Other",
      TRUE ~ NA_character_
    ),
    Ethnicity = factor(Ethnicity,
                       levels = c("American","Asian/Other","Hispanic")),
    Education = case_when(
      Education_raw %in% c("Less than 9th grade",
                           "9-11th grade (Includes 12th grade with no diploma)") ~
        "Less than High School",
      Education_raw %in% c("High school graduate/GED or equivalent") ~
        "HighSchoolGraduate/GED",
      Education_raw %in% c("Some college or AA degree",
                           "College graduate or above") ~
        "SomeCollegeorHigher",
      TRUE ~ NA_character_
    ),
    Education = factor(
      Education,
      levels = c("Less than High School",
                 "HighSchoolGraduate/GED",
                 "SomeCollegeorHigher")
    ),
    Gender = factor(Gender, levels = c("Male","Female"))
  )

# BMI & Waist
cycle_P_BMI <- get_nhanes_data("P_BMX", c("SEQN", "BMXBMI", "BMXWAIST")) %>%
  rename(BMI = BMXBMI, Waist_Circ = BMXWAIST)

# Blood pressure
cycle_P_BP <- nhanes("P_BPXO")[, c(
  "SEQN",
  "BPXOSY1","BPXOSY2","BPXOSY3",
  "BPXODI1","BPXODI2","BPXODI3"
)] %>%
  mutate(
    SYSTOLIC_BP  = rowMeans(select(., BPXOSY1, BPXOSY2, BPXOSY3), na.rm = TRUE),
    DIASTOLIC_BP = rowMeans(select(., BPXODI1, BPXODI2, BPXODI3), na.rm = TRUE)
  ) %>%
  select(SEQN, SYSTOLIC_BP, DIASTOLIC_BP)

# Smoking
cycle_P_Smoke <- get_nhanes_data("P_SMQ", c("SEQN", "SMQ020", "SMQ040")) %>%
  mutate(
    Smoking_status = case_when(
      SMQ020 != "Yes" ~ "Never Smoker",
      SMQ020 == "Yes" & SMQ040 %in% c("Every day","Some days") ~ "Current Smoker",
      SMQ020 == "Yes" & SMQ040 == "Not at all" ~ "Former Smoker",
      TRUE ~ NA_character_
    )
  ) %>%
  select(SEQN, Smoking_status)

# Alcohol
cycle_P_Alch <- get_nhanes_data("P_ALQ", c("SEQN", "ALQ121")) %>%
  rename(Alcohol_freq = ALQ121) %>%
  filter(!Alcohol_freq %in% c("Refused","Don't know")) %>%
  mutate(
    Alcohol = case_when(
      Alcohol_freq %in% c("Every day","Nearly every day",
                          "3 to 4 times a week","2 times a week","Once a week") ~
        "Frequent Drinker",
      Alcohol_freq %in% c("2 to 3 times a month","Once a month",
                          "7 to 11 times in the last year",
                          "3 to 6 times in the last year",
                          "1 to 2 times in the last year") ~
        "Occasional Drinker",
      Alcohol_freq %in% c("Never in the last year") ~ "Non-Drinker",
      TRUE ~ NA_character_
    )
  )

# Diabetes
cycle_P_Diabetes <- get_nhanes_data("P_DIQ", c("SEQN", "DIQ010")) %>%
  rename(Diabetes = DIQ010) %>%
  filter(Diabetes != "Don't know")

# CRP
cycle_P_CProt <- get_nhanes_data("P_HSCRP", c("SEQN","LBXHSCRP")) %>%
  rename(C_REATIVE_PROTEIN = LBXHSCRP)

# Total cholesterol
cycle_P_TChol <- get_nhanes_data("P_TCHOL", c("SEQN","LBXTC")) %>%
  rename(TOTAL_CHOLESTEROL = LBXTC)

# HDL cholesterol
cycle_P_HDL <- get_nhanes_data("P_HDL", c("SEQN","LBDHDD")) %>%
  rename(HDL_CHOLESTEROL = LBDHDD)

# CVD
cycle_P_CVD <- get_nhanes_data(
  "P_MCQ",
  c("SEQN","MCQ160B","MCQ160C","MCQ160D","MCQ160E","MCQ160F")
) %>%
  mutate(across(c(MCQ160B,MCQ160C,MCQ160D,MCQ160E,MCQ160F), as.character)) %>%
  mutate(across(c(MCQ160B,MCQ160C,MCQ160D,MCQ160E,MCQ160F),
                ~ ifelse(. == "Don't know", NA, .))) %>%
  mutate(
    CVD = ifelse(
      MCQ160B == "Yes" | MCQ160C == "Yes" |
        MCQ160D == "Yes" | MCQ160E == "Yes" |
        MCQ160F == "Yes", "Yes", "No"
    )
  ) %>%
  rename(
    HEART_FAILURE          = MCQ160B,
    CORONARY_HEART_DISEASE = MCQ160C,
    ANGINA                 = MCQ160D,
    HEART_ATTACK           = MCQ160E,
    STROKE                 = MCQ160F
  )

# Depression
cycle_P_Depres <- get_nhanes_data(
  "P_DPQ",
  c("SEQN","DPQ010","DPQ020","DPQ030","DPQ040","DPQ050",
    "DPQ060","DPQ070","DPQ080","DPQ090")
) %>%
  mutate(across(starts_with("DPQ0"), ~ ifelse(. %in% c(7,9), NA, .))) %>%
  mutate(
    PHQ9_Score = rowSums(select(., starts_with("DPQ0")), na.rm = TRUE),
    Depression = case_when(
      PHQ9_Score == 0 ~ NA_character_,
      PHQ9_Score >= 10 ~ "Yes",
      TRUE ~ "No"
    )
  ) %>%
  select(SEQN, Depression)

# CBC
cycle_P_CBC <- get_nhanes_data(
  "P_CBC",
  c("SEQN","LBDLYMNO","LBDNENO","LBDMONO","LBXPLTSI")
) %>%
  rename(
    Lymphocyte = LBDLYMNO,
    Neutrophil = LBDNENO,
    Monocyte   = LBDMONO,
    Platelet   = LBXPLTSI
  )

cycle_P_list <- list(
  cycle_P_demo, cycle_P_BMI, cycle_P_BP, cycle_P_Smoke, cycle_P_Alch,
  cycle_P_Diabetes, cycle_P_CProt, cycle_P_TChol, cycle_P_HDL,
  cycle_P_CVD, cycle_P_Depres, cycle_P_CBC
)

cycle_P_data <- reduce(cycle_P_list, full_join, by = "SEQN") %>%
  mutate(Cycle = "P")

# ===============================
# 5. Combine cycles
# ===============================
nhanes_combined <- bind_rows(cycle_L_data, cycle_P_data) %>%
  filter(!is.na(SEQN))

cat("\n[STEP 0] Initial combined n =", nrow(nhanes_combined), "\n")

# ===============================
# 6. Remove missing IBI inputs & compute IBI
# ===============================
need_ibi <- c("C_REATIVE_PROTEIN","Neutrophil","Lymphocyte")

df1 <- nhanes_combined %>%
  drop_na(any_of(need_ibi)) %>%
  mutate(IBI = C_REATIVE_PROTEIN * (Neutrophil / Lymphocyte))

cat("[STEP 1] After removing missing IBI inputs n =", nrow(df1), "\n")
cat("  (Removed =", nrow(nhanes_combined) - nrow(df1), ")\n")


ibi_quartiles <- df1 %>%
  filter(Cycle == "P") %>%
  pull(IBI) %>%
  quantile(probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

cat("\n[IBI QUARTILES based on Cycle P, all ages, IBI-complete]\n")
print(round(ibi_quartiles, 2))


# ===============================
# 7. Age restriction (>= 21) – before covariate listwise deletion
#    We'll define IBI quartile cuts here ON CYCLE P ONLY
# ===============================
df2 <- df1 %>%
  filter(Age >= 21)

cat("[STEP 2] Age >= 21 n =", nrow(df2), "\n")
cat("  (Age <21 removed =", nrow(df1) - nrow(df2), ")\n")


# ===============================
# 8. Derive Hypertension, binary Smoking/Alcohol, AgeGroup
# ===============================
df2 <- df2 %>%
  mutate(
    Hypertension = case_when(
      SYSTOLIC_BP >= 130 | DIASTOLIC_BP >= 80 ~ "Yes",
      is.na(SYSTOLIC_BP) | is.na(DIASTOLIC_BP) ~ NA_character_,
      TRUE ~ "No"
    ),
    Hypertension = factor(Hypertension, levels = c("No","Yes")),
    Smoking = case_when(
      Smoking_status == "Never Smoker" ~ "No",
      Smoking_status %in% c("Current Smoker","Former Smoker") ~ "Yes",
      TRUE ~ NA_character_
    ),
    Smoking = factor(Smoking, levels = c("No","Yes")),
    Alcohol_bin = case_when(
      Alcohol == "Non-Drinker" ~ "No",
      Alcohol %in% c("Occasional Drinker","Frequent Drinker") ~ "Yes",
      TRUE ~ NA_character_
    ),
    Alcohol_bin = factor(Alcohol_bin, levels = c("No","Yes")),
    AgeGroup = cut(
      Age,
      breaks = c(21, 40, 60, Inf),
      labels = c("21–39","40–59","60+"),
      right  = FALSE
    )
  )

# ===============================
# 9. Additional missing values across model covariates
# ===============================
reg_vars <- c(
  "CVD",
  "IBI","Age","Gender","Ethnicity","Education",
  "Smoking","Alcohol_bin","BMI","Diabetes",
  "TOTAL_CHOLESTEROL","HDL_CHOLESTEROL","Hypertension"
)

df3 <- df2 %>%
  filter(is.na(Diabetes) | Diabetes != "Don't know") %>%
  mutate(Diabetes = droplevels(as.factor(Diabetes))) %>%
  drop_na(any_of(reg_vars))

cat("[STEP 3] After additional missing values removed n =", nrow(df3), "\n")
cat("  (Removed =", nrow(df2) - nrow(df3), ")\n")

df3 %>% count(Cycle) %>% print()

# ===============================
# 10. Final factor baselines & IBI_Category using SAVED quartiles
# ===============================
df_final <- df3 %>%
  mutate(
    CVD          = factor(CVD, levels = c("No","Yes")),
    Gender       = relevel(factor(Gender),       ref = "Male"),
    Ethnicity    = relevel(factor(Ethnicity),    ref = "American"),
    Education    = relevel(factor(Education),    ref = "SomeCollegeorHigher"),
    Smoking      = relevel(factor(Smoking),      ref = "No"),
    Alcohol_bin  = relevel(factor(Alcohol_bin),  ref = "No"),
    Diabetes     = relevel(factor(Diabetes),     ref = "No"),
    Hypertension = relevel(factor(Hypertension), ref = "No"),
    Cycle        = factor(Cycle),
    IBI_Category = cut(
      IBI,
      breaks = ibi_quartiles,  # <-- PRE-SAVED CUTPOINTS (Cycle P adults)
      labels = c("Min-Q1","Q1-Q2","Q2-Q3","Q3-Max"),
      include.lowest = TRUE
    )
  )

cat("\n[FINAL] Analytic sample n =", nrow(df_final), "\n")
df_final %>% count(Cycle) %>% print()
df_final %>% count(IBI_Category) %>% print()

# ===============================
# 11. Bivariate table subset to Cycle == "P"
# ===============================
tab1 <- tableby(
  CVD ~ Gender + Age + AgeGroup + Ethnicity + Education +
    BMI + SYSTOLIC_BP + DIASTOLIC_BP +
    Smoking + Alcohol_bin + Diabetes +
    TOTAL_CHOLESTEROL + HDL_CHOLESTEROL +
    IBI + IBI_Category + Hypertension,
  data = df_final,
  numeric.stats = c("mean","sd"),
  numeric.test  = "kwt",  # Kruskal–Wallis
  cat.test      = "chisq"   # Pearson’s 
)

summary(tab1, text = TRUE, pfootnote = TRUE)


# Save for downstream use
saveRDS(df_final, file = "data_complete.rds")
