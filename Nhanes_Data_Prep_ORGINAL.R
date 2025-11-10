
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

# ===============================
# 2. Helper Function
# ===============================
get_nhanes_data <- function(nhanes_table, vars, translate = TRUE) {
  data <- nhanes(nhanes_table)[, vars, drop = FALSE]
  if (translate) data <- nhanesTranslate(nh_table = nhanes_table, names(data), data = data)
  data
}

# ===============================
# 3. Load Cycle L Data
# ===============================
cycle_L_demo <- get_nhanes_data("DEMO_L", c("SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH3", "DMDEDUC2", "INDFMPIR", 
                                            "WTINT2YR", "SDMVSTRA", "WTMEC2YR", "SDMVPSU")) %>%
  rename(Gender = RIAGENDR, Age = RIDAGEYR, Ethnicity = RIDRETH3, Education = DMDEDUC2,
         Poverty_Ratio = INDFMPIR, SAMPLEWEIGHT = WTINT2YR) %>%
  mutate(Ethnicity = case_when(
    Ethnicity %in% c("Mexican American", "Other Hispanic") ~ "Hispanic",
    Ethnicity %in% c("Non-Hispanic White") ~ "Non-Hispanic White",
    Ethnicity %in% c("Non-Hispanic Black") ~ "Non-Hispanic Black",
    Ethnicity %in% c("Non-Hispanic Asian", "Other Race - Including Multi-Racial") ~ "Other Muti-Racial",
    TRUE ~ NA_character_
  ) %>% as.factor(),
  Education = case_when(
    Education %in% c("Less than 9th grade", "9-11th grade (Includes 12th grade with no diploma)", "Don't know") ~ "Below High School",
    Education %in% c("High school graduate/GED or equivalent", "Some college or AA degree") ~ "High School",
    Education == "College graduate or above" ~ "Above High School"
  ) %>% as.factor())

cycle_L_BMI <- get_nhanes_data("BMX_L", c("SEQN", "BMXBMI", "BMXWAIST")) %>% rename(BMI = BMXBMI, Waist_Circ = BMXWAIST)

cycle_L_BP <- nhanes("BPXO_L")[, c("SEQN", "BPXOSY1", "BPXOSY2", "BPXOSY3", "BPXODI1", "BPXODI2", "BPXODI3")] %>%
  mutate(SYSTOLIC_BP = rowMeans(select(., BPXOSY1, BPXOSY2, BPXOSY3), na.rm = TRUE),
         DIASTOLIC_BP = rowMeans(select(., BPXODI1, BPXODI2, BPXODI3), na.rm = TRUE)) %>%
  select(SEQN, SYSTOLIC_BP, DIASTOLIC_BP)

cycle_L_Smoke <- get_nhanes_data("SMQ_L", c("SEQN", "SMQ020", "SMQ040")) %>%
  mutate(Smoking_status = case_when(
    SMQ020 != "Yes" ~ "Never Smoker",  # Has not smoked 100+ cigarettes
    SMQ020 == "Yes" & SMQ040 %in% c("Every day", "Some days") ~ "Current Smoker",
    SMQ020 == "Yes" & SMQ040 == "Not at all" ~ "Former Smoker",
    TRUE ~ NA_character_
  ))%>%
  select(SEQN, Smoking_status)


cycle_L_Alch <- get_nhanes_data("ALQ_L", c("SEQN", "ALQ121")) %>% rename(Alcohol = ALQ121) %>%
  filter(!Alcohol %in% c("Refused", "Don't know")) %>%
  mutate(Alcohol = case_when(
    Alcohol %in% c("Every day", "Nearly every day", "3 to 4 times a week", "2 times a week", "Once a week") ~ "Frequent Drinker",
    Alcohol %in% c("2 to 3 times a month", "Once a month", "7 to 11 times in the last year", "3 to 6 times in the last year", "1 to 2 times in the last year") ~ "Occasional Drinker",
    Alcohol %in% c("Never in the last year") ~ "Non-Drinker",
    TRUE ~ NA_character_
  ))

cycle_L_Diabetes <- get_nhanes_data("DIQ_L", c("SEQN", "DIQ010")) %>% rename(Diabetes = DIQ010) %>% filter(Diabetes != "Don't know")

cycle_L_VitD <- get_nhanes_data("VID_L", c("SEQN", "LBXVIDMS")) %>% rename(LBXVIDMS = LBXVIDMS) %>%
  mutate(VitaminD_Status = case_when(
    LBXVIDMS < 30 ~ "Deficient",
    LBXVIDMS >= 30 & LBXVIDMS < 50 ~ "Insufficient",
    LBXVIDMS >= 50 ~ "Sufficient",
    TRUE ~ NA_character_
  ))

cycle_L_CProt <- get_nhanes_data("HSCRP_L", c("SEQN", "LBXHSCRP")) %>% rename(C_REATIVE_PROTEIN = LBXHSCRP)

cycle_L_TChol <- get_nhanes_data("TCHOL_L", c("SEQN", "LBXTC")) %>% rename(TOTAL_CHOLESTEROL = LBXTC)

cycle_L_HDL <- get_nhanes_data("HDL_L", c("SEQN", "LBDHDD")) %>% rename(HDL_CHOLESTEROL = LBDHDD)

cycle_L_CVD <- get_nhanes_data("MCQ_L", c("SEQN", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F")) %>%
  mutate(across(c(MCQ160B, MCQ160C, MCQ160D, MCQ160E, MCQ160F), as.character)) %>%
  mutate(across(c(MCQ160B, MCQ160C, MCQ160D, MCQ160E, MCQ160F), ~ ifelse(. == "Don't know", NA, .))) %>%
  mutate(CVD = ifelse(MCQ160B == "Yes" | MCQ160C == "Yes" | MCQ160D == "Yes" | 
                        MCQ160E == "Yes" | MCQ160F == "Yes", "Yes", "No")) %>%
  rename(HEART_FAILURE = MCQ160B, CORONARY_HEART_DISEASE = MCQ160C,
         ANGINA = MCQ160D, HEART_ATTACK = MCQ160E, STROKE = MCQ160F)

cycle_L_Depres <- get_nhanes_data("DPQ_L", c(
  "SEQN", "DPQ010", "DPQ020", "DPQ030", "DPQ040", "DPQ050",
  "DPQ060", "DPQ070", "DPQ080", "DPQ090")) %>%
  mutate(across(starts_with("DPQ0"), ~ ifelse(. %in% c(7, 9), NA, .))) %>%
  mutate(PHQ9_Score = rowSums(select(., starts_with("DPQ0")), na.rm = TRUE)) %>%
  mutate(Depression = case_when(
    PHQ9_Score == 0 ~ NA_character_,
    PHQ9_Score >= 10 ~ "Yes",
    TRUE ~ "No"
  ))%>%
  select(SEQN, Depression)

cycle_L_CBC <- get_nhanes_data("CBC_L", c("SEQN", "LBDLYMNO", "LBDNENO", "LBDMONO", "LBXPLTSI")) %>%
  rename(
    Lymphocyte = LBDLYMNO,
    Neutrophil = LBDNENO,
    Monocyte = LBDMONO,
    Platelet = LBXPLTSI
  )


cycle_L_list <- list(cycle_L_demo, cycle_L_BMI, cycle_L_BP, cycle_L_Smoke, cycle_L_Alch,
                     cycle_L_Diabetes,  cycle_L_CProt, cycle_L_TChol, 
                     cycle_L_HDL, cycle_L_CVD, cycle_L_Depres, cycle_L_CBC)

cycle_L_data <- reduce(cycle_L_list, full_join, by = "SEQN") %>% mutate(Cycle = "L")

# ===============================
# 4. Load Cycle P Data (standardized using helper function)
# ===============================
cycle_P_demo <- get_nhanes_data("P_DEMO", c("SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH3", "DMDEDUC2", "INDFMPIR", "WTINTPRP", "WTMECPRP", "SDMVPSU", "SDMVSTRA")) %>%
  rename(Gender = RIAGENDR, Age = RIDAGEYR, Ethnicity = RIDRETH3, Education = DMDEDUC2,
         Poverty_Ratio = INDFMPIR, SAMPLEWEIGHT = WTINTPRP, WTMEC2YR = WTMECPRP) %>%
  mutate(Ethnicity = case_when(
    Ethnicity %in% c("Mexican American", "Other Hispanic") ~ "Hispanic",
    Ethnicity %in% c("Non-Hispanic White") ~ "Non-Hispanic White",
    Ethnicity %in% c("Non-Hispanic Black") ~ "Non-Hispanic Black",
    Ethnicity %in% c("Non-Hispanic Asian", "Other Race - Including Multi-Racial") ~ "Other Muti-Racial",
    TRUE ~ NA_character_
  ) %>% as.factor(),
  Education = case_when(
    Education %in% c("Less than 9th grade", "9-11th grade (Includes 12th grade with no diploma)", "Don't know") ~ "Below High School",
    Education %in% c("High school graduate/GED or equivalent", "Some college or AA degree") ~ "High School",
    Education == "College graduate or above" ~ "Above High School"
  ) %>% as.factor())

cycle_P_BMI <- get_nhanes_data("P_BMX", c("SEQN", "BMXBMI", "BMXWAIST")) %>% rename(BMI = BMXBMI, Waist_Circ = BMXWAIST)
cycle_P_BP <- nhanes("P_BPXO")[, c("SEQN", "BPXOSY1", "BPXOSY2", "BPXOSY3", "BPXODI1", "BPXODI2", "BPXODI3")] %>%
  mutate(SYSTOLIC_BP = rowMeans(select(., BPXOSY1, BPXOSY2, BPXOSY3), na.rm = TRUE),
         DIASTOLIC_BP = rowMeans(select(., BPXODI1, BPXODI2, BPXODI3), na.rm = TRUE)) %>%
  select(SEQN, SYSTOLIC_BP, DIASTOLIC_BP)


cycle_P_Smoke <- get_nhanes_data("P_SMQ", c("SEQN", "SMQ020", "SMQ040")) %>%
  mutate(Smoking_status = case_when(
    SMQ020 != "Yes" ~ "Never Smoker",  # Has not smoked 100+ cigarettes
    SMQ020 == "Yes" & SMQ040 %in% c("Every day", "Some days") ~ "Current Smoker",
    SMQ020 == "Yes" & SMQ040 == "Not at all" ~ "Former Smoker",
    TRUE ~ NA_character_
  ))%>%
  select(SEQN, Smoking_status)

cycle_P_Alch <- get_nhanes_data("P_ALQ", c("SEQN", "ALQ121")) %>% rename(Alcohol = ALQ121) %>%
  filter(!Alcohol %in% c("Refused", "Don't know")) %>%
  mutate(Alcohol = case_when(
    Alcohol %in% c("Every day", "Nearly every day", "3 to 4 times a week", "2 times a week", "Once a week") ~ "Frequent Drinker",
    Alcohol %in% c("2 to 3 times a month", "Once a month", "7 to 11 times in the last year", "3 to 6 times in the last year", "1 to 2 times in the last year") ~ "Occasional Drinker",
    Alcohol %in% c("Never in the last year") ~ "Non-Drinker",
    TRUE ~ NA_character_
  ))

cycle_P_Diabetes <- get_nhanes_data("P_DIQ", c("SEQN", "DIQ010")) %>% rename(Diabetes = DIQ010) %>% filter(Diabetes != "Don't know")
#cycle_P_VitD <- get_nhanes_data("P_VID", c("SEQN", "LBXVIDMS")) %>% rename(LBXVIDMS = LBXVIDMS) %>%
#  mutate(VitaminD_Status = case_when(
#    LBXVIDMS < 30 ~ "Deficient",
#    LBXVIDMS >= 30 & LBXVIDMS < 50 ~ "Insufficient",
#    LBXVIDMS >= 50 ~ "Sufficient",
#    TRUE ~ NA_character_
#  ))

cycle_P_CProt <- get_nhanes_data("P_HSCRP", c("SEQN", "LBXHSCRP")) %>% rename(C_REATIVE_PROTEIN = LBXHSCRP)
cycle_P_TChol <- get_nhanes_data("P_TCHOL", c("SEQN", "LBXTC")) %>% rename(TOTAL_CHOLESTEROL = LBXTC)
cycle_P_HDL <- get_nhanes_data("P_HDL", c("SEQN", "LBDHDD")) %>% rename(HDL_CHOLESTEROL = LBDHDD)
cycle_P_CVD <- get_nhanes_data("P_MCQ", c("SEQN", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F")) %>%
  mutate(across(c(MCQ160B, MCQ160C, MCQ160D, MCQ160E, MCQ160F), as.character)) %>%
  mutate(across(c(MCQ160B, MCQ160C, MCQ160D, MCQ160E, MCQ160F), ~ ifelse(. == "Don't know", NA, .))) %>%
  mutate(CVD = ifelse(MCQ160B == "Yes" | MCQ160C == "Yes" | MCQ160D == "Yes" | 
                        MCQ160E == "Yes" | MCQ160F == "Yes", "Yes", "No")) %>%
  rename(HEART_FAILURE = MCQ160B, CORONARY_HEART_DISEASE = MCQ160C,
         ANGINA = MCQ160D, HEART_ATTACK = MCQ160E, STROKE = MCQ160F)

cycle_P_Depres <- get_nhanes_data("P_DPQ", c(
  "SEQN", "DPQ010", "DPQ020", "DPQ030", "DPQ040", "DPQ050",
  "DPQ060", "DPQ070", "DPQ080", "DPQ090")) %>%
  mutate(across(starts_with("DPQ0"), ~ ifelse(. %in% c(7, 9), NA, .))) %>%
  mutate(PHQ9_Score = rowSums(select(., starts_with("DPQ0")), na.rm = TRUE)) %>%
  mutate(Depression = case_when(
    PHQ9_Score == 0 ~ NA_character_,
    PHQ9_Score >= 10 ~ "Yes",
    TRUE ~ "No"
  ))%>%
  select(SEQN, Depression)

cycle_P_CBC <- get_nhanes_data("P_CBC", c("SEQN", "LBDLYMNO", "LBDNENO", "LBDMONO", "LBXPLTSI")) %>%
  rename(
    Lymphocyte = LBDLYMNO,
    Neutrophil = LBDNENO,
    Monocyte = LBDMONO,
    Platelet = LBXPLTSI
  )

cycle_P_list <- list(cycle_P_demo, cycle_P_BMI, cycle_P_BP, cycle_P_Smoke, cycle_P_Alch,
                     cycle_P_Diabetes,  cycle_P_CProt, cycle_P_TChol, cycle_P_HDL, 
                     cycle_P_CVD, cycle_P_Depres, cycle_P_CBC)

cycle_P_data <- reduce(cycle_P_list, full_join, by = "SEQN") %>% mutate(Cycle = "P")


# Get all unique column names across both cycles
all_vars <- union(names(cycle_L_data), names(cycle_P_data))

# Add missing columns to each dataset
cycle_L_data <- cycle_L_data %>%
  mutate(across(setdiff(all_vars, names(.)), ~ NA)) %>%
  select(all_of(all_vars))

cycle_P_data <- cycle_P_data %>%
  mutate(across(setdiff(all_vars, names(.)), ~ NA)) %>%
  select(all_of(all_vars))

# Combine both cycles into one dataset
nhanes_combined <- bind_rows(cycle_L_data, cycle_P_data) %>%
  filter(!is.na(SEQN))
data.comb = nhanes_combined[nhanes_combined$Age > 21,]

#Removing missing data lines from Total Cholesterol,C-Reactive Protein & HDL-Cholesterol
MISSING = is.na(data.comb$CVD )

data.comb = subset(data.comb,subset=!MISSING)


# Calculate missing percentage for each column
missing_summary <- data.comb %>%
  summarise(across(everything(), ~ mean(is.na(.)) * 100)) %>%
  pivot_longer(cols = everything(), names_to = "Predictors", values_to = "Missing_Percent") %>%
  arrange(desc(Missing_Percent))  # Sort by highest missing percentage

# Print missing percentage table
print(missing_summary, n= 70)

data.complete = na.omit(data.comb)
dim(data.complete)

data.complete <- data.complete %>%
  filter(Diabetes != "Don't know") %>%     # remove rows with "Don't know"
  mutate(Diabetes = droplevels(Diabetes))  # drop the unused level


data.complete <- data.complete %>%
  mutate(IBI = C_REATIVE_PROTEIN * (Neutrophil / Lymphocyte))

data.complete <- data.complete %>%
  mutate(SII = Platelet * (Neutrophil / Lymphocyte))

data.complete <- data.complete %>%
  mutate(AISI = Neutrophil * Platelet * (Monocyte / Lymphocyte))

data.complete <- data.complete %>%
  mutate(SIRI = Neutrophil * (Monocyte / Lymphocyte))

data.complete <- data.complete %>%
  mutate(IBI_Q = cut(IBI,
                     breaks = quantile(IBI, probs = c(0, 0.25, 0.5, 0.75, 1), 
                                       na.rm = TRUE), labels = c("Low", "Medium-Low", "Medium-High", "High"),
                     include.lowest = TRUE))


data.complete <- data.complete %>%
  mutate(SII_Q = cut(SII,
                     breaks = quantile(SII, probs = c(0, 0.25, 0.5, 0.75, 1), 
                                       na.rm = TRUE), labels = c("Low", "Medium-Low", "Medium-High", "High"),
                     include.lowest = TRUE))

data.complete <- data.complete %>%
  mutate(AISI_Q = cut(AISI,
                      breaks = quantile(AISI, probs = c(0, 0.25, 0.5, 0.75, 1), 
                                        na.rm = TRUE), labels = c("Low", "Medium-Low", "Medium-High", "High"),
                      include.lowest = TRUE))

data.complete <- data.complete %>%
  mutate(SIRI_Q = cut(SIRI,
                      breaks = quantile(SIRI, probs = c(0, 0.25, 0.5, 0.75, 1), 
                                        na.rm = TRUE), labels = c("Low", "Medium-Low", "Medium-High", "High"),
                      include.lowest = TRUE))

data.complete <- data.complete %>%
  mutate(
    Smoking_status = as.factor(Smoking_status),
    Alcohol = as.factor(Alcohol),
    HEART_FAILURE = as.factor(HEART_FAILURE),
    CORONARY_HEART_DISEASE = as.factor(CORONARY_HEART_DISEASE),
    ANGINA = as.factor(ANGINA),
    HEART_ATTACK = as.factor(HEART_ATTACK),
    STROKE = as.factor(STROKE),
    CVD = as.factor(CVD),
    Depression = as.factor(Depression),
    Cycle = as.factor(Cycle)
  )

str(data.complete)

data.complete <- data.complete %>%
  mutate(WTMEC2YR = WTMEC2YR*2)  # For 2 cycles



# data.complete = readRDS("data_complete_cleaned.rds")


# Create Hypertension variable
# Common definition: Hypertension if Systolic BP >= 130 or Diastolic BP >= 80
data.complete$Hypertension <- ifelse(data.complete$SYSTOLIC_BP >= 130 | data.complete$DIASTOLIC_BP >= 80, 
                                     "Yes", "No")
data.complete$Hypertension <- factor(data.complete$Hypertension, levels = c("No", "Yes"))


#data.complete = data.complete[data.complete$Depression=="Yes",]
# Make Age categorical
data.complete$Age_Group <- cut(data.complete$Age,
                               breaks = c(-Inf, 40, 60, Inf),
                               labels = c("<40", "40-60", ">60"),
                               right = FALSE)
# Make BMI categorical

data.complete$BMI_Group <- cut(data.complete$BMI,
                               breaks = c(-Inf, 18.5, 25, 30, Inf),
                               labels = c("Underweight", "Normal", "Overweight", "Obese"),
                               right = FALSE)
# Make PIR categorical

data.complete$PIR_Group <- cut(data.complete$Poverty_Ratio,
                               breaks = c(-Inf, 1, 2, 4, Inf),
                               labels = c("Below Poverty", "Low", "Middle", "High"),
                               right = FALSE)



# Set reference levels
data.complete$Education <- relevel(data.complete$Education, ref = "Below High School")
data.complete$Smoking_status <- relevel(data.complete$Smoking_status, ref = "Never Smoker")
data.complete$Alcohol <- relevel(data.complete$Alcohol, ref = "Non-Drinker")
data.complete$Diabetes <- relevel(data.complete$Diabetes, ref = "No")
data.complete$Smoking_status <- relevel(data.complete$Smoking_status, ref = "Never Smoker")
data.complete$BMI_Group <- relevel(data.complete$BMI_Group, ref = "Normal")

saveRDS(data.complete, file = "data_complete_cleaned.rds")


library(arsenal)
# Use of "arsenal package" for table summary with p-values
tab1 <- tableby(CVD ~ ., data=data.complete, numeric.stats=c("mean","sd")
                , numeric.test="kwt")
summary(tab1, text=TRUE, pfootnote=TRUE)


library(survey)
data.design <- svydesign(id = ~SDMVPSU,
                         strata = ~SDMVSTRA,
                         weights = ~WTMEC2YR, 
                         nest    = TRUE,
                         survey.lonely.psu = "adjust",
                         data = data.complete)

summary(data.design)
