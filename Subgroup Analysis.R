setwd('C:/Users/danie/repos/595 Research/Code')


## ============================================================
## 0. Libraries
## ============================================================
library(survey)
library(dplyr)
library(purrr)
library(sjPlot)    # for get_model_data()
library(ggplot2)

## ============================================================
## 1. Load analytic data (Cycle P + L)
## ============================================================
data.complete <- readRDS("data_complete.rds")

cat("\n[INFO] Data loaded:\n")
cat(" - Rows:", nrow(data.complete), "  Cols:", ncol(data.complete), "\n")
cat(" - Cycle P Rows:", sum(data.complete$Cycle == "P"), "\n")
cat(" - Cycle L Rows:", sum(data.complete$Cycle == "L"), "\n")

## ============================================================
## 2. Restrict to TRAIN = Cycle P and create subgroup variables
## ============================================================
train_df <- subset(data.complete, Cycle == "P") %>%
  mutate(
    # Age subgroup: <60 vs ≥60
    Age60 = factor(ifelse(Age >= 60, "≥60", "<60"),
                   levels = c("≥60", "<60")),
    # BMI subgroup: <30 vs ≥30
    BMI30 = factor(ifelse(BMI >= 30, "≥30", "<30"),
                   levels = c("≥30", "<30")),
    # Diabetes subgroup: Yes vs No
    Diabetes_bin = factor(
      ifelse(Diabetes == "Yes", "Yes", "No"),
      levels = c("Yes", "No")
    )
  )

## ============================================================
## 3. Base formula and helpers
## ============================================================

# Main fully-adjusted model (your Model 3 structure)
base_formula <- CVD ~ IBI_Category + Age + Gender + Ethnicity + Education +
  Smoking + Alcohol_bin + BMI + Diabetes +
  TOTAL_CHOLESTEROL + HDL_CHOLESTEROL + Hypertension

# 3.1 Adjust formulas based on which subgroup we’re stratifying by
make_formulas_for_subgroup <- function(sg) {
  full   <- base_formula
  simple <- CVD ~ IBI_Category + Age + Gender  # fallback simpler model
  
  if (sg == "Gender") {
    
    full   <- update(full, . ~ . - Gender)
    simple <- CVD ~ IBI_Category + Age
  }
  
  if (sg == "Ethnicity") {
    full <- update(full, . ~ . - Ethnicity)
  }
  
  if (sg == "Hypertension") {
    full <- update(full, . ~ . - Hypertension)
  }
  
  if (sg == "Diabetes_bin") {
    
    full <- update(full, . ~ . - Diabetes)
  }
  
  # For Age60 and BMI30 we keep Age and BMI as continuous covariates – that’s fine
  
  list(full = full, simple = simple)
}

# 3.2 Placeholder for p for interaction (currently returns NA)
get_interaction_p <- function(data, subgroup_var, base_formula) {
  # Placeholder: returning NA keeps the column but leaves it blank
  NA_real_
}

# 3.3 Fit subgroup-specific survey-weighted model
fit_subgroup_model <- function(train_df,
                               subgroup_var,
                               subgroup_level,
                               full_formula,
                               simple_formula = CVD ~ IBI_Category + Age + Gender) {
  dat_sub <- train_df %>%
    dplyr::filter(.data[[subgroup_var]] == subgroup_level) %>%
    droplevels()
  
  # basic checks: enough data and variation in CVD
  if (nrow(dat_sub) < 50 || length(unique(dat_sub$CVD)) < 2) {
    message("Model failed for ", subgroup_var, " = ", subgroup_level,
            " (too few observations or no variation in CVD).")
    return(NULL)
  }
  
  des_sub <- svydesign(
    id     = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMEC2YR,
    nest   = TRUE,
    survey.lonely.psu = "adjust",
    data   = dat_sub
  )
  
  # try full model first
  fit_full <- try(
    svyglm(full_formula, design = des_sub, family = binomial("logit")),
    silent = TRUE
  )
  
  if (!inherits(fit_full, "try-error")) {
    return(fit_full)
  }
  
  message("Full model failed for ", subgroup_var, " = ", subgroup_level,
          ". Trying simpler model...")
  
  # fall back to simpler model
  fit_simple <- try(
    svyglm(simple_formula, design = des_sub, family = binomial("logit")),
    silent = TRUE
  )
  
  if (inherits(fit_simple, "try-error")) {
    message("Simple model ALSO failed for ", subgroup_var, " = ", subgroup_level,
            ". Skipping this subgroup.")
    return(NULL)
  }
  
  fit_simple
}

## ============================================================
## 4. Build summary table of ORs for Q3–Max vs Min–Q1 in each subgroup
## ============================================================

subgroup_vars <- c("Age60", "Gender", "Ethnicity", "BMI30", "Hypertension", "Diabetes_bin")

subgroup_results <- purrr::map_df(subgroup_vars, function(sg) {
  # force factor to control level order
  var_vec   <- factor(train_df[[sg]])
  levels_sg <- levels(droplevels(var_vec))
  
  # p for interaction
  p_int <- get_interaction_p(train_df, sg, base_formula)
  
  # get subgroup-specific formulas
  forms <- make_formulas_for_subgroup(sg)
  
  rows <- purrr::map_df(levels_sg, function(lv) {
    fit <- fit_subgroup_model(
      train_df,
      subgroup_var   = sg,
      subgroup_level = lv,
      full_formula   = forms$full,
      simple_formula = forms$simple
    )
    if (is.null(fit)) return(NULL)
    
    # get ORs (exp) and focus on Q3–Max vs Min–Q1
    md <- sjPlot::get_model_data(fit, type = "est", transform = "exp")
    row_ibihigh <- dplyr::filter(md, term == "IBI_CategoryQ3-Max")
    if (nrow(row_ibihigh) == 0L) return(NULL)
    
    n_lv <- sum(train_df[[sg]] == lv, na.rm = TRUE)
    
    row_ibihigh %>%
      dplyr::transmute(
        subgroup      = sg,
        level         = lv,
        number        = n_lv,
        OR            = estimate,
        CI_lower      = conf.low,
        CI_upper      = conf.high,
        p_value       = p.value,
        p_interaction = p_int
      )
  })
  
  rows
})

cat("\n[SUBGROUP RESULTS]\n")
print(subgroup_results)
cat("\nCounts by subgroup:\n")
print(dplyr::count(subgroup_results, subgroup))

## ============================================================
## 5. Forest plot for OR (Q3–Max vs Min–Q1 IBI) by subgroups
## ============================================================
if (!require(forestplot)) install.packages("forestplot")
library(forestplot)
library(stringr)


# Replace interaction p-values with subgroup p-values
subgroup_results <- subgroup_results %>%
  mutate(
    p_interaction = p_value   # use subgroup model p-value
  )

# 5.1 ordering and labels
subgroup_order <- c("Age60", "Gender", "Ethnicity", "BMI30", "Hypertension", "Diabetes_bin")

subgroup_labels <- c(
  Age60        = "Age",
  Gender       = "Gender",
  Ethnicity    = "Race/Ethnicity",
  BMI30        = "BMI",
  Hypertension = "Hypertension",
  Diabetes_bin = "Diabetes"
)

# Pretty level labels
pretty_level <- function(sg, lv) {
  lv_chr <- as.character(lv)
  if (sg == "Age60")        return(lv_chr)  # "≥60", "<60"
  if (sg == "BMI30")        return(lv_chr)  # "≥30", "<30"
  if (sg == "Diabetes_bin") return(lv_chr)  # "Yes", "No"
  lv_chr
}

# 5.2: Prepare df for plotting
df_plot <- subgroup_results %>%
  mutate(
    subgroup = factor(subgroup, levels = subgroup_order)
  ) %>%
  arrange(subgroup)

# interaction p-values by subgroup (currently all NA if placeholder)
p_int_df <- df_plot %>%
  group_by(subgroup) %>%
  summarise(
    p_int = suppressWarnings(first(na.omit(p_interaction))),
    .groups = "drop"
  )

# Build rows with group headers + level rows
rows_list <- list()

for (sg in subgroup_order) {
  tmp <- df_plot[df_plot$subgroup == sg, , drop = FALSE]
  if (nrow(tmp) == 0) next
  
  # header row
  p_int_val <- p_int_df$p_int[p_int_df$subgroup == sg]
  
  header_row <- data.frame(
    subgroup       = sg,
    subgroup_label = subgroup_labels[sg],
    level_label    = "",
    number         = NA_integer_,
    OR             = NA_real_,
    CI_lower       = NA_real_,
    CI_upper       = NA_real_,
    p_int_display  = ifelse(
      length(p_int_val) == 0 || is.na(p_int_val),
      "",
      ifelse(p_int_val < 0.05, "< 0.05", "≥ 0.05")
    ),
    stringsAsFactors = FALSE
  )
  
  # level rows
  level_rows <- tmp %>%
    dplyr::transmute(
      subgroup       = sg,
      subgroup_label = "",
      level_label    = pretty_level(sg, level),
      number         = number,
      OR             = OR,
      CI_lower       = CI_lower,
      CI_upper       = CI_upper,
      p_int_display  = ""
    )
  
  rows_list[[length(rows_list) + 1]] <- header_row
  rows_list[[length(rows_list) + 1]] <- level_rows
}

fp_df <- bind_rows(rows_list)

## 5.3: Build table text for forestplot
tabletext <- cbind(
  c("Variables", fp_df$subgroup_label),
  c("Label",     fp_df$level_label),
  c("Number",    ifelse(is.na(fp_df$number), "", format(fp_df$number, big.mark=","))),
  c("OR",        ifelse(is.na(fp_df$OR), "", sprintf("%.2f", fp_df$OR))),
  c("95% CI",    ifelse(
    is.na(fp_df$CI_lower),
    "",
    sprintf("(%.2f, %.2f)", fp_df$CI_lower, fp_df$CI_upper)
  )),
  c("P Value", fp_df$p_int_display) 
)

## 5.4: Prepare vectors for forestplot
mean_vals  <- c(NA, fp_df$OR)
lower_vals <- c(NA, fp_df$CI_lower)
upper_vals <- c(NA, fp_df$CI_upper)

is_summary <- c(TRUE, fp_df$subgroup_label != "")



# 5.5: Draw forest plot
forestplot::forestplot(
  labeltext  = tabletext,
  mean       = mean_vals,
  lower      = lower_vals,
  upper      = upper_vals,
  is.summary = is_summary,
  zero       = 1,
  clip       = c(0.2, 5),
  xlab       = "Odds ratio (IBI Q3–Max vs Min–Q1)"
)






