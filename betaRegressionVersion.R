## Setup ###################################################################################################

# To run at the command line, replace this with a function that calls setwd()
# to the local directory where you have the data and R files.
setwd_PROCOAST_pluralistic_ignorance_1()

# Results in packages and own-written analysis functions loaded into memory
source( "analysisFunctions.R", local=T )

# Results in a data frame dF loaded into memory
# dF only has data points that are complete cases
# (both personal and perceived attitude) for at least one of ESI and WTS
# dFWithIncompleteCases is also loaded, which has all data
source( "dataLoadAndPrep.R", local=T  )

# Load additional required packages
library(glmmTMB)

# Function to rescale ordinal data to (0,1) for beta regression
rescale_to_beta <- function(x, min_val = 1, max_val = 5) {
  # Rescale from [min, max] to (0, 1) using the transformation:
  # (y - min + c) / (max - min + 2c) where c = 0.5
  c <- 0.5
  (x - min_val + c) / (max_val - min_val + 2 * c)
}

# Function to summarize beta regression models with effect sizes
summarise_beta_model <- function(model, params) {
  print(summary(model))
  cat("\n")

  aic_value <- AIC(model)
  cat(sprintf("AIC: %.1f\n\n", aic_value))

  # For beta regression, we'll use conditional R-squared from performance package
  # or just skip if it causes issues
  tryCatch({
    r2 <- performance::r2(model)
    cat(sprintf("Conditional R²: %.3f\n", r2$R2_conditional))
    cat(sprintf("Marginal R²: %.3f\n\n", r2$R2_marginal))
  }, error = function(e) {
    cat("R² calculation skipped\n\n")
  })

  # Extract coefficients for specified parameters
  summ <- summary(model)
  coef_table <- coef(summ)$cond  # conditional model coefficients

  # For beta regression, coefficients are interpreted differently than odds ratios
  # We'll report exp(coef) as multiplicative effects on the odds
  effects <- exp(coef_table[params, "Estimate"])
  names(effects) <- params

  # Calculate confidence intervals
  ci_lower <- exp(coef_table[params, "Estimate"] - 1.96 * coef_table[params, "Std. Error"])
  ci_upper <- exp(coef_table[params, "Estimate"] + 1.96 * coef_table[params, "Std. Error"])
  ci <- cbind(ci_lower, ci_upper)
  rownames(ci) <- params

  p_values <- coef_table[params, "Pr(>|z|)"]
  names(p_values) <- params

  cat("Exponentiated Coefficients (multiplicative effects) and 95% Confidence Intervals:\n")
  for (param in params) {
    cat(sprintf("%s: Exp(B) = %.2f, 95%% CI [%.2f, %.2f], p = %s\n",
                param, effects[param], ci[param, 1], ci[param, 2],
                format.pval(p_values[param], digits = 3, eps = 0.001)))
  }
  cat("\n")

  result <- list(
    effects = effects,
    confidence_intervals = ci,
    p_values = p_values,
    aic = aic_value
  )

  invisible(result)
}

# Create comparison table for beta models
create_beta_comparison_table <- function(model_results, param) {
  n_models <- length(model_results)
  comparison_table <- data.frame(
    Model = paste("Model", 1:n_models),
    Effect_CI = sapply(1:n_models, function(i) {
      eff <- model_results[[i]]$effects[param]
      ci_lower <- model_results[[i]]$confidence_intervals[param, 1]
      ci_upper <- model_results[[i]]$confidence_intervals[param, 2]
      sprintf("%.2f [%.2f, %.2f]", eff, ci_lower, ci_upper)
    }),
    p = sapply(1:n_models, function(i) {
      sprintf("%.3f", model_results[[i]]$p_values[param])
    }),
    AIC = round(sapply(model_results, function(x) x$aic), 1)
  )
  return(comparison_table)
}

## ESI Beta Regression Models

cat("=== ESI BETA REGRESSION MODELS ===\n\n")

# Initialize list to store results
model_results_ESI_beta <- list()

# DEFINE THE ORIGINAL FORMULAS (these are needed for standardised_complete_cases)
formula_ESI_1 <- ESI_factor ~ ESI_norm
formula_ESI_2 <- ESI_factor ~ ESI_norm + stakeholderCategoryForESI
formula_ESI_3 <- ESI_factor ~ ESI_norm + stakeholderCategoryForESI + (1|caseStudy)
formula_ESI_4 <- ESI_factor ~ ESI_norm + stakeholderCategoryForESI + (1|caseStudy) + ESI_norm:stakeholderCategoryForESI
formula_ESI_5 <- ESI_factor ~ ESI_norm + stakeholderCategoryForESI + (1|caseStudy) + (0 + ESI_norm|caseStudy)

# Model 1: Simple model with only ESI_norm
cat("--- ESI Model 1 ---\n")
data_ESI_1 <- standardised_complete_cases(dF, formula_ESI_1)
# Convert factor to numeric properly: factor -> character -> numeric
data_ESI_1$ESI_factor_beta <- rescale_to_beta(as.numeric(as.character(data_ESI_1$ESI_factor)))

formula_ESI_1_beta <- ESI_factor_beta ~ ESI_norm
model_ESI_1_beta <- glmmTMB(
  formula_ESI_1_beta,
  data = data_ESI_1,
  family = beta_family()
)
model_results_ESI_beta[[1]] <- summarise_beta_model(model_ESI_1_beta, "ESI_norm")

# Model 2: Add stakeholder category
cat("--- ESI Model 2 ---\n")
data_ESI_2 <- standardised_complete_cases(dF, formula_ESI_2)
data_ESI_2$ESI_factor_beta <- rescale_to_beta(as.numeric(as.character(data_ESI_2$ESI_factor)))

formula_ESI_2_beta <- ESI_factor_beta ~ ESI_norm + stakeholderCategoryForESI
model_ESI_2_beta <- glmmTMB(
  formula_ESI_2_beta,
  data = data_ESI_2,
  family = beta_family()
)
model_results_ESI_beta[[2]] <- summarise_beta_model(model_ESI_2_beta, "ESI_norm")

# Model 3: Add random intercept for case study
cat("--- ESI Model 3 ---\n")
data_ESI_3 <- standardised_complete_cases(dF, formula_ESI_3)
data_ESI_3$ESI_factor_beta <- rescale_to_beta(as.numeric(as.character(data_ESI_3$ESI_factor)))

formula_ESI_3_beta <- ESI_factor_beta ~ ESI_norm + stakeholderCategoryForESI + (1|caseStudy)
model_ESI_3_beta <- glmmTMB(
  formula_ESI_3_beta,
  data = data_ESI_3,
  family = beta_family()
)
model_results_ESI_beta[[3]] <- summarise_beta_model(model_ESI_3_beta, "ESI_norm")

# Model 4: Add interaction
cat("--- ESI Model 4 ---\n")
data_ESI_4 <- standardised_complete_cases(dF, formula_ESI_4)
data_ESI_4$ESI_factor_beta <- rescale_to_beta(as.numeric(as.character(data_ESI_4$ESI_factor)))

formula_ESI_4_beta <- ESI_factor_beta ~ ESI_norm + stakeholderCategoryForESI + (1|caseStudy) +
  ESI_norm:stakeholderCategoryForESI
model_ESI_4_beta <- glmmTMB(
  formula_ESI_4_beta,
  data = data_ESI_4,
  family = beta_family()
)
model_results_ESI_beta[[4]] <- summarise_beta_model(model_ESI_4_beta, "ESI_norm")

# Model 5: Add random slope for ESI_norm
cat("--- ESI Model 5 ---\n")
data_ESI_5 <- standardised_complete_cases(dF, formula_ESI_5)
data_ESI_5$ESI_factor_beta <- rescale_to_beta(as.numeric(as.character(data_ESI_5$ESI_factor)))

formula_ESI_5_beta <- ESI_factor_beta ~ ESI_norm + stakeholderCategoryForESI + (1|caseStudy) +
  (0 + ESI_norm|caseStudy)
model_ESI_5_beta <- glmmTMB(
  formula_ESI_5_beta,
  data = data_ESI_5,
  family = beta_family()
)
model_results_ESI_beta[[5]] <- summarise_beta_model(model_ESI_5_beta, "ESI_norm")

# Create comparison table
comparison_table_ESI_beta <- create_beta_comparison_table(model_results_ESI_beta, "ESI_norm")

## WTS Beta Regression Models

cat("\n=== WTS BETA REGRESSION MODELS ===\n\n")

# Initialize list to store results
model_results_WTS_beta <- list()

# DEFINE THE ORIGINAL FORMULAS (these are needed for standardised_complete_cases)
formula_WTS_1 <- WTS_factor ~ WTS_norm
formula_WTS_2 <- WTS_factor ~ WTS_norm + stakeholderCategoryForWTS
formula_WTS_3 <- WTS_factor ~ WTS_norm + stakeholderCategoryForWTS + (1|caseStudy)
formula_WTS_4 <- WTS_factor ~ WTS_norm + stakeholderCategoryForWTS + (1|caseStudy) + WTS_norm:stakeholderCategoryForWTS
formula_WTS_5 <- WTS_factor ~ WTS_norm + stakeholderCategoryForWTS + (1|caseStudy) + (0 + WTS_norm|caseStudy)
formula_WTS_6 <- WTS_factor ~ WTS_norm + ESI + stakeholderCategoryForWTS + (1|caseStudy)

# Model 1: Simple model with only WTS_norm
cat("--- WTS Model 1 ---\n")
data_WTS_1 <- standardised_complete_cases(dF, formula_WTS_1)
data_WTS_1$WTS_factor_beta <- rescale_to_beta(as.numeric(as.character(data_WTS_1$WTS_factor)))

formula_WTS_1_beta <- WTS_factor_beta ~ WTS_norm
model_WTS_1_beta <- glmmTMB(
  formula_WTS_1_beta,
  data = data_WTS_1,
  family = beta_family()
)
model_results_WTS_beta[[1]] <- summarise_beta_model(model_WTS_1_beta, "WTS_norm")

# Model 2: Add stakeholder category
cat("--- WTS Model 2 ---\n")
data_WTS_2 <- standardised_complete_cases(dF, formula_WTS_2)
data_WTS_2$WTS_factor_beta <- rescale_to_beta(as.numeric(as.character(data_WTS_2$WTS_factor)))

formula_WTS_2_beta <- WTS_factor_beta ~ WTS_norm + stakeholderCategoryForWTS
model_WTS_2_beta <- glmmTMB(
  formula_WTS_2_beta,
  data = data_WTS_2,
  family = beta_family()
)
model_results_WTS_beta[[2]] <- summarise_beta_model(model_WTS_2_beta, "WTS_norm")

# Model 3: Add random intercept for case study
cat("--- WTS Model 3 ---\n")
data_WTS_3 <- standardised_complete_cases(dF, formula_WTS_3)
data_WTS_3$WTS_factor_beta <- rescale_to_beta(as.numeric(as.character(data_WTS_3$WTS_factor)))

formula_WTS_3_beta <- WTS_factor_beta ~ WTS_norm + stakeholderCategoryForWTS + (1|caseStudy)
model_WTS_3_beta <- glmmTMB(
  formula_WTS_3_beta,
  data = data_WTS_3,
  family = beta_family()
)
model_results_WTS_beta[[3]] <- summarise_beta_model(model_WTS_3_beta, "WTS_norm")

# Model 4: Add interaction
cat("--- WTS Model 4 ---\n")
data_WTS_4 <- standardised_complete_cases(dF, formula_WTS_4)
data_WTS_4$WTS_factor_beta <- rescale_to_beta(as.numeric(as.character(data_WTS_4$WTS_factor)))

formula_WTS_4_beta <- WTS_factor_beta ~ WTS_norm + stakeholderCategoryForWTS + (1|caseStudy) +
  WTS_norm:stakeholderCategoryForWTS
model_WTS_4_beta <- glmmTMB(
  formula_WTS_4_beta,
  data = data_WTS_4,
  family = beta_family()
)
model_results_WTS_beta[[4]] <- summarise_beta_model(model_WTS_4_beta, "WTS_norm")

# Model 5: Add random slope for WTS_norm
cat("--- WTS Model 5 ---\n")
data_WTS_5 <- standardised_complete_cases(dF, formula_WTS_5)
data_WTS_5$WTS_factor_beta <- rescale_to_beta(as.numeric(as.character(data_WTS_5$WTS_factor)))

formula_WTS_5_beta <- WTS_factor_beta ~ WTS_norm + stakeholderCategoryForWTS + (1|caseStudy) +
  (0 + WTS_norm|caseStudy)
model_WTS_5_beta <- glmmTMB(
  formula_WTS_5_beta,
  data = data_WTS_5,
  family = beta_family()
)
model_results_WTS_beta[[5]] <- summarise_beta_model(model_WTS_5_beta, "WTS_norm")

# Create comparison table
comparison_table_WTS_beta <- create_beta_comparison_table(model_results_WTS_beta, "WTS_norm")

# Model 6: WTS with ESI included (matching your final model)
cat("--- WTS Model 6 (with ESI) ---\n")
data_WTS_6 <- standardised_complete_cases(dF, formula_WTS_6)
data_WTS_6$WTS_factor_beta <- rescale_to_beta(as.numeric(as.character(data_WTS_6$WTS_factor)))

formula_WTS_6_beta <- WTS_factor_beta ~ WTS_norm + ESI + stakeholderCategoryForWTS + (1|caseStudy)
model_WTS_6_beta <- glmmTMB(
  formula_WTS_6_beta,
  data = data_WTS_6,
  family = beta_family()
)
summarise_beta_model(model_WTS_6_beta, c("WTS_norm", "ESI"))

# Print comparison tables
{
cat("\n=== COMPARISON TABLES ===\n\n")
cat("ESI Models:\n")
print(comparison_table_ESI_beta)
cat("\n")
cat("WTS Models:\n")
print(comparison_table_WTS_beta)
}