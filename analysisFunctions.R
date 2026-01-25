# Function to create results table

create_multi_var_sample_size_table <- function(df,
                                   row_var = "stakeholderCategory",
                                   col_var = "caseStudy",
                                   measure1_vars = c("ESI", "ESI_norm"),
                                   measure2_vars = c("WTS", "WTS_norm"),
                                   remove_columns = NULL
                                   ) {
  # Filter data to remove specified columns if requested
  if (!is.null(remove_columns)) {
    df <- df[!(df[[col_var]] %in% remove_columns), ]
  }

  # Create base table with all data
  table_all <- addmargins(table(df[, c(row_var, col_var)]))

  # Define complete case variables for each measure
  vars_measure1 <- c(measure1_vars, row_var, col_var)
  vars_measure2 <- c(measure2_vars, row_var, col_var)

  # Create tables for complete cases only
  table_measure1 <- addmargins(table(
    df[complete.cases(df[, vars_measure1]), c(row_var, col_var)]
  ))
  table_measure2 <- addmargins(table(
    df[complete.cases(df[, vars_measure2]), c(row_var, col_var)]
  ))

  # Create list of three separate tables plus measure names
  result <- list(
    table_all = table_all,
    table_measure1 = table_measure1,
    table_measure2 = table_measure2,
    measure1_name = measure1_vars[1],
    measure2_name = measure2_vars[1]
  )

  class(result) <- "multi_var_n_table"
  return(result)
}

# Print method for multi_var_n_table class
print.multi_var_n_table <- function(x, ...) {
  table_all <- x$table_all
  table_measure1 <- x$table_measure1
  table_measure2 <- x$table_measure2
  measure1_name <- x$measure1_name
  measure2_name <- x$measure2_name
  
  # Print header with measure names
  cat(sprintf("Total (complete for %s, complete for %s)\n\n",
              measure1_name, measure2_name))
  
  # Initialize combined table with same structure as table_all
  combined_table <- table_all
  dimnames(combined_table) <- dimnames(table_all)
  
  # Abbreviate row names (stakeholderCategory values) with minlength = 8
  rownames(combined_table) <- abbreviate(rownames(table_all), minlength = 8)
  
  # Abbreviate dimension name (stakeholderCategory header) with minlength = 8
  names(dimnames(combined_table))[1] <- abbreviate(names(dimnames(table_all))[1], minlength = 8)
  
  # Populate combined table
  for (i in seq_len(nrow(table_all))) {
    for (j in seq_len(ncol(table_all))) {
      row_name <- rownames(table_all)[i]
      col_name <- colnames(table_all)[j]
      
      val_measure1 <- if (row_name %in% rownames(table_measure1) &&
                          col_name %in% colnames(table_measure1)) {
        table_measure1[row_name, col_name]
      } else 0
      
      val_measure2 <- if (row_name %in% rownames(table_measure2) &&
                          col_name %in% colnames(table_measure2)) {
        table_measure2[row_name, col_name]
      } else 0
      
      combined_table[i, j] <- paste0(
        table_all[i, j],
        " (", val_measure1, ", ", val_measure2, ")"
      )
    }
  }
  
  # Print the combined table
  print(combined_table)
  invisible(x)
}

summarise_stakeholders <- function(dF) {
  # Create summary with counts
  summary_df <- dF %>%
    group_by(stakeholderCategory, caseStudy, stakeholderGroup) %>%
    summarise(n_participants = n(), .groups = 'drop') %>%
    arrange(stakeholderCategory, desc(n_participants))

  # Write to TSV file
  write.table(summary_df,
              file = "stakeholder_groups.tsv",
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)

  # Print formatted output
  for(i in 1:nrow(summary_df)) {
    cat(sprintf("%s, %s, %s: %d\n",
                summary_df$stakeholderCategory[i],
                summary_df$caseStudy[i],
                summary_df$stakeholderGroup[i],
                summary_df$n_participants[i]))
  }

  # Print summary statistics
  cat(sprintf("\nNumber of stakeholder groups: %d\n", nrow(summary_df)))
  cat(sprintf("Median size: %.1f\n", median(summary_df$n_participants)))
  cat(sprintf("Mean size: %.1f\n", mean(summary_df$n_participants)))
}

# Calculate grand mean across case studies
grand_mean <- function(values, complete_cases, case_study, transform = identity) {
  mean(tapply(transform(values[complete_cases]),
              case_study[complete_cases],
              function(x) mean(x, na.rm = TRUE)),
       na.rm = TRUE)
}

permutation_sign_test <- function(x, n_permutations = 99999, alternative = "two.sided", seed = NULL) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0) stop("No non-missing observations in x")
  if (!is.null(seed)) set.seed(seed)
  observed_mean <- mean(x)
  permuted_means <- numeric(n_permutations) # to store permuted means
  for (i in 1:n_permutations) {
    random_signs <- sample(c(-1, 1), size = n, replace = TRUE) # random sign permutation
    permuted_means[i] <- mean(x * random_signs)
  }
  p_value <- switch(alternative,
    "two.sided" = mean(abs(permuted_means) >= abs(observed_mean)),
    "greater" = mean(permuted_means >= observed_mean),
    "less" = mean(permuted_means <= observed_mean),
    stop("alternative must be 'two.sided', 'greater', or 'less'")
  )
  result <- list(
    observed_mean = observed_mean,
    p_value = p_value,
    permuted_means = permuted_means,
    n_permutations = n_permutations,
    n_observations = n,
    alternative = alternative
  )
  class(result) <- "permutation_sign_test"
  return(result)
}

# Print method for permutation_sign_test objects
print.permutation_sign_test <- function(x, digits = 4, ...) {
  cat("Permutation Sign Test\n")
  cat("=====================\n\n")
  cat("Observed mean:", round(x$observed_mean, digits), "\n")
  cat("Alternative hypothesis:", x$alternative, "\n")
  cat("Number of permutations:", x$n_permutations, "\n")
  cat("Number of observations:", x$n_observations, "\n")
  cat("P-value:", round(x$p_value, digits), "\n")
  sig <- ifelse(x$p_value < 0.001, "***",
         ifelse(x$p_value < 0.01, "**",
         ifelse(x$p_value < 0.05, "*",
         ifelse(x$p_value < 0.1, ".", ""))))
  if (sig != "") {
    cat("Significance:", sig, "\n")
    cat("---\n")
    cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }
  invisible(x)
}

calc_proportion_positive_groups <- function(model, random_var = "caseStudy") {
  # Get the design matrix for fixed effects
  X <- getME(model, "X")
  # Get fixed effects coefficients
  beta <- fixef(model)
  # Calculate fitted values from fixed effects only (excludes random effects)
  fixed_pred <- X %*% beta

  # Extract random intercept SD for specified random grouping variable
  vc <- as.data.frame(VarCorr(model))
  random_intercept_row <- vc$grp == random_var & vc$var1 == "(Intercept)"
  sigma_group <- vc$sdcor[random_intercept_row]

  # Get unique groups and their average fixed effect prediction
  model_data <- model.frame(model)
  group_ids <- model_data[[random_var]]

  # Calculate mean fixed effect for each group
  group_fixed_means <- tapply(fixed_pred, group_ids, mean)

  # For each case study, calculate probability that the total outcome is positive
  # where total outcome = fixed effect + random intercept
  group_probs <- pnorm(0, mean = group_fixed_means, sd = sigma_group, lower.tail = FALSE)

  # Return average probability across groups
  mean(group_probs)
}

create_pluralistic_ignorance_group_results_table <- function(
  dF, attitude_var, category_vars = c("caseStudy", "stakeholderCategory"), min_n = 3
) {
  perceived_norm_var <- paste0(attitude_var, "_norm")

  # Create empty data frame to store results with dynamic column names
  results <- data.frame(matrix(ncol = length(category_vars) + 3, nrow = 0))
  colnames(results) <- c(category_vars, "propAgree", "propPercNormAgree", "n")

  # Convert to appropriate data types
  for(i in 1:length(category_vars)) {
    results[[category_vars[i]]] <- character()
  }
  results$propAgree <- numeric()
  results$propPercNormAgree <- numeric()
  results$n <- numeric()

  # Get unique combinations of category variables
  category_combinations <- unique(dF[category_vars])
  category_combinations <- category_combinations[complete.cases(category_combinations), ]

  # Loop through each combination
  for (i in 1:nrow(category_combinations)) {
    # Create filter condition for this combination
    filter_condition <- rep(TRUE, nrow(dF))
    for (var in category_vars) {
      filter_condition <- filter_condition &
        (dF[[var]] == category_combinations[[var]][i]) &
        !is.na(dF[[var]])
    }

    # Subset data for this combination
    subset_data <- dF[filter_condition, ]

    # Only complete cases
    subset_data <- subset_data[!is.na(subset_data[[attitude_var]]) & !is.na(subset_data[[perceived_norm_var]]), ]
    n <- nrow(subset_data)

    # Skip if less than minimum participants
    if (n < min_n) next

    attitude_var_values <- as.numeric(subset_data[[attitude_var]])
    valid_main_values <- attitude_var_values[!is.na(attitude_var_values)]

    if (length(valid_main_values) == 0) {
      prop_agree <- NA
    } else {
      prop_agree <- sum(valid_main_values > 3) / length(valid_main_values)
    }

    # Calculate propPercNormAgree
    perceived_norm_values <- as.numeric(subset_data[[perceived_norm_var]])
    valid_other_values <- perceived_norm_values[!is.na(perceived_norm_values)]

    if (length(valid_other_values) == 0) {
      prop_believe_agree <- NA
    } else {
      prop_believe_agree <- mean(valid_other_values) / 100
    }

    # Create new row
    new_row <- category_combinations[i, ]
    new_row$propAgree <- prop_agree
    new_row$propPercNormAgree <- prop_believe_agree
    new_row$n <- n

    # Add row to results
    results <- rbind(results, new_row)
  }

  return(results)
}

# Faceted plot
create_pluralistic_ignorance_facet_plot <- function(
  data_list, attitude_vars, category_vars = c("caseStudy", "stakeholderCategory")
) {
  # Combine datasets with attitude variable labels
  combined_data <- bind_rows(
    data_list[[1]] %>% mutate(attitude = attitude_vars[1]),
    data_list[[2]] %>% mutate(attitude = attitude_vars[2])
  )

  # Define shape mapping for stakeholderCategory
  stakeholder_order <- c(
    "Commercial",
    "Tourists",
    "Government",
    "Local community",
    "Environment protectors",
    "Wild resource users"
  )

  # Define case study labels mapping
  case_study_labels <- c(
    "Estonia" = "Kloogaranna, Estonia",
    "Ireland" = "Streedagh, Ireland",
    "Italy" = "Frigole, Italy",
    "Malta" = "Comino, Malta",
    "Norway" = "Bergen, Norway",
    "Romania" = "Corbu, Romania",
    "Slovenia" = "Strunjan, Slovenia",
    "UK" = "Arne, UK"
  )

  # Set transparency level
  alpha_level <- 0.9

  # Determine shape scale based on whether we're using stakeholderCategory
  if (length(category_vars) >= 2 && category_vars[2] == "stakeholderCategory") {
    # Filter to only categories present in data and maintain order
    present_stakeholders <- stakeholder_order[stakeholder_order %in% unique(combined_data$stakeholderCategory)]

    shape_scale <- scale_shape_manual(
      name = "Stakeholder category",
      values = c(16, 17, 15, 18, 3, 4, 8, 9)[1:length(present_stakeholders)],
      breaks = present_stakeholders
    )
  } else {
    # Generic shape scale for other variables
    shape_scale <- scale_shape_manual(
      name = category_vars[2],
      values = c(16, 17, 15, 18, 3, 4, 8, 9)[1:length(unique(combined_data[[category_vars[2]]]))]
    )
  }

  # Create plot
  p <- ggplot(combined_data, aes(x = propAgree, y = propPercNormAgree)) +
    geom_point(aes(color = .data[[category_vars[1]]],
                   shape = .data[[category_vars[2]]],
                   size = n), alpha = alpha_level) +
    facet_wrap(~ attitude) +
    scale_color_brewer(palette = "Set2", name = "Case study", labels = case_study_labels) +
    shape_scale +
    scale_size_continuous(name = "Sample size", range = c(2, 8)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0.02, 0.02),
                       labels = function(x) paste0(x * 100)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.02, 0.02),
                       labels = function(x) paste0(x * 100)) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "solid") +
    labs(x = "Attitude (% agreeing)",
         y = "Perceived norm (% believed to agree)") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.box = "horizontal",
          legend.title.position = "top",
          legend.spacing.y = unit(0.1, "cm"),
          legend.key.height = unit(0.4, "cm"),
          strip.text = element_text(size = rel(1.0)),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(5, 5, 5, 5, "pt")) +
    coord_fixed(ratio = 1) +
    guides(
      color = guide_legend(override.aes = list(size = 4, alpha = alpha_level),
                          order = 1,
                          ncol = 1,
                          title.position = "top"),
      shape = guide_legend(override.aes = list(size = 4, alpha = alpha_level),
                          order = 2,
                          ncol = 1,
                          title.position = "top"),
      size = guide_legend(order = 3,
                         ncol = 1,
                         title.position = "top")
    )

  return(p)
}

standardised_complete_cases <- function(df, formula) {
  var_names <- all.vars(formula)
  selected_vars_df <- df[, var_names, drop = FALSE]
  complete_cases_df <- selected_vars_df[complete.cases(selected_vars_df), , drop = FALSE]
  for (col in names(complete_cases_df)) {
    if (is.numeric(complete_cases_df[[col]])) {
      complete_cases_df[[col]] <- scale(complete_cases_df[[col]])[, 1]
    }
  }
  return(complete_cases_df)
}

summarise_model_with_OR <- function(model, params) {
  print(summary(model))
  cat("\n")
  aic_value <- AIC(model)
  cat(sprintf("AIC: %.1f\n\n", aic_value))
  
  summ <- summary(model)
  coef_table <- summ$coefficients
  or <- exp(coef(model))[params]
  ci <- exp(confint(model))[params, , drop = FALSE]
  p_values <- setNames(
    coef_table[params, "Pr(>|z|)"],
    params
  )
  
  cat("Odds Ratios and 95% Confidence Intervals:\n")
  for (param in params) {
    cat(sprintf("%s: OR = %.2f, 95%% CI [%.2f, %.2f], p = %s\n",
                param, or[param], ci[param, 1], ci[param, 2],
                format.pval(p_values[param], digits = 3, eps = 0.001)))
  }
  cat("\n")
  
  result <- list(
    odds_ratios = or,
    confidence_intervals = ci,
    p_values = p_values,
    aic = aic_value
  )
  
  invisible(result)
}

create_comparison_table <- function(model_results, param) {
  n_models <- length(model_results)
  comparison_table <- data.frame(
    Model = paste("Model", 1:n_models),
    OR_CI = sapply(1:n_models, function(i) {
      or <- model_results[[i]]$odds_ratios[param]
      ci_lower <- model_results[[i]]$confidence_intervals[param, 1]
      ci_upper <- model_results[[i]]$confidence_intervals[param, 2]
      sprintf("%.2f [%.2f, %.2f]", or, ci_lower, ci_upper)
    }),
    p = sapply(1:n_models, function(i) {
      sprintf("%.3f", model_results[[i]]$p_values[param])
    }),
    AIC = round(sapply(model_results, function(x) x$aic), 1)
  )
  return(comparison_table)
}