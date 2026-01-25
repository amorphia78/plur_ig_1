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

## Sample description

# Number of survey responses
nrow(dFWithIncompleteCases)

# Number complete cases for at least one of ESI and WTS
nrow(dF)

# Proportion with university degrees
sum( table(dF$education)[c("6","8")] ) / sum( table(dF$education) )

# Proportions in age bands
dF %>%
  mutate(age_group = case_when(
    age <= 30 ~ "18-30",
    age <= 50 ~ "31-50",
    age <= 70 ~ "51-70",
    age > 70 ~ "71+"
  )) %>%
  count(age_group) %>%
  mutate(proportion = n / sum(n))

# Proportions of genders
table(dF$gender)/sum(table(dF$gender))

# Mean economic situation
mean(dF$econSituation,na.rm=T)

# Sample size table by group (Table S1)
# Montenegro is completely omitted because complete cases n = 0 for both ESI and and WTS
sampleCells <- create_multi_var_sample_size_table(
  df = dF, row_var = "stakeholderCategory", col_var = "caseStudy",
  measure1_vars = c("ESI", "ESI_norm"), measure2_vars = c("WTS", "WTS_norm"),
  remove_columns = c("Montenegro")
)
print( sampleCells )

# Describe stakeholder groups
# Note, stakeholder groups are specific to each case study,
# unlike stakeholder categories, which are generally applicable.
summarise_stakeholders(dF)

## Reliability

#ESI
cor.test(dF$envCit1,dF$envCit2)
sum(complete.cases(dF$envCit1, dF$envCit2))

#WTS
cor.test(dF$willToSac1,dF$willToSac2)
sum(complete.cases(dF$willToSac1, dF$willToSac2))

## Pluralistic analysis #############################################################

## Starting with agreement (H1)

## Variable descriptions

# Grand mean of percent agree with ESI item
grand_mean(dF$ESI >= 4, dF$complete_cases_ESI, dF$caseStudy, function(x) x * 100)

# Grand mean of estimated percent agree with ESI item
grand_mean(dF$ESI_norm, dF$complete_cases_ESI, dF$caseStudy)

# Grand mean of percent agree with WTS item
grand_mean(dF$WTS >= 4, dF$complete_cases_WTS, dF$caseStudy, function(x) x * 100)

# Grand mean of estimated percent agree with WTS item
grand_mean(dF$WTS_norm, dF$complete_cases_WTS, dF$caseStudy)

## Pre-registered tests

# H1 for ESI
H1_results_ESI <- permutation_sign_test(dF$ESI_EPI_score)
print(H1_results_ESI)

# H1 for WTS
H1_results_WTS <- permutation_sign_test(dF$WTS_EPI_score)
print(H1_results_WTS)

## Now disagreement (H2)

# Grand mean of percent disagree with ESI item
grand_mean(dF$ESI <= 2, dF$complete_cases_ESI, dF$caseStudy, function(x) x * 100)

# Grand mean of estimated percent disagree with ESI item
grand_mean(dF$ESI_norm_dis, dF$complete_cases_ESI, dF$caseStudy)

# Grand mean of percent disagree with WTS item
grand_mean(dF$WTS <= 2, dF$complete_cases_WTS, dF$caseStudy, function(x) x * 100)

# Grand mean of estimated percent disagree with WTS item
grand_mean(dF$WTS_norm_dis, dF$complete_cases_WTS, dF$caseStudy)

## Pre-registered tests

# H1 for ESI
H2_results_ESI <- permutation_sign_test(dF$ESI_EPI_score_dis)
print(H2_results_ESI)

# H1 for WTS
H2_results_WTS <- permutation_sign_test(dF$WTS_EPI_score_dis)
print(H2_results_WTS)


## Plotting pluralistic ignorance

table_for_plot_ESI <- create_pluralistic_ignorance_group_results_table(dF,"ESI")
table_for_plot_ESI
table_for_plot_WTS <- create_pluralistic_ignorance_group_results_table(dF,"WTS")
table_for_plot_WTS

combined_plot <- create_pluralistic_ignorance_facet_plot(
  data_list = list(table_for_plot_ESI, table_for_plot_WTS),
  attitude_vars = c("ESI", "WTS")
)

print(combined_plot)
ggsave("pluralistic_ignorance_plot.png",
       plot = combined_plot,
       width = 6.7,
       height = 6.2,
       units = "in",
       dpi = 300,
       bg = "white")


## Generality analysis

set.seed(123)
mod_EPI_ESI <- lmer( ESI_EPI_score ~ stakeholderCategory + ( 1 | caseStudy), data=dF )
mod_EPI_WTS <- lmer( WTS_EPI_score ~ stakeholderCategory + ( 1 | caseStudy), data=dF )
# Run the bootstrap
boot_results_ESI <- bootMer(mod_EPI_ESI, calc_proportion_positive_groups,
                            nsim = 99, use.u = FALSE, type = "parametric")
boot_results_WTS <- bootMer(mod_EPI_WTS, calc_proportion_positive_groups,
                            nsim = 99, use.u = FALSE, type = "parametric")
ci_ESI <- quantile(boot_results_ESI$t, c(0.025, 0.975))
cat(" ESI Point estimate:", calc_proportion_positive_groups(mod_EPI_ESI), "\n",
  "ESI Bootstrap mean:", mean(boot_results_ESI$t), "\n", 
  "ESI 95% CI: [", ci_ESI[1], ",", ci_ESI[2], "]\n")
ci_WTS <- quantile(boot_results_WTS$t, c(0.025, 0.975))
cat(" WTS Point estimate:", calc_proportion_positive_groups(mod_EPI_WTS), "\n",
  "WTS Bootstrap mean:", mean(boot_results_WTS$t), "\n", 
  "WTS 95% CI: [", ci_WTS[1], ",", ci_WTS[2], "]\n")

## Modelling of personal attitude ~ perceived norm ######################################################

# Linear models are just not OK with this heavily skewed ordinal data, as revealed by the residual plot.
# So as per prereg, we are going to do ordinal models instead.
formula_ESI_1 <- ESI ~ ESI_norm
model_ESI_1 <- lm(formula_ESI_1, data = standardised_complete_cases( dF, formula_ESI_1 ) )
plot(jitter(fitted(model_ESI_1)), jitter(residuals(model_ESI_1)), xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")

## ESI models

# Initialize list to store results of all models
model_results_ESI <- list()

formula_ESI_1 <- ESI_factor ~ ESI_norm
model_ESI_1_ordinal <- clm(formula_ESI_1, data = standardised_complete_cases(dF, formula_ESI_1))
model_results_ESI[[1]] <- summarise_model_with_OR(model_ESI_1_ordinal, "ESI_norm")

#Testing the proportional odds assumption for ESI
nominal_test(model_ESI_1_ordinal)

formula_ESI_2 <- ESI_factor ~ ESI_norm + stakeholderCategoryForESI
model_ESI_2_ordinal <- clm(formula_ESI_2, data = standardised_complete_cases(dF, formula_ESI_2))
model_results_ESI[[2]] <- summarise_model_with_OR(model_ESI_2_ordinal, "ESI_norm")

formula_ESI_3 <- ESI_factor ~ ESI_norm + stakeholderCategoryForESI + (1|caseStudy)
model_ESI_3_ordinal <- clmm(formula_ESI_3, data = standardised_complete_cases(dF, formula_ESI_3))
model_results_ESI[[3]] <- summarise_model_with_OR(model_ESI_3_ordinal, "ESI_norm")

formula_ESI_4 <- ESI_factor ~ ESI_norm + stakeholderCategoryForESI + (1|caseStudy) + ESI_norm:stakeholderCategoryForESI
model_ESI_4_ordinal <- clmm(formula_ESI_4, data = standardised_complete_cases(dF, formula_ESI_4))
model_results_ESI[[4]] <- summarise_model_with_OR(model_ESI_4_ordinal, "ESI_norm")

formula_ESI_5 <- ESI_factor ~ ESI_norm + stakeholderCategoryForESI + (1|caseStudy) + (0 + ESI_norm|caseStudy)
model_ESI_5_ordinal <- clmm(formula_ESI_5, data = standardised_complete_cases(dF, formula_ESI_5))
model_results_ESI[[5]] <- summarise_model_with_OR(model_ESI_5_ordinal, "ESI_norm")

comparison_table_ESI <- create_comparison_table(model_results_ESI, "ESI_norm")

## WTS models

# Initialize list to store results of all models
model_results_WTS <- list()

formula_WTS_1 <- WTS_factor ~ WTS_norm
model_WTS_1_ordinal <- clm(formula_WTS_1, data = standardised_complete_cases(dF, formula_WTS_1))
model_results_WTS[[1]] <- summarise_model_with_OR(model_WTS_1_ordinal, "WTS_norm")

#Testing the proportional odds assumption for WTS
nominal_test(model_WTS_1_ordinal)

formula_WTS_2 <- WTS_factor ~ WTS_norm + stakeholderCategoryForWTS
model_WTS_2_ordinal <- clm(formula_WTS_2, data = standardised_complete_cases(dF, formula_WTS_2))
model_results_WTS[[2]] <- summarise_model_with_OR(model_WTS_2_ordinal, "WTS_norm")

formula_WTS_3 <- WTS_factor ~ WTS_norm + stakeholderCategoryForWTS + (1|caseStudy)
model_WTS_3_ordinal <- clmm(formula_WTS_3, data = standardised_complete_cases(dF, formula_WTS_3))
model_results_WTS[[3]] <- summarise_model_with_OR(model_WTS_3_ordinal, "WTS_norm")

formula_WTS_4 <- WTS_factor ~ WTS_norm + stakeholderCategoryForWTS + (1|caseStudy) + WTS_norm:stakeholderCategoryForWTS
model_WTS_4_ordinal <- clmm(formula_WTS_4, data = standardised_complete_cases(dF, formula_WTS_4))
model_results_WTS[[4]] <- summarise_model_with_OR(model_WTS_4_ordinal, "WTS_norm")

formula_WTS_5 <- WTS_factor ~ WTS_norm + stakeholderCategoryForWTS + (1|caseStudy) + (0 + WTS_norm|caseStudy)
model_WTS_5_ordinal <- clmm(formula_WTS_5, data = standardised_complete_cases(dF, formula_WTS_5))
model_results_WTS[[5]] <- summarise_model_with_OR(model_WTS_5_ordinal, "WTS_norm")

comparison_table_WTS <- create_comparison_table(model_results_WTS, "WTS_norm")

print(comparison_table_ESI)
print(comparison_table_WTS)

formula_WTS_6 <- WTS_factor ~ WTS_norm + ESI + stakeholderCategoryForWTS + (1|caseStudy)
model_WTS_6_ordinal <- clmm(formula_WTS_6, data = standardised_complete_cases(dF, formula_WTS_6))
summarise_model_with_OR(model_WTS_6_ordinal, c( "WTS_norm", "ESI" ) )
exp(coef(model_WTS_6_ordinal))[c("WTS_norm", "ESI")]
exp(confint(model_WTS_6_ordinal))[c("WTS_norm", "ESI"),]

# To render:
# rmarkdown::render("YOUR_ABSOLUTE_PATH/analysis.R", output_format = rmarkdown::html_document( pandoc_args = c("--metadata", "author=Anon For Peer Review")))
