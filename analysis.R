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

# Describe stakeholder group sample sizes
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

## Confidence intervals for the actual-vs-perceived agreement gap ##########
## (case-study-stratified bootstrap of the grand-mean difference, i.e.
##  grand mean % actually agreeing minus grand mean % believed to agree)

library(boot)

set.seed(123)
B <- 99999

grand_mean_gap_ESI <- function(data, idx) {
  d <- data[idx, ]
  actual    <- grand_mean(d$ESI >= 4, d$complete_cases_ESI, d$caseStudy, function(x) x * 100)
  perceived <- grand_mean(d$ESI_norm, d$complete_cases_ESI, d$caseStudy)
  actual - perceived
}

grand_mean_gap_WTS <- function(data, idx) {
  d <- data[idx, ]
  actual    <- grand_mean(d$WTS >= 4, d$complete_cases_WTS, d$caseStudy, function(x) x * 100)
  perceived <- grand_mean(d$WTS_norm, d$complete_cases_WTS, d$caseStudy)
  actual - perceived
}

# Stratified by caseStudy so each replicate keeps the same per-case-study
# sample sizes as the observed data (case studies are fixed, not resampled).
boot_gap_ESI <- boot(dF, grand_mean_gap_ESI, R = B, strata = factor(dF$caseStudy))
boot_gap_WTS <- boot(dF, grand_mean_gap_WTS, R = B, strata = factor(dF$caseStudy))

ci_gap_ESI <- boot.ci(boot_gap_ESI, type = "perc")$percent[4:5]
ci_gap_WTS <- boot.ci(boot_gap_WTS, type = "perc")$percent[4:5]

cat(sprintf("ESI actual-perceived agreement gap: %.1f pp, 95%% CI [%.1f, %.1f]\n",
            boot_gap_ESI$t0, ci_gap_ESI[1], ci_gap_ESI[2]))
cat(sprintf("WTS actual-perceived agreement gap: %.1f pp, 95%% CI [%.1f, %.1f]\n",
            boot_gap_WTS$t0, ci_gap_WTS[1], ci_gap_WTS[2]))

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

# H2 for ESI
H2_results_ESI <- permutation_sign_test(dF$ESI_EPI_score_dis)
print(H2_results_ESI)

# H2 for WTS
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
                            nsim = 99999, use.u = FALSE, type = "parametric")
boot_results_WTS <- bootMer(mod_EPI_WTS, calc_proportion_positive_groups,
                            nsim = 99999, use.u = FALSE, type = "parametric")
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
# Note - it is actually unnecessary to standardise complete cases separately for each model,
# because complete cases are never different, as stakeholder category and case study is never missing.
# This is present in the code for historical reasons.

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






## ==========================================================================
##  Cross-construct specificity test  (Spearman, bootstrapped CIs and p)
## ==========================================================================
##
##  For each construct X, compare own X with its OWN perceived norm vs with the
##  OTHER construct's perceived norm:
##      matched    = rho(own X, perceived X-norm)
##      mismatched = rho(own X, perceived Y-norm)
##  A proposition-blind response disposition would inflate both equally, so a
##  positive (matched - mismatched) difference cannot arise from it: it is
##  evidence of proposition-specific own/perceived linkage.
##
##  Spearman because own-attitude items are 5-pt ordinals and norms are skewed;
##  case bootstrap because heavy ties break normal-theory tests for comparing
##  correlations, and because the two correlations share the `own X` variable.

library(boot)

set.seed(123)
B <- 10000
sp <- function(x, y) cor(x, y, method = "spearman")

## sample with all four quantities (the test needs every respondent to have both
## own attitudes and both perceived norms)
vars <- c("WTS", "WTS_norm", "ESI", "ESI_norm")
dX   <- dF[complete.cases(dF[, vars]), vars]
cat(sprintf("n with all four quantities: %d\n", nrow(dX)))

stat <- function(data, idx) {
  d <- data[idx, ]
  c(wts_match    = sp(d$WTS, d$WTS_norm),
    wts_mismatch = sp(d$WTS, d$ESI_norm),
    wts_diff     = sp(d$WTS, d$WTS_norm) - sp(d$WTS, d$ESI_norm),
    esi_match    = sp(d$ESI, d$ESI_norm),
    esi_mismatch = sp(d$ESI, d$WTS_norm),
    esi_diff     = sp(d$ESI, d$ESI_norm) - sp(d$ESI, d$WTS_norm))
}

bt <- boot(dX, stat, R = B)

## two-sided bootstrap p that the statistic differs from 0, with a floor
boot_p <- function(tstar) {
  tstar <- tstar[is.finite(tstar)]
  max(2 * min(mean(tstar <= 0), mean(tstar >= 0)), 1 / (length(tstar) + 1))
}

show <- function(label, idx) {
  est <- bt$t0[idx]
  ci  <- boot.ci(bt, index = idx, type = "perc")$percent[4:5]
  cat(sprintf("  %-36s %+.3f  95%% CI [%+.3f, %+.3f]  p = %.4f\n",
              label, est, ci[1], ci[2], boot_p(bt$t[, idx])))
}

cat("\n--- WTS ---\n")
show("rho(own WTS, WTS norm)  [matched]",    1)
show("rho(own WTS, ESI norm)  [mismatched]", 2)
show("difference (matched - mismatched)",    3)

cat("\n--- ESI ---\n")
show("rho(own ESI, ESI norm)  [matched]",    4)
show("rho(own ESI, WTS norm)  [mismatched]", 5)
show("difference (matched - mismatched)",    6)








# To render:
# rmarkdown::render("YOUR_ABSOLUTE_PATH/analysis.R", output_format = rmarkdown::html_document( pandoc_args = c("--metadata", "author=Anon For Peer Review")))
