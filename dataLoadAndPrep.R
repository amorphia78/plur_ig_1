library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(ordinal)

dF <- read.table("proCoastStakeholderSurveyWave1CompiledEsiWtsDemo.tsv",
                 header = TRUE,
                 sep = "\t",
                 stringsAsFactors = FALSE)

# Define easier variable names and calculate EPI scores for ESI and WTS
dF <- dF %>%
  group_by(caseStudy, stakeholderCategory) %>%
  mutate(
    # ESI EPI scores
    ESI = envCit1, # just renaming
    ESI_norm = envCit1OtherAgree * 10, # convert to percent (as on 10-pt scale). Note, this is perceived norm
    ESI_norm_real = mean(ESI >= 4, na.rm = T) * 100, # for caseStudy * stakeholderCat cell, because of grouping
    ESI_EPI_score = ESI_norm_real - ESI_norm, # Environmental Pluralistic Ignore score
    ESI_norm_dis = envCit1OtherDisagree * 10, # convert to percent (as on 10-pt scale). Note, this is perceived norm
    ESI_norm_real_dis = mean(ESI <= 2, na.rm = T) * 100, # for caseStudy * stakeholderCat cell, because of grouping
    ESI_EPI_score_dis = ESI_norm_real_dis - ESI_norm_dis, # Environmental Pluralistic Ignore score

    # WTS EPI scores
    WTS = willToSac2, # just renaming
    WTS_norm = willToSac2OtherAgree * 10, # convert to percent (as on 10-pt scale). Note, this is perceived norm
    WTS_norm_real = mean(WTS >= 4, na.rm = T) * 100, # for caseStudy * stakeholderCat cell, because of grouping
    WTS_EPI_score = WTS_norm_real - WTS_norm, # Environmental Pluralistic Ignore score
    WTS_norm_dis = willToSac2OtherDisagree * 10, # convert to percent (as on 10-pt scale). Note, this is perceived norm
    WTS_norm_real_dis = mean(WTS <= 2, na.rm = T) * 100, # for caseStudy * stakeholderCat cell, because of grouping
    WTS_EPI_score_dis = WTS_norm_real_dis - WTS_norm_dis, # Environmental Pluralistic Ignore score
  ) %>%
  ungroup()

#Define complete cases and keep only data that is complete case for at least one of ESI or WTS
dF$complete_cases_ESI <- complete.cases(dF[, c("ESI", "ESI_norm", "stakeholderCategory", "caseStudy")])
dF$complete_cases_WTS <- complete.cases(dF[, c("WTS", "WTS_norm", "stakeholderCategory", "caseStudy")])
dFWithIncompleteCases <- dF
dF <- dF[dF$complete_cases_ESI | dF$complete_cases_WTS, ]

## Model prep

# Make sure largest category is the reference level for stakeholderCategory
# Because of different patterns of missingness and major imbalance in this for WTS and ESI,
# we do it separately for the two analyses.

category_counts_ESI <- table(dF$stakeholderCategory[dF$complete_cases_ESI])
largest_category_ESI <- names(category_counts_ESI)[which.max(category_counts_ESI)]
dF$stakeholderCategoryForESI <- relevel(factor(dF$stakeholderCategory), ref = largest_category_ESI)

category_counts_WTS <- table(dF$stakeholderCategory[dF$complete_cases_WTS])
largest_category_WTS <- names(category_counts_WTS)[which.max(category_counts_WTS)]
dF$stakeholderCategoryForWTS <- relevel(factor(dF$stakeholderCategory), ref = largest_category_WTS)

# For ordinal modelling, we need the DVs as factors
dF$ESI_factor <- as.factor(dF$ESI)
dF$WTS_factor <- as.factor(dF$WTS)
