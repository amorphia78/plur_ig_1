## Analysis details

### Approach

All analyses are conducted with R. All data and analysis scripts (including R Markdown showing results) can be downloaded from our GitHub repository (Anon, 2026), which will be preserved with a snapshot on Zenodo following acceptance of the manuscript.

### Data preparation

Participants were included in the analyses if they provided complete data for at least one of the two primary attitude measures (ESI or WTS). Complete case requirements were present data for the attitude variable and the perceived norm variable.

### Variable computation

For each stakeholder group within each case study, we computed Environmental Pluralistic Ignorance (EPI) scores for individuals as follows:

- The actual proportion of individuals who agreed with the attitude (scored ≥4 on the 5-point scale), converted to a percentage (same value for all individuals within the group)
- The perceived norm (what percentage of others that respondents believed agreed), based on responses to a 0-10 scale converted to percentages
- The EPI score is the actual percentage minus the perceived norm percentage

We computed both agreement-based EPI scores (for testing H1) and disagreement-based EPI scores (for testing H2), where disagreement was defined as scores ≤2 on the 5-point scale.

For ordinal regression models, all continuous predictors (perceived norms and, where applicable, an attitude score) were standardized (mean-centred and scaled to unit variance) within the complete case dataset.

### Variable descriptions across case studies

Where variables are described with percentages, for example mean percentage agreeing with items, or mean percentage participants estimated would agree, these are grand means of case study means, with no weighting.

### Reliability analysis

We assessed reliability for both attitude measures using Pearson correlations between the two items administered for each construct, using all cases complete for the two items. For ESI, we correlated envCit1 ("Acting environmentally friendly is an important part of who I am") with envCit2 ("I see myself as an environmentally friendly person") (_n_ = 333). For WTS, we correlated willToSac1 ("I would be willing to accept cuts in my standard of living to protect the environment") with willToSac2 ("I would be willing to pay much higher taxes in order to protect the environment") (_n_ = 230).

### Permutation sign tests for pluralistic ignorance (H1 and H2)

We tested whether EPI scores differed significantly from zero using permutation-based sign tests (Oehlert, 2023). The test procedure:

- Calculate the observed mean EPI score
- Generate 99,999 random permutations by randomly multiplying each observation by +1 or -1
- Calculate the mean for each permutation
- To generate a two-sided test, compute the _p_\-value as the proportion of permuted absolute means that equal or exceed the observed absolute mean

### Parametric bootstrap for proportion of case studies with positive EPI

To quantify uncertainty in the proportion of case studies showing positive pluralistic ignorance (i.e., positive EPI scores), we conducted parametric bootstrap analyses for both ESI and WTS. This analysis addresses the question: across how many case studies, in the population of which these case studies are representative, do we expect to observe a positive pluralistic ignorance score? When the quantity of interest is a function of both fixed effects and random effects in a mixed model, parametric bootstrap is an appropriate resampling method (Austin & Leckie, 2020). The models have EPI scores as the outcome, stakeholder category as a fixed effect, and random intercepts for case study. They were fitted using the lmer() function from the lme4 package (Bates et al., 2015). For each attitude (ESI and WTS), we performed parametric bootstrap resampling with 99,999 iterations using the bootMer() function. Parametric bootstrapping simulates new datasets from the fitted model by:

- Drawing new random effects from their estimated distribution
- Drawing new residuals from their estimated distribution
- Constructing simulated response values based on the model structure
- Refitting the model to each simulated dataset

For each bootstrap iteration, we calculated the expected proportion of case studies that would show positive EPI scores. This was computed by:

- Extracting the mean fixed effect prediction for each case study
- For each case study, calculating the probability that a random intercept drawn from the estimated random effects distribution, when added to the fixed effect, would yield a positive value (using the normal CDF with the estimated standard deviation of random intercepts)
- Averaging these probabilities across all case studies to obtain the expected proportion

This quantity represents the expected proportion of case studies where pluralistic ignorance would be present, accounting for uncertainty in the random effects. We calculated 95% confidence intervals for this proportion using the 2.5th and 97.5th percentiles of the bootstrap distribution. We report the point estimate calculated from the original fitted model; the bootstrapped means showed little bias.

### Ordinal regression models for the associations between personal and perceived attitudes (H3 and beyond)

We used cumulative link models to examine whether perceived norms predicted individual attitudes while accounting for stakeholder category and case study. For both ESI and WTS, we fit a series of five models with increasing complexity:

**Model 1:** attitude ~ perceived norm

- Basic model with perceived norm as the sole predictor

**Model 2:** attitude ~ perceived norm + stakeholder category

- Adds stakeholder category as a fixed effect (the reference category is the stakeholder group with the largest sample size for that outcome)

**Model 3:** attitude ~ perceived norm + stakeholder category + (1 | case study)

- Adds random intercepts for case study

**Model 4:** attitude ~ perceived norm + stakeholder category + (1 | case study) + perceived norm : stakeholder category

- Adds interaction between perceived norm and stakeholder category, allowing the norm-attitude relationship to vary across stakeholder categories

**Model 5:** attitude ~ perceived norm + stakeholder category + (1 | case study) + (0 + perceived norm | case study)

- Replaces the stakeholder category interaction with random slopes for perceived norm by case study, allowing the norm-attitude relationship to vary across sites

Models were fitted using the cumulative link mixed model framework with a logit link function (proportional odds model). For Models 1 and 2 (fixed effects only), we used the clm() function from the ordinal package in R (Christensen, 2025). For Models 3-5 (including random effects), we used clmm() with Laplace approximation. The proportional odds assumption was tested by applying the nominal_test() function to Model 1.

For each model, we report regression coefficients for the perceived norm predictor converted to odds ratios through exponentiation, and 95% confidence intervals obtained from the profile method (Venzon & Moolgavkar, 1988).

### Cross-attitude model

We fit an additional model to examine whether how ESI and perceived WTS norms independently predicted WTS, which was as WTS Model 3 with the addition of ESI (which was treated as numeric and standardised):

**Model 6:** WTS ~ WTS perceived norm + ESI + stakeholder category + (1 | case study)