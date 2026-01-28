## Setup ###################################################################################################

# To run at the command line, replace this with a function that calls setwd()
# to the local directory where you have the data and R files.
setwd_PROCOAST_pluralistic_ignorance_1()

dF <- read.table("arneData.tsv",
                 header = TRUE,
                 sep = "\t",
                 stringsAsFactors = FALSE)

# We just need permutation_sign_test()
source("analysisFunctions.R")

# Sample size (there is no missing data)
nrow(dF)

# Percent support
percentSupport <- sum(dF$support>=4) / length(dF$support) * 100
percentSupport

# Perceived support (percent)
dF$perceivedSupportPercent <- dF$perceivedNormSupport * 10
mean(dF$perceivedSupportPercent)

# Calculating the EPI difference score
dF$diff <- percentSupport - dF$perceivedSupportPercent

# The EPI difference score is significantly different to zero
permutation_sign_test(dF$diff)

# There is a correlation between support and perceived support
cor.test(dF$support, dF$perceivedNormSupport)