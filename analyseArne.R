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

# Calculating the EPI difference score
meanSupport <- sum(dF$support>=4) / length(dF$support)
dF$diff <- meanSupport - ( dF$perceivedNormSupport / 10 )

# The EPI difference score is significantly different to zero
permutation_sign_test(dF$diff)

# There is a correlation between support and perceived support
cor.test(dF$support, dF$perceivedNormSupport)