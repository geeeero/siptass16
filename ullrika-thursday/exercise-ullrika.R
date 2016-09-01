###################################################
# code for RStan exercise Ullrika
###################################################

setwd("./ullrika-thursday")
##load the RSTAN library
library(rstan)
# recommended in startup message
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("stan_hazardassessment.R")
data_ssd <- list()
fitssd <- stan(model_code=code_ssd, data=data_ssd, iter=1000, chains=4)




#