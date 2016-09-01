####################################################################
## Stan code for hazard assessment using Species Sensitivity Distributions
## A demonstration of Bayesian Evidence Synthesis
## and Robust Bayesian Analysis using rough
## sensitivity analysis on priors or data
## Extra exercise - add a simple QSAR to fill data gap
## Ullrika Sahlin 2016
###################################################################

require('rstan')
require('shinystan')
###########################################################
## Load codes and some functions
source('stan_hazardassessment.R')
source('ssdcode.R')

###########################################

## Generate artificial toxicity data from a SSD with mu and sigma
ssd_data <- generate_data(mu = 2,sigma = 1,K = 20,s_sizes=c(1:3),seed = 1975)

## Extract data to learn the SSD based on estimates of toxicity and
## information on study specific uncertainty in these estimates
## Include hyper parameters for priors
p = 0.05 #is the largest fraction species we allow to be affected by the substance
## based on p we derive lambda which is the pth percentile of the normal distribution (which we use for the SSD)
dat = list(N=sum(ssd_data$s),K=ssd_data$K,y=unlist(ssd_data$y),s=ssd_data$s,
             sigma_y=unlist(ssd_data$sigma_y),#toxicity data
             mean_mu = -5, sig_mu = 10,upper_sigma = 10,#hyper priors
              w = unlist(sapply(ssd_data$s,FUN=function(x) rep(1/x,x))),
             lambda=qnorm(p))#quantile for calculating loss function

model = stan(model_name="model", model_code = code_ssd, data=dat,
               iter = 5000, chains = 4, verbose = FALSE)

print(model)
traceplot(model,c('mu','sigma'))
ppcheck_plot_toxicity(stanmodel=model,ssd_data)

ext = extract(model,c('LHCp','mu','sigma'))
pairs(ext)

#launch_shinystan(model)
out <- find_LHCp_hat(stanmodel=model,LINEX_alpha = 1.5,p = p,
                     dat=dat,makeplot=TRUE,how='cdf')
out

########################################
### Sensitivity analysis with respect to prior and weights on data


out <- robust_hazardassessment(p = 0.05,LINEX_alpha = 0.15,
                                   mean_mu = -5, sig_mu = 10,upper_sigma=10,
                                   w = unlist(sapply(ssd_data$s,FUN=function(x) rep(1/x,x))),
                                   ssd_data,
                                   iter=5000,chains=3)

out2 <- robust_hazardassessment(p = 0.05,LINEX_alpha = 0.15,
                               mean_mu = 0, sig_mu = 2,upper_sigma=3,
                               w = unlist(sapply(ssd_data$s,FUN=function(x) rep(1/x,x))),
                               ssd_data,
                               iter=5000,chains=3)

out
out2


################################################
## Learn the SSD based on point estimates of toxicity and
## information on study specific uncertainty
dat_1 = list(N=sum(ssd_data$s),K=ssd_data$K,y=unlist(ssd_data$y),s=ssd_data$s,
             sigma_y=unlist(ssd_data$sigma_y))
model_1 = stan(model_name="model_1", model_code = code_1, data=dat_1,
               iter = 10000, chains = 2, verbose = FALSE)

print(model_1)
ppcheck_plot_toxicity(stanmodel=model_1,ssd_data)

yobs = y=unlist(ssd_data$y)
launch_shinystan(model_1)
