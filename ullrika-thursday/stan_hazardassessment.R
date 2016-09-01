####################################################################
## Stan code for hazard assessment using Species Sensitivity Distributions
## A demonstration of Bayesian Evidence Synthesis
## and Robust Bayesian Analysis using rough
## sensitivity analysis on priors or data
## Extra exercise - add a simple QSAR to fill data gap
## Ullrika Sahlin 2016
###################################################################

# models included in this script:
# code_1 Learn a SSD based on toxcity data from K species
#         with info on estimate specific uncertainty
#         We use this one to get started

# code_ssd Learn a SSD based on toxcity data from K species
#         with info on estimate specific uncertainty
#         the code prepare for assessment of the hazardous concentration


# code_qsar_ssd Same as code 3 but fill data gap on a toxicity using 
#         a QSAR (simple regression)  

###########################################
## Learn a SSD based on toxicity data on K species
## With information on study specific uncertainty
## A simple version than the one used in hazard assessment
code_1 = "
data {
int<lower=0> N; // # observations
int<lower=0> K; // # of groups
vector[N] y; // observations
int s[K]; // group sizes
vector[N] sigma_y; // known study specific uncertainty
}
parameters {
real mu;
real<lower=0> sigma;
vector[K] alpha;
}
model {
int pos;
mu ~ uniform(-5,10);
sigma ~ uniform(0,10);

pos <- 1;
for (k in 1:K) {
alpha[k] ~ normal(0,sigma);
segment(y, pos, s[k]) ~ normal(mu+alpha[k],segment(sigma_y, pos, s[k]));
pos <- pos + s[k];
}
}
generated quantities { 
real y_rep[N];
int pos;
pos <- 1;
for (k in 1:K) {
for (r in 1:s[k]) {
y_rep[pos] <- normal_rng(mu+alpha[k],sigma_y[pos]);
pos <- pos + 1;
}
}
}
"
###########################################
## Learn a SSD based on toxicity data on K species
## With information on study specific uncertainty
## Prepare to assess hazardous concentration
code_ssd = "
functions { 
real calc_LHCp(real lambda1, real mu1, real sigma1){
return mu1+lambda1*sigma1; // function to calculate loss
}
}

data {
int<lower=0> N; // # observations
int<lower=0> K; // # of groups
vector[N] y; // observations
int s[K]; // group sizes
vector[N] sigma_y; // known group specific uncertainty

// hyper parameters for the priors
real mean_mu;//mean of mu
real<lower=0> sig_mu;// standard deviation of mu
real<lower=0> upper_sigma;// right value of the uniform for the sigma
vector[N] w;//weights on each study
// info to hazard assessment
real lambda; //quantile in the SSD - stan does not have a function to calculate this from p
}

transformed data {
vector[N] sigma_y_w;//weighted sigma_y
for(i in 1:N) 
sigma_y_w[i] = sigma_y[i] / w[i];
}

parameters {
real mu;
real<lower=0> sigma;
vector[K] alpha;
}

//transformed parameters {
//vector[K] t;
//t = mu + alpha; // we reparametrise from t to alpha to reduce correlated parameters
//}

model {
int pos;
mu ~ uniform(mean_mu,sig_mu);
sigma ~ uniform(0,upper_sigma);

pos <- 1;
for (k in 1:K) {
alpha[k] ~ normal(0,sigma);
segment(y, pos, s[k]) ~ normal(mu + alpha[k],segment(sigma_y_w, pos, s[k]));
pos <- pos + s[k];
}
}

generated quantities { 
real y_rep[N];
real LHCp;
int pos;
pos <- 1;
for (k in 1:K) {
for (r in 1:s[k]) {
y_rep[pos] <- normal_rng(mu + alpha[k],sigma_y_w[pos]);
pos <- pos + 1;
}
}
LHCp <- calc_LHCp(lambda,mu,sigma);
}

"



