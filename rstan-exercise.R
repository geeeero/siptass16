###################################################
# code for RStan exercise ISIPTA summer school 2016
###################################################

##load the RSTAN library
require(rstan)
# recommended in startup message
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

require(lattice) # for densityplot function
##set number of digits in output
options(digits = 3) 

model1 <- "
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real<lower=0> tau;
  real beta1;
  real beta2;
}
transformed parameters {
  real<lower=0> sigma;
  sigma=sqrt(1/tau);
}
model {
  tau ~ gamma(1e-3,1e-3);
  beta1 ~ normal(0,1e2);
  beta2 ~ normal(0,1e2);
  //for(i in 1:N) {
  //  y[i] ~ normal(beta1 + beta2*x[i], sigma);
  //}
  y ~ normal(beta1 + beta2*x, sigma);
}"

data <- list()
data$N <- 50
data$x <- rnorm(data$N)+30
data$y <- 3 + 5*data$x + rnorm(data$N,sd=1/10)

fit1 <- stan(model_code=model1, data=data, iter=1000, chains=4)
fit1a <- stan(model_code=model1, data=data, iter=1000, chains=4, thin=20)
plot(fit1)
print(fit1)
print(fit1a)

samples1 = extract(fit1, c("beta1", "beta2", "sigma"))
plot(samples1$sigma, type="l")
plot(samples1$beta1, type="l")
plot(samples1$beta2, type="l")

stan_trace(fit1, pars = c("beta1", "beta2", "sigma"), inc_warmup = TRUE, nrow = 3)
stan_trace(fit1, pars = c("beta1", "beta2", "sigma"), inc_warmup = FALSE, nrow = 3)

#densityplot(samples1$beta1)
#densityplot(samples1$beta2)
#densityplot(samples1$sigma)
stan_dens(fit1, pars = c("beta1", "beta2", "sigma"), inc_warmup = FALSE, ncol = 3)
stan_dens(fit1, pars = c("beta1", "beta2", "sigma"), inc_warmup = TRUE, ncol = 3)

# autocorrelation plot to look at thinning
stan_ac(fit1, pars = c("beta1", "beta2", "sigma"), ncol = 3)
samples1a <- As.mcmc.list(fit1, pars = c("beta1", "beta2", "sigma"))
coda::acfplot(samples1a, aspect="fill", lag.max=25)

# correlation problem is present
pairs(samples1, pch=19, cex=0.1)
pairs(fit1, pars = c("beta1", "beta2", "sigma"))
stan_scat(fit1, pars = c("beta1", "beta2"))
stan_scat(fit1, pars = c("beta1", "sigma"))

# pdfs for presentation
pdf("trace1w.pdf", width=8, height=5)
stan_trace(fit1, pars = c("beta1"), inc_warmup = TRUE)
dev.off()

pdf("trace1.pdf", width=8, height=5)
stan_trace(fit1, pars = c("beta1"), inc_warmup = FALSE)
dev.off()

pdf("acf1.pdf", width=4, height=2.5)
stan_ac(fit1, pars = c("beta1"))
dev.off()

# change the model to have a transformed data block
model2 <- "
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
transformed data {
  vector[N] xcentered;
  xcentered=x-mean(x);
}
parameters {
  real<lower=0> tau;
  real beta1c;
  real beta2;
}
transformed parameters {
  real<lower=0> sigma;
  real beta1;
  sigma=sqrt(1/tau);
  beta1=beta1c-beta2*mean(x);
}
model {
  tau ~ gamma(1e-3,1e-3);
  beta1c ~ normal(0,1e2);
  beta2 ~ normal(0,1e2);
  y ~ normal(beta1c + beta2*xcentered, sigma);
}"

fit2 <- stan(model_code=model2, data=data, iter=1000, chains=4)
plot(fit2)
print(fit2)

stan_trace(fit2, pars = c("beta1", "beta2", "sigma"), inc_warmup = TRUE, nrow = 3)
stan_trace(fit2, pars = c("beta1", "beta2", "sigma"), inc_warmup = FALSE, nrow = 3)

stan_dens(fit2, pars = c("beta1", "beta2", "sigma"), inc_warmup = FALSE, ncol = 3)
stan_ac(fit2, pars = c("beta1", "beta2", "sigma"), ncol = 3)    # much better now!
pairs(fit2, pars = c("beta1", "beta1c", "beta2"))               # much better now!

#