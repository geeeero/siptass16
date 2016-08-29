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

model0 <- "
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

fit0 <- stan(model_code=model0, data=data, iter=1000, chains=4)
fit0a <- stan(model_code=model0, data=data, iter=1000, chains=4, thin=20)
plot(fit0)
print(fit0)
print(fit0a)

samples0 = extract(fit0, c("beta1", "beta2", "sigma"))
plot(samples0$sigma, type="l")
plot(samples0$beta1, type="l")
plot(samples0$beta2, type="l")

stan_trace(fit0, pars = c("beta1", "beta2", "sigma"), inc_warmup = TRUE, nrow = 3)
stan_trace(fit0, pars = c("beta1", "beta2", "sigma"), inc_warmup = FALSE, nrow = 3)

#densityplot(samples0$beta1)
#densityplot(samples0$beta2)
#densityplot(samples0$sigma)
stan_dens(fit0, pars = c("beta1", "beta2", "sigma"), inc_warmup = FALSE, ncol = 3)
stan_dens(fit0, pars = c("beta1", "beta2", "sigma"), inc_warmup = TRUE, ncol = 3)

# autocorrelation plot to look at thinning
stan_ac(fit0, pars = c("beta1", "beta2", "sigma"), ncol = 3)
samples1 <- As.mcmc.list(fit_0, pars = c("beta1", "beta2", "sigma"))
coda::acfplot(samples1, aspect="fill", lag.max=100)

# correlation problem is present
pairs(samples0, pch=19, cex=0.1)
pairs(fit0, pars = c("beta1", "beta2", "sigma"))
stan_scat(fit0, pars = c("beta1", "beta2"))
stan_scat(fit0, pars = c("beta1", "sigma"))

# change the model to have a transformed data block


#### MATTHIAS TO GERO: STOPPED CONVERSION HERE!!! ####
quit()


###################################################
### code chunk number 11: IntroToBayes2.Rnw:252-258
###################################################
model <- jags.model(textConnection(model.txt), 
                    data=data, inits=function(){
                      list(beta1=rnorm(1,0,100),
                           beta2=rnorm(1,0,100),
                           tau=runif(1,0,100))
                    }, n.chains=3, quiet=TRUE)


###################################################
### code chunk number 12: IntroToBayes2.Rnw:260-262
###################################################
samples <- coda.samples(model, 
   c("beta1","beta2","tau"), n.iter=1e4)


###################################################
### code chunk number 13: IntroToBayes2.Rnw:267-269
###################################################
par(mar=c(3.3,3.3,1.5,1), mgp=c(2,1,0))
plot(samples)


###################################################
### code chunk number 14: IntroToBayes2.Rnw:285-287 (eval = FALSE)
###################################################
## model <- jags.model(...)
## update(model, n.iter=5e3)


###################################################
### code chunk number 15: IntroToBayes2.Rnw:290-292
###################################################
samples.BI <- coda.samples(model, 
   c("beta1","beta2","tau"), n.iter=1e4)


###################################################
### code chunk number 16: IntroToBayes2.Rnw:297-298
###################################################
xyplot(samples.BI)


###################################################
### code chunk number 17: IntroToBayes2.Rnw:303-304
###################################################
densityplot(samples.BI)


###################################################
### code chunk number 18: IntroToBayes2.Rnw:317-318
###################################################
effectiveSize(samples.BI)


###################################################
### code chunk number 19: IntroToBayes2.Rnw:321-322
###################################################
sapply(samples.BI, effectiveSize)


###################################################
### code chunk number 20: IntroToBayes2.Rnw:329-330
###################################################
acfplot(samples.BI, aspect="fill", lag.max=100)


###################################################
### code chunk number 21: IntroToBayes2.Rnw:339-341 (eval = FALSE)
###################################################
## samples.thin <- coda.samples(model, 
##    c("beta1","beta2","tau"), n.iter=1e4*1e2, thin=1e2)


###################################################
### code chunk number 22: IntroToBayes2.Rnw:345-346
###################################################
acfplot(samples.BI, aspect="fill", lag.max=100, thin=50)


###################################################
### code chunk number 23: IntroToBayes2.Rnw:359-360 (eval = FALSE)
###################################################
## model <-  jags.model(..., n.adapt=1e5)


###################################################
### code chunk number 24: IntroToBayes2.Rnw:368-369
###################################################
pairs(as.matrix(samples.BI), pch=19, cex=.1)


###################################################
### code chunk number 25: IntroToBayes2.Rnw:389-401
###################################################
model2.txt <- "model {
  #priors
  tau ~ dgamma(1e-3,1e-3)
  beta1 ~ dnorm(0,1e-4)
  beta2 ~ dnorm(0,1e-4)

  #likelihood
  x.bar <- mean(x)
  for(i in 1:N){
    y[i] ~ dnorm(beta1 + beta2*(x[i]-x.bar), tau)
  }
}"


###################################################
### code chunk number 26: IntroToBayes2.Rnw:409-415
###################################################
model2 <- jags.model(textConnection(model2.txt), 
                    data=data, inits=function(){
                      list(beta1=rnorm(1,0,100),
                           beta2=rnorm(1,0,100),
                           tau=runif(1,0,100))
                    }, n.chains=3, quiet=TRUE)


###################################################
### code chunk number 27: IntroToBayes2.Rnw:419-421
###################################################
samples2 <- coda.samples(model2, 
   c("beta1","beta2","tau"), n.iter=1e4)


###################################################
### code chunk number 28: IntroToBayes2.Rnw:427-428
###################################################
xyplot(samples2)


###################################################
### code chunk number 29: IntroToBayes2.Rnw:433-434
###################################################
densityplot(samples2)


###################################################
### code chunk number 30: IntroToBayes2.Rnw:439-440
###################################################
acfplot(samples2, aspect="fill", lag.max=100)


###################################################
### code chunk number 31: IntroToBayes2.Rnw:443-444
###################################################
effectiveSize(samples2)


###################################################
### code chunk number 32: IntroToBayes2.Rnw:454-456
###################################################
gelman.diag(samples)
gelman.diag(samples.BI)


###################################################
### code chunk number 33: IntroToBayes2.Rnw:468-470
###################################################
geweke.diag(samples.BI[[1]])
geweke.diag(samples2[[1]])


###################################################
### code chunk number 34: IntroToBayes2.Rnw:482-483
###################################################
heidel.diag(samples2[[1]])


###################################################
### code chunk number 35: IntroToBayes2.Rnw:496-498
###################################################
raftery.diag(samples.BI[[1]])
raftery.diag(samples2[[1]])


