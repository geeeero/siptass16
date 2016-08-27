###################################################
### code chunk number 1: IntroToBayes2.Rnw:84-90
###################################################
x = matrix(rnorm(1e5),nrow=100)
par(mfrow=c(1,2))
hist(colMeans(x), xlab="", main="mean(x)")
abline(v=0,col=2)
hist(apply(x,2,quantile,0.975), xlab="", main="quantile(x,0.975)")
abline(v=qnorm(0.975),col=2)

###################################################
### code chunk number 2: IntroToBayes2.Rnw:119-124
###################################################
##load the RSTAN library
require(rstan)
require(lattice) # for densityplot function
##set number of digits in output
options(digits = 3) 


###################################################
### code chunk number 3: IntroToBayes2.Rnw:183-194
###################################################
model_0 = "
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
model {
  tau ~ gamma(1e-3,1e-3);
  beta1 ~ normal(0,1e2);
  beta2 ~ normal(0,1e2);
  for(i in 1:N) {
    y[i] ~ normal(beta1 + beta2*x[i], sqrt(1/tau));
  }
}"


###################################################
### code chunk number 4: IntroToBayes2.Rnw:198-200
###################################################
data = list()
data$N = 50
data$x = rnorm(data$N)+30
data$y = 3 + 5*data$x + rnorm(data$N,sd=1/sqrt(100))

###################################################
### code chunk number 5: IntroToBayes2.Rnw:207-208
###################################################
# not needed for RSTAN?
#inits <- list(beta1=100, beta2=100, tau=1e4)


###################################################
### code chunk number 6: IntroToBayes2.Rnw:212-215
###################################################
fit_0 = stan(model_code=model_0, data=data, iter=1000, chains=4)

print(fit_0)
plot(fit_0)

###################################################
### code chunk number 7: IntroToBayes2.Rnw:219-221
###################################################
samples_0 = extract(fit_0, c("tau", "beta1", "beta2"))

###################################################
### code chunk number 8: IntroToBayes2.Rnw:227-229
###################################################
plot(samples_0$tau, type="l")
plot(samples_0$beta1, type="l")
plot(samples_0$beta2, type="l")
pairs(samples_0)

###################################################
### code chunk number 9: IntroToBayes2.Rnw:234-235
###################################################
### MATTHIAS TO GERO: not sure how to do
# xyplot(samples)


###################################################
### code chunk number 10: IntroToBayes2.Rnw:240-241
###################################################
densityplot(samples_0$tau)
densityplot(samples_0$beta1)
densityplot(samples_0$beta2)

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


