# for exercise after 'application in system reliability' part

#install.packages("devtools")
#library("devtools")
#install_github("louisaslett/ReliabilityTheory")
library("ReliabilityTheory")
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(xtable)

tlength <- 50
tvec <- seq(0, 10, length.out=tlength)

## ex 1: single component
data1 <- rgamma(10,5,1)
plot(function(x) 1-pgamma(x, 5, 1), xlim=c(0, 10), ylim=c(0, 1))
plot(function(x) 1-plnorm(x, 1.5, 0.25), xlim=c(0, 10), add=TRUE, col=2)

# create one-component system 
onecomp <- graph.formula(s -- 1 -- t)
onecomp <- setCompTypes(onecomp, list("C1" = 1))
onecompss <- computeSystemSurvivalSignature(onecomp)

# precise prior
y0 <- 1-plnorm(tvec, 1.5, 0.25)
n0 <- 5
onecomppost <- nonParBayesSystemInferencePriorSets(tvec, onecompss, list(C1 = data1), n0, n0, y0, y0)
plot(tvec, y0, ylim=c(0,1), type="l")
lines(tvec, onecomppost$lower, col=2)

# vacuous prior
onecomppost <- nonParBayesSystemInferencePriorSets(tvec, onecompss, list(C1 = data1), 1, 5, 1e-5, 1-1e-5)
plot(tvec, y0, ylim=c(0,1), type="l")
lines(tvec, onecomppost$lower, col=2)
lines(tvec, onecomppost$upper, col=2)

# semi-informative prior: y between 0.9 and 0.8 at time 2, and between 0.6 and 0.3 at time 5
y0u <- rep(1-1e-5, tlength)
y0u[tvec > 2] <- 0.9
y0u[tvec > 5] <- 0.6
y0l <- rep(0.8, tlength)
y0l[tvec > 2] <- 0.3
y0l[tvec > 5] <- 1e-5
plot(tvec, y0l, type="l", col=4, ylim=c(0,1))
lines(tvec, y0u)
onecomppost <- nonParBayesSystemInferencePriorSets(tvec, onecompss, list(C1 = data1), 1, 2, y0l, y0u)
lines(tvec, onecomppost$lower, col=4, lty=2)
lines(tvec, onecomppost$upper, lty=2)
onecomppost <- nonParBayesSystemInferencePriorSets(tvec, onecompss, list(C1 = data1), 10, 20, y0l, y0u)
lines(tvec, onecomppost$lower, col=4, lty=3)
lines(tvec, onecomppost$upper, lty=3)

# informative prior
y0l <- 1-pgamma(tvec, 3, 4)
y0u <- 1-pgamma(tvec, 3, 2)
plot(tvec, y0l, type="l", col=4)
lines(tvec, y0u)
onecomppost <- nonParBayesSystemInferencePriorSets(tvec, onecompss, list(C1 = data1), 1, 10, y0l, y0u)
lines(tvec, onecomppost$lower, col=4, lty=2)
lines(tvec, onecomppost$upper, lty=2)
# prior and posterior imprecision
plot(tvec, y0u-y0l, type="l")
lines(tvec, onecomppost$upper-onecomppost$lower, lty=2)

## ex 2: system with 2 or 3 component types
bridge <- graph.formula(s -- 2:3 -- 4 -- 5:6 -- 1 -- t, 2 -- 5, 3 -- 6)
bridge <- setCompTypes(bridge, list("T1"=c(2,3,5,6), "T2"=c(4), "T3"=c(1)))

#

bottomlegend <- theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank())
rightlegend <- theme(legend.title = element_blank())

#