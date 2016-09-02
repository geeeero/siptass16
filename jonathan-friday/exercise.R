####################################################################
# Exercise Friday morning Jonathan 
####################################################################

setwd("./jonathan-friday")
install.packages("linprog")
library(linprog)
library(ggplot2)

dataj <- read.csv2("data.txt", sep="", dec=".", header=F)
dataj <- data.frame(x=dataj[,1], y=dataj[,2])
qplot(data=dataj, x=x, y=y)

n <- dim(dataj)[1]
#pd <- 2 # polynomial degree of bounding functions
basisx <- cbind(dataj$x^2, dataj$x, 1)

ylow <- function(pl, pu, basis) # pl, pu are lower and upper vector of coefficients for polynomial basis
  (1/2 * (basis - abs(basis))) %*% pu + (1/2 * (basis + abs(basis))) %*% pl
yupp <- function(pl, pu, basis)
  (1/2 * (basis + abs(basis))) %*% pu + (1/2 * (basis - abs(basis))) %*% pl
dely <- function(pl, pu, basis)
  basis %*% (pu - pl) 

# phatl, phatu = argmin u,v of E_x[dely(basis, u, v): ylow() <= y_i <= yupp(), v-u >= 0]
solveLP()


#