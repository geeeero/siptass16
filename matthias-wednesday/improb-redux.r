# helper function to check that two matrices are equal within tolerance
.stopifnotalmostequal = function(a, b, tol=1e-10) {
  stopifnot(a + tol >= b)
  stopifnot(b + tol >= a)
}

###############################################################################
# calculating expectations and functions of expectations
###############################################################################

# Return a function that calculates the conditional expectation of
# random variables with respect to a set of probability mass
# functions. Giving it an event, this function returns another
# function which takes any set of random variables (specified as a
# simple vector), and will return a matrix, with one row per random
# variable, and one column per conditional probability mass function.
getconditionalexpectationsfunc = function(possibsize, pmfs, tol=1e-10) {
  pmfmatrix = matrix(pmfs, nrow=possibsize)
  # all columns sum to one?
  .stopifnotalmostequal(colSums(pmfmatrix), 1, tol=tol)
  # all entries are non-negative?
  stopifnot(pmfmatrix + tol >= 0)
  function(event) {
    stopifnot(length(event) == possibsize)
    cpmfmatrix = apply(
      # m[rows,,drop=FALSE] picks the rows but makes sure the result
      # remains a matrix
      pmfmatrix[which(event),,drop=FALSE], 2,
      function(col) { col / sum(col) })
    # return function that multiplies the matrices
    function(rvars) {
      matrix(rvars, ncol=sum(event), byrow=TRUE) %*% cpmfmatrix
    }
  }
}

# Convenience wrapper function to calculate unconditional expectations,
getexpectationsfunc = function(possibsize, pmfs, tol=1e-10) {
  event = rep(TRUE, possibsize)
  getconditionalexpectationsfunc(possibsize, pmfs, tol=tol)(event)
}

# Return a function which calls getexpectations, and then applies a
# vector function on each row (i.e. to calculate minimum expectation,
# maximum expectation, etc.).
getexpectationsapplyfunc = function(getexpectations, func) {
  function(rvarvalues) {
    apply(getexpectations(rvarvalues), 1, func)
  }
}

# Return a function which applies a pairwise vector function to
# expectations (e.g. to determine maximality, interval dominance,
# etc.). The result is a square matrix, with an entry for each pair of
# random variables.
getexpectationsapplyfunc2 = function(getexpectations, func) {
  function(rvarvalues) {
    expectations = getexpectations(rvarvalues)
    apply(expectations, 1, function(exp1) {
      apply(expectations, 1, function(exp2) { func(exp1, exp2) }) })
  }
}

# Return a function which evaluates the lower prevision.
getlowerprevisionsfunc = function(getexpectations) {
  getexpectationsapplyfunc(getexpectations, min)
}

# Return a function which evaluates the upper prevision.
getupperprevisionsfunc = function(getexpectations) {
  getexpectationsapplyfunc(getexpectations, max)
}

# Return a function which evaluates the "Hurwicz" prevision.
gethurwiczprevisionsfunc = function(getexpectations, optimism) {
  .hurwicz = function(expectations) {
    optimism * max(expectations) + (1 - optimism) * min(expectations)
  }
  getexpectationsapplyfunc(getexpectations, .hurwicz)
}

# Return a function which evaluates the precise prevision.
getpreciseprevisionsfunc = function(getexpectations) {
  .precise = function(expectations) {
    stopifnot(length(expectations) == 1)
    expectations[1]
  }
  getexpectationsapplyfunc(getexpectations, .precise)
}

###############################################################################
# optimality criteria based on the full distribution (no data)
###############################################################################

# Return a function which tells you which random variables are
# Gamma-maxi-"something", where "something" is any function from random
# variables to values (for example a function created with
# getlowerprevisionsfunc, getupperprevisionsfunc, or
# gethurwiczprevisionsfunc).
isgammamaxisomethingfunc = function(getsomethings, tol=1e-10) {
  function(rvarvalues) {
    somethings = getsomethings(rvarvalues)
    somethings >= (max(somethings) - tol)
  }
}

# Partial orders based on expectation values.
rbayescompare = function(exp1, exp2, tol=1e-10) {
  all(exp1 > exp2 + tol)
}

waldcompare = function(exp1, exp2, tol=1e-10) {
  all(exp1 + tol >= exp2) && any(exp1 > exp2 + tol)
}

intervalcompare = function(exp1, exp2, tol=1e-10) {
  min(exp1) > max(exp2) + tol
}

# Return maximal elements based on partial order.
# This is *not* the most efficient implementation, but it is simple.
ismaximalfunc = function(getexpectations, compare) {
  function(rvarvalues) {
    getcomparisonmatrix = getexpectationsapplyfunc2(
      getexpectations, compare)
    apply(getcomparisonmatrix(rvarvalues), 1, function(row) { !any(row) })
  }
}

# Return a function which tells you which random variables are Bayes
# maximal (i.e. Bayes maximal with respect to some element of the
# credal set specified).
isrbayesadmissiblefunc = function(getexpectations, tol=1e-10) {
  function(rvarvalues) {
    expectations = getexpectations(rvarvalues)
    # for each probability mass function (column), identify the random
    # variables (rows) that have maximal expectation
    maxexpectations = apply(expectations, 2, function(col) { col >= (max(col) - tol) })
    # for each random variable (row), check if it had maximal
    # expectation with respect to any probability mass function
    apply(maxexpectations, 1, function(row) { any(row) })
  }
}

###############################################################################
# optimality criteria that also take data into account (strategies)
###############################################################################

################################################################################
# Bayes strategies
################################################################################

# Apply Bayes's theorem.
getposteriorpmf = function(paramsize, priorpmf, likelihoodpmf) {
  likelihoodmat = matrix(likelihoodpmf, nrow=paramsize, byrow=TRUE)
  c(apply(likelihoodmat, 2, function(col) {
    prod = col * priorpmf
    prod / sum(prod)
  }))
}

isbayesstrategy = function(paramsize, posteriorpmf, utility) {
  # Repeat Bayes analysis for each posterior.
  posteriormat = matrix(posteriorpmf, ncol=paramsize, byrow=TRUE)
  # utility matrix lists outcome per decision, so it is already in the
  # right format to be seen as a set of random variables
  rvars = utility
  apply(posteriormat, 1, function(pmf) {
    getposteriorexpectation = getexpectationsfunc(paramsize, pmf)
    getpreciseprevision = getpreciseprevisionsfunc(getposteriorexpectation)
    isbayes = isgammamaxisomethingfunc(getpreciseprevision)
    isbayes(rvars)
    })
}

################################################################################
# Wald's optimality
################################################################################

# list all strategies in a matrix (one row per strategy, enumerating a
# decision for each possible data value)
getstrategies = function(datasize, decisionsize) {
  args = lapply(1:datasize, function(x) { 1:decisionsize })
  result = as.matrix(expand.grid(args))
  dimnames(result) <- NULL
  result
}

# Get Wald's expected utility for any given strategy.
getwaldutilityfunc = function(datasize, likelihoodpmf, utility) {
  paramsize = length(likelihoodpmf) %/% datasize
  decisionsize = length(utility) %/% paramsize
  utilitymat = matrix(utility, nrow=decisionsize, byrow=TRUE)
  likelihoodmat = matrix(likelihoodpmf, nrow=paramsize, byrow=TRUE)
  getwaldexpectations = sapply(
    1:paramsize, function(param) {
      getexpectationsfunc(datasize, likelihoodmat[param,])
    })

  function(strategy) {
    stopifnot(datasize == length(strategy))
    sapply(
      1:paramsize,
      function(param) {
        rvar = sapply(strategy, function(decision) { utilitymat[decision,param] })
        getwaldexpectations[[param]](rvar)[1,1]
      })
  }
}

# Get all of Wald's expected utilities.
getwaldutilities = function(datasize, likelihoodpmf, utility) {
  paramsize = length(likelihoodpmf) %/% datasize
  decisionsize = length(utility) %/% paramsize
  getwaldutility = getwaldutilityfunc(datasize, likelihoodpmf, utility)
  strategies = getstrategies(datasize, decisionsize)
  c(apply(strategies, 1, getwaldutility))
}

# Find all Wald admissible strategies.
iswaldstrategy = function(datasize, likelihoodpmf, utility) {
  paramsize = length(likelihoodpmf) %/% datasize
  decisionsize = length(utility) %/% paramsize
  getexpectations = function(rvarvalues) {
    matrix(rvarvalues, ncol=paramsize, byrow=TRUE)
  }
  ismaximal = ismaximalfunc(getexpectations, waldcompare)
  ismaximal(getwaldutilities(datasize, likelihoodpmf, utility))
}
