# author: Matthias C. M. Troffaes
# date: 8 Aug 2015
# license: GPLv3

# install.packages('e1071', dependencies=TRUE)
library(e1071)

# helper functions for type checking
###############################################################################

.stopifnotdataframe = function(x) {
  if (!is.data.frame(x))
    stop("expected a data frame")
}

.stopifnottable = function(x) {
  if (!is.table(x))
    stop("expected a table")
}
 
.stopifnotlistoftables = function(xs) {
  if (!is.list(xs)) stop("expected a list of tables")
  lapply(xs, function(x) { .stopifnottable(x) })
  NULL
}
 
.stopifnotnumericvalue = function(x) {
  if (!is.numeric(x)) stop("expected a single numeric value")
  if (length(x) != 1) stop("expected a single numeric value")
}
 
.stopifnotnumericvector = function(x) {
  if (!is.numeric(x)) stop("expected a numeric vector")
}

# functions for cross fold validation
# https://gist.github.com/mcmtroffaes/709908
###############################################################################

kfcv.sizes = function(n, k=10) {
  # generate sample sizes for k-fold cross validation on a data set of
  # size n
  
  # usage:
  #
  #   kfcv.sizes(n, k=...)
  #

  sizes = c()
  for (i in 1:k) {
    first = 1 + (((i - 1) * n) %/% k)
    last = ((i * n) %/% k)
    sizes = append(sizes, last - first + 1)
  }
  sizes
}

kfcv.testing = function(n, k=10) {
  # generate testing sample indices for k-fold cross validation on a
  # data set of size n
  
  # usage:
  #
  #   kfcv.testing(n, k=...)
  #

  indices = list()
  sizes = kfcv.sizes(n, k=k)
  values = 1:n
  for (i in 1:k) {
    # take a random sample of given size
    s = sample(values, sizes[i])
    # append random sample to list of indices
    indices[[i]] = s
    # remove sample from values
    values = setdiff(values, s)
  }
  indices
}

kfcv.classifier = function(data, attribs, class, classifier, k=10) {
  # run k-fold cross validation with an arbitrary classifier

  # usage:
  #
  #   kfcv.classifier(data, class, classifier, k=...)
  #
  # where data is the data frame (each column is an attribute, and
  # each row is an observation), class is the column index for the
  # attribute to be predicted

  # classifier$trainer is a function which takes a training set,
  # attribute column indices, and a class column index; it returns a
  # model (a data structure) that can be used by the tester function
  # described next

  # classifier$tester is a function which takes a trained model and a
  # single row from a test set; it returns list of numerical test
  # results (e.g. whether the classifier predicted correctly, utility
  # for misclassification, ...)

  colMeans(
    do.call(
      rbind.data.frame,
      lapply(
        kfcv.testing(dim(data)[1]),
        function(testingindices) {
          model = classifier$trainer(data[-testingindices,], attribs, class)
          do.call(
            rbind.data.frame,
            lapply(
              testingindices,
              function(rowid) { classifier$tester(model, data[rowid,]) }
              ))
        })
      ),
    na.rm=TRUE)
}

# compose classifiers from a list of classifiers
###############################################################################

classifier.composed = function(classifiers) {
  list(
    trainer=function(train, attribs, class) {
      # the model is simply a list of models, one for each classifier
      lapply(classifiers, function(cls) cls$trainer(train, attribs, class))
    },
    tester=function(model, testrow) {
      # apply each tester to each model
      do.call(
        c, lapply(
          1:length(model),
          function(i) (classifiers[[i]]$tester)(model[[i]], testrow)))
    }
    )
}

# trainer and tester functions for the e1071 naiveBayes classifier
###############################################################################

classifier.naive = function(laplace=0) {
  list(
    trainer=function(train, attribs, class) {
      list(
        model=naiveBayes(train[,attribs], train[,class], laplace=laplace),
        attribs=attribs,
        class=class)
    },
    tester=function(model, testrow) {
      accuracy =
        predict(model$model, testrow[,model$attribs]) == testrow[,model$class]
      list(acc=accuracy)
    }
    )
}

# our own implementation of the naive Bayes classifier
###############################################################################

# calculate empirical probability mass function
# rows = the data frame; each column is a variable, each row is a joint
#        observation
# colid = column of the variable of which to calculate probability mass function
# s1, s2, s3 = values to change the counts (see code)
.helper.empirical.prob.table = function(s1, s2, s3) {
  function(rows, colid) {
    .stopifnotdataframe(rows)
    counts = table(rows[colid])
    prob = (counts + s1) / (sum(counts + s2) + s3)
    stopifnot(all(levels(factor(rows[,colid])) %in% names(prob)))
    prob
  }
}

.empirical.prob.table = .helper.empirical.prob.table(0, 0, 0)

.laplace.empirical.prob.table = function(alpha=1)
  .helper.empirical.prob.table(alpha, alpha, 0)
 
.lower.empirical.prob.table = function(s=2)
  .helper.empirical.prob.table(0, 0, s)

.upper.empirical.prob.table = function(s=2)
  .helper.empirical.prob.table(s, 0, s)

# calculate empirical conditional probability mass function
# rows = the data frame; each column is a variable, each row is a joint
#        observation
# colid1 = column with values of the conditioning variable
# colid2 = column with values of the non-conditioning random variable
# returns a table where each row is a probability mass function
.helper.empirical.conditional.prob.table = function(s1, s2, s3) {
  function(rows, colid1, colid2) {
    .stopifnotdataframe(rows)
    counts = table(rows[c(colid1, colid2)])
    (counts + s1) / (apply((counts + s2), 1, sum) + s3)
  }
}

.empirical.conditional.prob.table =
  .helper.empirical.conditional.prob.table(0, 0, 0)

.laplace.empirical.conditional.prob.table = function(alpha=1)
  .helper.empirical.conditional.prob.table(alpha, alpha, 0)
 
.lower.empirical.conditional.prob.table = function(s=2)
  .helper.empirical.conditional.prob.table(0, 0, s)

.upper.empirical.conditional.prob.table = function(s=2)
  .helper.empirical.conditional.prob.table(s, 0, s)

# calculate joint probability p(c,a) = p(c) * p(a1|c) * ... * p(ak|c)
.classifier.naive.joint.prob = function(c.prob.table, a.prob.tables, a.colids, c.level, testrow) {
  # check input types and values
  .stopifnottable(c.prob.table)
  .stopifnotlistoftables(a.prob.tables)
  .stopifnotnumericvector(a.colids)
  .stopifnotdataframe(testrow)
  # debug: we should have one p(c|a) per attribute
  stopifnot(length(a.colids) == length(a.prob.tables))
  # calculate p(a1|c), p(a2|c), ..., p(ak|c)
  p.acs = sapply(
    1:length(a.colids),
    function(i) {
      colid = a.colids[i]
      a.level = testrow[1, colid]
      prob.table.a.given.c = a.prob.tables[[i]]
      prob.table.a.given.c[c.level, a.level]
    }
    )
  # calculate p(c)
  p.c = c.prob.table[c.level]
  # joint probability under naive assumption
  # p(c,a)=p(c) * p(a1|c) * ... * p(ak|c)
  p.c * prod(p.acs)
}

# calculate all class probabilities from a test row
.classifier.naive.joint.probs = function(c.prob.table, a.prob.tables, a.colids, testrow) {
  .stopifnottable(c.prob.table)
  .stopifnotlistoftables(a.prob.tables)
  .stopifnotnumericvector(a.colids)
  .stopifnotdataframe(testrow)
  # calculate p(c,a) for each class level c
  # (the a_i are stored in testrow, the class levels are c.levels)
  lapply(
      names(c.prob.table),
      function(c.level) {
        .classifier.naive.joint.prob(c.prob.table, a.prob.tables, a.colids, c.level, testrow)
      }
      )
}

# predict the class from a test row
.classifier.naive.predict = function(c.prob.table, a.prob.tables, a.colids, testrow) {
  .stopifnottable(c.prob.table)
  .stopifnotlistoftables(a.prob.tables)
  .stopifnotnumericvector(a.colids)
  .stopifnotdataframe(testrow)
  # calculate p(c,a) for each class level c
  # (the a_i are stored in testrow, the class levels are c.levels)
  c.levels = names(c.prob.table)
  p.cas = .classifier.naive.joint.probs(c.prob.table, a.prob.tables, a.colids, testrow)
  c.levels[which.max(p.cas)]
}

classifier.naive2 = function(alpha=0) {
    list(
      trainer=function(train, attribs, class) {
          list(
            attribs=attribs,
            class=class,
            c.prob.table=.laplace.empirical.prob.table(alpha)(train, class),
            a.prob.tables=lapply(
              attribs,
              function(attrib) {
                  .laplace.empirical.conditional.prob.table(alpha)(train, class, attrib)
              })
            )
      },
      tester=function(model, testrow) {
          predictedclass = .classifier.naive.predict(
            model$c.prob.table, model$a.prob.tables, model$attribs, testrow)
          list(acc=(predictedclass == testrow[1, model$class]))
      }
      )
}

# a simple credal naive Bayes classifier based on interval dominance
###############################################################################

classifier.credal = function(s=2) {
    list(
      trainer=function(train, attribs, class) {
          list(
            attribs=attribs,
            class=class,
            c.lprob.table=.lower.empirical.prob.table(s)(train, class),
            c.uprob.table=.upper.empirical.prob.table(s)(train, class),
            a.lprob.tables=lapply(
              attribs,
              function(attrib) {
                  .lower.empirical.conditional.prob.table(s)(train, class, attrib)
              }),
            a.uprob.tables=lapply(
              attribs,
              function(attrib) {
                  .upper.empirical.conditional.prob.table(s)(train, class, attrib)
              })
            )
      },
      tester=function(model, testrow) {
          lprobs = unlist(.classifier.naive.joint.probs(model$c.lprob.table, model$a.lprob.tables, model$attribs, testrow))
          uprobs = unlist(.classifier.naive.joint.probs(model$c.uprob.table, model$a.uprob.tables, model$attribs, testrow))
          c.levels = names(model$c.lprob.table)
          result = c.levels[uprobs >= max(lprobs)]
          accuracy = (testrow[1, model$class] %in% result)
          list(
            acc=accuracy,
            deter=(length(result) == 1),
            singleacc=if (length(result) == 1) accuracy else NA,
            setsize=if (length(result) != 1) length(result) else NA,
            setacc=if (length(result) != 1) accuracy else NA)
      }
      )
}

# reading the data
###############################################################################

getdata = function() {
  mydata = read.csv(
    "data.txt",
    na.strings="?",
    col.names=c("BIRADS","Age","Shape","Margin","Density","Severity"))
  # remove rows with missing items
  mydata = na.omit(mydata)
  # fix issues in the BIRADS column (all values meant to be between 1 and 5)
  mydata$BIRADS[mydata$BIRADS==0] = 1
  mydata$BIRADS[mydata$BIRADS==55] = 5
  # discretise age
  mydata$Age = cut(mydata$Age, breaks=c(0,45,55,75,Inf),labels=FALSE)
  # turn all columns into factors
  for (name in names(mydata)) mydata[,name] = factor(mydata[,name])
  # return selected columns in order
  mydata[,c("BIRADS", "Shape", "Margin", "Density", "Age", "Severity")]
}

# examples
###############################################################################

test.naive = function() {
    mammo = getdata()
    myclassifier = classifier.naive2(0)
    model = myclassifier$trainer(mammo, 1:5, 6)
    # correct case
    testrow = mammo[6,]
    print(testrow)
    print(myclassifier$tester(model, testrow))
    # incorrect case
    testrow = mammo[5,]
    print(testrow)
    print(myclassifier$tester(model, testrow))
}

test.credal = function() {
    # set up the classifier
    myclassifier = classifier.composed(
      list(classifier.naive2(0),
           classifier.naive2(1),  # naive with laplace correction
           classifier.credal(2)))

    # test with just first 30 observations to get some interesting effects
    mammo = getdata()[1:30,]
    # predict severity (column 6) from all attributes
    print(kfcv.classifier(mammo, 1:5, 6, myclassifier))
    # predict severity from all attributes except BIRADS (column 1)
    print(kfcv.classifier(mammo, 2:5, 6, myclassifier))
    # predict severity from BIRADS only
    print(kfcv.classifier(mammo, 1, 6, myclassifier))

    # test with full data
    mammo = getdata()
    # predict BIRADS from other attributes
    print(kfcv.classifier(mammo, 2:5, 1, myclassifier))
}
