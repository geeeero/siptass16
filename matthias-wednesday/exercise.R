# Matthias' exercise Wednesdaypmfs <- c(0.5, 0.5, 0.8, 0.2)
pmfs <- c(0.5, 0.5, 0.8, 0.2)
rvars <- c(440, 260, 420, 300, 370, 370)
getexpectations <- getexpectationsfunc(2, pmfs) # produces a function
getexpectations(rvars)

getlowerprevisions <- getlowerprevisionsfunc(getexpectations)
getlowerprevisions(rvars)

isgammamaximin <- isgammamaxisomethingfunc(getlowerprevisions)
isgammamaximin(rvars)

getupperprevisions <- getupperprevisionsfunc(getexpectations)
isgammamaximax <- isgammamaxisomethingfunc(getupperprevisions)
isgammamaximax(rvars)

isintervalmaximal <- ismaximalfunc(getexpectations, intervalcompare)
isintervalmaximal(rvars)

isrbayesmaximal <- ismaximalfunc(getexpectations, rbayescompare)
isrbayesmaximal(rvars)

isrbayesadmissible <- isrbayesadmissiblefunc(getexpectations)
isrbayesadmissible(rvars)

# adding a third distribution -> second option is now also optimal
pmfs <- c(0.5, 0.5, 0.8, 0.2, 0.65, 0.35)
getexpectations <- getexpectationsfunc(2, pmfs) # produces a function
isrbayesadmissible <- isrbayesadmissiblefunc(getexpectations)
isrbayesadmissible(rvars)
getexpectations(rvars)

# --------------------------------------------------------------------

setwd("./matthias-wednesday")
mammo <- getdata ()
myclassifier <- classifier.naive2(0) # argument = s parameter
model <- myclassifier$trainer(mammo, 1:5, 6) # train the classifier
testrow = mammo[6,]
print(testrow)
print(myclassifier$tester(model, testrow))
testrow = mammo[5,]
print(testrow )
print(myclassifier$tester(model, testrow))

myclassifier <- classifier.credal(2) # redal classifyer with s = 2
myclassifier$tester(model, testrow)

myclassifier = classifier.composed(
  list(classifier.naive2(0),
       classifier.naive2(1),
       classifier.credal(2)))

model <- myclassifier$trainer(mammo, 1:5, 6) # train the classifier
myclassifier$tester(model, testrow)

mammo = getdata()[1:30,]
print(kfcv.classifier(mammo, 1:5, 6, myclassifier))
print(kfcv.classifier(mammo, 2:5, 6, myclassifier))
print(kfcv.classifier(mammo, 1, 6, myclassifier))

mammo = getdata()
print(kfcv.classifier(mammo, 1:5, 6, myclassifier)) # predicting the expert score by the x-ray

# higher s
myclassifier <- classifier.credal(10) # redal classifyer with s = 2
print(kfcv.classifier(mammo, 1:5, 6, myclassifier))



#