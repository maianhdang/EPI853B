###Homework 3
(due Monday 9/26/2016, submit va e-mail to gustavoc@msu.edu )

1. Write a function that implements a backfitting algorithm for OLS. Name your function `LS.backfit`, your function must: have, 'y' and 'X' 
arguments for the response and incidence matrix, and options for adding and intercept and centering (by default `int=TRUE` and `center=TRUE`).
The return value must be eual to that of `lm`. To get full credit your function must pass all the tests listed below.


Hint: if you center, the intercept changes, you can get the original intercept (the one you would obtain if you do not center) by adding the 
intercept of the centered model plus `b[-1]'xBar`, where `b` are the regression coefficients you obtained and `xBar=colMeans(X)` are the means for
all the predictors.

