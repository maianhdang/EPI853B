###Homework 3
(due Monday 9/26/2016, submit va e-mail to gustavoc@msu.edu )

1. Write a function that implements a backfitting algorithm for OLS. Name your function `LS.backfit`, your function must: have, `y` and `X`
arguments for the response and incidence matrix, and options for tolerance for convergengence (`tol`) and for adding an intercept (`int`). To get full credit your function must pass all the tests listed below.


Hint: if you center, the intercept changes, you can get the original intercept (the one you would obtain if you do not center) by adding the 
intercept of the centered model plus `b[-1]'xBar`, where `b` are the regression coefficients you obtained and `xBar=colMeans(X)` are the means for
all the predictors.

**Simulation**

```R
 n=1000
 p=20
 Z=matrix(nrow=n,ncol=p,data=runif(n*p))
 b=rgamma(shape=1,rate=1,n=p)*sample(c(-1,1),size=p,replace=T)
 signal=Z%*%b+25
 error=rnorm(n=n,sd=sd(signal))
 y=error+signal
```

**Tests**
```R
#Test 1
 TEST=rep(NA,3)
 bHat_1=coef(lm(y~Z))
 bHat_2=LS.backfit(y=y,X=Z)
 TEST[1]=max(abs(bHat_1-bHat_2))<1e-5

#Test 2
 bHat_1=coef(lm(y~Z))
 bHat_2=LS.backfit(y=y,X=Z,int=T)
 TEST[2]=max(abs(bHat_1-bHat_2))<1e-5 

#Test 3
 bHat_1=coef(lm(y~Z-1))
 bHat_2=LS.backfit(y=y,X=Z,int=F,tol=1e-7)
 TEST[3]=max(abs(bHat_1-bHat_2))<1e-5  

```
