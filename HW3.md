###Homework 3
(due Monday 9/26/2016, submit va e-mail to gustavoc@msu.edu )

1. Write a function that implements a backfitting algorithm for OLS. Name your function `LS.backfit`, your function must: have, 'y' and 'X' 
arguments for the response and incidence matrix, and options for adding and intercept and centering (by default `int=TRUE` and `center=TRUE`).
The return value must be eual to that of `lm`. To get full credit your function must pass all the tests listed below.


Hint: if you center, the intercept changes, you can get the original intercept (the one you would obtain if you do not center) by adding the 
intercept of the centered model plus `b[-1]'xBar`, where `b` are the regression coefficients you obtained and `xBar=colMeans(X)` are the means for
all the predictors.


**Tests**. Your function must pass all these tests. Your grade will be the sum of tests that passed.

*Simulation*

```R
 n=150
 p=20
 X=matrix(nrow=n,ncol=p,data=runif(n*p))
 b=rgamma(shape=1,rate=1,n=p)*sample(c(-1,1),size=p,replace=T)
 signal=X%*%b+25
 error=rnorm(n=n,sd=sd(signal))
 y=error+signal
 TEST=rep(NA,8)  
```

*Test 1*.
```R
 bHat_1=coef(lm(y~X))
 bHat_2=LS.backfit(y=y,X=X)
 TEST[1]=max(abs(bHat_1-bHat_2))<1e-5
````

*Test 2*.
```R
 bHat_1=coef(lm(y~X))
 bHat_2=LS.backfit(y=y,X=X,int=T)
 TEST[2]=max(abs(bHat_1-bHat_2))<1e-5  
```

*Test 3*.
```R
 bHat_1=coef(lm(y~X))
 bHat_2=LS.backfit(y=y,X=X,int=F)
 TEST[3]=max(abs(bHat_1-bHat_2))<1e-5  
```

*Test 4*.
```R
 bHat_1=coef(lm(y~X-1))
 bHat_2=LS.backfit(y=y,X=X,int=F,center=F)
 TEST[4]=max(abs(bHat_1-bHat_2))<1e-5  
```


*Test 5*.
```R
 bHat_1=coef(lm(y~X))[-1]
 bHat_2=LS.backfit(y=y,X=X,int=F,center=T)[-1]
 TEST[5]=max(abs(bHat_1-bHat_2))<1e-5  
```

*Test 6*.
```R
 yHat_1=predict(lm(y~X))
 yHat_2=cbind(1,X)%*%LS.backfit(y=y,X=X)
 TEST[6]=max(abs(yHat_1-yHat_2))<1e-5  
```

