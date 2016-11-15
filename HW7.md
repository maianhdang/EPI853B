# HW 7: Maximum likelihood with censored data

Due Tuesday Nov. 22nd in class (bring your printed report to class).


**1.** For each of the following problem estimate the mean and variance using maximum likelihood using: (a) the `survreg()` function, and
(b) `optim()`. For survreg present your code and the results of `summary(fm)`. For (b) reproduce the results of `survreg()` using `optim()`.


*1.1. Right censored data*

Use only the data in `Y` to fit your models.

```R
 mu=10
 SD=2
 n=1000
 y=rnorm(mean=mu,sd=SD,n=n)
 d=rbinom(n=n,prob=.6,size=1)
 z=runif(min=.02,max=1,n=n) 
 time=ifelse(d==1,y,y-z)
 Y=data.frame(d=d,time=time)
```

*1.2. Left censored data*

Use only the data in `Y` to fit your models.

```R
 mu=10
 SD=2
 n=1000
 y=rnorm(mean=mu,sd=SD,n=n)
 d=rbinom(n=n,prob=.6,size=1)
 z=runif(min=.02,max=1,n=n) 
 time=ifelse(d==1,y,y+z)
 Y=data.frame(d=d,time=time)
```

**2. Estimating the bias induced by ignoring censoring**

Using the code in 1.1 estimate, using 5000 MC replicates, the bias of  the estimate of the intercept and of the variance  of each of these methods:
  - Maximum Likelihood accounting for censoring (user `survreg()` for that),
  - The naieve estimators `muHat=mean(Y$time)` and `vHat=var(Y$time)`.
  
Report a table with estimate bias for each of the parameters (intercept and variance) and each of the methods and sumarize your finding in no
more than two sentences.


