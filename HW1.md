###  HW1: Ordinary least squares


**Due: Thursday, Sept. 14, in class**.

#### (1) Derivation of the OLS estimator

1.a. Consider a `ray-regression` of the form **y**=**Xb**+**e**. For this regression derive the OLS estimator analythically.

1.b. Expand the model by including an intercept, that is, consider y=a+xb+e. Derive the OLS analythically.


#### (2) OLS: comparison of `lm`, `lsfit` and computation using matrix operations.

Develop a function that will take as arguments a response vector (`y`) and a matrix of covariates (`X`), computes internally OLS estimates and return them.


**2.1**. Using the following simulation compute OLS estimates using `lm`, `lsfit` and your function and report the coefficients you got with
each of these functions.


```r
  n=1000
  b=c(130,12,-4,3)
  p=length(b)
  X=matrix(nrow=n,ncol=length(b),data=rnorm(n*p))
  signal=X%*%b
  error=rnorm(n=n,sd=sd(signal))
  y=signal+error
```

**2.2.** Computational performance

Estimate the computational time needed to obtain 10,000 times OLS estimates for the above simulation using each of the three functions considered above.
Note: carry out the simulation only once and fit OLS estimates 10,000 times.


[Home](https://github.com/gdlc/EPI853B)

