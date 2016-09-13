## Computing OLS esitmates

In this examples we show alternative ways of computing OLS estimates. We begin with R functions (`lm` and `lsfit`) and then include
alternative ways of computing estimates using matrix operations, factorizations and with iterative procedures.

**A simple simulated data set**

```R
 p=5
 n=5000
 b=rnorm(p)
 X=matrix(nrow=n,ncol=p,data=runif(n*p))
 signal=150+X%*%b
 error=rnorm(sd=sd(signal),n=n)
 y=signal+error

```

**The `lm` function**
```R
 fm<-lm(y~X)
 bHat=coef(fm)
 summary(fm)
```
**The `lsfit` function**
```R

 fm2<-lsfit(y=y,x=X)
 bHat=coef(fm)
 ls.print(fm2) # the summary method is not very useful with lsfit
```
**Our own lm using matrix operations**

**Inversion using cholesky decomposition**

**OLS using the QR-decomposition**

**OLS using the singular-value decomposition**

**Inversion using iterative procedures (Gauss-Seidel method)**


