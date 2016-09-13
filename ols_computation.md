## Computing OLS esitmates

In this examples we show alternative ways of computing OLS estimates. We begin with R functions (`lm` and `lsfit`) and then include
alternative ways of computing estimates using matrix operations, factorizations and with iterative procedures.

**A simple simulated data set**

```R
 p=100
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


myLS.solve=function(y,X,int=TRUE){
  if(int){
         X=cbind(1,X)
  }
  C=crossprod(X)
  rhs=crossprod(X,y)
  CInv=solve(C)  
  sol=crossprod(CInv,rhs)
  return(sol)
}

timeIn=proc.time()
for(i in 1:1000){
  fm=lm(y~X)
}
timeOut=proc.time()

timeIn=proc.time()
for(i in 1:1000){
  fm=myLS.solve(y,X)
}
timeOut=proc.time()

myLS.chol=function(y,X,int=TRUE){
  if(int){ 
         X=cbind(1,X)
  }
  C=crossprod(X)
  rhs=crossprod(X,y)
  CInv=chol2inv(chol(C))
  sol=crossprod(CInv,rhs)
  return(sol)
}

myLS.qr=function(y,X,int=TRUE){
  if(int){
      X=cbind(1,X)
  }
  QR=qr(X)
  Q=qr.Q(QR)
  gHat=crossprod(Q,y)
  R=qr.R(QR)
  RInv=solve(R)
  sol=RInv%*%gHat
  return(sol)
}
**Inversion using iterative procedures (Gauss-Seidel method)**


