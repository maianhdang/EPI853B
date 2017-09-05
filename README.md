

# Statistical Computing  (EPI-853B)

In this course we will cover computational methods commonly used in statistics, including algorithms used for fitting and non-linear regressions, maximum likelihood estimation, simulation of random variables, bootstrap, cross-validation and algorithms for implementing high dimensional regressions.

**Software**: The course will be mostly based on [R](https://www.r-project.org/). If time permits we will also work with [Julia](http://julialang.org/).

**Approach**: Although the focus of the course is on computational methods, for each topic we will first describe the problem from a statistical perspective. If they exist, exact analytical solutions will be discussed and implemented. Otherwise numerical methods will be presented. Derivations will be presented in class and students are expected to take their own notes. Scripts for computations will be developed in class and a summary will be posted in this repository. Students are expected to bring their own laptops. If you do not have access to a laptop, please check with the instructor to get access to one.

**Evaluation**: The evaluation will be based on 4-5 HW and two in-class exams.

**Textbook**: [An Introduction to Statistical Learning](http://www-bcf.usc.edu/~gareth/ISL/index.html). This book covers many of the topics we will discuss. For topics not covered in the book we will provide aditional materials.

**Instructor**: Gustavo de los Campos (gustavoc@msu.edu)
<div id="Outline" />

[**Syllabus**](https://github.com/gdlc/EPI853B/blob/master/EPI853B_Syllabus.docx)

## Course Content

Note: this is a tentative list of topics, if time permits we will try to cover all of them; however, the list of topics is ambitious and we may not cover all the topics listed.

  * The R-software
  * Matrix Algebra in R
  * Ordinary Least Squares I: Estimation
  * Ordinary Least Squares II: Inference
  * Non-linear regression using Splines
  * Multi-core computing in R
  * Maximum Likelihood
  * Simulation of random variables
  * Monte Carlo Markov Chain Methods
  * Cross-Validation
  * Penalised Regressions

## Homework
  * [HW1](https://github.com/gdlc/EPI853B/blob/master/HW1.md)
  * [HW2](https://github.com/gdlc/EPI853B/blob/master/HW2.md)
  * [HW3](https://github.com/gdlc/EPI853B/blob/master/HW3.md)
  * [HW4](https://github.com/gdlc/EPI853B/blob/master/HW4.md)
  * [HW5](https://github.com/gdlc/EPI853B/blob/master/HW5.md)

<<<<<<< HEAD
<div id="intro" />
## Introduction 

   * Goals
   * Rules
   * Approach
   * Evaluation and grading

[Back to Outline](#Outline)
___

<div id="R" />
## (1) The [R](https://www.r-project.org/)-software 

**Topics you should be familiar with**
   * Installation [CRAN](https://cran.r-project.org/mirrors.html)
   * Types (boolean, integer, numeric, character, factors)
   * Variables and simple operations
   * Arrays 
   	* Vectors
   	* Matrices
   	* Multi-dimensional arrays
	* Lists
	* Data Frames	
   * Conditional statements
   * For and while loops
   * I/O
	* `read.table`  and `write.table`
	* `load` and `save`
	* `dput` and `dget`
	* `scan` and `wirte`
	* `readBin` and `writeBin`
   * Functions 

**Materials**
   * Section 2.3 of the book
   * [R website](https://www.r-project.org/) (see entry for manuals)
   * [R for Data Science](http://r4ds.had.co.nz/)

[Back to Outline](#Outline)
___

<div id="Matrix" />
## (2) Matrix Algebra (definitions and computational methods)

**A few matrices and vectors for demonstration**
```R
  x1=c(1,2,0)
  x2=c(1,-1,3)
  A=cbind(x1,x2)
  x3=c(3,4,1)
  x4=c(2,0,1)
  B=rbind(x3,x4)
  dim(A)
  dim(B)
```

**Matrix addition, subtraction and other cell-by-cell operations**
```R
  # cell-by-cell operations
  B+B
  B-B
  B^2
  exp(B) # note, this is very different than the matrix exponential, see below
  B*B    # note, this is very different than the matrix product
```

**Transpose**
```R
  Bt=t(B)
  all.equal(t(Bt),B)
```
Let's write our own function to transopose a matrix

```R
 ## A very inneficient function that illustrates matrix transposition and the creation of functions in R
 myT=function(X){
   nCol.out=nrow(X)
   nRow.out=ncol(X)
   OUT=matrix(nrow=nRow.out,ncol=nCol.out,NA)
   
   for(i in 1:nRow.out){
   	for(j in 1:nCol.out){
   		OUT[i,j]=X[j,i] # simply flip the indexes
   	}
   }
   return(OUT)
 }
 
 all(t(B)==myT(B))

```

**Matrix product** (`%*%`, `crossprod` and `tcrossprod`)

```R
  x=matrix(nrow=5,ncol=3,data=rnorm(15))
  t(X)%*%X   # X'X
  crossprod(X) # same but a bit faster
  X%*%t(X)   # XX'
  tcrossprod(X) # same but a bit faster
```

**QR-decomposition**

The QR-decomposition can be used to factorize a matrix into the product of an orthonomal basis (Q) and a triangular matrix (R) so that X=QR with Q'Q=I.

```R
 TMP=qr(X) # computes QR decomposition, returns a list
 qr.Q(TMP) # recovers Q from a QR decomposition
 qr.R(TMP) # recovers R from a QR decomposition
 all(qr.X(TMP)==X) # recovers X
 round(crossprod(qr.Q),5) # Q is orthonormal
```
**Rank**

```R
  qr(X)$rank
```
**Determinant**

```R
 X=matrix(nrow=3,ncol=3,data=rnorm(9))
 det(X)
 XtX=crossprod(cbind(X,X[,1])) # det of a singular (or rank defficient) matrix is zero
 det(XtX)
 qr(XtX)$rank
 
```
**Matrix Inversion**
 
```R
   XInv=solve(X)
   round(XInv%*%X,5)
   round(X%*%XInv,5) 
   solve(XtX) # rank-deficient matrix does not have an inverse
   
```
**Singular Value Decomposition**
 Any matrix X (nxp) can be decomposed as  X=UDV', where U(nxq), is an orthonormal basis for the row-space of X, V(qxn) is an orthonormal basis for the column-space of X and  D(qxq) is a diagnoal matrix with singular values in the diagnoal. Here q=min(n,p).
  
```R
   SVD=svd(X)
   str(SVD)
   round(crossprod(SVD$u),5) #$u gives an orthonormal basis for the row-space of X
   round(crossprod(SVD$v),5) #$v gives an orthonormal basis for the row-space of X
   SVD$d
   
   ## Now for a rank-deficient matrix
   SVD=svd(XtX) 
   sum(SVD$d>1e-10) # the rank is the sum of positive singular values
```

**Materials**
   * Book: pages 9-12
   * A Review of [Matrix Algebra ](http://cs229.stanford.edu/section/cs229-linalg.pdf)

[Back to Outline](#Outline)
___

<div id="OLS-I" />
## (3) Ordinary-least squares [Chapter 3, plus materials provided below]

####The OLS Problem

  Consider a regression problem of the form
  
      y=Xb+e
      
 where: `y` (nx1) is a response vector, X (nxp) is a design matrix of effects, `b` (px1) is a vector of regression coefficients and `e` (nx1) is a vector of model residuals.
 
 The ordinary least-squares estimate of b is obtained by miniminzing the residual sum of squares, `RSS=e'e=(y-Xb)'(y-Xb)`, with respect to `b`.
 
**Derivation**
The steps to find the analythical solution arre:
   - Differentiate RSS(y,X,b) with respect to the jth coefficient (bj)
   - Set all the queations equal to zero (this gives a stationary point).
   - Because e'e is a quadratic form, the matrix of second deriviatives is PSD, so we know that the stationary point is *a* solution.

**Analythical Solutions**
 
 We show in class that the solution is given by the vector bHat that solves the follwing systems of equations
 
 ```R     
 (X'X)bHat=X'y 
 
 ```
      
 If X is full-column rank, then we have


 ```R    
 bHat=solve(X'X)%*%X'y 
 
 ```
     
**Notation**: We usually call `C=X'X` the matrix of coefficients of the system of equations, `rhs=X'y`, the *right-and-side* of the sytem and `sol=bHat.` The short notation for the OLS system is

   ```R 
   C%*%sol=rhs 
   ```
    
## Computing OLS esitmates

In this examples we show alternative ways of computing OLS estimates. We begin with R functions `lm` and `lsfit` and then include
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

This method will work only if X has a rank equal to the number of columsn of X.

```R
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
```
**Inversion using cholesky decomposition**

This method will work only if X has a rank equal to the number of columsn of X. This guarantees that `X'X` can be factorized using the [Cholesky decomposition](https://en.wikipedia.org/wiki/Cholesky_decomposition). Sinche the cholesky factor is triangular (upper-triangular for R, we have `X'X=U'U`, where `U` is the upper-triangular Cholesky), we can easily compute the inverse of `X'X` from its Cholesky. 

Let's check that the inverse obtained via Cholesky is equivalent than the one computed using `solve`
 ```R
  Z=matrix(nrow=10000,ncol=1000,rnorm(1000*10000))
  C=crossprod(cbind(Z))
  system.time(CInv.solve<-solve(C))
  system.time(CInv.chol<-chol2inv(chol(C)))
  all(round(CInv.solve,8)==round(CInv.chol,8))

 ```

**OLS via Cholesky**

```R
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
```
**OLS using the QR-decomposition**

```R
z

#or

myLS.qr=function(y,X,int=TRUE){
  if(int){
      X=cbind(1,X)
  }
  QR=qr(X)
  sol=qr.coef(QR,y)
  return(sol)
}
```

**OLS using the singular-value decomposition**

```R
myLS.svd=function(y,X,int=TRUE){
  if(int){
      X=cbind(1,X)
  }
  p=min(dim(X))
  SVD=svd(X,nu=p,nv=p)
  gHat=crossprod(SVD$u,y)
  sol=tcrossprod(SVD$v,t(gHat/SVD$d))
  return(sol)
}
```

**Computing OLS estimates iterative procedures (Gauss-Seidel method)**

We can also find OLS estimates using iterative procedures. The following function implements the Gauss-Seidel method.
This method is equivalent to the back-fitting algorithm and it is also a coordinate descent gradient algorithm.

```R
 GS.solver=function(X,y,tol=1e-5,int=TRUE,center=TRUE){
  if(center){
       X=scale(X,center=TRUE,scale=FALSE) 
  }
  if(int){
   X=cbind(1,X)
  }
  C=crossprod(X)
  rhs=crossprod(X,y)
  p=ncol(X)
  b=rep(0,p)
  if(center&int){ b[1]=mean(y) }
  ready=FALSE
  iter=0
   while(!ready){
     iter=iter+1
     b0=b
     for(j in 1:p){
      tmp=sum(C[j,-j]*b[-j])
      b[j]=(rhs[j]- tmp )/C[j,j]
     }
     tmp=abs(b-b0)
     ready= all(tmp<tol)
     print(paste0('Iter=',iter, ' Max abs dif=',(max(tmp))))
   }
   return(b)
 }
```

Let's test it with a simple simulation.

```R
 n=5000
 p=100
 X=matrix(nrow=n,ncol=p,data=rnorm(n*p))
 bTRUE=c(runif(p))
 signal=X%*%bTRUE+100
 error=rnorm(n=length(signal),mean=0,sd=sd(signal)/2)
 y=signal+error

 system.time(bOLS<-coef(lm(y~X)))
 system.time(bOLS.GS<-GS.solver(y=y,X=X,center=TRUE,int=TRUE))

 max(abs(bOLS[-1]-bOLS.GS[-1])) # note: since we are centering for GS, the est. intercepts are different
```
  

## Regression with categorical predictors (`model.matrix`)


The `model.matrix()` function can be used to generate incidence matrices for regression models.
The following example illustrates how to produce contrats for a means-model.


#### Contrasts for categorical predictors

The following [webpage](http://www.ats.ucla.edu/stat/r/library/contrast_coding.htm) has a summary of some of the contrasts commonly used.

**Contrasts for Means Model**

```r
  sex=sample(c('M','F'),size=10,replace=T)
  Z=model.matrix(~sex-1) # by default R uses treatment contrasts, '-1' means removes the intercept
  crossprod(Z)
  table(sex)
```

R offers various types of contrasts, the following functions generate contrasts for factors with n levels:
  * `contr.treatment(n,...)` is the default option and generates contrasts for differences with a reference category.
  * `contr.helmert(n,...)` sequencial differences (group 2-group1, group3-mean of groups 1 and 2,...etc.
  * `contr.poly(n,...)`
  * `contr.sum(n,...)`
  * `contr.SAS(n,...)`
 
```R
  contr.treatment(4)
  contr.helmert(4)
  contr.sum(4)
```
We can tell `model.matrix` whant contrasts to use using the `contrasts.arg` argument.

```R
x1=c('a','a','b','b','b','c','c','c')
x2=c('1','1','2','2','1','3','3','1')
X=model.matrix(~x1+x2,contrasts.arg=list(x1=contr.treatment,x2=contr.helmert))
```

We can also recover the contrasts used by model matrix using the `attr()` function.

```R
attr(X,'contrasts')
```

We can also pass these arguments into lm.

```R
 y=rnorm(8)
 
 fm1=lm(y~x1+x2) # uses standars dummy coding (contr.tratment)
 fm2=lm(y~x1+x2,contrasts=list(x1=contr.helmert,x2=contr.sum))
 fm1$contrasts 
```
We can control what level will be the reference for the treatment contrasts by changhing the order of the levls of a factor.

```R
 x=factor(c('a','a','b','b','b','c','c','c'))
 levels(x)
 y=rnorm(length(x))
 summary(lm(y~x))
 z=factor(x,levels=c('c','a','b'))
 summary(lm(y~z))
```

[Back to Outline](#Outline)
___

<div id="OLS-II" />
## (4) Inference in the OLS regression [Chapters 3 and 5, plus materials provided below]

   
**Bias** Under relatively minimal conditions, OLS estimates are unbiased. What are these conditions? How do we show that OLS are unbiased? What conditions can lead to systematic bias?  
   
   
Consider a linear regession of the form  `y=Xb+e`, if `X` is a full-rank matrix and `b` is a vector of effects. The OLS estimate of `b` is given by 
   
       `bHat=Inv(X'X)X'y` 
       
The expected value of this estimator, for fixed X is
   
   `E[bHat|X]=E[Inv(X'X)X'y|X]=Inv(X'X)X'E[y|X]`
   
*Case 1:* If the linear model holds, that is if `E[y|X]=Xb`, then we have `E[bHat|X]==Inv(X'X)X'Xb=b`; therfore, we conclude that OLS estimates are unbiased. This result involves two assumptions: (i) that `X` is full rank and therfore `X'X` can be inverted and (ii) that the expected value of the error terms are zero. However, other assumptions that are sometimes made for inferential purposues, e.g., that the errors to be independen and/or normally distributed or that they have equal variance (homoskedasticity) are not required for the OLS estimates to be un-biased. Only (i) and (ii) are needed.

*Case 2:* Even if the model is non linear, the error terms from projection of `y` on `X` are, by construction, orthognal to `X`. So, even if the model is non-linear we have E[X'e]=0 and this makes OLS unbiased estimates of the projection of `y` on `X`, provided that `X` is full rank.

*Rank deficient case* If `X` is rank deficient, then we need to use a generalized-inverse, in this case, `Ginv(X'X)X'X` is not equal to an identity matrix and therfore, the expected value of the OLS estimate is not equal to the true value of the regression. However, even in this case we can still get un-biased estimates of 'estimable functions' for instance, contrasts between treatments.

**Variance** The conditional varinace of the OLS estimates is given by:  

`Var(bHat|X)=Var(Inv(X'X)X'y|X)=Inv(X'X)X'Var(y)XInv(X'X)=Inv(X'X)X'GXInv(X'X)`

where 'g' is the variance-covariance matrix of the error terms. If the errors are indpeendent and homoskedastic, then we have `Var(e)=I*vE`, and therefore,


`Var(bHat|X)=Inv(X'X)X'GXInv(X'X)=Inv(X'X)X'X(X'X)*vE=Inv(X'X)*vE`

Therofore, under (iii) independence and (iv) homoskedasticity of the error terms, the variance covaraince matrix of the OLS estimates is the product of the inverse of the coefficient matrix, `Inv(X'X)` times the error variance.

How do we estimate the error variance? An un-biased estimator can be obtained as: `eHat'eHat/(n-p)` where `eHat` are the OLS residuals (`eHat=y-XbHat`), `n` is sample size and `p` is the rank of X.

**Standard errors and p-values**. The square-root of the diagonal elements of the variance co-variance matrix of the OLS estimates are the SE of the estimates. Using this and the estimated effects we can obtain t-statistics and the corresponding p-values (`see pt()`).


**Omitted Variable Bias**. Suppose that the true model is given by `y=Xa+Zb+e`, and that instead we regress `y` on `X` only, are the OLS estimates still unbiased? The answer is yes, the OLS estimates from the regression of `y` on `X` are unbiased estimates of the projection of `y` on `X`, however, these are not un-biased estimates of `a`. To see this is useful to consider two models:


Full Model, or long regression:  `y=Xa+Zb+e`
Short regression:                 `y=Xc+d`

The OLS estimate of `c` is `Inv(X'X)X'y`. The expected value of this estimator is 

`E[cHat]=Inv(X'X)X'E[y]=Xa+Zb=a+Inv(X'X)X'Zb=a+Tb`

where `T` is a matrix contianing the regressions of the columns of Z on X. Therfore, the OLS estimates of the regression coefficients in the short regression can be bias with respect to the 'true' effects of X (`a`), however, they are still un-biased with respect to `c`, the projection of y on X.

## (5) Monte Carlo Methods
   * Introduction
   * Computing mean, variance and bias of estimates using MC methods
   * Monte Carlo Error
   
<div id="splines" />
## (7) Non-linear regression using splines [Chapter 7]	

**Approximating a conditional expectation function using bins**

```R
stepFunction=function(y,x,nW){
   if(nW<3){stop("The minimum # of windows is 3") }
   thresholds=quantile(x,prob= seq(from=1/nW,to=1-1/nW,by=1/nW))

   X=matrix(nrow=length(y),ncol=length(thresholds)+1)
   X[,1]=x<=thresholds[1]
   for(i in 2:(length(thresholds))){
      X[,i]=as.integer(x>thresholds[i-1] & x<=thresholds[i] )
   }
   X[,ncol(X)]<-x>max(thresholds)
   fm=lm(y~X-1)
   return(predict(fm))
 }
```

```R
 n=100
 x=seq(from=0,to=4*pi,length=n)
 signal=sin(x)+sin(x/2)
 error=rnorm(n)
 y=signal+error
 plot(y~x,col=4)
 lines(x=x,y=signal,col=2,lwd=2)
 lines(x=x,stepFunction(y,x,20),col=4,lty=2,lwd=2)
```

As we increase the number of windows, the bias of the estimator is reduced (we can approximate the true conditional expectation function better) but the variance increases.  How could we chose an optimal number of windows? One possiblity is to use cross-validation (see HW5).

While the step function is very flexible, it does not render a 'smooth' approximation; for instnace, the function is not continous. 


```R
 # This function evaluates (x-tau)+^degree for parametes x, tau and degree
 bf=function(x,tau,degree){ 
   z=x-tau	
   ifelse(z>0,z^degree,0)
 }
 

 nKnots=5
 degree=3
 
 knots=quantile(x,prob= seq(from=1/nKnots,length=nKnots,by=1/nKnots))
 Z=matrix(ncol=degree+1+length(knots),nrow=length(x),NA)
 
 # First let's through in our incidence matrix the polynomials
 Z[,1]=1
 for(i in 1:degree){ Z[,i+1]=x^i }
 
 # Now the local basis functions
 for(i in 1:nKnots){
   Z[,i+degree+1]=bf(x,tau=knots[i],degree=degree)
 }
 
  fm=lm(y~Z-1)
  plot(y~x)
  lines(x=x,y=signal,col=2,lwd=2)
  lines(x=x,y=predict(fm),col=4)
  
# Lets' compare with bs()
  W=bs(x,knots=knots,degree=degree) # note: by default bs does not include an intercept
  fm2=lm(y~W)
  points(x=x,y=predict(lm(y~W)),col='red',cex=.5)
  
```

The same approach can be used with higher order polynomials.


There are many different ways of creating basis functions for a spline. The `spline` package offeres functions for the so-called B-Splines (`bs`) and the Natural Spline (`ns`) which is a cubic spline with interior and boundary knots.

```R
 library(splines)
```
[Back to Outline](#Outline)

## (8) Bootstrap
   * Chapter 5, Section 5.2 
   * The Bootstrap method  [Efron & Gong, Am. Stat., 2012](http://www.tandfonline.com/doi/pdf/10.1080/00031305.1983.10483087?needAccess=true)
[Back to Outline](#Outline)
   

   
<div id="parallel" />
## (9) Parallel computing in R (Alexander Grueneberg) 

   * Introduction to multi-core computing
   * The [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf) package
   * Apply-type operations in parallel
   * Matrix computations in parallel
   * Monte Carlo Simulations using parallele computing

### Introduction

The [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf) package has been included in R since version 2.14.0. It combines the work of the two packages `snow` and `multicore` that both pursued different approaches in bringing parallel computing capabilities to R (as it lacks those features by design). The easiest way to use the approach of the `snow` package is by using the `parLapply` function. For the approach of the `multicore` package it is the `mclapply` function.

Both functions have important caveats: among those, `mclapply` will only work with UNIX-like operating systems (R, macOS) while `parLapply` will work with all operating systems but takes several manual steps to use. Both functions act similarly to the `lapply` function that comes with R that runs a function on each element of a list, returning the results of those functions calls as a list of the same length as `X`:

* `lapply(X, FUN)`
* `mclapply(X, FUN, mc.cores = detectCores())`
* `parLapply(cluster, X, FUN)`

The `mc.cores` argument of `mclapply` defines how many CPUs/cores should be used. `parLapply` takes a `cluster` argument that needs to be defined using `cluster <- makeCluster(detectCores())`.


### Example: Bootstrap in parallel

Task: Compute 90% bootstrap confidence interval of the bootstrap distribution of the sample mean for a small sample.

Serial as a `for` loop:

```R
start <- proc.time()
x <- c(8, 3, 2, 5, 6, -1, 2, 3, 5, 6, 3, 9)
nRep <- 1000000
means <- numeric(nRep)
seed <- 1000
for (i in 1:nRep) {
  set.seed(seed + i)
  z <- sample(x, replace = T, size = length(x))
  means[i] <- mean(z)
}
quantile(means, c(0.05, 0.95))
proc.time() - start
```

Serial using `apply`:

```R
start <- proc.time()
x <- c(8, 3, 2, 5, 6, -1, 2, 3, 5, 6, 3, 9)
nRep <- 1000000
seed <- 1000
meanFunction <- function(i) {
  set.seed(seed + i)
  z <- sample(x, replace = T, size = length(x))
  return(mean(z))
}
means <- sapply(1:nRep, meanFunction) # or lapply followed by unlist
quantile(means, c(0.05, 0.95))
proc.time() - start
```

In parallel using `mclapply` (this will be serial on Windows):

```R
start <- proc.time()
library(parallel)
x <- c(8, 3, 2, 5, 6, -1, 2, 3, 5, 6, 3, 9)
nRep <- 1000000
seed <- 1000
meanFunction <- function(i) {
  set.seed(seed + i)
  z <- sample(x, replace = T, size = length(x))
  return(mean(z))
}
means <- mclapply(1:nRep, meanFunction, mc.cores = detectCores())
means <- unlist(means)
quantile(means, c(0.05, 0.95))
proc.time() - start
```

In parallel using `parLapply` (this should work with all operating systems):

```R
start <- proc.time()
library(parallel)
x <- c(8, 3, 2, 5, 6, -1, 2, 3, 5, 6, 3, 9)
nRep <- 1000000
seed <- 1000
meanFunction <- function(i) {
  set.seed(seed + i)
  z <- sample(x, replace = T, size = length(x))
  return(mean(z))
}
cluster <- makeCluster(detectCores())
clusterExport(cluster, c("x", "seed"))
means <- parLapply(cluster, 1:nRep, meanFunction)
means <- unlist(means)
quantile(means, c(0.05, 0.95))
stopCluster(cluster)
proc.time() - start
```

Benchmarks on *my* machine (Core i5, 4 cores);

```
		      min        lq      mean    median        uq       max neval
for		12.662996 12.747717 13.041865 12.845498 13.070518 13.882597     5
apply		14.074011 14.250181 14.840325 14.264057 14.683107 16.930270     5
mclapply	 6.591227  6.643596  7.054739  6.790404  6.893870  8.354599     5
parLapply	 8.253781  8.363489  8.408053  8.389866  8.413008  8.620121     5
```

[Back to Outline](#Outline)
___

<div id="ML" />
## (10) Maximum Likelihood estimation [Chapter 4, sections 4.1-4.3]

For a paper describing the origins of the principle of Maximum Lilkelihood see [Aldrich, 1997](https://projecteuclid.org/download/pdf_1/euclid.ss/1030037906)

** 1. Closed Form Solutions**
The liklelihood function is the probability of the data given the parameters viewed as a function of the parameters with data fixed.
Maximum likelihood Estimates (MLEs) are the values of the parameters that maximize the likelihood function. 

In simple models (e.g., Poisson, Bernoully, multiple-linear regression with Gaussian IID error terms) we can find a closed-form for the ML estimates. In these cases MLEs are obtained by:

  - Taking the 1st derivative of the log-likelihood with respect to each of the parameters
  - (First Order Conditions, FOC) Set all the derivatives equal to zero, this gives a stationary point and renders a system with as many equations as parameters.
  - Solve for the parameter values that satisfy the FOC.
  - Check, based on the sign of the 2nd derivatives that the function is concave.
  
  
However, in the vast majority of the cases we cannot find a closed-form for the ML estimates. In these cases we use numerical methods. There are many numerical methdos for optimizing a function, we will briefly consider three approaches: grid search, Newton-Rapson and general purpouse optimization algorithm (functions `optim` and `optimize` in R).


**1. Grid search**

In it's simplest form a grid search requires:
 - Implementing a function for evaluating the log-likelihood
 - Defining a grid of values for the parameters
 - Evaluating the log-likelihood for every value in the grid
 - Finding the value within the grid that yields the largest value of the likelihood

The following R-code illustrates this for the case of Poisson data

NOTE: In this case the MLE has a closed form (it is simply the sample mean)

 ```R
 ## Simulating data
  trueLambda=10
  y=rpois(lambda=trueLambda,n=100)

 ## A function to evaluate the Poisson log-likelihood
  logLikPois=function(y,lambda){     
       logLik=log(lambda)*sum(y)-length(y)*lambda
       return(logLik)
  }

## A grid
 myGrid=seq(from=.01,to=100,length=10000)

 ## Search
 logLik=rep(NA,length(myGrid))
 for(i in 1:length(myGrid)){
  logLik[i]=logLikPois(y=y,lambda=myGrid[i])
  print(i)
 }
 plot(logLik~myGrid,type='l',col=4)
 abline(v=mean(y),lty=2,col=2)
```

If we have a parameter vector instead of a single parameter our grid will be multi-dimensional (as many dimensions as parameters) but the princple of the grid search is the same.

Grid Search is computationally intensive, most of the values of the grid could be ruled out y looking at the derivatives of the function.

The above approach very crude. If we know that our function is concave, the we can set a very sparse grid, then find three points where the interior point gives higher likelihood than the flanking points. This interval must contain the maximum. We can then refine our grid by making a denser grid within that interval, and interate in this fashion untlil our estimate does not change more than a pre-established precision.


**2 Newton Rapson (NR) Method**

This method can be used to find a stationary point (either a minima or maxima) of a twice differentiable function. Starting from an intial guess `theta_0`, the NR method moves from `theta_0` to `tehta_1=theta_0+delta_0` where `delta_0` is the negative of the ratio of the first and seconde derivatives of the objective function both evaluated at `theta_0`. The method can be motivated using a 2nd order Taylor expansion to the objective function. We will discuss this further in class.

The following example illustrates the aplictaion of the NR methods for estimating the parameter of a Poisson distribution.

```R
## Data
 ## Simulating data
  trueLambda=10
  y=rpois(lambda=trueLambda,n=100)


 ## Derivatives
  firstD=function(lambda,y){  
     sum(y)/lambda -length(y)
  }
  
  secondD=function(lambda,y){
     -sum(y)/(lambda^2)
  }
  
  
  theta=3
  tol=1e-5
  ready=FALSE
  counter=0
  while(!ready){
  	print(theta)
  	delta=- firstD(theta,y)/secondD(theta,y)
	new_theta=theta+delta
	ready=(abs(theta-new_theta)<tol)
	counter=counter+1
	print(paste0(counter,'  ', new_theta))
	theta=new_theta
  }
```


**3 General Purpose Optimization Functions**

In this section we illustrate how to obtain MLEs using `optim()`, this is a general purpose optimizaton function that can be used to minimize an objective function with respect to a set of parameters. The inputs for `optim` are: (i) the function to be minimized (e.g., negative-log-likelihood), (ii) a parameter vector (all the parameters to be estimated must be provided in a single vector), and (iii) additional arguments to the function that are to be kept fixed during the optimization procedure (e.g., data, covariates, etc.). We illustrate the use `optim` in the context of logistic regression and censored regression.


**3.1. Logistic Regression**


*Toy simulation*

```R 
 x=rnorm(1000)
 risk=2+x*2
 prob=exp(risk)/(1+exp(risk))
 y=rbinom(prob=prob,n=length(x),size=1)
```

*Logistic Regression using GLM*
```R
 fm0=glm(y~x,family=binomial)
 fm0$coef
```
*Logistic regression using optim*

```R
 fm1=optim(fn=negLogLikLogistic,y=y,X=cbind(1,x),par=c(log(mean(y)/(1-mean(y))),0))
 fm1$par
```

**Inference**

Under [regularity conditons](http://www.stat.yale.edu/~pollard/Books/Asymptopia/Classical.pdf), the large-sample distribution of  ML estimates is Multivariate Normal with mean equal to the true parmeter values (i.e., MLs are asymptotically un-biased) and (co)variance matrix equal to the inverse of [Fisher's Information matrix](https://en.wikipedia.org/wiki/Fisher_information) which is the expected value of the matrix of second derivatieves of the log-likelihood function with respect to the parameters. This matrix often depends on the true parameter values which are uknown. Thus, in practice, Expected Information is replaced with the *Observed Information* which is the inverse of the matrix of the second derivatives of the log-likelihood evaluated at the ML estimates. We can obtain this matrix from `optim()` by setting `hessian=TRUE`. The followin example illsutrates this.

```R
 fm1=optim(fn=negLogLikLogistic,y=y,X=cbind(1,x),par=c(log(mean(y)/(1-mean(y))),0),hessian=T)
 OI=solve(fm1$hessian)
 cbind(fm1$par,sqrt(diag(OI)))
 summary(fm0)$coef# GLM
```

**3.2. Censored Regression**

Biomedical and engeneering data is usually censored. Censored means that for some subjects we do not know the actual value of the response but we have information (bounds) about it. There are three types of censoring:
 - Right censored, we know that y is greater than some constant, `y>c`.
 - Left censored, we know that y is smaler than some constant, `y<c`.
 - Interval censored, means that we know that y lies in an interval, `c1<y<c2`.
 
 As always, the likelihood function is the probability of the data given the parameters viewed as a function of the parameters with data fixed. For un-censored data points this is just the density function. For right-censored data point the likelihood is the Survival function (i.e., `S(c)=p(y>c)=1-CDF(c)`). For left-censored data points the likelihood is the CDF, `p(y<c)=CDF(c)`, and for interval censored the likelihood is the difference in the CDF function evaluated at the right and left-bounds of the interval, `p(c1<y<c2)=CDF(c2)-CDF(c1)`.  The following example illustrates how to compute ML estimates with censored data assuming that data follows a Gamma distribution.
 
 
*Suppose we have gamma data*

```R
  x=rgamma(shape=8,rate=.5,n=100)
```


*ML estimation with the complete data*

```R
 negLogLikGamma=function(x,theta){
	logLik=dgamma(x=x,shape=exp(theta[1]),rate=exp(theta[2]),log=TRUE)
	return(-sum(logLik))
 }

 fm=optim(fn=negLogLikGamma,x=x,par=log(c(2,2)))
 exp(fm$par)
 
```

*Suppose now that a proportion of the data are right-censored*

```R
 event=runif(length(x))<.7
 y=x
 y[!event]=x[!event]-runif(n=sum(!event),min=1,max=3)
```


*Log-likelihood function for right-censored data*

```R
 negLogLikGammaRightCen=function(time,theta){
 	logLik=pgamma(q=time, shape=exp(theta[1]), rate = exp(theta[2]), lower.tail = FALSE, log.p = TRUE)
 	 return(-sum(logLik))
 }

 logLikGammaCen=function(time,event,theta){
 	
 		time_observed=time[event]
 		time_right_cen=time[!event]
 		
 		negLogLik_1=negLogLikGamma(x=time_observed,theta)
 		negLogLik_2=negLogLikGammaRightCen(time=time_right_cen,theta=theta)
 		return((negLogLik_1+negLogLik_2))
 }
 
 fm2=optim(fn=logLikGammaCen,time=x,event=event==T,par=log(c(2,2)))
 
  exp(fm2$par)
```

**The `survaival` R-package**



[Back to Outline](#Outline)
___




<div id="RV" />
## (8) Simulation of Random Variables

   * Random numbers
   * Transformations of RV
   * Generating dras from Multivariate Normal
   * Inverse CDF

[Back to Outline](#Outline)
___

<div id="MCMC" />
## (9) Monte Carlo Markov Chain Methods [References provided below]

   * Introduction
   * The Monte-Carlo Error
   * Basic sampling methods
   	* Rejection Sampling
   	* Gibbs Sampler
   	* Metropolis Hastings

[Back to Outline](#Outline)
___

<div id="CV" />
## (10) Cross-validation Methods [Chapter 5]

   * Why we need it?
   * Different types of prediction errors
   * Validation methods
   	* Training-testing
   	* Replicated training-testint
   	* q-fold CV
   * Leave-one-out CV in OLS
   * Choosing the optimal number of knots in a spline using CV.

[Back to Outline](#Outline)
___

<div id="penalised" />

## (11) Penalised Regressions [Chapter 6]

   * Why are penalized regressions needed? The variance-bias trade off
   * Penalized RSS 
   * Standard penalty functions and the solutions they induce
   * Ridge Regression
	* Lasso
	* Bridge Regression
	* The coordinate descent gradient algorithm 
=======
>>>>>>> origin/master
