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

  * [Introduction](#intro)
  * [The R-software](#R)
  * [Matrix Algebra in R](#Matrix)   
  * [Ordinary Least Squares I: Estimation](#OLS-I)
  * [Ordinary Least Squares II: Inference](#OLS-II)
  * [Maximum Likelihood](#ML)
  * [Non-linear regression using Splines](#splines)
  * [Multi-core computing in R](#parallel)
  * [Simulation of random variables](#RV)
  * [Monte Carlo Markov Chain Methods](#MCMC)
  * [Cross-Validation](#CV)
  * [Penalised Regressions](#penalised)

## Homework
  * [HW1](https://github.com/gdlc/EPI853B/blob/master/HW1.md)
  * [HW2](https://github.com/gdlc/EPI853B/blob/master/HW2.md)
  * [HW3](https://github.com/gdlc/EPI853B/blob/master/HW3.md)
  * [HW4](https://github.com/gdlc/EPI853B/blob/master/HW4.md)
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
myLS.qr=function(y,X,int=TRUE){
  if(int){
      X=cbind(1,X)
  }
  QR=qr(X)
  Q=qr.Q(QR)
  R=qr.R(QR)
  gHat=crossprod(Q,y)
  sol=backsolve(R,gHat)
  return(sol)
}
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
 bf=function(x,tau,degree){ 
   z=x-tau	
   ifelse(z>0,z^degree,0)
 }
 DF=5
 thresholds=quantile(x,prob= seq(from=1/DF,to=1-1/DF,by=1/DF))
 Z=matrix(ncol=DF+1,nrow=length(x),NA)
 Z[,1]=1
 Z[,2]=x
 for(i in 3:ncol(Z)){
   Z[,i]=bf(x,thresholds[i-2],1)
 }
 
  fm=lm(y~Z-1)
  plot(y~x)
  lines(x=x,y=signal,col=2,lwd=2)
  lines(x=x,y=predict(fm),col=4)

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
   * Bootstrap in parallel
   
[Back to Outline](#Outline)
___

<div id="ML" />
## (10) Maximum Likelihood estimation [Chapter 4, sections 4.1-4.3]

   * The likelihood function
   * Analytical solution in the Gaussian and Bernoulli models
   * Numerical optimisation (application to GLM)
   * The Newton-Rapson Method
   * Data Augmentation and the EM-algoritm

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
