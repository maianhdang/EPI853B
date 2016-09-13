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
      
 where: y (nx1) is a response vector, X (nxp) is a design matrix of effects, b (px1) is a vector of regression coefficients and e (nx1) is a vector of model residuals.
 
 The ordinary least-squares estimate of b is obtained by miniminzing the residual sum of squares, RSS=e'e=(y-Xb)'(y-Xb), with respect to b.
 
**Derivation**
The steps to find the analythical solution arre:
   - Differentiate RSS(y,X,b) with respect to the jth coefficient (bj)
   - Set all the queations equal to zero (this gives a stationary point).
   - Because e'e is a quadratic form, the matrix of second deriviatives is PSD, so we know that the stationary point is *a* solution.

**Analythical Solutions**
 
 We show in class that the solution is given by the vector bHat that solves the follwing systems of equations
 
 ```R     (X'X)bHat=X'y ```
      
 
 If X is full-column rank, then we have
 
 ```R    bHat=solve(X'X)%*%X'y ```
     
**Notation**: We usually call ```R C=X'X ``` the matrix of coefficients of the system of equations, ```R rhs=X'y ```, the *right-and-side* of the sytem and `sol=bHat.` The short notation for the OLS system is

   ```R C%*%sol=rhs ```
    
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

The coefficient matrix (X'X) is symetric. If X is full-rank, then X'X is postivie definite (pd), meaning that a'X'Xa>0 for all vectors a different than zero. If X is rank deficient X'X is positive semi-definite (psd), meaning that  a'X'Xa>=0 for all vectors a different than zero. Next we will show how to obtain the inverse of the coefficient matrix using the Cholesky decomposition of X'X.

Symmetric PD matrices can be decomposed using the [Cholesky]() decomposition.

Lower-triangular Cholesky of PD matrix A:    A=LL'
Upper-triangular Cholesky of PD matrix A: U=L', A=U'U.

The `chol()` function of R returns the upper-triangular Cholesky factor of a matrix. 

Using solve(A%*%B)=solve(B)%*%solve(A), we get that 

solve(C)=solve(U)%*%t(solve(U))

where U is the Cholesky of C. Inverting a triangular matrix is much simpler than inverting any general non-singular sqaure matrix. The function

    chol2inv()
    
 takes as an argument the cholseky factor of a PD matrix and returns the inverse of the same matrix, but the inversion is much faster than using solve when p is large.
 
 ```R
  Z=matrix(nrow=10000,ncol=1000,rnorm(1000*10000))
  C=crossprod(cbind(Z))
  system.time(CInv.solve<-solve(C))
  system.time(CInv.chol<-chol2inv(chol(C)))
  all(round(CInv.solve,8)==round(CInv.chol,8))

 ```
**OLS using the QR-decomposition**

```R
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
```
**OLS using the singular-value decomposition**




**Inversion using iterative procedures (Gauss-Seidel method)**


* Regression with categorical predictors (`model.matrix`)


[Back to Outline](#Outline)
___

<div id="OLS-II" />
## (4) Inference in the OLS regression [Chapters 3 and 5, plus materials provided below]

   * Bias and variance of OLS estimates
   * Omitted variable bias
   * Evaluation of Bias and Variance Using Monte Carlo Methods
   * The Bootstrap method [Chapter 5, Section 5.2]

<div id="ML" />
## (5) Maximum Likelihood estimation [Chapter 4, sections 4.1-4.3]

   * The likelihood function
   * Analytical solution in the Gaussian and Bernoulli models
   * Numerical optimisation (application to GLM)
   * The Newton-Rapson Method
   * Data Augmentation and the EM-algoritm

[Back to Outline](#Outline)
___

<div id="splines" />
## (6) Non-linear regression using splines [Chapter 7]	

   * Basis functions
   * Non-linear regression using splines

[Back to Outline](#Outline)
___

<div id="parallel" />
## (7)  Multi-core computing in R [the parrallel R-package](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf)

   * Introduction to multi-core computing
   * The `parallel` package

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
