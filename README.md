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

**Rank** (`qr(X)$rank`)

**Determinant**

**Matrix Inversion**

**Singular Value Decomposition**

**Generalized Inverse**

## Homework 1

Create a function to carry out:
	- Matrix products (`matProd(x,y)`)
	- crossprod (`crossprod_b(x,y)`, allow for the function to have only argument x)
	- tcrossprod (`tcrossprod_b(x,y)`, allow for the function to have only argument x)
	- Test your function using a randomnly genreated matrix with 5 rows and 10 columns, conduct 100 random tests.


**Materials**
   * Book: pages 9-12
   * A Review of [Matrix Algebra ](http://cs229.stanford.edu/section/cs229-linalg.pdf)


[Back to Outline](#Outline)
___

<div id="OLS-I" />
## (3) Ordinary-least squares [Chapter 3, plus materials provided below]

   * The problem
   * Analytical solution
   * The lm and ls.fit functions
   * Regression with categorical predictors (`model.matrix`)
   * Computing OLS estimates using matrix operations
   * OLS using the QR and SVD decompositions
   * Iterative Procedures: Gauss-Seidel method

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
