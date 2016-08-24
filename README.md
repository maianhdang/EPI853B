# Statistical Computing  (EPI-853B)

In this course we will cover computational methods commonly used in statistics, including algorithms used for fitting and non-linear regressions, maximum likelihood esitmation, simulation of random variables, bootstrap, cross-validation and algorithms for implementing high dimensional regressions.

**Software**: The course will be mostly based on [R](https://www.r-project.org/). If time permits we will also work with [Julia](http://julialang.org/).

**Approach**: Althought the focus of the coruse is on computational methods, for each topic  we will first describe the problem from an statistical perspective. If they exhist, exact analythical solutions will be discuss and implemented. Otherwise numerical methods will be presented. Derivations will be presented in class and students are expected to take their own notes. Scripts for computations will be develped in class and a summary will be posted in this repository. Students are expected to bring their own laptops. If you do not have access to a laptop, please check with the instructor to get access to one.

**Evaluation**: The evaluation will be based on HW (number to be determined) and two in-class exam.

**Textbook**: [An Introduction to Statistical Learning](http://www-bcf.usc.edu/~gareth/ISL/index.html). This book covers many of the topics we will discuss. For topics not covered in the book we will provide aditional materials.

**Instructor**: Gustavo de los Campos (gustavoc@msu.edu)

## Outline
<div id="Outline" />

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
  * [Penalised Regressions](#penalized)



## Introduction 
<div id="intro" />
   * Goals
   * Rules
   * Approach
   * Evaluation and grading

[Back to Outline](#intro)
___

## (1) The R-software  [Section 2.3, plus materials provided below]
<div id="R" />
   * Installation [CRAN]()
   * Variables and simple operations
   * Arrays 
   	* Vectors
   		* Creation
   		* Indexing and replacment (three methods)
   		* Attributes (length, names)
   		 Binding vectors
   	* Matrices
   		* Creation
   		* Indexing and replacment (three methods)
   		* Atrributes (rownames, colnames, dim, nrow, ncol)
   	* Multi-dimensional arrays
   		* Creation
   		* Indexing and replacment (three methods)
   		* Attributes
	* Lists
		* Why we need lits
		* Creation
		* Indexing and replacement (three methods)
	* Data Frames	
		* Why are these needed
		* Creating data frames (using read table, using data.frame and using as.data.frame)
		* Indexing and replacement (four methods)
	* Conditional statements
	* For and while loops
	* I/O
		* `read.table`  and `write.table`
		* `load` and `save`
		* `dput` and `dget`
		* `scan` and `wirte`
		* `readBin` and `writeBin`

[Back to Outline](#intro)
___

## (2) Matrix Algebra (definitions and computational methods)
<div id="Matrix" />
   * Basic definitions
   * Matrix addition and subtraction
   * Transpose
   * Matrix product
   * Inversion
   * Singular Value Decomposition
   * Generalized Inverse

[Back to Outline](#intro)
___

## (3) Ordinary-least squares [Chapter 3, plus materials provided below]
<div id="OLS-I" />
   * The problem
   * Analytical solution
   * The lm and ls.fit functions
   * Regression with categorical predictors (`model.matrix`)
   * Computing OLS estimates using matrix operations
   * OLS using the QR and SVD decompositions
   * Iterative Procedures: Gauss-Seidel method

[Back to Outline](#intro)
___

## (4) Inference in the OLS regression [Chapters 3 and 5, plus materials provided below]
<div id="OLS-II" />
   * Bias and variance of OLS estimates
   * Omitted variable bias
   * Evaluation of Bias and Variance Using Monte Carlo Methods
   * The Bootstrap method [Chapter 5, Section 5.2]

## (5) Maximum Likelihood estimation [Chapter 4, sections 4.1-4.3]
<div id="ML" />
   * The likelihood function
   * Analytical solution in the Gaussian and Bernoulli models
   * Numerical optimisation (application to GLM)
   * The Newton-Rapson Method
   * Data Augmentation and the EM-algoritm

[Back to Outline](#intro)
___

## (6) Non-linear regression using splines [Chapter 7]	
<div id="splines" />
   * Basis functions
   * Non-linear regression using splines

[Back to Outline](#intro)
___

## (7)  Multi-core computing in R [the parrallel R-package](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf)
<div id="parallel" />
   * Introduction to multi-core computing
   * The `parallel` package

[Back to Outline](#intro)
___

## (8) Simulation of Random Variables
<div id="RV" />
   * Random numbers
   * Transformations of RV
   * Generating dras from Multivariate Normal
   * Inverse CDF

[Back to Outline](#intro)
___

## (9) Monte Carlo Markov Chain Methods [References provided below]
<div id="MCMC" />
   * Introduction
   * The Monte-Carlo Error
   * Basic sampling methods
   	* Rejection Sampling
   	* Gibbs Sampler
   	* Metropolis Hastings

[Back to Outline](#intro)
___

## (10) Cross-validation Methods [Chapter 5]
<div id="CV" />
   * Why we need it?
   * Different types of prediction errors
   * Validation methods
   	* Training-testing
   	* Replicated training-testint
   	* q-fold CV
   * Leave-one-out CV in OLS
   * Choosing the optimal number of knots in a spline using CV.

[Back to Outline](#intro)
___

## (11) Penalised Regressions [Chapter 6]
<div id="penalised" />
   * Why are penalized regressions needed? The variance-bias trade off
   * Penalized RSS 
   * Standard penalty functions and the solutions they induce
   * Ridge Regression
	* Lasso
	* Bridge Regression
	* The coordinate descent gradient algorithm 
