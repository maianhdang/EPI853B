# EPI-809b


## Introduction 
   * Goals
   * Rules
   * Approach
   * Evaluation and grading

## (1) The R-software
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

## (2) Matrix Algebra (definitions and computational methods)
	* Basic definitions
	* Matrix addition and subtraction
	* Transpose
	* Matrix product
	* Inversion
	* Singular Value Decomposition
	* Generalized Inverse
	
## (3) Ordinary-least squares
	* The problem
	* Analytical solution
	* The lm and ls.fit functions
	* Regression with categorical predictors (`model.matrix`)
	* Computing OLS estimates using matrix operations
	* OLS using the QR and SVD decompositions
	* Iterative Procedures: Gauss-Seidel method
	
## (4) Inference in the OLS regression
	* Bias and variance of OLS estimates
	* Evaluation of Bias and Variance Using Monte Carlo Methods
	* The Bootstrap method
	
## (5) Maximum Likelihood estimation
	* The likelihood function
	* Analytical solution in the Gaussian and Bernoulli models
	* Numerical optimisation (application to GLM)
	* The Newton-Rapson Method
	* Data Augmentation and the EM-algoritm

## (6)  Multi-core computing in R
	* Introduction to multi-core computing
	* The `parallel` package

## (7) Simulation of Random Variables
	* Random numbers
	* Transformations of RV
	* Generating dras from Multivariate Normal
	* Inverse CDF
	
## (8) Monte Carlo Markov Chain
	* Introduction
	* Monte-Carlo Errors
	* Rejection Sampling
	* Gibbs Sampler
	* Metropolis Hastings
	
## (9) Cross-validation Methods
	* Why we need it?
	* Different types of prediction errors
	* Training-testing
	* Replicated training-testint
	* q-fold CV
	* Leave-one-out CV in OLS
	
## (10) Penalised Regressions
	* Why we needed it
	* Penalized RSS 
	* Standard penalty functions and the solutions they induce
	* Ridge Regression
	* Lasso
	* Bridge Regression
	* The coordinate descent gradient algorithm 
	
	
	
	
		   		
