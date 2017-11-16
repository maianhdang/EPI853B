

# Statistical Computing  (EPI-853B)

In this course we will cover computational methods commonly used in statistics, including algorithms used for fitting and non-linear regressions, maximum likelihood estimation, simulation of random variables, bootstrap, cross-validation and algorithms for implementing high dimensional regressions.

**Software**: The course will be mostly based on [R](https://www.r-project.org/). If time permits we will also work with [Julia](http://julialang.org/).

**Approach**: Although the focus of the course is on computational methods, for each topic we will first describe the problem from a statistical perspective. If they exist, exact analytical solutions will be discussed and implemented. Otherwise numerical methods will be presented. Derivations will be presented in class and students are expected to take their own notes. Scripts for computations will be developed in class and a summary will be posted in this repository. Students are expected to bring their own laptops. If you do not have access to a laptop, please check with the instructor to get access to one.

**Evaluation**: The evaluation will be based on 4-5 HW and two in-class exams.

**Textbook**: [An Introduction to Statistical Learning](http://www-bcf.usc.edu/~gareth/ISL/index.html). This book covers many of the topics we will discuss. For topics not covered in the book we will provide aditional materials.

**Instructor**: Gustavo de los Campos (gustavoc@msu.edu)


[**Syllabus**](https://github.com/gdlc/EPI853B/blob/master/EPI_863B_Syllabus.pdf)


<div id="Outline" />

## Course Content

Note: this is a tentative list of topics, if time permits we will try to cover all of them; however, the list of topics is ambitious and we may not cover all the topics listed.

  * [The R-software](#R)
  * [Matrix Algebra in R](#matrix)
  * Ordinary Least Squares I: Estimation & Inference
    * Derivation of closed-form solution
    * Computation using `lm`, `lsfit` [OLS](https://github.com/gdlc/EPI853B/blob/master/OLS.md)
    * Computing estimates, SEs and p-values using matrix operations [OLS](https://github.com/gdlc/EPI853B/blob/master/OLS.md)
    * Approximating the sampling distribution of estimates using [Bootstrap](https://github.com/gdlc/EPI853B/blob/master/Bootstrap.md)
  * Matrix factorizations:  [singular-value decomposition](https://github.com/gdlc/EPI853B/blob/master/matrixFactor.md), [eigen-decomposition](https://github.com/gdlc/EPI853B/blob/master/matrixFactor.mdF), [QR-decomposition](https://github.com/gdlc/EPI853B/blob/master/matrixFactor.md) and [Cholesky](https://github.com/gdlc/EPI853B/blob/master/matrixFactor.md).
  * Maximum Likelihood
    * Reiew of maximum likelihood
    * Estimation using `optim`
    * Maximum likelihood estimation and inference in the logistic regression.
    * Maximum likelihood estimation and inference in parametric survival regression.
  * Permutations 
    * [A gentle introduction to permutation tests](http://www.tandfonline.com/doi/abs/10.1198/000313008X269576)
    * [Multiple testing and permutations in GWAS](https://www.nature.com/articles/nrg3706)
    * [Examples](https://github.com/gdlc/EPI853B/blob/master/permutations.md)
  * Multi-core computing in R
  * Assesment of prediction accuracy using cross-validation methods
  * Introduction to linux and distributed computing in High-performance Computing Clusters
  
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

<div id="R" />

## (1) The [R](https://www.r-project.org/)-software 

   * [R-Intro](https://github.com/quantgen/RIntro)   
   * [R website](https://www.r-project.org/) (see entry for manuals)
   * [R for Data Science](http://r4ds.had.co.nz/)

[Back to Outline](#Outline)


<div id="matrix" />

## (2) Matrix Algebra (definitions and computational methods)

**Definition**. A matrix is a 2-dimensional array of values of the same type.

**Internal structure**. By default R sorts matrices by column.
```r
  X=matrix(nrow=3,ncol=2,data=1:6)
  X
  # you can use the `byrow` argument to tell R that you are providing the data to fill the matrix sorted by rows.
  X=matrix(nrow=3,ncol=2,data=1:6,byrow=TRUE)
  X
```

**Column and rownames**. We can append names to rows and columns.

```r
  rownames(X)=c('a','b','c')
  colnames(X)=c('c1,'c2')
  colnames(X)
  rownames(X)
  dimnames(X)
  
  X=X[3:1,]
  rownames(X)
  
  X=X[,2:1]
  colnames(X)
```

**Basic matrix operators**

```r
  X=matrix(nrow=3,ncol=2,data=1:6)
  Y=X
  
  X+Y # matrix addition, cell-by-cell
  X-Y # substraction
  log(X) # any function when called on a matrix it is applied to each of its cells
  X^2
  
  X*X # Note, this is the Haddamard (i.e. cell by cell product)
  
  # To obtain the matrix product use `%*%` instead of `*`
  A=X
  B=matrix(nrow=ncol(X),ncol=10,data=runif(ncol(X)*10)
  
  C=A%*%B
  
  # to verify that the result is correct
  sum(A[3,]*B[,2])==C[3,2]
```

**Apply function**

There is a family of functions `lapply`, `tapply`, `applay`, etc. that can be used to apply operations to dimensions of an array (of different kinds). For matrices we use `apply`, this function can be used to apply functions to rows or comumns of a matrix.

```r
n=1000;p=500
X=matrix(nrow=n,ncol=p,runif(n*p))

cSums<-apply(X=X,FUN=sum,MARGIN=2)
rSums<-apply(X=X,FUN=sum,MARGIN=1)

## passing your own function
  sumsOfLogs=function(x){ sum(log(x)) }
  tmp=apply(X=X,FUN=sumsOfLogs,MARGIN=2)

## column and row suma are already build in
  cSums2=colSums(X)
  rSums2=rowSums(X)

```

If we have a vector and an index set (e.g., male/female) we can apply a function to the vector for every level of the index using `tapply`.

```r
 x=rnorm(100)
 id=sample(c("M","F"),size=100,replace=T)
 tapply(X=x,INDEX=id,FUN=sum)
 sum(x[id=='M'])
 sum(x[id=='F'])
 

```
[Back to Outline](#Outline)

