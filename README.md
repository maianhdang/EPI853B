

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

  * [The R-software](https://github.com/gdlc/EPI853B/blob/master/RIntro.md)
  * [Matrix Algebra in R](https://github.com/gdlc/EPI853B/blob/master/matrixAlgebraR.md)
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


