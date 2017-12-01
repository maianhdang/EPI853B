## HW5 EPI 853b

**Due**: Thursday, Dec. 7th (in class)


### Question 1: power analysis


A researcher wants to determine wether race has an effect on the likelihood of
developing a disease after accounting for differences due to sex and age.

The plan is to test this by using a logistic regression of the disease outcome on sex, age and race. Significance will be assessed
based on the p-value for the race effect in this logistic regression.

The researcher wants to determine the type-I error rate and power of such a study as a function
of sample size and the size of the effect of race.


**Task**: conduct a simulation study (with sufficient number of Monte Carlo replicates) to estimate
the type-I error rate and power for  each of the scenarios that can be formed by
combining sample size of N=50,100, 200, 300 and 500 and race effects 0 (null), .01,.05,.1,.2,...,.8).


**Parameters**: In your simulation assume equal proportion of male/female, equal proportion of black/white
and assume that age~N(50,20). Assume that these three predictors are independent.

**(i) Report power curves** (in a single plot) with estimated rejection rate in the vertical axis versus sample size, by size of effect.

**Interpret your results**

**Report your code in an appendix**



### Question 2: multiple testing

Use the following simulation to estimate the power and false discovery rate single-prediction linear regression for
following decision rules:

  - Reject based on Bonferroni-adjusted p-values at a target family-wise error rate of 0.05
  - Reject based on FDR-adjusted p-values with target FDR of 5 and 10%.


Usue the following code to generate each of your Monte Carlo replicates.

```r

         X=matrix(nrow=n,ncol=p,rbinom(n=n*p,size=2,p=.25))
	 b=rep(0,p)
	 b[sample(1:p,size=q)]=runif(min=.1,max=.3,n=q) # q coefficients come from Ha, the rest from H0
		
	 signal=X%*%b
	 error=rnorm(n=n,sd=sd(signal)*sqrt(2/3))
	 y=signal+error

```
