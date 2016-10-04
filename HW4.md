### HW4: Inferences in OLS


###Homework 3
(due Monday 10/10/2016, submit va e-mail to gustavoc@msu.edu )

**Tasks**: 

  -1- Create a function that fits OLS estimates and retunrs a table equivalent to the summary table of lm. Your function `getOLS(y,X)` should take as arguments
  and response vector and an incidence matrix for effects. Internally, the function should produce OLS estimates, SE, t-value and p-value, just like the
  table you will get from `summary(lm(y~X))$coef`. Internally you are only allowed to use matrix operations, you cannot use `lm` or `lsfit` or similar.
  
  -2- Using the code provided below carry out 5000 MC simulations and estimate the bias, variance and MSE of OLS estimates using MC methods.
      Embed your code in a function `MC.OLS(X,beta,nRep,R2)` should take as arguments a design matrix, 'true effects', the # of mC replicates and the model
      R2. Your function must return a table with colums corresponding to the columsn of X and in rows the estimated bias, variance and MSE of OLS estimates.
      
**Simulation**

```R
 X1=sample(1:3,size=N,replace=T)
 X2=sample(1:2,size=N,replace=T)
 b1=c(100,80,120)
 b2=c(10,20)
 signal=b1[X1]+b2[X2]
 error=rnorm(N,sd=sd(signal)*sqrt(1-R2))
 y=error+signal
```

**Tests**
```R
#Test 1
 all(round(summary(lm(y~X-1))$coef,5)==round(getOLS(y,X),5))
 
#Test 2
 all(round(ourFunctionForMCSimulations(X,beta,nRep,R2),5)==round(MC.OLS(X,beta,nRep,R2),5))
```
