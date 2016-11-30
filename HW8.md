##HW8: Estimation of power and Type-I error rate.

Due: Tuesday Dec. 6. Printed copy at the class.

**1.** Consider a binary outcome with success probability (theta) being 

           theta_i=exp(mu+x_i*b)/(1+exp(mu+x_i*b))
           
Above: theta_i is the success probability of the ith individual, mu=.2 is an intercept,
x_i is a covariate and b is an effect. Assume x_i~norm(mean=0,sd=1).

Develop a MC simulation for estimating type-I error rate assuming a sample size of N=100
and a significance level of .05.

Provide your code, your results and a 1-sentence conclusion.  
 
**2.** Estimate the power of the experiment for b=c(.1,.2,.4) and N=30,50,100.

Report the estimated power and summarize your conclusions. 
  
**3.** Extend the setting of problem 1 by adding a second covariate, z_i~rnorm(0,1), 
independent of x_i.

The goal is to test the effect of x_i while accounting for z_i. Assume 

           theta_i=exp(mu+x_i*b1+z_i*b2 )/(1+exp(mu+x_i*b1+ z_i*b2 ))
with b2=.2

Develop a MC simulation for estimating power and type-I error rate. For power analysis
consider the same grid of values of effect (b1=c(.1,.2,.4) and sample size (N=30,50,100).

**4.** Multiple testing

Consider the setting of problem 2, but now modified such that x_i and z_i are correlated (see below).

Develop a MC study to estimate the family-wise error rate for an experiment with
sample size N=100 when rejection is done using a nominal alpha-level of 0.05. 

Report your code and the estimated family-wise error rate.

To simulate two normal random variables with null mean, unit variance and correlation 
of 0.8 you can use the following code

```R
 Z=matrix( nrow=N,ncol=2,rnorm(2*N)) # 2 iid normal random variables  
 COV=matrix(nrow=ncol(Z),ncol=ncol(Z),.9)
 diag(COV)  
 X=Z%*%chol(COV)
```

**5.** Estimate the nominal alpha-level that you should use to achieve a family-wise error rate of
0.05.

Report your code and the estimated alpha.
