## Maximum Likelihood

The likelihood function is the probability of the data viewed as a function of the parameters.


For instance, consider IID bernoulli data with success probability p(xi=1)=theta, p(xi= 0)=1-theta.

The likelihood function is: `L=dbinom(x=x,size=1,prob= ?)`. In Maximum Likelihood, we find the value of theta that 
maximizes this function, with `x=[1,0,0,1,0,...]`, your data, fixed. 

The solution to an opimization problem is invariant under monotonic transformations (addition of a constant, multiplying by a postivie constant, logarithm, square-root) will
of the objective function (the likelihood). Futhermore, the argument that maximizes a function, `f(theta)` is equal to the value that minimizes
the negative of the function, `-f(theta)`. Thus, in practice, instead of maximizing the likelihood we minmize the negative of the log-likelihood.

In the IID bernoulli example the log-likelihood has a very simple form (discussed in class)

```r
  logLik=function(x,theta){
    logLik=mean(x)*theta + (1-mean(x))*(1-theta)
    return( logLik)
  }
```

#### Obtaining ML estimates using grid search

In the above example we showed in class that the ML estimate of the success probability is the sample mean. However, in most
cases the ML estimator does not have a closed form. In such cases we can obtain ML estimates using grid search or iterative algorithms.
Grid search simply evaluates the objective function of a grid of value of the parameter. This of course works when the number of parameters
to be estimated is small (one in the bernoulli case).


```r
## Simulating IID Bernulli data
 theta0=.2 # True value
 n=20
 x=rbinom(size=1,n=n,prob=theta0)

## Grid search

  # A function to evaluate the logarith of the likelhiood
	myF=function(x,theta){
       logLik=mean(x)*log(theta) + (1-mean(x))*log(1-theta)
        return(logLik)
    }
myF(x,.01)
myF(x,.05)

## Grid of values of theta
 theta=seq(from=0,to=1,length=1000)

## Evaluate the log-likelihood over values of theta
 logLik=rep(NA,length(theta))
 for(i in 1:length(theta)){
    logLik[i]=myF(x,theta[i])
 }

## Plot
plot(logLik~theta,col=2,type='l')
abline(v=which(logLik==max(logLik)),col=4)
```

**Suggested excercise**: Repeat the above problem for 100 samples, each of size 20. For every sample you will get a different log-likelihood function. Plot all the log-likelihood functions and add vertical lines denoting the ML estimate for each sample. Repeat this with N=1000. Observe how the concavity of the log-likelihood (and consequently the variance of the ML estimates) increases (decreases) with sample size.


#### Obtaining ML estimates using general-purpouses optimizers

```r
 ## Write a function to evaluate the negative of the log-likelihood
  negLogLik=function(x,theta){
   logLik=mean(x)*log(theta) + (1-mean(x))*log(1-theta)
   return(-logLik)
  }
  
 ## Pass it to optmize, together with the data an an interval for the search.
   optimize(f=negLogLik,x=x,interval=c(0,1))

```
