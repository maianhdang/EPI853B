### Maximum Likelihood in the Binomial Model


**Due**: Tuesday Nov 1st. Bring a printed copy to the class.



In this HW we will use the [Gout](https://github.com/gdlc/EPI853B/blob/master/gout.md) dataset.

1. Bernoulli model with homogeneous success probability.

Let Y={y1,..,yn} be a random sample of n IID Bernoulli random variables with success probability theta.

1.1. Write the likelihood function

1.2. Derive the Maximum Likelihood estimator

1.3. Using the Gout data set and Y$Gout as data compute and report the ML estimate of the probability of developing Gout and an approximate 
95% CI.


Note: to derive a 95% CI you can use the fact that ML estimates follows, asymptotically, a normal distribution with mean equal to the true parameter
value. 

2. Logistic Regression

   Consider a logistic regression of Gout on Sex, Age, Ethnicity and Serum Urate.
   
   2.1. Write the likelihood function
   2.2. Create an R-function that evaluates the likleihood function
   2.3. Estimate parameters using optim
