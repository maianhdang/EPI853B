### HW5: Splines

 (due Monday 10/17/2016, submit via e-mail to gustavoc@msu.edu )

**Tasks**: 

  -1- Evaluate, using MC simulations (at least 1000 replicates) the bias, variance and MSE for a cubic natural spline 
  as a function of the number of knots used. Evaluate, variance, bias and MSE at `x=c(2,5,11)` for 2,3,5,20,20,50 knots. 
  Report also the MSE, variance and bias at the same points for a quadratic polynomial. Report one plot per value of x (2,5,11)
  with the variance, bias and MSE versus # knots. Add in the same plot with a horizontal line the estimated variance, bias and MSE
  of the quadratic polynomial.
  
  -2- Use cross-validation methods (100 training-testing partitions with 30 points in testing and 70 in training) to choose the optimal number of knots of a natural spline as a function of sample size.
  For this, use `n=c(10,20,50,100)`, `knots=seq(from=2,to=30)`; for each sample size estimate the prediction MSE and plot
  prediction MSE versus number of knots. Indicate the optimal number of knots.

Hint: you can get the basis functions for a natural spline using `ns()` from the `splines` library.
