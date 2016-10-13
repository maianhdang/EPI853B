### HW5: Splines

 (due Monday 10/17/2016, submit via e-mail to gustavoc@msu.edu )

**Problem 1**: Variance, Bias and MSE of spline and of cuadratic polynomials

Using the following data generating process

```R
 n=100
 x=seq(from=0,to=4*pi,length=n)
 signal=sin(x)+sin(x/2)
 error=rnorm(n)
 y=signal+error
 plot(y~x,col=4)
 lines(x=x,y=signal,col=2,lwd=2)
```

Generate 5000 MC samples and estimate the Variance, Bias and MSE of predictions at `x=5` for a quadratic spline with 2,3,5,10,20 knots. Do the same for a regular quadratic polynomial. **Report** a plot per quantity of interest (variance, bias and MSE) with # of knots in the horizontal axis and the estimated quantity (either variance, bias or MSE) in the vertical axis for the spline. Add, with a horizontal line the same quantity for the quadratic polynomial.

Summarize your conclusions.


**Problem 2** (Bootstrap)

Using [Galton’s height data]( http://www.math.uah.edu/stat/data/Galton.html), regress height on the average height of the two parents using a natural spline with 5 DF `ns(x=x,df=5`). Using 10,000 Bootstrap samples estimate a 95% confidence band. **Report** a plot with average parental height in the horizontal axis and height of the offspring in the vertical axis. Add as lines your estimated regression function and 95% confidence bands.

