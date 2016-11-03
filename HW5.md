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



## Examples of solutions

**Problem 1**

```R
library(splines)

 MC.REP=5000
 n=100
 x=seq(from=0,to=4*pi,length=n)
 signal=sin(x)+sin(x/2)
 
 
 
 nKnots=c(2,3,5,10,20)
 x0=5 
 

 YHat=matrix(nrow=MC.REP,ncol=length(DF))
 
 for(i in 1:MC.REP){
 		error=rnorm(n)
 		y=signal+error
 		for(j in 1:length(DF)){
 		  Z=bs(x=x,df=nKnots[j]+2,degree=2)
 		  bHat=coef(lm(y~Z))
 		  W=bs(x=x0,degree=attr(Z, "degree"), 
 		  		   knots=attr(Z, "knots"), 
 		           Boundary.knots=attr(Z, "Boundary.knots"),
 		           intercept=attr(Z, "intercept"))
 		           
 		   YHat[i,j]=sum(c(1,W)*bHat)
 		   
 		           
 	    }
}

BIAS=colMeans(YHat)-sin(x0)-sin(x0/2)	  
VAR=apply(X=YHat,MARGIN=2,FUN=var)

MSE=VAR+BIAS^2
plot(BIAS~nKnots)
plot(VAR~nKnots)
plot(MSE~nKots)

```

**Problem 2**

```R

Y=read.table('~/Desktop/galton.txt',header=T)

library(splines) 

Y$PA=(Y$Mother+Y$Father)/2	  
Y=Y[order(Y$PA),] # to get nice lines in the plot sort by x

fm0=lm(Height~ns(PA,df=5),data=Y) 

plot(Height~PA,data=Y,col=8)
lines(x=Y$PA,y=predict(fm0),col=2,lwd=2)
	
nRep=10000
YHat=matrix(nrow=nrow(Y),ncol=nRep)

for(i in 1:nRep){
	
	## 1st get a bootstrap sample
	tmp=sample(1:nrow(Y),size=nrow(Y),replace=T)
	BOOTS.DATA=Y[tmp,]
	
	
	## Now fit the model to that data
	Z=ns(BOOTS.DATA$PA,df=5)
	fm=lm(BOOTS.DATA$Height~Z)
	
	# Now we predict over a grid of values of x (in this case the original values for PA
	W=ns(x=Y$PA,  
 		  		   knots=attr(Z, "knots"), 
 		           Boundary.knots=attr(Z, "Boundary.knots"),
 		           intercept=attr(Z, "intercept"))
	YHat[,i]=cbind(1,W)%*%coef(fm)
	print(i)
}


TMP=apply(FUN=quantile,X=YHat,MARGIN=1,prob=c(.025,.975))

 lines(x=Y$PA,y=TMP[1,],col=4)
lines(x=Y$PA,y=TMP[2,],col=4)

```
