
### Family-Wise Error Reate, False Discovery Rate (FDR) and power 


Consider testing N (N=N1+N2+N3+N4) hypotheses. For each hypothesis H0 may hold or not (rows in the table below),
and we may reject the null or fail to reject it (columns).

|           | Do not reject H0  | Reject H0          |
|-----------|-------------------|---------------------|
| H0 holds  | True Negative (N1) | False Positive (N2)|
| Ha holds  | False Negative (N3)| True positive (N4) |


From the above table we can compute several important probabilities:

   - Type-I error rate: P(reject|H0 holds)=N2/(N1+N2)
   - False-discovery rate (i.e., proportion of cases for which H0 holds among all rejections): FDR=N2/(N2+N4)
   - Power: P(Reject|Ha)=N4/(N3+N4)
   


#### Estimating Type-I error rate, False Discover Rate and Power using Monte Carlo Simulations

##### (1) Single-test case

The following example shows how to estimate type-I error rate, power and FDR for a single test. 

We simulate under diffeerent values 

```r
  R2=c(0,.01,.03,.05,.1) # Model R-sq.
  N=50 # sample size
  nRep=10000 # number of Monte Carlo replicates
  significance=0.05 # significance for rejection
   
  countRejections=rep(0, length(R2)) # We count rejections for every scenario
  for(i in 1:nRep){
      x=rnorm(N)
      for(j in 1:length(R2)){
        signal=x*sqrt(R2[j])
        error=rnorm(sd=sqrt(1-R2[j]),n=N) 
        y=signal+error
        fm=lsfit(y=y,x=x)     
        countRejections[j]=countRejections[j]+(ls.print(fm,print.it = F)$coef[[1]][2,4]<significance)
      }
      if(i%%100==0){print(i)}
  }
  plot(y=countRejections/nRep,type='o',col=2,x=R2,ylab='Power',xlab='R2',ylim=c(0,1))
  abline(h=significance,col=4,lty=2,main='Power Curve',ylim=c(0,1))
```


**Task**: Modify the code to estimate power as a function of R-sq and sample size, for N=30, 50, 100,500. Produce a plot of power versus sample size, by R-sq (i.e., different power curves per R-sq. level).

#### (2) Multiple Testing

#### Distribution of p-values under H0 and under Ha


**A toy-simulation**

```r
  N=100 
  p=50000
  X=matrix(nrow=N,ncol=p,data=rbinom(size=2,p=.3,N*p))
  hasEffect=rep(F,p) 
  hasEffect[sample(1:p,size=500)]=T
  pValues=rep(NA,p)
  
  for(i in 1:p){ 
 	x=X[,i]
 	if(hasEffect[i]){ 
 	  b=rnorm(n=1,mean=.05,sd=1)
  	  signal=x*b
      error=rnorm(n=N,sd=1)
      y=signal+error
    }else{
      y=rnorm(N)
    }
    fm=lsfit(y=y,x=x)     
    pValues[i]=ls.print(fm,print.it = F)$coef[[1]][2,4]
    if(i%%10000==0){print(i)}
  }
```

**Distribution plots**

```r
  par(mfrow=c(2,2))

  # p-values under the null follow a uniform distribution
  hist(pValues[!hasEffect],100,main='H0')

  # under Ha the distribution is not uniform
  hist(pValues[hasEffect],100,main='Ha')
  
  # overall the distrubtion is a mixture
  hist(pValues,100,main='Overall')


```  


