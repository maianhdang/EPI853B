
### Family-wise error reate, False Discovery Rate (FDR) and power 


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

