
### Family-wise error reate, False Discovery Rate (FDR) and power 


Considerb testing N=N1+N2+N3+N4 hypotheses. For each hypothesis H0 may hold or not (rows in the table below),
and we may reject/fail to reject H0 (columns).

+-----------+--------------------+--------------------+
|           | Do not reject H0  | Reject H0          |
+===========+====================+====================+
| H0 holds  | True Negative (N1) | False Positive (N2)|
+-----------+--------------------+--------------------+
| Ha holds  | False Negative (N3)| True positive (N4) |
+-----------+--------------------+--------------------+

From the above table we can compute several important probabilities:

   - Type-I error rate: P(reject|H0 holds)=N2/(N1+N2)
   - False-discovery rate (i.e., proportion of cases for which H0 holds among all rejections): FDR=N2/(N2+N4)
   - Power: P(Reject|Ha)=N4/(N3+N4)
   


#### Estimating Type-I error rate, False Discover Rate and Power using Monte Carlo Simulations

##### Single-test case

We can simulate data under H0 and under various Ha, apply our decision rule, and cout how many cases fall in each of the cells of the table above-described. The following example illustrates this.


```r
  h2=c(0,.001,.005,.01,.02,.05,.1)# heritability parameter
  N=100 # sample size
  nRep=1000 # number of Monte Carlo replicates
  significance=0.01 # significance for rejection
  
  countRejections=rep(0, length(h2)) # We count rejections for every scenario
  for(i in 1:nRep){
      x=rbinom(size=2,n=N,p=.2) # we assume effect=1, and scale errors to get the desired h2
      Vg=var(x)
      Ve=Vg/h2*(1-h2)
      for(j in 1:length(h2)){
        if(j==1){ y=rnorm(N) }# simulating under the null
        if(j>1){ y=x+rnorm(sd=sqrt(Ve[j]),n=N) }
        fm=lsfit(y=y,x=x)     
        countRejections[j]=countRejections[j]+(ls.print(fm,print.it = F)$coef[[1]][2,4]<significance)
      }
      #print(i)
  }
  plot(y=countRejections/nRep,type='o',col=2,x=h2);abline(h=significance,col=4,lty=2,main='Power Curve',ylim=c(0,1))
```
