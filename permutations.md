
## Permutation test
A simple example: testing sex differences using a permutation test.

**Simulating a small data set**

```r
n=100
sex=sample(0:1,size=n,replace=T)
y=10+sex*.2+rnorm(n)
```


**Fitting the model to the original data**

```r
 fm0=lm(y~sex)

```

**Estimating a p-value using a permutation test**

```r
 nSamples=10000
 testStat=rep(NA,nSamples)
 n=length(y)

 for(i in 1:nSamples){
   z=sex[order(runif(n))]
   fm=lm(y~z)
   testStat[i]=summary(fm)$coef[2,3]
   print(i)
 }

 isBelow=testStat<  -summary(fm0)$coef[2,3]
 isAbove=testStat>  summary(fm0)$coef[2,3]

 mean(isBelow+isAbove) # compare this with summary(fm0)$coef[2,3]

```

