### Type-I error rate

Suppose we have a test statistic `T(y)` and a decision rule: reject if `T(Y)>k`. 

The type-I error rate such decision rule is the probability of rejecting given that H0 holds, that is `P(T(Y)>k|H0)`. Here `T(Y)`, the test-statistic
is a function of the data. From a Frequentist perspective data, `Y`, is random, therfore `T(Y)` is a random variable. 
What threshold (`k`) should we choose? The use of p-values circunvents this problem: we set `k` to be equal to the obseved value of the test statistic, 
`k=T(y)` (small-letter `y` is used to denote the realized sample). The p-value of the test is `P(T(Y)>T(y)|H0)`.  If we reject at 0.05 we expect that 
the corresponding test would have a 5% Type-I error rate.

**Example 1:** Estimating type-I error rate for a single test

```r
 n=100 # sample size
 nRep=10000
 
 reject<-rep(NA,nRep)
 x=rbinom(n=n,size=2,p=.3)
 
 for(i in 1:nRep){
    y=rnorm(n) # simulating under the null
    fm=lm(y~x)
    reject[i]=summary(lm(y~x))$coef[2,4]<.05
 }
 mean(reject)
```

### Multiple testing and experiment-(or family-)wise error rate

Suppse we conduct 3 tests simultaneously. The experiment-wise error rate is the probability of making at least 1 mistake, that is:
`p(T1(Y)> T1(y) or T2(Y)> T1(y) or T3(Y)> T1(y) | H0)`. If we reject at each test with 0.05 significance and the tests are independent, then 
the probability of making at least one mistke is 3x0.05=0.15. The following example illustrates this with a simulation. Because the columns of `X`
are independent, the test are also independent.

```r
 n=100 # sample size
 nRep=10000
 
 reject<-matrix(NA,nrow=nRep,ncol=3)
 X=matrix(nrow=n,ncol=3,data=rbinom(n=3*n,size=2,p=.3))
 
 for(i in 1:nRep){
    y=rnorm(n) # simulating under the null
    X=matrix(nrow=n,ncol=3,data=rbinom(n=3*n,size=2,p=.3))
    pValues=ls.print(fm,print.it=F)$coef[[1]][-1,4]
    fm=lsfit(y,x=X)
    reject[i,]= pValues <.05
    print(i)
 }
 colMeans(reject) # each test has 0.05 type-I error rate
 anyRejection= reject[,1]|reject[,2]|reject[,3]
 mean(anyRejection)
```
 

