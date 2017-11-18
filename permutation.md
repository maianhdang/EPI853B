
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

 mean(abs(testStat)>abs(summary(fm0)$coef[2,3]))
 summary(fm0)$coef[2,4]
 
```

### An example with 3 covariates

By permuting one co-variate at a time we can get p-values for hypotheses involving a single coefficient.


```r
DATA=read.table('~/Desktop/gout.txt',header=F)[1:50,]
X=as.matrix(model.matrix(~V1+V2+V3,data=DATA))[,-1]
y=DATA$V4
fm0=lm(y~X)
summary(fm0)

TSTAT=matrix(nrow=100000,ncol=3)

X1=X[,-1]
x1=X[,1]

X2=X[,-2]
x2=X[,2]

X3=X[,-3]
x3=X[,3]

n=nrow(X)

for(i in 1:nrow(TSTAT)){
	x1=sample(x1,replace=F,size=n)
	x2=sample(x2,replace=F,size=n)
	x3=sample(x3,replace=F,size=n)
	
	fm1=lm(y~x1+X1)
	fm2=lm(y~x2+X2)
	fm3=lm(y~x3+X3)
	
	TSTAT[i,1]=summary(fm1)$coef[2,3]
	TSTAT[i,2]=summary(fm2)$coef[2,3]
	TSTAT[i,3]=summary(fm3)$coef[2,3]
	print(i)
}

tStat=summary(fm0)$coef[2:4,3]

cbind(c(mean(abs(TSTAT[,1])>abs(tStat[1])),
		mean(abs(TSTAT[,2])>abs(tStat[2])),
		mean(abs(TSTAT[,3])>abs(tStat[3])) ), summary(fm0)$coef[2:4,4])
		

```

