### Parallel and GPU Computing in R

Most computers are equiped with with multiple processors and have a graphical proccesing unit. So far we have carried out computations in a single processor.
Here we explore ways to run multiple tasks in parallel. Parallelism (multiple processes running a the same time) can happen within a computer
(e.g., multi-core computing) or between computers in a network (e.g., in a High Performance Computing Cluster). We focus on the first case.


Multiple packages provide functons for parallel computing, a summary is presented in the following 
[CRAN-view](https://cran.r-project.org/web/views/HighPerformanceComputing.html). 
The [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf) package provides various option for multi-core computing.


**Detecting the number of cores available**

```r
 library(parallel)
 detectCores()
```

**Using `mclapply`**


The following example (adapted from the mclapply help) illustartes how to use `mclapply`. In the example we use the function `rnorm`, the 
only argument to `rnorm` without default value is `n`, the numbrer of samples to be collected. The example calls `rnorm(n)` with
values of `n=C(2,3,5,4)` in parallel at 4 cores.
```r
  n<-c(2,3,5,4)
  mclapply(X=n,rnorm,mc.cores=4)
```

We can use the same idea to implement any task that can be splitted into chunks. For instance, we can carry out the sum of the entries
of a vector into chunks, apply the sum to each chunk and then aggregate results (we discuss aggregation in the next example).

```r
 x=rnorm(4e6)
 sum(x) 
 
 index=rep(1:4,each=1e6)  # we will use this to split our task in 8 chunks
   
 sumChunk=function(x,index,j){
  	sum(x[index==j])
 }
 
 sumChunk(x,index,1)
 sum(x[index==1])
 
 mclapply(FUN=sumChunk,x=x,index=index,X=1:4,mc.cores=4)
  
```

Generating a function that calls `mclapply` and aggregates results

```r
  pSum<-function(x,nTasks,ncores){
    n<-length(x)
    index=rep(1:nTasks,each=ceiling(n/nTasks))[1:n]
    sums<- mclapply(FUN=sumChunk,x=x,index=index,X=1:nTasks,mc.cores=ncores)
    sum(unlist(sums))
  }
  
  sum(x)
  pSum(x,nTasks=8,ncores=4)
  
```

In the previous example the computations done with `mclapply` took more time than the serial computation. Carrying out computations in parallel has an overhead. Parallel computing is convinient when each task takes a substantial amount of time. 
We illustrate this implementing a function that computes crossproducts in parallel.

Consider the operation `XtY<-crossprod(X,Y)` for the case where `X` has a large number of columns. The dimensions of `XtY` are `dim(XtY)=c(ncol(X),ncol(Y)`. Furthermore, the ith row of `XtY`  corespond to `crossprod(X[,i],Y)`. More in general, if `cols` represent a set of columns of `X` we have that `XtY[cols,]=crossprod(X[,cols],Y)`. The following code verifies this numerically.


```r
  n=10000
  p=500
  q=10
  X=matrix(nrow=n,ncol=p,rnorm(n*p))
  Y=matrix(nrow=n,ncol=q,rnorm(n*q))
  cols=sample(1:p,10)
  
  XtYa<-crossprod(X,Y)[cols,]
  XtYb<-crossprod(X[,cols],Y)
  all(XtYa==XtYb)
```

We can then compute crossproducts for subsets of columns of X in parallel and then use `rbind` to put together all the results. This idea is implemented in the following `pCrossprod` function.

```r
  crossprod.chunk<-function(W,Y,index,chunk){
    W<-W[,index==chunk]
    crossprod(W,Y)
  }
  
  # Then we produce the wrapper
  pCrossprod<-function(X,Y,nTasks,nCores){ 
     n=ncol(X)
     index=rep(1:nTasks,each=ceiling(n/nTasks))[1:n]
     results<- mclapply(FUN=crossprod.chunk,W=X,Y=Y,index=index,X=1:nTasks,mc.cores=nCores)
     XY=matrix(nrow=ncol(X),ncol=ncol(Y))
     lastRow<-0
     for(i in 1:length(results)){
        firstRow<-lastRow+1
        lastRow<-firstRow+nrow(results[[i]])-1
        XY[firstRow:lastRow,]<-results[[i]]
     }
     
     return(XY)
  }

```

The following example illustrates the use of the `pCrossprod` function. The following example does not show great speed-up, the reason being that the number of columns of Y is relatively small, hence, `crossprod` is very fast and the overhead of parallel computing is too high. 

```r
  library(parallel)
  n=20000
  p=5000
  q=10
  X=matrix(nrow=n,ncol=p,rnorm(n*p))
  Y=matrix(nrow=n,ncol=q,rnorm(n*q))
  
  # Example 1: X'Y, no speed-up with parallel computing because q is small.
   system.time(XY1<-crossprod(X,Y) )
   system.time(XY2<-pCrossprod(X,Y,nTasks=12,nCores=4) )
 
  # Example 1: X'X, significant speed-up, in my computer by a factor of ~2.
   system.time(XX1<-crossprod(X) )
   system.time(XX2<-pCrossprod(X,X,nTasks=12,nCores=4) )
```

In the follwoing example parallel computing can get a significant speedup. We use the `crossprod_parallel` of the `BGData` package.

```r

system.time(XX<-crossprod(X) )
system.time(XX2<-crossprod_parallel(X,nCores=4,nTasks=12))

```

The `BGData` package offers `tcrossprod_parallel(), `crossprod_parallel()` and parallel versions of apply.


### GPU-computing

The `gpuR` package offer functions for matrix operations using GPUs. In this example GPU computing gives a speed-up by a factor of 15+!!!


```r
 library(gpuR)
 Z<-gpuMatrix(X,'float')
 system.time(ZZ<-crossprod(Z))
 
```
