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

```
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

