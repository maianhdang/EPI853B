
###Homework 2 
(due Tuesday 9/20/2016, submit va e-mail to gustavoc@msu.edu )

1. Solving a "triangular" system of linear equations.

  Goal: to solve `Cx=r` for x, assuming that C is: (i) lower-triangular and (ii) upper-triangular.
  Tasks:
    - Create two R functions ( `backsolve_2(C,r)` and `forwardsolve_2(C,r)` ) that solve upper and lower-triangular systems using a loop.
    - Show with a small example that your function delivers the right solution (compare with either `solve(C,r)` or `backsolve()` or `forwardsolve()`.

2. Benchmark of the different strategies we discussed for computing OLS estimates
  - Goal: to estimate the computing time needed to obtain OLS using the different methods we discuss in clase ('lm`,`lsfit`,`solve`,`chol`,`qr`,`svd`).
  - Task: compute for an array of values of `n` and `p` the time it takes to compute OLS for each of the mehtods listed.
  - Report a table per method containing the average time, relative to lm, for each of the methods.
  
  Note: for time estimation it is usually good to run each process multiple times, especially if the time of the task is small.
  You can time your examples as follows:
  
```R
  timeIn=proc.time()
    # A LOOP HERE WITH THE EVALUATION OF YOUR TASK
  timeOut=proc.time()
  timeOut-timeIn
```
 
 
```R
  nRep=c(1e3,100,50,10)
  n=c(100,1e3,1e4,1e5)
  p=c(10,100,1000,1e4)
  GRID=cbind(nRep,n,p)
```
[Back to Outline](https://github.com/gdlc/EPI853B/#Outline)
