
### EPI-853b

### HW 5

**Due:** Tuesday Nov. 21st, in class (hard copy)


**Data**: Ovarian data set, can be obtained as follows


```r
 library(survival)
 data(ovarian)
 help(ovarian) # describes the variables
```


#### Q1. Kaplan-Meier (KM)

Develop an R-function that computes and plots KM estimates of survival curves. Your function should take as input a time variable and a status indicator (event/right censored) and produce results comparable to those of survfit()   and of plot(survfit()).

Test your function using ovarian$futime and ovarian$fustat as your time and event variables, respectively. Compare your results with those of survfit.


#### Q2. Standard errors in parametric survival regression

**2.1.** Use survreg to fit a parametric regression using a log-normal distribution, use ovarian$futime and ovarian$fustat as your time and event variables and ovarian$age, ovarian$resid.ds and ovarian$rx as predictor. Note, the last two variables are grouping factor. Report estimates, SE and p-values.


**2.2.** Conduct a bootstrap analyses (e.g., N=10K bootstrap samples) to estimate SEs, compare results with those of 2.1.
