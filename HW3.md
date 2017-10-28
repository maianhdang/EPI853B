## Homework 3

**Due**: Tuesday Oct 31st in class (hard copy).


**Goal**: to assess finite sample properties of ML estimates in logistic regression using Bootstrap.

**Data**: the gout data set.

**Model**: logistic regression of gout on age, sex, race and serum urate.

**Part 1**. Fit the logistic regression model above described using the whole data set. 

  1.1. Report parameter estimates, SEs and p-values.
  
  1.2. Sumarize your conclusions.
  
  1.3. Repot a table with the predicted risk of developing gout for a white male, age 65 and serum urate levels of 4, 6 and 9. 

**Part 2**. Using the model above-described and 10,000 Bootstrap samples estimate the expected value and the SE of the estimated effect 
of age for data sets with sample size 100, 200, 300, 500 and 1000.
  
  Note: with sample size 100, and also with sample size 200, you may have bootstrap samples that lead to non-existen ML estimates (i.e., flat likelihood). In this case the algorithm will not converge. You should discard results from models that did not converge. This could be done, for instance, using a structure like the following inside the loop that you use to generate bootstrap estimates.
  
  ```r
    ready=FALSE
    while(!ready){
      sample your data
      fit the model
      ready<- a check for converengence, e.g. fm$converged in glm  or $convergence if you use optim.
    }
  
  ```
  2.1. Report a plot with the estimated expected value vs. sample size and a plot of the estimated SE versus sample size. 
  
  2.2. Summarize your conclusions.
  
  2.3. Do we have evidence of finite sample bias?

**Part 3**. Using the whole data set, estimate, using 10,000 bootstrap samples, the sampling distribution of the 
estimated probability of developing gout for the cases of question 1.3. 

3.1. Report, in the same plot, estimated densities for each of the cases.

3.2. Summarize your conclusions.

**To read the data you can use the following code**

```r
rm(list=ls())
library(BGData)
load('../../genos2015/output/genosZ.RData')

# Reading the data
DATA=read.table('~/Dropbox/EPI853B/gout.txt',header=T,stringsAsFactors=F)

# Sample sizes we will consider
n=nrow(DATA)
N=c(100,200,300,500,1000)
DATA$gout=ifelse(DATA$gout=='Y',1,0)
```
