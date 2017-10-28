## Homework 3

**Due**: Tuesday Oct 31st in class (hard copy, please no e-mails).


**Goal**: to assess finite sample properties of ML estimates in logistic regression using Bootstrap.

**Data**: the gout data set.

**Model**: logistic regression of gout on age, sex, race and serum urate.

**Part 1**. Fit the logistic regression model above described using the whole data set. 

  1.1. Report parameter estimates, SEs and p-values.
  
  1.2. Sumarize your conclusions.
  
  1.3. Repot a table with the predicted risk of developing gout for a white male, age 65 and serum urate levels of 4, 6 and 9. 

**Part 2**. Using the model above-described and 10,000 Bootstrap samples estimate the expected value and the SE of the estimated effect 
of age for data sets with sample size 100, 200, 300, 500 and 1000.

Report a plot with the estimated expected value vs. sample size and a plot of the estimated SE versus sample size. 

**Part 3**. For a sample size of 1000, estimate, using 10,000 bootstrap samples, the sampling distribution of the 
estimated probability of developing gout for the cases of question 1.3. Report, in the same plot, estimated densities for each of the cases.

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
