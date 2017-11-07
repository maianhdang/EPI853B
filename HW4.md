### HW4: Survival Regression

Due Tuesday, Nov 14 in class.



Using the ovarian data set

(1) Fit and report Kaplan-Meier estimates of survival curves for individuals under treatment and control (hint, use `survfit`, consider `summary(survit(Y~...))` and `plot(survfit(...))`).


(2) Fit a parametric survival regression assuming a Gaussian likelihood using survreg. Include all the predictors available in the data set


(3) Obtain the results of (2), parameter estimates, SE, p-values, using optim.


