### Ordnary Least Squares (OLS)

Consider a linear regression of the form **y**=**Xb**+**e**. The OLS estimates of the vector of regression coefficients is given by
<br />
<br />
      **bHat**=argmin{  RSS(**y**,**X**,**b**) }
<br />
<br />

Where:   RSS(**y**,**X**,**b**)=(**y**-**Xb**)'(**y**-**Xb**) is the residual sum of squares.
<br />
<br />

The solution of the above problem can be obtained from the following systems of equatrions
<br />
<br />
(**X**'**X**)**b**=**X**'**y**
<br />
<br />

####(1) Estimation

**Computation of OLS estimates using `lm`, `lsfit` and with matrix operations.**


**Computation of OLS estimates using various decompositions (`chol`, `qr`,`svd`).**


####(2) Inference

**SE, CIs, p-values, F and t-test from asymptotic normality**


**SE and CIs using Bootstrap**
   
