Matrix factorisations are used often for (i) obtaining OLS estimates, (ii) inverting matrices, (iii) determining the rank of a matrix, 
(iii) describing features of the data and (iv) in simulations. 


The following table described often used matrix factorisations. 

|  Name   | Description  | R-function |
|-----------|-----------|----------|
| Singular Value Decomposition           |  X=UDV' applies to any matrix        |  `svd`        |
| Eigen Decomposition | X=UDU' applies to symmetric matrices  |  `eigen` |
| QR | X=QR applies to any marrix  |  `qr` |
| Cholesky| X=U'U applies to symmetric matrices  |  `chol` |

### (1) Singular-value decomposition

Finds orthonormal-basis for the row (**U**) and column (**V**) linear spaces spanned by **X**.

**1.1. Computation**

```r 
  X=matrix(nrow=20,ncol=4,data=runif(80))
  SVD=svd(X)
  str(X) 
  dim(SVD$u)    # left-singular vectors
  length(SVD$d) # the singular values
  dim(SVD$v)    # righ-singular vectors
```
**1.2. Veriffying properties **

```r
  # Recovering X from the factorization
   Xnew=SVD$u%*%diag(SVD$d)%*%t(SVD$v) 
   all(round(Xnew,8)==round(X,8))
  
  # Rank of the matrix is equal to the number of positive singular values
   sum(SVD$>1e-8)
   
  # U and V have orthonormal columns
  round(crossprod(SVD$u),8)
  round(crossprod(SVD$v),8)
  
  # For full-column rank matrices (as in this example) the rows or V are also orthonormal 
  # (for full-row rank matrices UU'=I, but this does not happen in this case because rank(X)<nrow(X)
  round(tcrossprod(SVD$v),8)
```

The left-singular vectors (U) can be used to describe features of the rows (subjects) and the right-singular vectors (V)
can be used to describe patterns of the columns of X. Thus, the first two left-singular vectors are often used to describe
data structure, you can use `plot(SVD$u[,1:2])` to display this. Likewise `plot(SVD$v[,1:2])` can be used to uncover
relationships among columns. Biplots mixes these two features in a single plot.

**1.2. Solving systems of equations and obtaining OLS estimates**

Consider `y=Xb+e`, replace `X` with it's SVD decomposition, that is, write: `y=UDV'b+e`. Let `d=DV'b`, then `y=Ud+e`. 
This is a regression of `y` on an orthonormal basis for X. Because `U'U=I` the OLS estimate of `d` is `dHat=U'y`. This orthogonal projection can always be pefrormed, even if `X` is rank-deficient. In the full-column-rank case, we have `bHat=VD^-1dHat`. The following example verifies this.

```r
  # OLS estimates
   y=rnorm(nrow(X))
   XX=crossprod(X)
   rhs=crossprod(X,y)
   bHat=solve(XX,rhs)

  # Now the same using SVD
   dHat=crossprod(SVD$u,y)
   bHat2=SVD$v%*%diag(1/SVD$d)%*%dHat
   cbind(bHat,bHat2)
```


### (1) QR-decomposition

This factorization is commonly used to obtain OLS estimates.  The following example shows how to obtain the QR decomposition.

```r
  QR=qr(X)
  Q=qr.Q(QR) # gives you the Q matrix
  R=qr.R(QR) # gives you the R matrix
  
  ## Verifying orthogonality
  round(crossprod(Q), 8)
  
  ## Recovering X
  XNew=Q%*%R
  
  all(round(XNew,8)==round(X,8))
  
```

Following the same ideas we discussed before, substitute in the linear model `y=Xb+e` `X` with the QR-decomposition, `y=QRb+e`, let `d=Rb`, to get `y=Qd+e`. Since `Q` has orthogonal columns, the OLS estimate of d is `dHat=Q'y`. And, in the full-colmun rank case, `bHat=RInv*dHat`. Because `R` is triangular, the inverse can be computed easily using a recursive algorithm (discussed in class).

```r
 dHat=crossprod(Q,y)
 bHat3= backsolve(R,dHat) # solves a system of equations Rb=dHat, for the case where the coefficient matrix, R is triangular.
 cbind(bHat,bHat2, bHat3)

```


