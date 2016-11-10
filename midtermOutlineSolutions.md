## Midterm, outline of solutions

```R
## 1.1.
 load('~/Desktop/gout.RData')
 fm=lm(SBP~Sex+Race+Age,data=Y)
 summary(fm)

## 1.2. 
 age=c(rep(50,4),rep(65,4))
 sex=rep(c(1,1,0,0),2)
 race=rep(c(1,0),4)
 W=cbind(1,sex,race,age)
 yHat=W%*%coef(fm)

## 2
 X=as.matrix(model.matrix(~Sex+Race+Age,data=Y))
 varE=sum(residuals(fm)^2)/(nrow(X)-ncol(X))

 bHat=solve(crossprod(X),crossprod(X,Y$SBP))
 residuals=Y$SBP-X%*%bHat
 varE=sum(residuals^2)/(nrow(X)-ncol(X))
 SE=sqrt(diag(solve(crossprod(X)))*varE)
 cbind(bHat,SE)


## 3.1
 y=ifelse(Y$Gout=="Y",1,0)
 fm=glm(y~Sex+Race+UricAcid,data=Y,family=binomial)
 summary(fm)
 
## 3.2

negLogLik=function(y,X,beta){
	logOdds=X%*%beta
	theta=exp(logOdds)/(1+exp(logOdds))	
	logLik=sum(ifelse(y==1,log(theta),log(1-theta)))
	return(-logLik)
}


X=as.matrix(model.matrix(~Sex+Race+UricAcid,data=Y))

mu=log(mean(y)/(1-mean(y)))
theta=c(mu,rep(0,ncol(X)-1))

fm2=optim(y=y,X=X,fn=negLogLik,par=theta)
round(cbind(coef(fm),fm2$par),4)

## 3.3


fm2$convergence

## 3.4

BHat=matrix(nrow=3000,ncol=ncol(X))

for(i in 1:nrow(BHat)){
	tmp=sample(1:nrow(X),size=nrow(X),replace=T)
	
	tmpY=y[tmp]
	tmpX=X[tmp,]
	fit=glm(tmpY~tmpX-1,family=binomial)
	BHat[i,]=coef(fit)
}
SE=sqrt(apply(FUN=var,X=BHat,MARGIN=2))
round(cbind(summary(fm)$coef[,2] ,SE),4)


## 1.3. bonus
 X=as.matrix(model.matrix(~Sex+Race+Age,data=Y))
 varE=sum(residuals(fm)^2)/(nrow(X)-ncol(X))
 V=W%*%solve(crossprod(X))%*%t(W)*varE
 SE=sqrt(diag(V))
 round(cbind(yHat,yHat-1.96*SE,yHat+1.96*SE),2)



## 3.bonus
  
  Y=Y[order(Y$UricAcid),]   
  y=ifelse(Y$Gout=="Y",1,0)
  fm=glm(y~Sex+Race+UricAcid,data=Y,family=binomial)
  bHat=coef(fm)
  
  x=seq(from=3,to=15,by=1/100)
  
  ETA=cbind(sum(bHat[1:3])+bHat[4]*x,
  			sum(bHat[c(1,3)])+bHat[4]*x,
  			sum(bHat[c(1,2)])+bHat[4]*x,
  			sum(bHat[1])+bHat[4]*x)
  colnames(ETA)=c('M-W','F-W','M-B','F-B')
  THETA=exp(ETA)/(1+exp(ETA))
  
 plot(THETA[,1]~x,xlab='Serum Urate',ylab='Predicted Risk',ylim=range(THETA),type='l')
 for(i in 2:4){
 	lines(x=x,y=THETA[,i],col=i)
 }
 
 text(x=rep(4,4),y=c(.9,.85,.8,.75),col=1:4,label=colnames(ETA))
 
 ```
