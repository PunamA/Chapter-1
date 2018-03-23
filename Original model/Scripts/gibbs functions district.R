tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}
#----------------------------------------------------------------------------------------------
acceptMH <- function(p0,p1,x0,x1,BLOCK){   #accept for M, M-H
  # if BLOCK, then accept as a block,
  # otherwise, accept individually
  
  nz           <- length(x0)  #no. to accept
  if(BLOCK) nz <- 1
  
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)
  
  if(BLOCK & length(keep) > 0) x0 <- x1
  if(!BLOCK)                   x0[keep] <- x1[keep]           
  accept <- length(keep)        
  
  list(x = x0, accept = accept)
}
#-------------------------------
rmvnorm1=function (n, sigma, pre0.9_9994 = FALSE) 
{
  #   retval <- chol(sigma, pivot = TRUE)
  #   o <- order(attr(retval, "pivot"))
  #   retval <- retval[, o]
  s. <- svd(sigma)
  if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
    warning("sigma is numerically not positive definite")
  }
  R = t(s.$v %*% (t(s.$u) * sqrt(s.$d)))
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*% R
  retval
}
#----------------------------
create.dmat=function(dat){
  nloc=length(unique(dat$village))
  dmat=matrix(0,nrow(dat),nloc)
  for (i in 1:nloc){
    cond=dat$village==i
    dmat[cond,i]=1
  }
  dmat
  
#   teste=apply(dmat,2,sum)
#   plot(teste,table(dat$village))
#   unique(teste-table(dat$village))
}
#--------------------------
update.alpha=function(param){
  prec=dtd+(1/param$sig2)*param$invSigma
  var1=solve(prec)
  err=param$z-xmat%*%param$betas
  pmedia=t(dmat)%*%err+(1/param$sig2)*param$invSigma%*%wmat%*%param$gamma
  t(rmvnorm(1,var1%*%pmedia,var1))
}
#--------------------------
update.rho=function(param){
  err=param$alpha-wmat%*%param$gamma
  k=1/2

  res=rep(NA,length(rho.vals))
  for (i in 1:length(rho.vals)){
    res[i]=t(err)%*%invSigma.vals[[i]]%*%err
  }
  res=-(1/2)*ldet.vals-((nvil-1)/2)*log(res)
  res=exp(res-max(res))
  res1=res/sum(res)
  tmp=rmultinom(1,size=1,res1)
  ind=which(tmp==1)
  return(list(rho=rho.vals[ind],invSigma=invSigma.vals[[ind]]))
}
#--------------------------
update.gamma=function(param){
  k=c(10,rep(1,ncol(wmat)-1))
  invTmat=diag(1/k)
  invSigma=(1/param$sig2)*param$invSigma
  prec=t(wmat)%*%invSigma%*%wmat+invTmat
  var1=solve(prec)
  pmedia=t(wmat)%*%invSigma%*%param$alpha
  t(rmvnorm(1,var1%*%pmedia,var1))
}
#--------------------------
update.betas=function(param){
  prec=xtx+diag(1,ncol(xmat))
  var1=solve(prec)
  err=param$z-param$alpha[dat$village]
  pmedia=t(xmat)%*%err
  t(rmvnorm(1,var1%*%pmedia,var1))
}
#--------------------------
update.sig2=function(param){
  a=(nvil-1)/2
  err=param$alpha-wmat%*%param$gamma
  b=t(err)%*%param$invSigma%*%err/2  
  1/rgamma(1,a,b)
}
#--------------------------
update.z=function(param){
  media=param$alpha[dat$village]+xmat%*%param$betas

  z=rep(NA,nrow(dat))
  cond=!is.na(dat$Mal) & dat$Mal==1
  z[cond]=tnorm(sum(cond),lo=0,hi=Inf ,mu=media[cond],sig=1)
  cond=!is.na(dat$Mal) & dat$Mal==0
  z[cond]=tnorm(sum(cond),lo=-Inf,hi=0,mu=media[cond],sig=1)
  cond= is.na(dat$Mal)
  z[cond]=rnorm(sum(cond),mean=media[cond],sd=1)
  z
}
#-----------------------------------------
update.y=function(param){
  y=ifelse(param$z > 0, 1,0)
  y
}
#-----------------------------------------
update.misc=function(param){ 
  misc=as.data.frame(cbind(mal=dat$Mal,dat$Mal,param$y)[which(is.na(dat$Mal)),])
  mean(misc$mal!=misc$V3)
} 