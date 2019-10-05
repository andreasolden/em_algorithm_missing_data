# LOG-LIKELIHOOD
logl=function(th,x,y)
{
  x=as.matrix(x)
  px=ncol(x)
  a=th[1:2]
  b=th[3:(2+px)]
  p=matrix(1,n,3)
  #for(k in 1:2) p[,k+1]=exp(a[k]+x%*%b) #Original
  p[,2]=exp(a[1]+x%*%b) #Uglier, but 127 times faster
  p[,3]=exp(a[2]+x%*%b)
  p=p/apply(p,1,sum)
  ll=sum(y*log(p))
  mll=-ll
  return(mll)
}


# EM-ALGORITHM TO ESTIMATE PARAMETERS
em_algo=function(y1,x,niter=100,conv=1e-6, onlyparam=FALSE) #y1=caught, x=matrix(10*rnorm(px*n),n,px)
{
  x=as.matrix(x) #ensure matrix
  n=nrow(x)
  them=c(a,b)
  themold=0*them # Vector (0,0,0,0,0)
  iter=1 #number of iterations
  while(sum((them-themold)^2/(themold+0.0001)^2)>conv && iter<niter) #&& saves time, no eval of second expression if first breaks
  {
    themold=them #save old to get criteria
    em_a=them[1:2] #TWO first numbers of them
    em_b=them[-(1:2)] # NOT 1 and 2 of them
    ystar=matrix(1,n,2) #matrix full of ones
    ystar[,2]=exp(em_a[1]+x%*%em_b) #second column of Ystar
    ystar=ystar/apply(ystar,1,sum)
    yy=cbind(ystar,y1)
    yy[yy[,3]==1,1:2]=0
    mlogitem=optim(par = them,fn = logl,x=x,y=yy, hessian = FALSE, method = "BFGS") 
    them=mlogitem$par ###THE LOGL FORCES col 2 and 3 to have same coeff. 
    iter=iter+1
  }
  
  if(onlyparam==TRUE){
    return( list(em_a=em_a, em_b=em_b))
  } else{
    return( list(em_a=em_a,em_b=em_b, yy=yy, iter=iter) )
  }
}

# DATA GENERATING FUNCTIONS

data_generation = function(n, px, a, b, xsd, xcor=99) {
  
  if(px==1) {
    x=matrix(rnorm(px*n),n,px)
    x[,1]=xsd[1]*x[,1]
    
    p=matrix(1,n,3)
    p[,2]=exp(a[1]+x%*%b) 
    p[,3]=exp(a[2]+x%*%b)
    p=p/apply(p,1,sum)
    y=t(apply(p, 1, function(x) rmultinom(1,1,x))) 
    y1=y[,3] # OBSERVED FRAUDULANT CLAIMS
    return(list(p=p,y=y,y1=y1, x=x)) #return list of elements 
  } else{
    
    covarmat = xcor * tcrossprod(xsd)
    
    x=(mvrnorm(n=n, mu=c(rep(0,px)), Sigma=covarmat))
    
    p=matrix(1,n,3)
    p[,2]=exp(a[1]+x%*%b) 
    p[,3]=exp(a[2]+x%*%b)
    p=p/apply(p,1,sum)
    y=t(apply(p, 1, function(x) rmultinom(1,1,x))) 
    y1=y[,3] # OBSERVED FRAUDULANT CLAIMS
    return(list(p=p,y=y,y1=y1, x=x)) #return list of elements 
  }
}

















