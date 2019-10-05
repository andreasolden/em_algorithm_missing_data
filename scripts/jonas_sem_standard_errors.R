
# LOG-LIKELIHOOD
logl=function(th,x,y)
{
  x=as.matrix(x)
  px=ncol(x)
  a=th[1:2]
  b=th[3:(2+px)]
  p=matrix(1,n,3)
  for(k in 1:2) p[,k+1]=exp(a[k]+x%*%b)
  p=p/apply(p,1,sum)
  ll=sum(y*log(p))
  mll=-ll
  return(mll)
}

# EM-ALGORITHM TO ESTIMATE PARAMETERS
em_alg2=function(y1,x,niter=100,conv=1e-5,verbose=TRUE) #y1=caught, x=matrix(10*rnorm(px*n),n,px)
{
  x=as.matrix(x) #re-lagrer x som matries
  logit=glm(y1 ~ x,family=binomial())
  binomial_coeff = logit$coefficients
  bin_se = sqrt(diag(vcov(logit)))
  them= c(logit$coefficients[[1]], binomial_coeff)# theta em, vector (1,1,1,0,0) starting values
  themold=0*them # Vector (0,0,0,0,0)
  iter=1 #number of iterations
  while(sum((them-themold)^2)>1e-5 & iter<niter) #Double convergence criteria
  {
    if(iter==1 & verbose) print(c("iteration","a2","a3",paste("b",1:px,sep="")))
    themold=them #save old to get criteria
    em_a=them[1:2] #TWO first numbers of them
    em_b=them[-(1:2)] # NOT 1 and 2 of them
    ystar=matrix(1,n,2) #matrix full of ones
    ystar[,2]=exp(em_a[1]+x%*%em_b) #second column of Ystar
    ystar=ystar/apply(ystar,1,sum)
    yy=cbind(ystar,y1)
    yy[yy[,3]==1,1:2]=0
    mlogitem=optim(par = them,fn = logl,x=x,y=yy,hessian=TRUE, method = "BFGS")
    them=mlogitem$par ###THE LOGL FORCES col 2 and 3 to have same coeff. 
    if(verbose) print(c(iter,them))
    iter=iter+1
  }
  return(list(em_a=em_a,em_b=em_b,hess=mlogitem$hess))
}

# ONE ITERATION OF IN THE EM-ALGORITHM
em_one=function(th0,y1,x)
{
  a=th0[1:2]
  b=th0[-(1:2)]
  ystar=matrix(1,n,2)
  ystar[,2]=exp(a[1]+x%*%b)
  ystar=ystar/apply(ystar,1,sum)
  yy=cbind(ystar,y1)
  yy[yy[,3]==1,1:2]=0
  mlogitem=optim(par = th0,fn = logl,x=x,y=yy,hessian=TRUE, method = "BFGS")
  th1=mlogitem$par
  return(th1)
}

# COMPUTE THE DM (CORRECTION TERM TO COMPUTE ASYMTOTIC VARIANCE-COVARIANCE MATRIX)
dm=function(thopt,th0,y1,x,niter=10,verbose=FALSE)
{
  d=length(thopt)
  r=matrix(0,d,d)
  for(iter in 1:niter)
  {
    th0i=thopt
    th1=em_one(th0,y1,x)
    rold=r
    for(i in 1:d)
    {
      th0i[i]=th0[i]
      th1i=em_one(th0i,y1,x)
      for(j in 1:d)
      {
        r[i,j]=(th1i[j]-thopt[j])/(th0i[i]-thopt[i])
      }
      th0=th1
    }
    if(verbose) print(r-rold)
  }
  return(r)
}
  
sem=function(thopt,th0,y1,x,niter=10,verbose=FALSE)
{
  thopt=c(mod1$em_a,mod1$em_b)
  DM=dm(thopt,th0,y1,x,niter=niter,verbose=verbose)
  Ioc=mod1$hess
  V=-solve(Ioc%*%solve(diag(4)-DM))
  se=sqrt(diag(V))
  return(se)
}

# SIMULATE DATA FROM 3NOMIAL LOGISTIC REGRESSION
# P(Z=k)=exp(a[k]+b[k]*x)/(1+sum(exp(a[k]+b[k]*x),k=1,2))
# P(Z=3)=1/(1+sum(exp(a[3]+b[3]*x)))
set.seed(123)
n=1000
px=2
x=matrix(0, nrow = n, ncol = px)
x[,1] = rnorm(n, mean = 0, sd = 12.3)
x[,2] = rnorm(n, mean = 0, sd = 1.8)

a=c(-1.8,-1.5)
b=c(-0.02, 0.2)
thtrue=c(a,b)
p=matrix(1,n,3)
for(k in 1:2) p[,k+1]=exp(a[k]+x%*%b)
p=p/apply(p,1,sum)
y=matrix(0,n,3)
for(i in 1:n) y[i,]=rmultinom(1,1,prob=p[i,])
y1=y[,3] # OBSERVED FRAUDULANT CLAIMS

# ESTIMATE PARAMETERS ONLY OBSERVING y1=1 IF FRAUDULANT, 0 OTHERWISE
mod1=em_alg2(y1,x)

# COMPUTE STANDARD ERRORS
logit=glm(y1 ~ x,family=binomial())
th0 = logit$coefficients
th0=c(rep(th0[1],ncol(x)),th0[-1])
#th0=c(th0[1],0,th0[-1])
thopt=c(mod1$em_a,mod1$em_b)
se=sem(thopt,th0,y1,x,verbose=TRUE,niter=10)
thopt
se
