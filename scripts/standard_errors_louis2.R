
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

logf=function(th,X,y)
{
  px=length(X)
  a=th[1:2]
  b=th[3:(2+px)]
  p=matrix(1,3,1)
  for(k in 1:2) p[k+1]=exp(a[k]+X%*%b)
  p=p/sum(p)
  lf=sum(y*log(p))
  mlf=-lf
  return(mlf)
}

# EM-ALGORITHM TO ESTIMATE PARAMETERS
em_alg2=function(y1,x,niter=500,conv=1e-10,verbose=FALSE) #y1=caught, x=matrix(10*rnorm(px*n),n,px)
{
  x=as.matrix(x) #re-lagrer x som matries
  n=nrow(x)
  #logit=glm(y1 ~ x,family=binomial())
  #binomial_coeff = logit$coefficients
  #bin_se = sqrt(diag(vcov(logit)))
  #them= c(logit$coefficients[[1]], binomial_coeff)# theta em, vector (1,1,1,0,0) starting values
  them=c(-1.8,-1.5,-0.02,-0.2)
  themold=0*them # Vector (0,0,0,0,0)
  iter=1 #number of iterations
  while(sum((them-themold)^2/(themold+0.0001)^2)>conv & iter<niter) #Double convergence criteria
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
  require(numDeriv)
  SS=matrix(0,length(them),length(them))
  B=matrix(0,length(them),length(them))
  for(i in 1:n) 
    {
      g=grad(logf,them,X=x[i,],y=yy[i,])
      SS=SS+g%*%t(g)
      h=hessian(logf,them,X=x[i,],y=yy[i,])
      B=B+h
  }
  se=sqrt(diag(solve(B-SS)))
  return(list(em_a=em_a,em_b=em_b,hess=mlogitem$hess,SS=SS,B=B,se=se))
}


# SIMULATE DATA FROM 3NOMIAL LOGISTIC REGRESSION
# P(Z=k)=exp(a[k]+b[k]*x)/(1+sum(exp(a[k]+b[k]*x),k=1,2))
# P(Z=3)=1/(1+sum(exp(a[3]+b[3]*x)))

nrepl=500
pars=matrix(0,nrepl,4)
ses=pars

for(irepl in 1:nrepl)
{

n=1000
px=2
#px=1
x=matrix(rnorm(px*n),n,px)
x[,1]=12.3*x[,1]
x[,2]=1.8*x[,2]
#x=12.3*rnorm(n)
a=c(-1.8,-1.5)
b=c(-0.02,0.2)
#b=as.matrix(-0.02)

thtrue=c(a,b)
p=matrix(1,n,3)
for(k in 1:2) p[,k+1]=exp(a[k]+x%*%b)
p=p/apply(p,1,sum)
y=matrix(0,n,3)
for(i in 1:n) y[i,]=rmultinom(1,1,prob=p[i,])
y1=y[,3] # OBSERVED FRAUDULANT CLAIMS

# ESTIMATE PARAMETERS ONLY OBSERVING y1=1 IF FRAUDULANT, 0 OTHERWISE
mod1=em_alg2(y1,x,conv=1e-6,niter=200)
thpr=c(mod1$em_a,mod1$em_b)
pars[irepl,]=thpr
ses[irepl,]=mod1$se

}

pars
ses

apply(pars, 2, function(x) sum(is.na(x)))
apply(ses, 2, function(x) sum(is.na(x)))


thtrue
colMeans(pars, na.rm = TRUE)
colMeans(ses, na.rm = TRUE)



















