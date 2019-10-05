# Parallel EM-algorithm with bootstrap

#Load log-likelihood and EM-algorithm
source("functions.r")

#Load relevant packages
library(foreach); library(doParallel); library(doRNG)

# SIMULATE DATA FROM 3NOMIAL LOGISTIC REGRESSION
start_time = Sys.time()
nrepl=1000
bootrepl=1000
n=1000
px=2
a=c(-1.8,-1.5)
b=c(-0.02,0.2)
thtrue=c(a,b)

cl <- makeCluster(detectCores())
registerDoParallel(cl)
registerDoRNG(seed = 3108)

monteboot <- 
  foreach(irepl=1:nrepl, .combine=rbind) %dopar% {
    
    x=matrix(rnorm(px*n),n,px)
    x[,1]=12.3*x[,1]
    x[,2]=1.8*x[,2]
    
    p=matrix(1,n,3)
    p[,2]=exp(a[1]+x%*%b) 
    p[,3]=exp(a[2]+x%*%b)
    p=p/apply(p,1,sum)
    y=t(apply(p, 1, function(x) rmultinom(1,1,x))) 
    y1=y[,3] # OBSERVED FRAUDULANT CLAIMS
    
    # EM-ALG
    emres=em_algo(y1,x,conv=1e-6,niter=100)
    thpr=c(emres$em_a,emres$em_b) #Parameters
    
    # Bootstrapped standard errors
    df_boot = cbind(y1, x)
    prop_sigma=matrix(0,bootrepl,4)
    
    for(ibot in 1:bootrepl)
    {
      xboot = df_boot[sample(nrow(df_boot), size = n, replace = TRUE),]
      modboot = em_algo(xboot[,1], xboot[,2:(px+1)], conv=1e-6, niter=100, onlyparam = TRUE)
      prop_sigma[ibot,] = c(modboot$em_a,modboot$em_b)
    }
    
    #Bootstrapped standard errors
    ses=as.vector(apply(prop_sigma, 2, sd))
    rm(prop_sigma)
    
    #Perfect information trinomial logit
    mlogitpi = optim(thtrue,logl,x=x,y=y, method = "BFGS", hessian = TRUE)
    mlogitpises = sqrt(diag(solve(mlogitpi$hessian)))
    #Binomial logit
    mlogitbi=summary(glm(y1 ~ x,family=binomial()))
    
    #OUT
    c(thpr, ses, emres$iter, sum(y[,1]), sum(y[,2]), sum(y[,3]), mlogitpi$par, mlogitpises, mlogitbi$coefficients[,1], mlogitbi$coefficients[,2])
  }

stopCluster(cl)

end_time = Sys.time()

tot_time = start_time - end_time

colnames(monteboot) <- c("EMa1", "EMa2", "EMb1", "EMb2", "EMa1SE", "EMa2SE", "EMb1SE", "EMb2SE", 
                         "iter", "Y1", "Y2", "Y3",  
                         "PIa1", "PIa2", "PIb1", "PIb2", "PIa1SE", "PIa2SE", "PIb1SE", "PIb2SE",
                         "BIa1", "BIb1", "BIb2", "BIa1SE", "BIb1SE", "BIb2SE"
)



saveRDS(monteboot, file = "monteboot.rds")
saveRDS(tot_time, file = "tot_time_monteboot.rds")

#saveRDS(monteboot, file = "monteboot03092019.rds")
#saveRDS(tot_time, file = "tot_time_monteboot03092019.rds")

