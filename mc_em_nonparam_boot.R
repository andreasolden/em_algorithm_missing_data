# Parallel EM-algorithm with bootstrap

#Load log-likelihood and EM-algorithm
source("functions.r")

#Load relevant packages
library(foreach); library(doParallel); library(doRNG); library(MASS); library(here)

# SIMULATE DATA FROM 3NOMIAL LOGISTIC REGRESSION
start_time = Sys.time()
nrepl=4
bootrepl=10
n=300
px=2
a=c(-1.8,-1.5)
b=c(-0.02,0.2)
thtrue=c(a,b)
xsd=c(12.3, 1.8)
xcor=matrix(c(1,0.0,0.0,1), nrow = px, ncol = px)

cl <- makeCluster(detectCores())
registerDoParallel(cl)
registerDoRNG(seed = 3108)

monteboot <- 
  foreach(irepl=1:nrepl, .combine=rbind, .packages=("MASS")) %dopar% {
    
    # Generate data
    data=data_generation(n=n, px=px, a=a, b=b, xsd=xsd, xcor=xcor)
    
    # EM-ALG
    emres=em_algo(data$y1,data$x,conv=1e-6,niter=100)

    # Bootstrapped standard errors
    df_boot = cbind(data$y1, data$x)
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
    mlogitpi = optim(thtrue,logl,x=data$x,y=data$y, method = "BFGS", hessian = TRUE)
    mlogitpises = sqrt(diag(solve(mlogitpi$hessian)))
    #Binomial logit
    mlogitbi=summary(glm(data$y1 ~ data$x,family=binomial()))
    
    #OUT
    c(EMa1=emres$em_a[1], EMa2=emres$em_a[2], EMb1=emres$em_b[1], EMb2=emres$em_b[2], #EM coeff
      EMa1SE=ses[1], EMa2SE=ses[2], EMb1SE=ses[3], EMb2SE=ses[4], #EM ses
      PIa1=mlogitpi$par[1], PIa2=mlogitpi$par[2], PIb1=mlogitpi$par[3], PIb2=mlogitpi$par[4], #PI coeff
      PIa1SE=mlogitpises[1], PIa2SE=mlogitpises[2], PIb1SE=mlogitpises[3], PIb2SE=mlogitpises[4], #PI ses
      BIa=mlogitbi$coefficients[1,1], BIb1=mlogitbi$coefficients[2,1], BIb2=mlogitbi$coefficients[3,1], #BI coeff
      BIaSE=mlogitbi$coefficients[1,2], BIb1SE=mlogitbi$coefficients[2,2], BIb2SE=mlogitbi$coefficients[3,2], #BI ses
      iter=emres$iter, Y1=sum(data$y[,1]), Y2=sum(data$y[,2]), Y3=sum(data$y[,3]) #Other info
      )

  }

stopCluster(cl)

end_time = Sys.time()
tot_time = start_time - end_time

#Save time 
saveRDS(monteboot, file = here("simulation_results", "monteboot.rds"))
saveRDS(tot_time, file = here("simulation_results", "tot_time_monteboot.rds"))

tot_time
apply(monteboot, 2, function(x) sum(is.na(x)))
colMeans(monteboot, na.rm = TRUE)



