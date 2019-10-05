# Parallel EM-algorithm with bootstrap

#Load log-likelihood and EM-algorithm 
source("functions.r")

#Load relevant packages ####
library(foreach); library(doParallel); library(doRNG); library(MASS); library(here)

# SIMULATE- CONFIGURATIONS ####
# Note
# If px=1: need: nrepl, bootrepl, n, px, a, b, c(a,b), xsd

# ONE VARIABLE

# Simulation 1 ####
  # One variable, alfa(-1.8,-1.5), beta(-0.02), xsd(12.3) ####
  # START- SIMULATE
  start_time = Sys.time()
  nrepl=1000
  bootrepl=200
  n=1000
  px=1
  a=c(-1.8,-1.5)
  b=c(-0.02)
  xsd=c(12.3)

  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  registerDoRNG(seed = 3108)
  
  monteboot <- 
    foreach(irepl=1:nrepl, .combine=rbind, .packages=("MASS")) %dopar% {
      
      # Generate data
      data=data_generation(n=n, px=px, a=a, b=b, xsd=xsd)
      
      # EM-ALG
      emres=em_algo(data$y1,data$x,conv=1e-6,niter=100)
      
      # Bootstrapped standard errors
      df_boot = cbind(data$y1, data$x)
      prop_sigma=matrix(0,bootrepl,(px+2))
      
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
      mlogitpi = optim(c(a,b),logl,x=data$x,y=data$y, method = "BFGS", hessian = TRUE)
      mlogitpises = sqrt(diag(solve(mlogitpi$hessian)))
      #Binomial logit
      mlogitbi=summary(glm(data$y1 ~ data$x,family=binomial()))
      
      #OUT
      c(EMa1=emres$em_a[1], EMa2=emres$em_a[2], EMb=emres$em_b[1],  #EM coeff
        EMa1SE=ses[1], EMa2SE=ses[2], EMbSE=ses[3], #EM ses
        PIa1=mlogitpi$par[1], PIa2=mlogitpi$par[2], PIb=mlogitpi$par[3], #PI coeff
        PIa1SE=mlogitpises[1], PIa2SE=mlogitpises[2], PIbSE=mlogitpises[3], #PI ses
        BIa=mlogitbi$coefficients[1,1], BIb=mlogitbi$coefficients[2,1], #BI coeff
        BIaSE=mlogitbi$coefficients[1,2], BIbSE=mlogitbi$coefficients[2,2], #BI ses
        iter=emres$iter, Y1=sum(data$y[,1]), Y2=sum(data$y[,2]), Y3=sum(data$y[,3]) #Other info
      )
    }
  
  stopCluster(cl)
  
  end_time = Sys.time()
  tot_time = start_time - end_time
  # Save ####
  saveRDS(monteboot, file = here("simulation_results", "one_var_base.rds"))
  saveRDS(tot_time, file = here("simulation_results", "one_var_base_time.rds"))

# Simulation 2 ####
  # One variable, alfa(-3,-1.5), beta(-0.02), xsd(12.3) ####
  # This is 200 % instead of 20% alpha-deviation
  # START- SIMULATE
  start_time = Sys.time()
  nrepl=1000
  bootrepl=200
  n=1000
  px=1
  a=c(-3,-1.5)
  b=c(-0.02)
  xsd=c(12.3)

  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  registerDoRNG(seed = 3108)
  
  monteboot <- 
    foreach(irepl=1:nrepl, .combine=rbind, .packages=("MASS")) %dopar% {
      
      # Generate data
      data=data_generation(n=n, px=px, a=a, b=b, xsd=xsd)
      
      # EM-ALG
      emres=em_algo(data$y1,data$x,conv=1e-6,niter=100)
      
      # Bootstrapped standard errors
      df_boot = cbind(data$y1, data$x)
      prop_sigma=matrix(0,bootrepl,(px+2))
      
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
      mlogitpi = optim(c(a,b),logl,x=data$x,y=data$y, method = "BFGS", hessian = TRUE)
      mlogitpises = sqrt(diag(solve(mlogitpi$hessian)))
      #Binomial logit
      mlogitbi=summary(glm(data$y1 ~ data$x,family=binomial()))
      
      #OUT
      c(EMa1=emres$em_a[1], EMa2=emres$em_a[2], EMb=emres$em_b[1],  #EM coeff
        EMa1SE=ses[1], EMa2SE=ses[2], EMbSE=ses[3], #EM ses
        PIa1=mlogitpi$par[1], PIa2=mlogitpi$par[2], PIb=mlogitpi$par[3], #PI coeff
        PIa1SE=mlogitpises[1], PIa2SE=mlogitpises[2], PIbSE=mlogitpises[3], #PI ses
        BIa=mlogitbi$coefficients[1,1], BIb=mlogitbi$coefficients[2,1], #BI coeff
        BIaSE=mlogitbi$coefficients[1,2], BIbSE=mlogitbi$coefficients[2,2], #BI ses
        iter=emres$iter, Y1=sum(data$y[,1]), Y2=sum(data$y[,2]), Y3=sum(data$y[,3]) #Other info
      )
    }
  
  stopCluster(cl)
  
  end_time = Sys.time()
  tot_time = start_time - end_time
  #Save ####
  saveRDS(monteboot, file = here("simulation_results", "one_var_a_up.rds"))
  saveRDS(tot_time, file = here("simulation_results", "one_var_a_up_time.rds"))

# Simulation 3 ####
  # One variable, alfa(-1.8,-1.5), beta(-0.02), xsd(123) ####
  # Ten times the SD
  # START- SIMULATE
  start_time = Sys.time()
  nrepl=1000
  bootrepl=200
  n=1000
  px=1
  a=c(-1.8,-1.5)
  b=c(-0.02)
  xsd=c(123)

  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  registerDoRNG(seed = 3108)
  
  monteboot <- 
    foreach(irepl=1:nrepl, .combine=rbind, .packages=("MASS")) %dopar% {
      
      # Generate data
      data=data_generation(n=n, px=px, a=a, b=b, xsd=xsd)
      
      # EM-ALG
      emres=em_algo(data$y1,data$x,conv=1e-6,niter=100)
      
      # Bootstrapped standard errors
      df_boot = cbind(data$y1, data$x)
      prop_sigma=matrix(0,bootrepl,(px+2))
      
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
      mlogitpi = optim(c(a,b),logl,x=data$x,y=data$y, method = "BFGS", hessian = TRUE)
      mlogitpises = sqrt(diag(solve(mlogitpi$hessian)))
      #Binomial logit
      mlogitbi=summary(glm(data$y1 ~ data$x,family=binomial()))
      
      #OUT
      c(EMa1=emres$em_a[1], EMa2=emres$em_a[2], EMb=emres$em_b[1],  #EM coeff
        EMa1SE=ses[1], EMa2SE=ses[2], EMbSE=ses[3], #EM ses
        PIa1=mlogitpi$par[1], PIa2=mlogitpi$par[2], PIb=mlogitpi$par[3], #PI coeff
        PIa1SE=mlogitpises[1], PIa2SE=mlogitpises[2], PIbSE=mlogitpises[3], #PI ses
        BIa=mlogitbi$coefficients[1,1], BIb=mlogitbi$coefficients[2,1], #BI coeff
        BIaSE=mlogitbi$coefficients[1,2], BIbSE=mlogitbi$coefficients[2,2], #BI ses
        iter=emres$iter, Y1=sum(data$y[,1]), Y2=sum(data$y[,2]), Y3=sum(data$y[,3]) #Other info
      )
    }
  
  stopCluster(cl)
  
  end_time = Sys.time()
  tot_time = start_time - end_time
  # Save ####
  saveRDS(monteboot, file = here("simulation_results", "one_var_xsd_up.rds"))
  saveRDS(tot_time, file = here("simulation_results", "one_var_xsd_up_time.rds"))

# TWO VARIABLES 

# Simulation 4 ####
  # Two variables, alfa(-1.8,-1.5), beta(-0.02, 0.2), xsd(12.3, 1.8), corr(0) ####
  # New base
  # START- SIMULATE
  start_time = Sys.time()
  nrepl=1000
  bootrepl=200
  n=1000
  px=2
  a=c(-1.8,-1.5)
  b=c(-0.02,0.2)
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
      prop_sigma=matrix(0,bootrepl,(px+2))
      
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
      mlogitpi = optim(c(a,b),logl,x=data$x,y=data$y, method = "BFGS", hessian = TRUE)
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
  # Save ####
  saveRDS(monteboot, file = here("simulation_results", "two_var_base.rds"))
  saveRDS(tot_time, file = here("simulation_results", "two_var_base_time.rds"))

# Simulation 5 ####
  # Two variables, alfa(-1.8,-1.5), beta(-0.02, 0.2), xsd(12.3, 1.8), corr(0.5) ####
  # Introduce correlation
  # START- SIMULATE
  start_time = Sys.time()
  nrepl=1000
  bootrepl=200
  n=1000
  px=2
  a=c(-1.8,-1.5)
  b=c(-0.02,0.2)
  xsd=c(12.3, 1.8)
  xcor=matrix(c(1,0.5,0.5,1), nrow = px, ncol = px)
  
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
      prop_sigma=matrix(0,bootrepl,(px+2))
      
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
      mlogitpi = optim(c(a,b),logl,x=data$x,y=data$y, method = "BFGS", hessian = TRUE)
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
  # Save ####
  saveRDS(monteboot, file = here("simulation_results", "two_var_05corr.rds"))
  saveRDS(tot_time, file = here("simulation_results", "two_var_05corr_time.rds"))

# Simulation 6 ####
  # Two variables, alfa(-1.8,-1.5), beta(-0.02, 0.2), xsd(12.3, 1.8), corr(0.9) ####
  # Introduce correlation
  # START- SIMULATE
  start_time = Sys.time()
  nrepl=1000
  bootrepl=200
  n=1000
  px=2
  a=c(-1.8,-1.5)
  b=c(-0.02,0.2)
  xsd=c(12.3, 1.8)
  xcor=matrix(c(1,0.9,0.9,1), nrow = px, ncol = px)
  
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
      prop_sigma=matrix(0,bootrepl,(px+2))
      
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
      mlogitpi = optim(c(a,b),logl,x=data$x,y=data$y, method = "BFGS", hessian = TRUE)
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
  # Save ####
  saveRDS(monteboot, file = here("simulation_results", "two_var_09corr.rds"))
  saveRDS(tot_time, file = here("simulation_results", "two_var_09corr_time.rds"))

  # Simulation 7 ####
  # Two variables, alfa(-1.8,-1.5), beta(-0.02, 0.2), xsd(123, 1.8), corr(0.0) ####
  # Increase sd for x1
  # START- SIMULATE
  start_time = Sys.time()
  nrepl=1000
  bootrepl=200
  n=1000
  px=2
  a=c(-1.8,-1.5)
  b=c(-0.02,0.2)
  xsd=c(123, 1.8)
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
      prop_sigma=matrix(0,bootrepl,(px+2))
      
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
      mlogitpi = optim(c(a,b),logl,x=data$x,y=data$y, method = "BFGS", hessian = TRUE)
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
  # Save ####
  saveRDS(monteboot, file = here("simulation_results", "two_var_sdx1_up.rds"))
  saveRDS(tot_time, file = here("simulation_results", "two_var_sdx1_up_time.rds"))
# Simulation 8 ####
  # Two variables, alfa(-1.8,-1.5), beta(-0.02, 0.2), xsd(123, 18), corr(0.0) ####
  # Increase (ten times) sd for x1 and x2
  # START- SIMULATE
  start_time = Sys.time()
  nrepl=1000
  bootrepl=200
  n=1000
  px=2
  a=c(-1.8,-1.5)
  b=c(-0.02,0.2)
  xsd=c(123, 18)
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
      prop_sigma=matrix(0,bootrepl,(px+2))
      
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
      mlogitpi = optim(c(a,b),logl,x=data$x,y=data$y, method = "BFGS", hessian = TRUE)
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
  # Save ####
  saveRDS(monteboot, file = here("simulation_results", "two_var_sdx_up.rds"))
  saveRDS(tot_time, file = here("simulation_results", "two_var_sdx_up_time.rds"))

  # Simulation 9 ####
  # Two variables, alfa(-1.8,-1.5), beta(-0.02, 0.2), xsd(123, 18), corr(0.9) ####
  # Increase (ten times) sd for x1 and high correlation
  # START- SIMULATE
  start_time = Sys.time()
  nrepl=1000
  bootrepl=200
  n=1000
  px=2
  a=c(-1.8,-1.5)
  b=c(-0.02,0.2)
  xsd=c(123, 18)
  xcor=matrix(c(1,0.9,0.9,1), nrow = px, ncol = px)
  
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
      prop_sigma=matrix(0,bootrepl,(px+2))
      
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
      mlogitpi = optim(c(a,b),logl,x=data$x,y=data$y, method = "BFGS", hessian = TRUE)
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
  # Save ####
  saveRDS(monteboot, file = here("simulation_results", "two_var_sdx_and_cor_up.rds"))
  saveRDS(tot_time, file = here("simulation_results", "two_var_sdx_and_cor_up_time.rds"))

  # Simulation 10 ####
  # Two variables, alfa(-1.8,-1.5), beta(-0.02, 0.2), xsd(12.3, 1.8), corr(0) change sseed ####
  # START- SIMULATE
  start_time = Sys.time()
  nrepl=1000
  bootrepl=200
  n=1000
  px=2
  a=c(-1.8,-1.5)
  b=c(-0.02,0.2)
  xsd=c(12.3, 1.8)
  xcor=matrix(c(1,0.0,0.0,1), nrow = px, ncol = px)
  
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  registerDoRNG(seed = 1703)
  
  monteboot <- 
    foreach(irepl=1:nrepl, .combine=rbind, .packages=("MASS")) %dopar% {
      
      # Generate data
      data=data_generation(n=n, px=px, a=a, b=b, xsd=xsd, xcor=xcor)
      
      # EM-ALG
      emres=em_algo(data$y1,data$x,conv=1e-6,niter=100)
      
      # Bootstrapped standard errors
      df_boot = cbind(data$y1, data$x)
      prop_sigma=matrix(0,bootrepl,(px+2))
      
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
      mlogitpi = optim(c(a,b),logl,x=data$x,y=data$y, method = "BFGS", hessian = TRUE)
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
  
  # Save ####
  saveRDS(monteboot, file = here("simulation_results", "two_var_base_new_seed.rds"))
  saveRDS(tot_time, file = here("simulation_results", "two_var_base_time_new_seed.rds"))

# Simulation 11 #### 
  # Two variables, alfa(-1.8,-1.5), beta(-0.02, 0.2), xsd(12.3, 1.8), corr(0) ####
  # 1000 bootrepl
  # START- SIMULATE
  start_time = Sys.time()
  nrepl=1000
  bootrepl=1000
  n=1000
  px=2
  a=c(-1.8,-1.5)
  b=c(-0.02,0.2)
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
      prop_sigma=matrix(0,bootrepl,(px+2))
      
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
      mlogitpi = optim(c(a,b),logl,x=data$x,y=data$y, method = "BFGS", hessian = TRUE)
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
  # Save ####
  saveRDS(monteboot, file = here("simulation_results", "two_var_1000_bootrepl.rds"))
  saveRDS(tot_time, file = here("simulation_results", "two_var_1000_bootrepl_time.rds"))
