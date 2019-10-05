#This script runs the Monte Carlo simulations and saves the results in folder results

#       Make sure you have installed checkpoint 
#       Make sure you opened the Rproj file first or you might not get results
#       Note that doSNOW might need tweaking and custom cores. Ask.

#Load necessary packages and funtions ####
library(checkpoint) # Run initiate_project before this. Makes sure we use same package versions
checkpoint("2019-01-21") #All packages will load as they were on CRAN this day
library(here) #Relative file paths. here("folderA", "file.extension") or  here("folderA/file.extension")
library(foreach) #For running simulation
library(doSNOW) #For parallelization
library(MASS) # For mvrnorm simulating from a multivariate distribution

#Source script where we have stored all functions
source(here("scripts", "define_functions.R"))


#Begin simulations ####


#Simulations 1 ####
#One variable, alfa(-1.8,-1.5), beta(-0.02), xval(0), xsd(12.3)

#Define parameters ####
set.seed(123)
n=1000
nrepl=1000
niter=100
#conv=1e-5
px=1
alfas=c(-1.8,-1.5)
betas=c(-0.02)
xmeans=c(0)
xsd=c(12.3)
xcor <- matrix(c(1), nrow = px, ncol = px)
test_simparam(n, px, alfas, betas, xmeans, xsd, xcor) #Simple check on lengths

# actual simulation ####
cl<-makeCluster(4); registerDoSNOW(cl)
system.time({
        dfs=as.data.frame( 
                foreach(irepl=1:nrepl, .combine=rbind, .packages=("MASS")) %dopar%      
                {
                        #Generate data
                        data=datagen(n, alfas, betas, xmeans, xsd, xcor)
                        
                        #EM algorithm. Note binomial model embedded in EM
                        em_res = em_alg(data$y1,data$x)
                        
                        #Perfect info trinomial logit with start true*2
                        mlogit=optim(c(alfas, betas)*2,logl,x=data$x,y=data$y) 
                        
                        #Create list of output
                        c(n, px, em_res$iter, nrepl, sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]),
                          alfas, betas, em_res$em_a, em_res$em_b, em_res$binomial_coeff, mlogit$par)
                } )
})
stopCluster(cl) 
colnames(dfs) <- givecolumnnames()

#Save results ####
saveRDS(dfs, file = here("results", "simres1.rds"))

#End simulation 1 ####


#Simulations 2 ####
#One variable, changed alfa1: alfa(-3,-1.5). This is 200% diff instead of 20%

#Define parameters ####
set.seed(123)
n=1000
nrepl=1000
niter=100
#conv=1e-5
px=1
alfas=c(-3,-1.5)
betas=c(-0.02)
xmeans=c(0)
xsd=c(12.3)
xcor <- matrix(c(1), nrow = px, ncol = px)
test_simparam(n, px, alfas, betas, xmeans, xsd, xcor) #Simple check on lengths

# actual simulation ####
cl<-makeCluster(4); registerDoSNOW(cl)
system.time({
        dfs=as.data.frame( 
                foreach(irepl=1:nrepl, .combine=rbind, .packages=("MASS")) %dopar%      
                {
                        #Generate data
                        data=datagen(n, alfas, betas, xmeans, xsd, xcor)
                        
                        #EM algorithm. Note binomial model embedded in EM
                        em_res = em_alg(data$y1,data$x)
                        
                        #Perfect info trinomial logit with start true*2
                        mlogit=optim(c(alfas, betas)*2,logl,x=data$x,y=data$y) 
                        
                        #Create list of output
                        c(n, px, em_res$iter, nrepl, sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]),
                          alfas, betas, em_res$em_a, em_res$em_b, em_res$binomial_coeff, mlogit$par)
                } )
})
stopCluster(cl) 
colnames(dfs) <- givecolumnnames()

#Save results ####
saveRDS(dfs, file = here("results", "simres2.rds"))

#End simulation 2 ####


#Simulations 3 ####
#One variable, changed standard dev of x-var (123 instead of 12.3)

#Define parameters ####
set.seed(123)
n=1000
nrepl=1000
niter=100
#conv=1e-5
px=1
alfas=c(-1.8,-1.5)
betas=c(-0.02)
xmeans=c(0)
xsd=c(123)
xcor <- matrix(c(1), nrow = px, ncol = px)
test_simparam(n, px, alfas, betas, xmeans, xsd, xcor) #Simple check on lengths

# actual simulation ####
cl<-makeCluster(4); registerDoSNOW(cl)
system.time({
        dfs=as.data.frame( 
                foreach(irepl=1:nrepl, .combine=rbind, .packages=("MASS")) %dopar%      
                {
                        #Generate data
                        data=datagen(n, alfas, betas, xmeans, xsd, xcor)
                        
                        #EM algorithm. Note binomial model embedded in EM
                        em_res = em_alg(data$y1,data$x)
                        
                        #Perfect info trinomial logit with start true*2
                        mlogit=optim(c(alfas, betas)*2,logl,x=data$x,y=data$y) 
                        
                        #Create list of output
                        c(n, px, em_res$iter, nrepl, sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]),
                          alfas, betas, em_res$em_a, em_res$em_b, em_res$binomial_coeff, mlogit$par)
                } )
})
stopCluster(cl) 
colnames(dfs) <- givecolumnnames()

#Save results ####
saveRDS(dfs, file = here("results", "simres3.rds"))

#End simulation 3 ####


#Simulations 4 ####
#Two variables, alfa(-1.8,-1.5), beta(-0.02, 0.2), xval(0), xsd(12.3, 1.8), corr(0)

#Define parameters ####
set.seed(123)
n=1000
nrepl=1000
niter=100
#conv=1e-5
px=2
alfas=c(-1.8,-1.5)
betas=c(-0.02, 0.2)
xmeans=c(0,0)
xsd=c(12.3, 1.8)
xcor <- matrix(c(1,0,0,1), nrow = px, ncol = px)
test_simparam(n, px, alfas, betas, xmeans, xsd, xcor) #Simple check on lengths

# actual simulation ####
cl<-makeCluster(4); registerDoSNOW(cl)
system.time({
        dfs=as.data.frame( 
                foreach(irepl=1:nrepl, .combine=rbind, .packages=("MASS")) %dopar%      
                {
                        #Generate data
                        data=datagen(n, alfas, betas, xmeans, xsd, xcor)
                        
                        #EM algorithm. Note binomial model embedded in EM
                        em_res = em_alg(data$y1,data$x)
                        
                        #Perfect info trinomial logit with start true*2
                        mlogit=optim(c(alfas, betas)*2,logl,x=data$x,y=data$y) 
                        
                        #Create list of output
                        c(n, px, em_res$iter, nrepl, sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]),
                          alfas, betas, em_res$em_a, em_res$em_b, em_res$binomial_coeff, mlogit$par)
                } )
})
stopCluster(cl) 
colnames(dfs) <- givecolumnnames()

#Save results ####
saveRDS(dfs, file = here("results", "simres4.rds"))

#End simulation 4 ####


#Simulations 5 ####
#Two variables, introduce 0.5 correlation between x's

#Define parameters ####
set.seed(123)
n=1000
nrepl=1000
niter=100
#conv=1e-5
px=2
alfas=c(-1.8,-1.5)
betas=c(-0.02, 0.2)
xmeans=c(0,0)
xsd=c(12.3, 1.8)
xcor <- matrix(c(1,0.5,0.5,1), nrow = px, ncol = px)
test_simparam(n, px, alfas, betas, xmeans, xsd, xcor) #Simple check on lengths

# actual simulation ####
cl<-makeCluster(4); registerDoSNOW(cl)
system.time({
        dfs=as.data.frame( 
                foreach(irepl=1:nrepl, .combine=rbind, .packages=("MASS")) %dopar%      
                {
                        #Generate data
                        data=datagen(n, alfas, betas, xmeans, xsd, xcor)
                        
                        #EM algorithm. Note binomial model embedded in EM
                        em_res = em_alg(data$y1,data$x)
                        
                        #Perfect info trinomial logit with start true*2
                        mlogit=optim(c(alfas, betas)*2,logl,x=data$x,y=data$y) 
                        
                        #Create list of output
                        c(n, px, em_res$iter, nrepl, sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]),
                          alfas, betas, em_res$em_a, em_res$em_b, em_res$binomial_coeff, mlogit$par)
                } )
})
stopCluster(cl) 
colnames(dfs) <- givecolumnnames()

#Save results ####
saveRDS(dfs, file = here("results", "simres5.rds"))

#End simulation 5 ####


#Simulations 6 ####
#Two variables, introduce 0.9 correlation between x's

#Define parameters ####
set.seed(123)
n=1000
nrepl=1000
niter=100
#conv=1e-5
px=2
alfas=c(-1.8,-1.5)
betas=c(-0.02, 0.2)
xmeans=c(0,0)
xsd=c(12.3, 1.8)
xcor <- matrix(c(1,0.9,0.9,1), nrow = px, ncol = px)
test_simparam(n, px, alfas, betas, xmeans, xsd, xcor) #Simple check on lengths

# actual simulation ####
cl<-makeCluster(4); registerDoSNOW(cl)
system.time({
        dfs=as.data.frame( 
                foreach(irepl=1:nrepl, .combine=rbind, .packages=("MASS")) %dopar%      
                {
                        #Generate data
                        data=datagen(n, alfas, betas, xmeans, xsd, xcor)
                        
                        #EM algorithm. Note binomial model embedded in EM
                        em_res = em_alg(data$y1,data$x)
                        
                        #Perfect info trinomial logit with start true*2
                        mlogit=optim(c(alfas, betas)*2,logl,x=data$x,y=data$y) 
                        
                        #Create list of output
                        c(n, px, em_res$iter, nrepl, sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]),
                          alfas, betas, em_res$em_a, em_res$em_b, em_res$binomial_coeff, mlogit$par)
                } )
})
stopCluster(cl) 
colnames(dfs) <- givecolumnnames()

#Save results ####
saveRDS(dfs, file = here("results", "simres6.rds"))

#End simulation 6 ####


#Simulations 7 ####
#Two variables, same as base but with sd*10 for X1 only

#Define parameters ####
set.seed(123)
n=1000
nrepl=1000
niter=100
#conv=1e-5
px=2
alfas=c(-1.8,-1.5)
betas=c(-0.02, 0.2)
xmeans=c(0,0)
xsd=c(123, 1.8)
xcor <- matrix(c(1,0,0,1), nrow = px, ncol = px)
test_simparam(n, px, alfas, betas, xmeans, xsd, xcor) #Simple check on lengths

# actual simulation ####
cl<-makeCluster(4); registerDoSNOW(cl)
system.time({
        dfs=as.data.frame( 
                foreach(irepl=1:nrepl, .combine=rbind, .packages=("MASS")) %dopar%      
                {
                        #Generate data
                        data=datagen(n, alfas, betas, xmeans, xsd, xcor)
                        
                        #EM algorithm. Note binomial model embedded in EM
                        em_res = em_alg(data$y1,data$x)
                        
                        #Perfect info trinomial logit with start true*2
                        mlogit=optim(c(alfas, betas)*2,logl,x=data$x,y=data$y) 
                        
                        #Create list of output
                        c(n, px, em_res$iter, nrepl, sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]),
                          alfas, betas, em_res$em_a, em_res$em_b, em_res$binomial_coeff, mlogit$par)
                } )
})
stopCluster(cl) 
colnames(dfs) <- givecolumnnames()

#Save results ####
saveRDS(dfs, file = here("results", "simres7.rds"))

#End simulation 7 ####


#Simulations 8 ####
#Two variables, same as base but with sd*10 for X1 AND X2

#Define parameters ####
set.seed(123)
n=1000
nrepl=1000
niter=100
#conv=1e-5
px=2
alfas=c(-1.8,-1.5)
betas=c(-0.02, 0.2)
xmeans=c(0,0)
xsd=c(123, 18)
xcor <- matrix(c(1,0,0,1), nrow = px, ncol = px)
test_simparam(n, px, alfas, betas, xmeans, xsd, xcor) #Simple check on lengths


# actual simulation ####
cl<-makeCluster(4); registerDoSNOW(cl)
system.time({
        dfs=as.data.frame( 
                foreach(irepl=1:nrepl, .combine=rbind, .packages=("MASS")) %dopar%      
                {
                        #Generate data
                        data=datagen(n, alfas, betas, xmeans, xsd, xcor)
                        
                        #EM algorithm. Note binomial model embedded in EM
                        em_res = em_alg(data$y1,data$x)
                        
                        #Perfect info trinomial logit with start true*2
                        mlogit=optim(c(alfas, betas)*2,logl,x=data$x,y=data$y) 
                        
                        #Create list of output
                        c(n, px, em_res$iter, nrepl, sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]),
                          alfas, betas, em_res$em_a, em_res$em_b, em_res$binomial_coeff, mlogit$par)
                } )
})
stopCluster(cl) 
colnames(dfs) <- givecolumnnames()

#Save results ####
saveRDS(dfs, file = here("results", "simres8.rds"))

#End simulation 8 ####


#Simulations 9 ####
#Two variables, same as base but with sd*10 for X1 AND X2 AND CORR(0.9)

#Define parameters ####
set.seed(123)
n=1000
nrepl=1000
niter=100
#conv=1e-5
px=2
alfas=c(-1.8,-1.5)
betas=c(-0.02, 0.2)
xmeans=c(0,0)
xsd=c(123, 18)
xcor <- matrix(c(1,0.9,0.9,1), nrow = px, ncol = px)
test_simparam(n, px, alfas, betas, xmeans, xsd, xcor) #Simple check on lengths

# actual simulation ####
cl<-makeCluster(4); registerDoSNOW(cl)
system.time({
        dfs=as.data.frame( 
                foreach(irepl=1:nrepl, .combine=rbind, .packages=("MASS")) %dopar%      
                {
                        #Generate data
                        data=datagen(n, alfas, betas, xmeans, xsd, xcor)
                        
                        #EM algorithm. Note binomial model embedded in EM
                        em_res = em_alg(data$y1,data$x)
                        
                        #Perfect info trinomial logit with start true*2
                        mlogit=optim(c(alfas, betas)*2,logl,x=data$x,y=data$y) 
                        
                        #Create list of output
                        c(n, px, em_res$iter, nrepl, sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]),
                          alfas, betas, em_res$em_a, em_res$em_b, em_res$binomial_coeff, mlogit$par)
                } )
})
stopCluster(cl) 
colnames(dfs) <- givecolumnnames()

#Save results ####
saveRDS(dfs, file = here("results", "simres9.rds"))

#End simulation 9 ####

#Simulations 10 ####
#Three variables, as age license and records

#Define parameters ####
set.seed(123)
n=1000
nrepl=1000
niter=100
#conv=1e-5
px=3
alfas=c(-1.8,-1.5)
betas=c(0.023, -0.005, -0.200)
xmeans=c(38.0, 14.2, 1.4)
xsd=c(12.3, 9.1, 1.8)
xcor <- matrix(c(1, 0.4 ,0.6,
                 0.4, 1, 0.5,
                 0.6, 0.5, 1), nrow = px, ncol = px)
test_simparam(n, px, alfas, betas, xmeans, xsd, xcor) #Simple check on lengths

# actual simulation ####
cl<-makeCluster(4); registerDoSNOW(cl)
system.time({
        dfs=as.data.frame( 
                foreach(irepl=1:nrepl, .combine=rbind, .packages=("MASS")) %dopar%      
                {
                        #Generate data
                        data=datagen(n, alfas, betas, xmeans, xsd, xcor)
                        
                        #EM algorithm. Note binomial model embedded in EM
                        em_res = em_alg(data$y1,data$x)
                        
                        #Perfect info trinomial logit with start true*2
                        mlogit=optim(c(alfas, betas)*2,logl,x=data$x,y=data$y) 
                        
                        #Create list of output
                        c(n, px, em_res$iter, nrepl, sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]),
                          alfas, betas, em_res$em_a, em_res$em_b, em_res$binomial_coeff, mlogit$par)
                } )
})
stopCluster(cl) 
colnames(dfs) <- givecolumnnames()

#Save results ####
saveRDS(dfs, file = here("results", "simres10.rds"))

#End simulation 10 ####

#Sim1 with standard error tweak ####
#One variable, alfa(-1.8,-1.5), beta(-0.02), xval(0), xsd(12.3)

#Define parameters ####
set.seed(123)
n=1000
nrepl=1000
niter=100
#conv=1e-5
px=1
alfas=c(-1.8,-1.5)
betas=c(-0.02)
xmeans=c(0)
xsd=c(12.3)
xcor <- matrix(c(1), nrow = px, ncol = px)
test_simparam(n, px, alfas, betas, xmeans, xsd, xcor) #Simple check on lengths

# actual simulation ####
cl<-makeCluster(4); registerDoSNOW(cl)
system.time({
        dfs=as.data.frame( 
                foreach(irepl=1:nrepl, .combine=rbind, .packages=("MASS")) %dopar%      
                {
                        #Generate data
                        data=datagen(n, alfas, betas, xmeans, xsd, xcor)
                        
                        #EM algorithm. Note binomial model embedded in EM
                        em_res = em_alg(data$y1,data$x)
                        
                        #Perfect info trinomial logit with start true*2
                        mlogit=optim(c(alfas, betas)*2,logl,x=data$x,y=data$y) 
                        
                        #Create list of output
                        c(n, px, em_res$iter, nrepl, sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]),
                          alfas, betas, em_res$em_a, em_res$em_b, em_res$binomial_coeff, mlogit$par)
                } )
})
stopCluster(cl) 
colnames(dfs) <- givecolumnnames()

#End simulation 10 ####