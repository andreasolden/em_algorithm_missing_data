# Andreas bootstrap parallel

myinstall <- function(x){
        for( i in x ){
                #  require returns TRUE invisibly if it was able to load package
                if( ! require( i , character.only = TRUE ) ){
                        #  If package was not able to be loaded then re-install
                        install.packages( i , dependencies = TRUE )
                        #  Load package after installing
                        require( i , character.only = TRUE )
                }
        }
}

#  Then try/install packages...
myinstall( c("foreach" , "doParallel", "tidyverse") )

#Andreas bootstrap 
library(foreach) #For running simulation
library(doParallel) #For parallelization
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
em_alg2=function(y1,x,niter=100,conv=1e-6,verbose=FALSE) #y1=caught, x=matrix(10*rnorm(px*n),n,px)
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
        
        return(list(em_a=em_a,em_b=em_b, yy = yy))
}

# SIMULATE DATA FROM 3NOMIAL LOGISTIC REGRESSION

nrepl=24
bootrepl=1000


start_time <- Sys.time()

cl <- makeCluster(detectCores())
registerDoParallel(cl)

monteboot <- 
        foreach(irepl=1:nrepl, .combine=rbind) %dopar% {
                
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
                mod1=em_alg2(y1,x,conv=1e-6,niter=100)
                
                thpr=c(mod1$em_a,mod1$em_b)
                
                # Create loop that draws n samples with replacement from x
                # Estimates same model
                # Saves all estimates? 
                #
                
                df_boot = cbind(y1, x)
                
                prop_sigma=matrix(0,bootrepl,4)
                
                for(ibot in 1:bootrepl)
                {
                        xboot = df_boot[sample(nrow(df_boot), size = n, replace = TRUE),]
                        
                        y1boot = xboot[, 1]
                        xboot = xboot[, 2:3]
                        
                        modboot = em_alg2(y1boot, xboot, conv=1e-6, niter=100)
                        
                        prop_sigma[ibot,] = c(modboot$em_a,modboot$em_b)
                }
                
                
                #Out
                
                
                
                thpr=as.vector(thpr)
                ses=as.vector(apply(prop_sigma, 2, sd))
                
                c(thpr, ses)
                
                #        c(n, px, em_res$iter, nrepl, sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]),
                #alfas, betas, em_res$em_a, em_res$em_b, em_res$binomial_coeff, mlogit$par)
                
        }

stopCluster(cl)

end_time <- Sys.time()

tot_time <- end_time - start_time

#head(pars)
#head(ses)

head(monteboot)

tot_time
apply(monteboot, 2, function(x) sum(is.na(x)))


#thtrue
colMeans(monteboot, na.rm = TRUE)


saveRDS(monteboot, file = "monteboot.rds")
saveRDS(tot_time, file = "tot_time_monteboot.rds")














