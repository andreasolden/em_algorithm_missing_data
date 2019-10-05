# Andreas bootstrap parallel

library(foreach); library(doParallel); library(doRNG)

#library(profvis); library(rbenchmark)


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
em_alg2=function(y1,x,niter=100,conv=1e-6, onlyparam=FALSE) #y1=caught, x=matrix(10*rnorm(px*n),n,px)
{
        x=as.matrix(x) #ensure matrix
        n=nrow(x)
        them=c(-1.8,-1.5,-0.02,-0.2)
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
                #for(k in 1:2) p[,k+1]=exp(a[k]+x%*%b) # Original
                p[,2]=exp(a[1]+x%*%b) #Uglier, but a lot faster (127x)
                p[,3]=exp(a[2]+x%*%b)
                p=p/apply(p,1,sum)
                #for(i in 1:n) y[i,]=rmultinom(1,1,prob=p[i,]) # Original
                y=t(apply(p, 1, function(x) rmultinom(1,1,x))) #ALso faster than original
                y1=y[,3] # OBSERVED FRAUDULANT CLAIMS
                
                # EM-ALG
                emres=em_alg2(y1,x,conv=1e-6,niter=100)
                
                #Param
                thpr=c(emres$em_a,emres$em_b)
                
                # Bootstrap
                
                aboot = emres$em_a
                bboot = emres$em_b
                
                prop_sigma=matrix(0,bootrepl,4)
                
                for(ibot in 1:bootrepl)
                {
                        
                        xboot = rnorm( length(x), mean=mean(x), sd=sd(x) )
                        
                        xboot[,1]=sd(x)*xboot[,1]
                        xboot[,2]=sd(x)*xboot[,2]
                        
                        pboot=matrix(1,n,3)

                        pboot[,2]=exp(aboot[1]+xboot%*%bboot) #Uglier, but a lot faster (127x)
                        pboot[,3]=exp(aboot[2]+xboot%*%bboot)
                        pboot=pboot/apply(pboot,1,sum)

                        yboot=t(apply(pboot, 1, function(x) rmultinom(1,1,xboot))) #ALso faster than original
                        y1boot=yboot[,3] # OBSERVED FRAUDULANT CLAIMS
                        
                        xboot = cbind(y1boot, xboot)

                        modboot = em_alg2(xboot[,1], xboot[,2:(px+1)], conv=1e-6, niter=100, onlyparam = TRUE)
                        
                        prop_sigma[ibot,] = c(modboot$em_a,modboot$em_b)
                }
                
                #SES
                ses=as.vector(apply(prop_sigma, 2, sd))
                rm(prop_sigma)
                
                #PI Trinomial
                mlogitpi = optim(thtrue,logl,x=x,y=y, method = "BFGS", hessian = TRUE)
                mlogitpises = sqrt(diag(solve(mlogitpi$hessian)))
                #Binomial logit
                mlogitbi=summary(glm(y1 ~ x,family=binomial()))
                
                ############################################################################################################
                
                #OUT
                c(thpr, ses, emres$iter, sum(y[,1]), sum(y[,2]), sum(y[,3]), mlogitpi$par, mlogitpises, mlogitbi$coefficients[,1], mlogitbi$coefficients[,2])
                
        }

#Ac(n, px, em_res$iter, nrepl, sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]),
#alfas, betas, em_res$em_a, em_res$em_b, em_res$binomial_coeff, mlogit$par)
#} )

stopCluster(cl)

end_time = Sys.time()

tot_time = start_time - end_time
tot_time

colnames(monteboot) <- c("EMa1", "EMa2", "EMb1", "EMb2", "EMa1SE", "EMa2SE", "EMb1SE", "EMb2SE", 
                         "iter", "Y1", "Y2", "Y3",  
                         "PIa1", "PIa2", "PIb1", "PIb2", "PIa1SE", "PIa2SE", "PIb1SE", "PIb2SE",
                         "BIa1", "BIb1", "BIb2", "BIa1SE", "BIb1SE", "BIb2SE"
)

#Controls missing and mean
apply(monteboot, 2, function(x) sum(is.na(x)))
colMeans(monteboot, na.rm = TRUE)

saveRDS(monteboot, file = "monteboot.rds")
saveRDS(tot_time, file = "tot_time_monteboot.rds")












