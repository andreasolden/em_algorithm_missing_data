library(rbenchmark)

# LOG-LIKELIHOOD
logl=function(th,x,y)
{
        x=as.matrix(x)
        px=ncol(x)
        a=th[1:2]
        b=th[3:(2+px)]
        p=matrix(1,n,3)
        p[,2]=exp(a[1]+x%*%b) #Styggere, men 127 ganger kjappere
        p[,3]=exp(a[2]+x%*%b)
        p=p/apply(p,1,sum)
        ll=sum(y*log(p))
        mll=-ll
        return(mll)
}



# EM-ALGORITHM TO ESTIMATE PARAMETERS
em_alg2=function(y1,x,niter=100,conv=1e-6) #y1=caught, x=matrix(10*rnorm(px*n),n,px)
{
        x=as.matrix(x) #re-lagrer x som matries
        n=nrow(x)
        them=c(-1.8,-1.5,-0.02,-0.2)
        themold=0*them # Vector (0,0,0,0,0)
        iter=1 #number of iterations
        while(sum((them-themold)^2/(themold+0.0001)^2)>conv && iter<niter) #&& sparer tid da den ikke evaluerer kriterie 2 når 1 bryter
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
        
        return(list(em_a=em_a,em_b=em_b, yy = yy))
}


nrepl=3
bootrepl=10
n=300

px=2
x=matrix(rnorm(px*n),n,px)
x[,1]=12.3*x[,1]
x[,2]=1.8*x[,2]

a=c(-1.8,-1.5)
b=c(-0.02,0.2)

thtrue=c(a,b)

p=matrix(1,n,3)
p[,2]=exp(a[1]+x%*%b)
p[,3]=exp(a[2]+x%*%b)
p=p/apply(p,1,sum)


y=t(apply(p, 1, function(x) rmultinom(1,1,x)))
y1=y[,3] # OBSERVED FRAUDULANT CLAIMS
df_boot = cbind(y1, x)

my_boot_function <- 
        function(x) {
                xboot = df_boot[sample(nrow(x), size = n, replace = TRUE),]
                modboot = em_alg2(xboot[,1], xboot[,2:3], conv=1e-6, niter=20)
                return(c(modboot$em_a,modboot$em_b))
        }

benchmark(
        test1 = {
                prop_sigma=matrix(0,bootrepl,4)
                for(ibot in 1:bootrepl) {
                prop_sigma[ibot,] = my_boot_function(df_boot)
                }
                
        },
        test2 = {
                prop_sigma=matrix(0,bootrepl,4)
                
                for(ibot in 1:bootrepl)
                {
                        xboot = df_boot[sample(nrow(df_boot), size = n, replace = TRUE),]
                        
                        #y1boot = xboot[, 1]
                        #xboot = xboot[, 2:3]
                        
                        modboot = em_alg2(xboot[,1], xboot[,2:3], conv=1e-6, niter=20)
                        
                        prop_sigma[ibot,] = c(modboot$em_a,modboot$em_b)
                } 
        },
        replications = 100
)




#Alternative code



prop_sigma=matrix(0,bootrepl,4)
for(ibot in 1:bootrepl) {
        prop_sigma[ibot,] = my_boot_function(df_boot)
}

testtest <- my_boot_function(df_boot)

# Original code 
prop_sigma=matrix(0,bootrepl,4)

for(ibot in 1:bootrepl)
{
        xboot = df_boot[sample(nrow(df_boot), size = n, replace = TRUE),]
        
        y1boot = xboot[, 1]
        xboot = xboot[, 2:3]
        
        modboot = em_alg2(xboot[,1], xboot[,2:3], conv=1e-6, niter=20)
        
        prop_sigma[ibot,] = c(modboot$em_a,modboot$em_b)
}
