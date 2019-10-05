#Testing functions 2

#Testing functions
library(here)# allows us to use relative file paths
library(doSNOW); library(foreach)
library(MASS)

################################################################################
#This script defines our functions and nothing more 
################################################################################

#Log-likelihood function
################################################################################
logl=function(th,x,y) #
{
        x=as.matrix(x)
        px=ncol(x)
        a=th[1:2]     
        b=th[3:(2+px)]
        p=matrix(1,n,3) #Note, used to be zeroes and not ones
        for(k in 1:2) p[,k+1]=exp(a[k]+x%*%b)
        p=p/apply(p,1,sum)
        ll=sum(y*log(p))
        mll=-ll
        return(mll)
}

# Function to create column names that takes into account px ####
givecolumnnames <- function(){
        sim_names <- c("N", "px", "iter", "nrepl", "y1", "y2", "y3")
        trueval_names <- c("Alfa1", "Alfa2", paste0("Beta", 1:px))
        em_names <- c("EMa1", "EMa2", paste0("EMb", 1:px))
        bi_names <- c("BIa1", paste0("BIb", 1:px))
        tri_names <- c("Tri_a1", "Tri_a2", paste0("Trib", 1:px))
        columnnames <- c(sim_names, trueval_names, em_names, bi_names, tri_names)
        return(columnnames)
}

# Function that tests if the parameters of a simulation seems consistent ####
test_simparam <- function(n, px, alfas, betas, xmeans, xsd, xcor){
        if(n>0 && alfas ==2 && length(betas)==px && length(xmeans)==px && 
           length(xsd)==px && nrow(xcor)==px && ncol(xcor)==px && sum(diag(xcor))==px)
        {print("Error. Input inconsistent") } else {print("Input seems ok")}
}

#EM-algorithm ####
em_alg=function(y1,x,niter=100,conv=1e-5,verbose=TRUE) #y1=caught, x=matrix(10*rnorm(px*n),n,px)
{
        x=as.matrix(x) #re-lagrer x som matries
        logit=glm(data$y1 ~ data$x,family=binomial())
        binomial_coeff = logit$coefficients
        
        them= c(logit$coefficients[[1]], binomial_coeff)# theta em, vector (1,1,1,0,0) starting values
        themold=0*them # Vector (0,0,0,0,0)
        iter=1 #number of iterations
        while(sum((them-themold)^2)>1e-5 & iter<niter) #Double convergence criteria
        {
                themold=them #save old to get criteria
                em_a=them[1:2] #TWO first numbers of them
                em_b=them[-(1:2)] # NOT 1 and 2 of them
                ystar=matrix(1,n,2) #matrix full of ones
                ystar[,2]=exp(em_a[1]+x%*%em_b) #second column of Ystar
                ystar=ystar/apply(ystar,1,sum)
                yy=cbind(ystar,y1)
                yy[yy[,3]==1,1:2]=0
                mlogitem=optim(par = them,fn = logl,x=x,y=yy,hessian=TRUE)
                them=mlogitem$par ###THE LOGL FORCES col 2 and 3 to have same coeff. 
                iter=iter+1
        }
        return(list(em_a=em_a,em_b=em_b,yy=yy,iter=iter, binomial_coeff=binomial_coeff))
}

# define a data-generating process ####

datagen <- function(n, alfas, betas, xmeans, xsd, xcor){
        
        #Create covariance matrix from correlation matrix
        intermediate <- xsd %*% t(xsd) #i*j t is transposed
        covar_matrix <- intermediate * xcor
        #Use mvrnorm to draw correlated data based on covariance matrix
        x=as.matrix(mvrnorm(n=n, mu=xmeans, Sigma=covar_matrix))
        
        #Create outcome variables
        p=matrix(1,n,3) #probability matrix, filled with 1's to start
        for(k in 1:2) p[,k+1]=exp(alfas[k]+x%*%betas) # fill p's by formula
        p=p/apply(p,1,sum) # force sum to one
        y=matrix(0,n,3) # Creates outcome matrix, filled with 0's to start
        for(i in 1:n) y[i,]=rmultinom(1,1,prob=p[i,]) # y's based on p's
        y1=y[,3] # Separate observed fraudulent claims for easy operations later
        return(list(p=p,y=y,y1=y1, x=x)) #return list of elements 
} 


#Simulations 1 ####
n=200
px=3
alfas=c(-1.8,-1.5)
betas=c(-0.02, 0.2, 0.3)
xmeans=c(0,0,0)
xsd=c(12.3, 4, 6)
xcor <- matrix(c(1, 0.5, 0.4,
                 0.5, 1, 0.6,
                 0.4, 0.6, 1),
               nrow = px, ncol = px)
nrepl=16
niter=50 #Default 100
n=300
# conv=1e-5

test_simparam(n, px, alfas, betas, xmeans, xsd, xcor)

set.seed(123)
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

#saveRDS(dfs, file = here("results", "dfs1.rds"))



