#Simulations v2

#This script runs the Monte Carlo simulations and saves the results in folder results

#Initiate project
library(checkpoint) #Checkpoint assures us that we use the same package versions
checkpoint("2019-01-21") #As they were at the date set in this function. i.e "2019-01-21"

library(here)# allows us to use relative file paths
#Since this is an Rproject you should always open Rproj to initiate session
#This exports project options and sets WD to folder of Rproj file so here works
#Example: here("folderA", "file.extension") or  here("folderA/file.extension")

#Load packages for the simulations
library(doSNOW); library(foreach)

#Source script with functions 
source(here("scripts", "define_functions.R"))


#Begin simulations


#Simulations 1 ####
#One variable, alfa(-1.8,-1.5), beta(-0.02), xval(0), xsd(12.3)

#Define parameters ####
set.seed(123)
a=c((-1.5*1.2),-1.5)
b=c(-0.02, 0.2)
xval=c(0,0)
sd=c(12.3, 1.8)
px=1
xcor=0.0
nrepl=1000
niter=100
n=1000
cl<-makeCluster(4); registerDoSNOW(cl)

# actual simulation ####
system.time({
        dfs=foreach(irepl=1:nrepl,.combine=rbind, .packages=("MASS")) %dopar%      
        {
                data=datagen2(px=px, n=n, a=a, b=b, xval=xval, sd=sd, xcor=xcor)
                
                #EM and BI
                
                # CHECK IF MODEL ABLE TO SEPARATE 00 FROM 01 CORRECTLY
                emmod=emhid2(data$y1,data$x)
                
                #OUT
                commonout = c(n, px, emmod$iter, nrepl, ((a[1]-a[2])/a[1]),
                              sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]))
                
                BIcoeffout = emmod$BIcoeffout# c(emmod$BIcoeffout[[1]], emmod$BIcoeffout[[2:(1+px)]])
                
                Realcoeffout = data$thtrue
                
                EMcoeffout = c(emmod$a,emmod$b)
                
                #TRINOMIAL LOGIT w/ PI
                thpi=data$thtrue*2 #Need some starting values. [2:(3+px)]?
                mlogit=optim(thpi,logl,x=data$x,y=data$y) #logl function with three cat. 
                alogit=mlogit$par[1:2]
                blogit=mlogit$par[3:(2+px)]
                
                PIcoeffout = mlogit$par
                
                #c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
                c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
        } 
})
stopCluster(cl) 

dfs=as.data.frame(dfs)
if(px==1){
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb", "EMa1",
                           "EMa2" , "EMb", "BIa", "BIb", "PIa1", "PIa2", "PIb")
} else{
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb1", "tb2", 
                           "EMa1", "EMa2" , "EMb1", "EMb2", "BIa", "BIb1", "BIb2", "PIa1", "PIa2", "PIb1", 
                           "PIb2") 
}

#Save results ####
saveRDS(dfs, file = here("results", "dfs1.rds"))

#End simulation 1 ####


#Simulations 2 ####
#One variable, changed alfa1: alfa(-3,-1.5). This is 200% diff instead of 20%

#Define parameters ####
set.seed(123)
a=c((-1.5*2),-1.5)
b=c(-0.02, 0.2)
xval=c(0,0)
sd=c(12.3, 1.8)
px=1
xcor=0.0
nrepl=1000
niter=100
n=1000
cl<-makeCluster(4); registerDoSNOW(cl)

# actual simulation ####
system.time({
        dfs=foreach(irepl=1:nrepl,.combine=rbind, .packages=("MASS")) %dopar%      
        {
                data=datagen2(px=px, n=n, a=a, b=b, xval=xval, sd=sd, xcor=xcor)
                
                #EM and BI
                
                # CHECK IF MODEL ABLE TO SEPARATE 00 FROM 01 CORRECTLY
                emmod=emhid2(data$y1,data$x)
                
                #OUT
                commonout = c(n, px, emmod$iter, nrepl, ((a[1]-a[2])/a[1]),
                              sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]))
                
                BIcoeffout = emmod$BIcoeffout# c(emmod$BIcoeffout[[1]], emmod$BIcoeffout[[2:(1+px)]])
                
                Realcoeffout = data$thtrue
                
                EMcoeffout = c(emmod$a,emmod$b)
                
                #TRINOMIAL LOGIT w/ PI
                thpi=data$thtrue*2 #Need some starting values. [2:(3+px)]?
                mlogit=optim(thpi,logl,x=data$x,y=data$y) #logl function with three cat. 
                alogit=mlogit$par[1:2]
                blogit=mlogit$par[3:(2+px)]
                
                PIcoeffout = mlogit$par
                
                #c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
                c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
        } 
})
stopCluster(cl) 

dfs=as.data.frame(dfs)
if(px==1){
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb", "EMa1",
                           "EMa2" , "EMb", "BIa", "BIb", "PIa1", "PIa2", "PIb")
} else{
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb1", "tb2", 
                           "EMa1", "EMa2" , "EMb1", "EMb2", "BIa", "BIb1", "BIb2", "PIa1", "PIa2", "PIb1", 
                           "PIb2") 
}

#Save results ####
saveRDS(dfs, file = here("results", "dfs2.rds"))

#End simulation 2 ####


#Simulations 3 ####
#One variable, changed standard dev of x-var (123 instead of 12.3)

#Define parameters ####
set.seed(123)
a=c((-1.5*1.2),-1.5)
b=c(-0.02, 0.2)
xval=c(0,0)
sd=c(123, 1.8)
px=1
xcor=0.0
nrepl=1000
niter=100
n=1000
cl<-makeCluster(4); registerDoSNOW(cl)

# actual simulation ####
system.time({
        dfs=foreach(irepl=1:nrepl,.combine=rbind, .packages=("MASS")) %dopar%      
        {
                data=datagen2(px=px, n=n, a=a, b=b, xval=xval, sd=sd, xcor=xcor)
                
                #EM and BI
                
                # CHECK IF MODEL ABLE TO SEPARATE 00 FROM 01 CORRECTLY
                emmod=emhid2(data$y1,data$x)
                
                #OUT
                commonout = c(n, px, emmod$iter, nrepl, ((a[1]-a[2])/a[1]),
                              sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]))
                
                BIcoeffout = emmod$BIcoeffout# c(emmod$BIcoeffout[[1]], emmod$BIcoeffout[[2:(1+px)]])
                
                Realcoeffout = data$thtrue
                
                EMcoeffout = c(emmod$a,emmod$b)
                
                #TRINOMIAL LOGIT w/ PI
                thpi=data$thtrue*2 #Need some starting values. [2:(3+px)]?
                mlogit=optim(thpi,logl,x=data$x,y=data$y) #logl function with three cat. 
                alogit=mlogit$par[1:2]
                blogit=mlogit$par[3:(2+px)]
                
                PIcoeffout = mlogit$par
                
                #c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
                c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
        } 
})
stopCluster(cl) 

dfs=as.data.frame(dfs)
if(px==1){
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb", "EMa1",
                           "EMa2" , "EMb", "BIa", "BIb", "PIa1", "PIa2", "PIb")
} else{
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb1", "tb2", 
                           "EMa1", "EMa2" , "EMb1", "EMb2", "BIa", "BIb1", "BIb2", "PIa1", "PIa2", "PIb1", 
                           "PIb2") 
}

#Save results ####
saveRDS(dfs, file = here("results", "dfs3.rds"))

#End simulation 3 ####


#Simulations 4 ####
#Two variables, alfa(-1.8,-1.5), beta(-0.02, 0.2), xval(0), xsd(12.3, 1.8), corr(0)

#Define parameters ####
set.seed(123)
a=c((-1.5*1.2),-1.5)
b=c(-0.02, 0.2)
xval=c(0,0)
sd=c(12.3, 1.8)
px=2
xcor=0.0
nrepl=1000
niter=100
n=1000
cl<-makeCluster(4); registerDoSNOW(cl)

# actual simulation ####
system.time({
        dfs=foreach(irepl=1:nrepl,.combine=rbind, .packages=("MASS")) %dopar%      
        {
                data=datagen2(px=px, n=n, a=a, b=b, xval=xval, sd=sd, xcor=xcor)
                
                #EM and BI
                
                # CHECK IF MODEL ABLE TO SEPARATE 00 FROM 01 CORRECTLY
                emmod=emhid2(data$y1,data$x)
                
                #OUT
                commonout = c(n, px, emmod$iter, nrepl, ((a[1]-a[2])/a[1]),
                              sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]))
                
                BIcoeffout = emmod$BIcoeffout# c(emmod$BIcoeffout[[1]], emmod$BIcoeffout[[2:(1+px)]])
                
                Realcoeffout = data$thtrue
                
                EMcoeffout = c(emmod$a,emmod$b)
                
                #TRINOMIAL LOGIT w/ PI
                thpi=data$thtrue*2 #Need some starting values. [2:(3+px)]?
                mlogit=optim(thpi,logl,x=data$x,y=data$y) #logl function with three cat. 
                alogit=mlogit$par[1:2]
                blogit=mlogit$par[3:(2+px)]
                
                PIcoeffout = mlogit$par
                
                #c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
                c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
        } 
})
stopCluster(cl) 

dfs=as.data.frame(dfs)
if(px==1){
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb", "EMa1",
                           "EMa2" , "EMb", "BIa", "BIb", "PIa1", "PIa2", "PIb")
} else{
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb1", "tb2", 
                           "EMa1", "EMa2" , "EMb1", "EMb2", "BIa", "BIb1", "BIb2", "PIa1", "PIa2", "PIb1", 
                           "PIb2") 
}

#Save results ####
saveRDS(dfs, file = here("results", "dfs4.rds"))

#End simulation 4 ####


#Simulations 5 ####
#Two variables, introduce 0.5 correlation between x's

#Define parameters ####
set.seed(123)
a=c((-1.5*1.2),-1.5)
b=c(-0.02, 0.2)
xval=c(0,0)
sd=c(12.3, 1.8)
px=2
xcor=0.5
nrepl=1000
niter=100
n=1000
cl<-makeCluster(4); registerDoSNOW(cl)

# actual simulation ####
system.time({
        dfs=foreach(irepl=1:nrepl,.combine=rbind, .packages=("MASS")) %dopar%      
        {
                data=datagen2(px=px, n=n, a=a, b=b, xval=xval, sd=sd, xcor=xcor)
                
                #EM and BI
                
                # CHECK IF MODEL ABLE TO SEPARATE 00 FROM 01 CORRECTLY
                emmod=emhid2(data$y1,data$x)
                
                #OUT
                commonout = c(n, px, emmod$iter, nrepl, ((a[1]-a[2])/a[1]),
                              sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]))
                
                BIcoeffout = emmod$BIcoeffout# c(emmod$BIcoeffout[[1]], emmod$BIcoeffout[[2:(1+px)]])
                
                Realcoeffout = data$thtrue
                
                EMcoeffout = c(emmod$a,emmod$b)
                
                #TRINOMIAL LOGIT w/ PI
                thpi=data$thtrue*2 #Need some starting values. [2:(3+px)]?
                mlogit=optim(thpi,logl,x=data$x,y=data$y) #logl function with three cat. 
                alogit=mlogit$par[1:2]
                blogit=mlogit$par[3:(2+px)]
                
                PIcoeffout = mlogit$par
                
                #c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
                c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
        } 
})
stopCluster(cl) 

dfs=as.data.frame(dfs)
if(px==1){
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb", "EMa1",
                           "EMa2" , "EMb", "BIa", "BIb", "PIa1", "PIa2", "PIb")
} else{
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb1", "tb2", 
                           "EMa1", "EMa2" , "EMb1", "EMb2", "BIa", "BIb1", "BIb2", "PIa1", "PIa2", "PIb1", 
                           "PIb2") 
}

#Save results ####
saveRDS(dfs, file = here("results", "dfs5.rds"))

#End simulation 5 ####


#Simulations 6 ####
#Two variables, introduce 0.9 correlation between x's

#Define parameters ####
set.seed(123)
a=c((-1.5*1.2),-1.5)
b=c(-0.02, 0.2)
xval=c(0,0)
sd=c(12.3, 1.8)
px=2
xcor=0.9
nrepl=1000
niter=100
n=1000
cl<-makeCluster(4); registerDoSNOW(cl)

# actual simulation ####
system.time({
        dfs=foreach(irepl=1:nrepl,.combine=rbind, .packages=("MASS")) %dopar%      
        {
                data=datagen2(px=px, n=n, a=a, b=b, xval=xval, sd=sd, xcor=xcor)
                
                #EM and BI
                
                # CHECK IF MODEL ABLE TO SEPARATE 00 FROM 01 CORRECTLY
                emmod=emhid2(data$y1,data$x)
                
                #OUT
                commonout = c(n, px, emmod$iter, nrepl, ((a[1]-a[2])/a[1]),
                              sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]))
                
                BIcoeffout = emmod$BIcoeffout# c(emmod$BIcoeffout[[1]], emmod$BIcoeffout[[2:(1+px)]])
                
                Realcoeffout = data$thtrue
                
                EMcoeffout = c(emmod$a,emmod$b)
                
                #TRINOMIAL LOGIT w/ PI
                thpi=data$thtrue*2 #Need some starting values. [2:(3+px)]?
                mlogit=optim(thpi,logl,x=data$x,y=data$y) #logl function with three cat. 
                alogit=mlogit$par[1:2]
                blogit=mlogit$par[3:(2+px)]
                
                PIcoeffout = mlogit$par
                
                #c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
                c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
        } 
})
stopCluster(cl) 

dfs=as.data.frame(dfs)
if(px==1){
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb", "EMa1",
                           "EMa2" , "EMb", "BIa", "BIb", "PIa1", "PIa2", "PIb")
} else{
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb1", "tb2", 
                           "EMa1", "EMa2" , "EMb1", "EMb2", "BIa", "BIb1", "BIb2", "PIa1", "PIa2", "PIb1", 
                           "PIb2") 
}

#Save results ####
saveRDS(dfs, file = here("results", "dfs6.rds"))

#End simulation 6 ####


#Simulations 7 ####
#Two variables, same as base but with sd*10 for X1 only

#Define parameters ####
set.seed(123)
a=c((-1.5*1.2),-1.5)
b=c(-0.02, 0.2)
xval=c(0,0)
sd=c(123, 1.8)
px=2
xcor=0.0
nrepl=1000
niter=100
n=1000
cl<-makeCluster(4); registerDoSNOW(cl)

# actual simulation ####
system.time({
        dfs=foreach(irepl=1:nrepl,.combine=rbind, .packages=("MASS")) %dopar%      
        {
                data=datagen2(px=px, n=n, a=a, b=b, xval=xval, sd=sd, xcor=xcor)
                
                #EM and BI
                
                # CHECK IF MODEL ABLE TO SEPARATE 00 FROM 01 CORRECTLY
                emmod=emhid2(data$y1,data$x)
                
                #OUT
                commonout = c(n, px, emmod$iter, nrepl, ((a[1]-a[2])/a[1]),
                              sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]))
                
                BIcoeffout = emmod$BIcoeffout# c(emmod$BIcoeffout[[1]], emmod$BIcoeffout[[2:(1+px)]])
                
                Realcoeffout = data$thtrue
                
                EMcoeffout = c(emmod$a,emmod$b)
                
                #TRINOMIAL LOGIT w/ PI
                thpi=data$thtrue*2 #Need some starting values. [2:(3+px)]?
                mlogit=optim(thpi,logl,x=data$x,y=data$y) #logl function with three cat. 
                alogit=mlogit$par[1:2]
                blogit=mlogit$par[3:(2+px)]
                
                PIcoeffout = mlogit$par
                
                #c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
                c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
        } 
})
stopCluster(cl) 

dfs=as.data.frame(dfs)
if(px==1){
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb", "EMa1",
                           "EMa2" , "EMb", "BIa", "BIb", "PIa1", "PIa2", "PIb")
} else{
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb1", "tb2", 
                           "EMa1", "EMa2" , "EMb1", "EMb2", "BIa", "BIb1", "BIb2", "PIa1", "PIa2", "PIb1", 
                           "PIb2") 
}

#Save results ####
saveRDS(dfs, file = here("results", "dfs7.rds"))

#End simulation 7 ####


#Simulations 8 ####
#Two variables, same as base but with sd*10 for X1 AND X2

#Define parameters ####
set.seed(123)
a=c((-1.5*1.2),-1.5)
b=c(-0.02, 0.2)
xval=c(0,0)
sd=c(123, 18)
px=2
xcor=0.0
nrepl=1000
niter=100
n=1000
cl<-makeCluster(4); registerDoSNOW(cl)

# actual simulation ####
system.time({
        dfs=foreach(irepl=1:nrepl,.combine=rbind, .packages=("MASS")) %dopar%      
        {
                data=datagen2(px=px, n=n, a=a, b=b, xval=xval, sd=sd, xcor=xcor)
                
                #EM and BI
                
                # CHECK IF MODEL ABLE TO SEPARATE 00 FROM 01 CORRECTLY
                emmod=emhid2(data$y1,data$x)
                
                #OUT
                commonout = c(n, px, emmod$iter, nrepl, ((a[1]-a[2])/a[1]),
                              sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]))
                
                BIcoeffout = emmod$BIcoeffout# c(emmod$BIcoeffout[[1]], emmod$BIcoeffout[[2:(1+px)]])
                
                Realcoeffout = data$thtrue
                
                EMcoeffout = c(emmod$a,emmod$b)
                
                #TRINOMIAL LOGIT w/ PI
                thpi=data$thtrue*2 #Need some starting values. [2:(3+px)]?
                mlogit=optim(thpi,logl,x=data$x,y=data$y) #logl function with three cat. 
                alogit=mlogit$par[1:2]
                blogit=mlogit$par[3:(2+px)]
                
                PIcoeffout = mlogit$par
                
                #c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
                c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
        } 
})
stopCluster(cl) 

dfs=as.data.frame(dfs)
if(px==1){
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb", "EMa1",
                           "EMa2" , "EMb", "BIa", "BIb", "PIa1", "PIa2", "PIb")
} else{
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb1", "tb2", 
                           "EMa1", "EMa2" , "EMb1", "EMb2", "BIa", "BIb1", "BIb2", "PIa1", "PIa2", "PIb1", 
                           "PIb2") 
}

#Save results ####
saveRDS(dfs, file = here("results", "dfs8.rds"))

#End simulation 8 ####


#Simulations 9 ####
#Two variables, same as base but with sd*10 for X1 AND X2 AND CORR(0.9)

#Define parameters ####
set.seed(123)
a=c((-1.5*1.2),-1.5)
b=c(-0.02, 0.2)
xval=c(0,0)
sd=c(123, 18)
px=2
xcor=0.9
nrepl=1000
niter=100
n=1000
cl<-makeCluster(4); registerDoSNOW(cl)

# actual simulation ####
system.time({
        dfs=foreach(irepl=1:nrepl,.combine=rbind, .packages=("MASS")) %dopar%      
        {
                data=datagen2(px=px, n=n, a=a, b=b, xval=xval, sd=sd, xcor=xcor)
                
                #EM and BI
                
                # CHECK IF MODEL ABLE TO SEPARATE 00 FROM 01 CORRECTLY
                emmod=emhid2(data$y1,data$x)
                
                #OUT
                commonout = c(n, px, emmod$iter, nrepl, ((a[1]-a[2])/a[1]),
                              sum(data$y[,1]), sum(data$y[,2]), sum(data$y[,3]))
                
                BIcoeffout = emmod$BIcoeffout# c(emmod$BIcoeffout[[1]], emmod$BIcoeffout[[2:(1+px)]])
                
                Realcoeffout = data$thtrue
                
                EMcoeffout = c(emmod$a,emmod$b)
                
                #TRINOMIAL LOGIT w/ PI
                thpi=data$thtrue*2 #Need some starting values. [2:(3+px)]?
                mlogit=optim(thpi,logl,x=data$x,y=data$y) #logl function with three cat. 
                alogit=mlogit$par[1:2]
                blogit=mlogit$par[3:(2+px)]
                
                PIcoeffout = mlogit$par
                
                #c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
                c(commonout,Realcoeffout,EMcoeffout,BIcoeffout,PIcoeffout)
        } 
})
stopCluster(cl) 

dfs=as.data.frame(dfs)
if(px==1){
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb", "EMa1",
                           "EMa2" , "EMb", "BIa", "BIb", "PIa1", "PIa2", "PIb")
} else{
        colnames(dfs) <- c("n", "px", "iter", "nrepl", "ad", "y1", "y2", "y3", "ta1", "ta2", "tb1", "tb2", 
                           "EMa1", "EMa2" , "EMb1", "EMb2", "BIa", "BIb1", "BIb2", "PIa1", "PIa2", "PIb1", 
                           "PIb2") 
}

#Save results ####
saveRDS(dfs, file = here("results", "dfs9.rds"))

#End simulation 9 ####

