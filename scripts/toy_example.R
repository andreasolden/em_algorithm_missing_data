#Toy example 

set.seed(123)

# LOG-LIKELIHOOD
logl=function(th,x,y)
{
        x=as.matrix(x)
        px=ncol(x)
        a=th[1:2]
        b=th[3:(2+px)]
        p=matrix(1,n,3)
        for(k in 1:2) p[,k+1]=exp(a[k]+x%*%b) #SHOULD THIS BE a[i+1]? NO! 
        p=p/apply(p,1,sum)
        ll=sum(y*log(p))
        mll=-ll
        return(mll)
}

# FIT MODEL WITH UNOBSERVED STATE
emhidden=function(y1,x,niter=100,conv=1e-5,verbose=TRUE)
{
        
        x=as.matrix(x)
        px=ncol(x)
        them=c(rep(1,2),rep(0,px))
        themold=c(1.0, 2.0, 0.5, 0.8)
        #themold=0*them
        iter=1
        x=as.matrix(x)
        while(sum((them-themold)^2)>1e-5 & iter<niter)
        {
                if(iter==1) print(c("iteration","a2","a3",paste("b",1:px,sep="")))
                themold=them
                a=them[1:2]
                b=them[-(1:2)]
                ystar=matrix(1,n,2)
                #ystar[,1]=exp(a[1]+x%*%b)
                ystar[,2]=exp(a[1]+x%*%b)
                ystar=ystar/apply(ystar,1,sum)
                yy=cbind(ystar,y1)
                yy[yy[,3]==1,1:2]=0
                mlogitem=optim(par = them,fn = logl,x=x,y=yy,hessian=TRUE)
                them=mlogitem$par
                if(verbose) print(c(iter,them))
                iter=iter+1
        }
        return(list(yy=yy,x=x,a=a,b=b))
        #return(yy)
}

# PREDICT
#This funtion takes whatever is saved as mod, and extracts a, b and x. 
#Sets n=nrow, then makes probs. EM gives estimates, this gives probs?
probs=function(mod)
{
        a=mod$a
        b=mod$b
        x=mod$x
        n=nrow(x)
        for(i in 1:2) p[,i+1]=exp(a[i]+x%*%b) 
        p=p/apply(p,1,sum)
        return(p)
}


# # CONVERT INDICATORS TO FACTOR
# ind2fac=function(x)
# {
#     for(i in 1:length(x))
#     {
#         if(x[i]) k=i
#     }
#     return(k)
# }



# SIMULATE DATA FROM 3NOMIAL LOGISTIC REGRESSION
# P(Z=k)=exp(a[k]+b[k]*x)/(1+sum(exp(a[k]+b[k]*x),k=1,2))
# P(Z=3)=1/(1+sum(exp(a[k]+b[k]*x),k=1,2))
n=1000
px=2
#x=matrix(10*rnorm(px*n),n,px) #DEtte er Jonas sin original, den under er min for å ødelegge. 
x=matrix(rnorm(px*n,10,5),n,px)
a=c(1,2)
b=c(0.5,0.8)
thtrue=c(a,b)
p=matrix(1,n,3)
for(k in 1:2) p[,k+1]=exp(a[k]+x%*%b) #
p=p/apply(p,1,sum)
y=matrix(0,n,3)
for(i in 1:n) y[i,]=rmultinom(1,1,prob=p[i,])
y1=y[,3] # OBSERVED FRAUDULANT CLAIMS

# ESTIMATE PARAMETERS ONLY OBSERVING y1=1 IF FRAUDULANT, 0 OTHERWISE
mod1=emhidden(y1,x)
thtrue # Compare estimates with true values

# CHECK IF MODEL ABLE TO SEPARATE 00 FROM 01 CORRECTLY
# This essentially looks only at the honest fraudulent, puts in values and see
#The probabilities. Keep larger than 0.5 compare wiht true. 
y01=y[,2] # 1 if miss-classified fraudulent
y01pred=mod1$yy[,2]>0.5 # 1 if estimated probability of miss-classified fraudulent (class 01) >0.5
table(y01,y01pred) # Element (2,2) of the table says how many of the real 01-obs are correctly classified

# CHECK DETECTION OF FRAUD (CLASSES 2 AND 3) -  CAUDILL ET AL METHOD
#fobs: count # in each col. Get prob from EM. Count max. table. 
fobs=apply(y,1,which.max) # true y, 1 implies row, fun=which.max. eg: teller antall i hver kolonne
table(fobs) #contingency table of counts at each combination of factor levels 
ppred=probs(mod1) #Fun that gets prob based on a,b and x. 
fpred=apply(ppred,1,which.max)# counts positition as fobs. 

#table(fobs,fpred)
#fobs ikke 1, fpred ikke 1. Hvor skjer det samtidig? TRUE/TRUE
#fobs=sanne y'er, fpred. get probs based on function (use a,b,x)
class.table=table(fobs!=1,fpred!=1) #fobs on left, fpred on top. 
class.table

#ER DETTE BARE RIKTIG PLASSERING GENERELT????

# CHECK DETECTION OF FRAUD (CLASSES 2 AND 3) - BINOMIAL LOGIT MODEL
logit=glm(y1~x,family=binomial())
plogit=predict(logit,type="response") #type = "response" gives the predicted probabilities
lpred=plogit>0.5 #True false table. fra 1:n
ltable=table(fobs!=1,lpred) #Sammenlikne med fobs igjen.
#fobs!1 betyr ikke kol 1, betyr ikke hoenst honest, betyr fraudulent
#SÃ¥ sammenlikner fraudulent med fraudulent. lpred true/false. 
#Kunne gjort dette med lpred direkte, bare ppred>0.5 i dette uttrykket. 
ltable

# CHECK DETECTION OF FRAUD (CLASSES 2 AND 3) - TRINOMIAL LOGIT MODEL
#If all observed case.
#estimer parametere, bruk parametere for Ã¥ predikere, stÃ¸rst ss betyr 
#
th=thtrue*2 #Need some starting values. 
mlogit=optim(th,logl,x=x,y=y) #logl function with three cat. 
alogit=mlogit$par[1:2]
blogit=mlogit$par[3:4]
pml=matrix(1,n,3)
for(i in 1:2) pml[,i+1]=exp(a[i]+x%*%b)
pml=pml/apply(pml,1,sum)
fml=apply(pml,1,which.max)
#table(fobs,fml)
ml.table=table(fobs!=1,fml!=1)
ml.table

#Notes, still not out of sample.
#what is most interesting? 
#Does sum mis-placed make sense. 
