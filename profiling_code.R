library(microbenchmark)
library(rbenchmark)
library(MASS)


mu <- rep(0,2)
Sigma <- matrix(0.0, nrow=2, ncol=2) + diag(2)*1

set.seed(123)
var1 <- mvrnorm(n=1000, mu=mu, Sigma=Sigma)
set.seed(123)
var2 <- cbind(rnorm(n=1000, mean = 0, sd = 1), rnorm(n=1000, mean = 0, sd = 1))

summary(var1)
summary(var2)

mbm <- microbenchmark(
  "newway" = {
    var1 <- mvrnorm(n=1000, mu=mu, Sigma=Sigma)
  },
  "oldway" = {
    var2 <- cbind(rnorm(n=1000, mean = 0, sd = 1), rnorm(n=1000, mean = 0, sd = 1))
    
  }
)

mbm

#data gen
nrepl=1000
bootrepl=1000
n=1000
px=2
a=c(-1.8,-1.5)
b=c(-0.02,0.2)
thtrue=c(a,b)
xsd=c(12.3, 1.8)
xcor=matrix(c(1,0.5,0.5,1), nrow = px, ncol = px)
covarmat = xcor * tcrossprod(xsd)

mbm <- benchmark(
  "newway" = {
    var1 = mvrnorm(n=1000, mu=c(rep(0,px)), Sigma=covarmat)
  },
  "oldway" = {
    var2 = matrix(nrow=n, ncol=px)
    var2 = mvrnorm(n=1000, mu=c(rep(0,px)), Sigma=covarmat)
    
  }
)

mbm


##### Naming objects

