# Package Demo
library(mvtnorm)
library(ScoreTest)
#library(glmnet)

#data generation
d = 100
a = 0
rho = 0.25
signal = 0.45
n = 1000
s  = 3
mean = rep(0,d)
Sgm  = diag(d) 
coi   = 4 

# Generate covariance matrix
for (i in 1:d)
{
  for (j in 1:d)
  {
    Sgm[i,j] = rho^(abs(i-j))
  }
}

X  = rmvnorm(n=n,mean=mean,sigma = Sgm) 
#X = matrix(rnorm(n*d),n)
bt = c(rep(signal,s),rep(0,d-s)) 

#gaussian
Y  = X%*%bt + rnorm(n)
#true alternative
dScore(Y, X, 1,family = 'gaussian')$pval
#true null
dScore(Y, X, 4,family = 'gaussian')$pval

#binomial
Y = sign(X%*%bt + rlogis(n)) 
#Y = (Y+1)/2
#true alternative
dScore(Y, X, 1,family = 'binomial')$pval
#true null
dScore(Y, X, 5,family = 'binomial')$pval

#poisson
Y = rpois(n,exp(X%*%bt))
#true alternative
dScore(Y, X, 1,family = 'poisson')$pval
#true null
dScore(Y, X, 4,family = 'poisson')$pval
