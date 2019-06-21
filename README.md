# ScoreTest
This is an R package containing the implementation of the score test function from the paper ''A General Theory of Hypothesis Tests and Confidence Regions for Sparse High Dimensional Models''
by Yang Ning and Han Liu. 

To install this package, please follow the following steps:

1. Install package devtools from CRAN. 
```{R}
install.packages("devtools")
```

2. Load the devtools package.

```{R}
library(devtools)
```

3. Install ScoreTest.

```{R}
install_github("Maosbrother/ScoreTest")
```

To test the package, feel free to play with the following toy example (more examples in Demo.R):

```{R}
#data generation
d  <-  200
a  <-  0
rho  <-  0.25
signal  <-  0.45
n  <-  100
s   <-  3
mean  <-  rep(0,d)
Sgm   <-  diag(d) 
coi    <-  4 

# Generate covariance matrix
for (i in 1:d)
{
  for (j in 1:d)
  {
    Sgm[i,j]  <-  rho^(abs(i-j))
  }
}

X   <-  rmvnorm(n=n,mean=mean,sigma = Sgm) 
#X = matrix(rnorm(n*d),n)
bt  <-  c(rep(signal,s),rep(0,d-s)) 

#gaussian
Y   <-  X%*%bt + rnorm(n)
#true alternative
dScore(Y, X, 1,family = 'gaussian')$pval
#true null
dScore(Y, X, 4,family = 'gaussian',refit = T)$pval
```
