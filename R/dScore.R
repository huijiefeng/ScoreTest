#' Decorrelated Score Test.
#' 
#' @param Y The response vector.
#' @param X The data matrix.
#' @param coi Coefficient of interest. It should be a single integer.
#' @param beta The fitted value specified by the user. The default is NULL and the function will call glmnet package for Lasso autometically.
#' @param maxit Maximal number of iterations. The default is 100,000.
#' @param nfolds Number of folds for cross validation. The default is 5.
#' @param lambda Tuning parameter for fitting coefficients. The default is NULL that uses default grid from glmnet.
#' @param family 'gaussian', 'binomial' or 'poisson' for linear, logistic and poisson regression, respectively.
#' @param refit Logical parameter. If refit is TRUE then for coefficient estimates, the function will first use glmnet for variable selection, and then redo the regression using only the selected variables. The default is FALSE.
#' @return A list containing:
#' \itemize{ 
#' \item fitted_value - The fitted coefficients using LASSO.
#' \item test_stats - The Score test statistic.
#' \item pval - P-value.
#' }
#' @import glmnet
#' @export
#' @examples
#' X <- matrix(rnorm(10000),1000)
#' bt <- c(0.5,0.5,0.5,0)
#' Y <- X%*%bt + rnorm(n)
#' dScore(Y, X, 1,family = 'gaussian')$pval

dScore <- function(Y, X, coi, beta = NULL,maxit = 100000, nfolds = 5, lambda = NULL,family = 'gaussian',
                   refit = FALSE){
  
  # Lots of error messages
  if(!is.numeric(Y)){
    stop('Y must be a numeric vector')
  }
  
  if(!is.numeric(X)){
    stop('X must be a numeric vector')
  }
  
  if(!is.numeric(coi) || (coi%%1!=0)){
    stop('Coeffcient of interest must be integer')
  }
  
  Y <-  as.matrix(Y)
  
  if(ncol(Y) != 1){
    stop('Y must be one dimensional')
  }
  
  if(nrow(Y) != nrow(X)){
    stop('Dimensions of X and Y do not match')
  }
  
  d = ncol(X) # Get dimensions
  n = nrow(X)
  
  if(ncol(X) < coi){
    stop('Coefficient of interest exceed dimensions')
  }
  
  if (family == 'binomial'){
    if (length(unique(Y)) != 2){
      stop('Response Y does not have 2 and only 2 levels')
    }
    temp_lev = unique(Y)
    idx_1 <-  Y == temp_lev[1]
    idx_2 <-  Y == temp_lev[2]
    Y[idx_1] = 0
    Y[idx_2] = 1
  }
  
  if (is.null(beta) == FALSE){
    if(length(beta) != ncol(X)){
      stop('Dimensions of coefficient vector and data matrix do not match')
    }
  } else{
      cv.fit <- cv.glmnet(X,Y,maxit = maxit,nfolds = nfolds, family = family, lambda = lambda)
      tmp <- which(cv.fit$glmnet.fit$lambda == cv.fit$lambda.min)
      beta <- cv.fit$glmnet.fit$beta[,tmp]
      if (refit == TRUE){
        idx_nonzero = which(beta != 0)
        beta_nonzero <- glm(Y~X[,idx_nonzero]-1,family = family)$coefficients
        beta[idx_nonzero] <- beta_nonzero
      }
    } 
  
  
  test_func <- get(paste('.',family,'_helper',sep = ''))
  
  U <- test_func(beta,Y,X,coi,nfolds)
  pval <- 2*pnorm(-abs(U))
  
  results <- list(fitted_value = beta,
                  test_stats = U,
                  pval = as.vector(pval))
  
  return(results)
}