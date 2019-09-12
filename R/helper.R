#' Lasso fit wrapper

.fit_lasso_wrapper <- function(X,Y,maxit = 10000,nfolds = 5,family = 'gaussian',lambda = NULL,refit = FALSE,refit_ratio = 1){
  cv.fit <- cv.glmnet(X,Y,maxit = maxit,nfolds = nfolds, family = family, lambda = lambda)
  tmp <- which(cv.fit$glmnet.fit$lambda == cv.fit$lambda.min)
  beta <- cv.fit$glmnet.fit$beta[,tmp]
  if (sum(beta!=0) == 0){
    tmp = which(cv.fit$nzero >= 1)[1]
    beta <- cv.fit$glmnet.fit$beta[,tmp]
  }
  if (refit == TRUE){
    idx_nonzero = which(beta != 0)
    if(length(idx_nonzero) > floor(d*refit_ratio)){
      tmp = which(cv.fit$nzero >= floor(d*refit_ratio))[1]
      beta <- cv.fit$glmnet.fit$beta[,tmp]
      idx_nonzero = which(beta != 0)
    }
    
    beta_nonzero <- glm(Y~X[,idx_nonzero]-1,family = family)$coefficients
    beta[idx_nonzero] <- beta_nonzero
  }
  return(beta)
}


#' Gaussian helper function
#' @param beta Fitted value.
#' @param Y Response vector.
#' @param X Data matrix.
#' @param coi Coefficient of interest.
#' @param nfolds Number of folds for cross-validation.
#' @keywords internal
#' 
.gaussian_helper <- function(beta,Y,X,coi,nfolds){
  ahat <- beta[coi] 
  ghat <- beta[-coi]  
  Hs <- t(X)%*%X 

  la <- X[,coi] # Z_i in paper. Where Q_i = (Z_i,X_i) is all covariate data
  lg <- X[, -coi] # X_i in paper (note).
  
  sig2 = mean((Y - X%*%beta)^2)
  
  what <- .fit_lasso_wrapper(lg,la,refit = F)
  
  res.n <- as.vector(Y-lg%*%ghat)
  S <- mean(res.n*la) - t(what)%*%(colMeans(res.n*lg))
  
  var <- (Hs[coi,coi] - t(what)%*%Hs[-coi,coi])/n
  
  U <- sqrt(n)*S/sqrt(var*sig2)
  return(U)
}

#' Logistic/binomial helper function
#' @param beta Fitted value.
#' @param Y Response vector.
#' @param X Data matrix.
#' @param coi Coefficient of interest.
#' @param nfolds Number of folds for cross-validation.
#' @keywords internal
#' 
.binomial_helper <-  function(beta,Y,X,coi,nfolds){
  la <- X[,coi] # Z_i in paper. Where Q_i = (Z_i,X_i) is all covariate data
  lg <- X[,-coi] # X_i in paper (note).
  weights <- exp(X%*%beta)
  weights <- sqrt(weights)/(1+weights)
  la_w <- la*weights
  lg_w <- sapply(1:ncol(lg), function(i){
    weights*lg[,i]
  })
  
  what <- .fit_lasso_wrapper(lg_w,la_w,refit = F)
  
  sig2 <- mean(la_w * (la_w - lg_w%*%what))
  
  weights <- exp(lg%*%beta[-coi])
  weights <- weights/(1+weights)
  S <- -mean((Y- weights)*(la - lg%*%what) )
  
  U <- sqrt(length(Y))*S/sqrt(sig2)
  return(U)
}

#' Poisson helper function
#' @param beta Fitted value.
#' @param Y Response vector.
#' @param X Data matrix.
#' @param coi Coefficient of interest.
#' @param nfolds Number of folds for cross-validation.
#' @keywords internal
#' 
.poisson_helper <- function(beta,Y,X,coi,nfolds){
  la <- X[,coi] # Z_i in paper. Where Q_i = (Z_i,X_i) is all covariate data
  lg <- X[, -coi] # X_i in paper (note).
  weights <- sqrt(exp(X%*%beta))
  la_w <- la*weights
  lg_w <- sapply(1:ncol(lg), function(i){
    weights*lg[,i]
  })
  
  what <- .fit_lasso_wrapper(lg_w,la_w,refit = F)
  
  sig2 <- mean(la_w * (la_w - lg_w%*%what))
  
  weights <- exp(lg%*%beta[-coi])
  
  S <- -mean((Y- weights)*(la - lg%*%what) )
  
  U <- sqrt(length(Y))*S/sqrt(sig2)
  return(U)
}
