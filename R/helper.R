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
  
  cv.fitw <- cv.glmnet(lg,la,nfolds = nfolds)
  fitw <- cv.fitw$glmnet.fit
  tmp <- which(fitw$lambda == cv.fitw$lambda.min)
  what <- fitw$beta[,tmp]
  
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
  #fit w_hat
  cv.fitw <- cv.glmnet(lg_w,la_w,nfolds = nfolds)
  tmp <- which(cv.fitw$glmnet.fit$lambda == cv.fitw$lambda.min)
  what <- cv.fitw$glmnet.fit$beta[,tmp]
  
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
  #fit w_hat
  cv.fitw <- cv.glmnet(lg_w,la_w,nfolds = nfolds)
  tmp <- which(cv.fitw$glmnet.fit$lambda == cv.fitw$lambda.min)
  what <- cv.fitw$glmnet.fit$beta[,tmp]
  
  sig2 <- mean(la_w * (la_w - lg_w%*%what))
  
  weights <- exp(lg%*%beta[-coi])
  
  S <- -mean((Y- weights)*(la - lg%*%what) )
  
  U <- sqrt(length(Y))*S/sqrt(sig2)
  return(U)
}
