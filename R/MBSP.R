################################################
# Function to implement Gibbs sampler for MBSP #
# with the three parameter beta-normal prior   #          
################################################


######################
# FUNCTION ARGUMENTS #
######################
# Y = nxq matrix of response variables.
# X = nxp design matrix of p covariates. We apply global-local shrinkage to these covariates
#     in order to facilitate variable selection from them.
# confounders = optional nxr matrix of confounding variables Z that should *always* be included
#               in the model. If there are confounders, then a uniform prior is placed on their
#               regression coefficients.
# u = first parameter in the TPBN density (default is u=0.5 for horseshoe)
# a = second parameter in the TPBN density (default is a=0.5 for horseshoe) 
# tau = global tuning parameter. In general, it is recommended that tau be between 0 and 1, 
#       but the user may specify any tau greater than 0. A smaller tau gives greater shrinkage.
#       If not specified by user, defaults to 1/(p*sqrt(n log(n)))
# max_steps = number of total iterations to run MCMC
# burnin = # of samples in burn-in
# save_samples = Boolean variable for whether or not to return all posterior samples
# model_criteria = Boolean variable for whether or not to return DIC and WAIC

##################
# RETURNS A LIST #
##################
# B_est = pxq posterior median point estimator for covariates' regression coefficients B.

# B_CI_lower = 2.5th percentile of posterior for the B matrix
# B_CI_upper = 97.5th percentile of posterior for the B matrix
# active_predictors = row indices of the active predictors
# B_samples = all posterior samples of B
# C_est = rxq posterior median point estimator for confounders' regression coefficients C.
#                   This is not returned if there are no confounders.
# C_CI_lower = 2.5th percentile of posterior for the C matrix. This is not returned if
#              there are no confounders.
# C_CI_upper = 97.5th percentile of the posterior for the C matrix. This is not returned
#              if there are no confounders.
# C_samples = all posterior samples of C. This is not returned if there are no confounders.
# Sigma_est = posterior median point estimator for covariance matrix Sigma
# Sigma_CI_lower = 2.5th percentile of the posterior for the covariance matrix Sigma
# Sigma_CI_upper = 97.th percentile of the posterior for the covariance matrix Sigma
# DIC = Deviance Information Criterion (DIC) for the MBSP model. Can be used for model selection.
#       If model_criteria==FALSE, then this is returned as NA.
# WAIC =  Widely Applicable Information Criterion (WAIC) for the MBSP model. Can be used for model selection.
#         If model_criteria==FALSE, then this is returned as NA.


MBSP = function(Y, X, confounders=NULL, u=0.5, a=0.5, tau=NA, 
                max_steps = 6000, burnin=1000, save_samples=TRUE,
                model_criteria=FALSE) {
  
  # Burnin time should not exceed the total number of iterations.
  if (burnin > max_steps){
    stop("ERROR: Burn-in cannot be greater than # of iterations. \n")
  }
  
  if(nrow(Y) != nrow(Y))
    stop("Y and X must have the same number of rows.")
  
  # Extract dimensions n, p, and q
  n <- nrow(X)
	p <- ncol(X)
	q <- ncol(Y)
	
	# Center X if not already done
	X <- scale(X, center=TRUE, scale=FALSE)
	# if tau = NA, set it equal to 1/(p*sqrt(n*log(n)))
	if(is.na(tau)){
	  tau <- 1/(p*sqrt(n*log(n)))	  
	}
	# if user specified tau, must be strictly positive
	if (tau<0){
	  stop("ERROR: tau should be strictly positive. \n")
	} 
	# If there are confounders, check that number of confounders is not
	# greater than sample size
	if(!is.null(confounders)){
	  if(nrow(confounders) != n)
	    stop("The data matrix for the confounders should have the same number of rows as Y and X.")
	  
	  r <- dim(confounders)[2]
	  if(r > n){
	    stop("The number of confounders should not exceed sample size.")
	  }
	  # Center the confounders matrix as well
	  Z <- scale(confounders, center=TRUE, scale=FALSE)
	 }
	
	# Time saving
	XtX <- t(X) %*% X	
	XtY <- t(X) %*% Y
	# List to hold the draws of B
	B_samples <- rep(list(matrix(0,p,q)), max_steps)
	# If there are confounders, also hold the draws of C
	if(!is.null(confounders)){
	  ZtZ_inv <- chol2inv(chol(t(Z) %*% Z))
	  C_samples <- rep(list(matrix(0,r,q)), max_steps)
	}
	# List to hold the draws of Sigma
	Sigma_samples <- rep(list(matrix(0,q,q)), max_steps)
	if(model_criteria==TRUE){
	  # List to hold log-likelihood. Needed to compute DIC and WAIC
	  loglik_samples <- matrix(0, n, max_steps)
	}
	
	#########################
	# Initial guesses for B #
	#########################
	if(is.null(confounders)){
	  lambda <- 0.01
	  # Ridge estimator
	  if(p <= n){
	    B <- chol2inv(chol(XtX + lambda*diag(p))) %*% XtY
	  } else if(p > n) { 
	    # Using Woodbury matrix identity
	    term1 <- (1/lambda)*XtY
	    term2 <- (1/lambda^2)*t(X) %*% chol2inv(chol((1/lambda)*X%*%t(X) + diag(n))) %*% X %*% XtY
	    B <- term1 - term2
	  }
	}
	
	# If there are confounders, need an initial guess for C too
	if(!is.null(confounders)){
	  XZ <- cbind(X,Z)
	  XZtY <- t(XZ) %*% Y
	  XZ_cp <- t(XZ) %*% XZ
	  lambda <- 0.01
	  
	  if((p+r) <= n){
	    BC <- chol2inv(chol(XZ_cp + lambda*diag(p+r))) %*% XZtY
	  } else if((p+r) > n) { 
	    # Using Woodbury matrix identity
	    term1 <- (1/lambda)*XZtY
	    term2 <- (1/lambda^2)*t(XZ) %*% chol2inv(chol((1/lambda)*XZ%*%t(XZ) + diag(n))) %*% XZ %*% XZtY
	    BC <- term1 - term2
	  }
	  # Initial guesses for B and C
	  B <- BC[1:p,]
	  C <- BC[(p+1):(p+r),]
	  # Free memory
	  rm(BC)
	}
	
	###########################
	# Initial guess for Sigma #
	###########################
  if(is.null(confounders)){
	  resid <- Y - X%*%B
	  Sigma <- (n-1)/n * stats::cov(resid) #  Initial guess for Sigma
  } else {
    resid <- Y - X%*%B - Z%*%C
    Sigma <- (n-1)/n * stats::cov(resid)
  }
	#####################
	# Initial guess for #
	# xi_1, ..., xi_p,  # 
	# nu_1, ..., nu_p   #
	#####################
	zeta <- rep(a*tau,p) # Initial guesses for zeta_i's
	nu <- u*zeta # Initial guesses for nu_i's
	
	###########################################
  # Set hyperparameters for inverse Wishart #
	###########################################
  d <- q+2
	resid_vec <- as.vector(resid)
  k <- stats::var(resid_vec) 
  
  ###########################
	# Start the Gibbs sampler #
  ###########################
	j <- 0
	while (j < max_steps) {
		j <- j + 1

		if (j %% 1000 == 0) {
			cat("Iteration:", j, "\n")
		}
		
		############
		# Sample B #
		############
		if (p <= n){
		  ridge <- XtX + diag(1/zeta)
		  inv_ridge <- chol2inv(chol(ridge))
		  
		  # Conditional posterior mean
		  if(is.null(confounders)){
	      post_M <- inv_ridge %*% XtY 
		  } else {
		    post_M <- inv_ridge %*% t(X) %*% (Y-Z%*%C)    
		  }
		    
		  # Draw from matrix normal density
		  B <- matrix_normal(post_M, inv_ridge, Sigma)
		
		} else if (p>n){
		  # Use the more efficient sampling algorithm if p>n
		  
		  # Draw U ~ MN(O, D, Sigma)
		  U <- matrix(0, nrow=p, ncol=q)
		  for(i in 1:p){
		    U[i, ] <- mvtnorm::rmvnorm(1, rep(0,q), zeta[i]*Sigma)
		  }
		  # Draw M ~ MN(O, I_n, Sigma)
		  M <- matrix(0, nrow=n, ncol=q)
		  for(i in 1:n){
		    M[i, ] <- mvtnorm::rmvnorm(1, rep(0,q), Sigma)
		  }
		  # Set V = X%*%U + M
		  V <- X%*%U + M
		  # Solve for W
		  if(is.null(confounders)){
		    W <- chol2inv(chol(X%*%(zeta*t(X))+diag(n)))%*%(Y-V)
		  } else {
		    W <- chol2inv(chol(X%*%(zeta*t(X))+diag(n)))%*%(Y-Z%*%C-V)
		  }
		  # Draw B from conditional distribution and return Theta 
		  B <- U + (zeta*t(X))%*%W
		}
		
		#############################
		# If there are confounders, #
		# then sample from C        #
		#############################
		if(!is.null(confounders)){
		  post_C <- ZtZ_inv %*% t(Z) %*% (Y-X%*%B)
		  C <- matrix_normal(post_C, ZtZ_inv, Sigma)
		}
		
		# Calculate square root of Sigma^-1
		inv_Sigma <- solve(Sigma)
		invSigma_eig <- eigen(inv_Sigma)
		invSigma_half <- invSigma_eig$vectors %*% (sqrt(invSigma_eig$values)*t(invSigma_eig$vectors)) 
		
		##############################
		# Sample zeta_i's and nu_i's #
		##############################
		for (i in 1:p){
		  norm_term <- sum((t(B[i,])%*%invSigma_half)^2) 
		  v <- max(norm_term, .Machine$double.eps) # to prevent chi parameter from collapsing to 0
		  zeta[i] <- GIGrvg::rgig(n=1, lambda=u-q/2, chi=v, psi=2*nu[i])
	    nu[i] <- stats::rgamma(n=1, shape=a, scale=1/(tau+zeta[i]))
		}
		
		#####################################		
		# Sample Sigma from inverse-Wishart #
		#####################################
		if(is.null(confounders)){
		  resid <- Y - X %*% B
		} else{
		  resid <- Y - X%*% B - Z %*% C  
		}
		sum1 <- t(resid)%*%resid
		sum2 <- t(B)%*%(B/zeta)
		Sigma <- MCMCpack::riwish(n+d+p, sum1+sum2+k*diag(q))
		
		# Save the most recent estimate of B and Sigma
		B_samples[[j]] <- B
		Sigma_samples[[j]] <- Sigma
		
		# If there are confounders, also save recent estimate of C to the list
		if(!is.null(confounders)){
		  C_samples[[j]] <- C
		}
		
		if(model_criteria==TRUE){
		  # Compute log-likelihood
		  if(is.null(confounders)){
		    # If no confounders
		    mu = X%*%B
		    for(i in 1:n){
		      loglik_samples[i,j] <- mvtnorm::dmvnorm(Y[i,], mean=mu[i,], sigma=Sigma, log=TRUE)
		    }
		  } else {
		    # If there are confounders
		    mu = X%*%B + Z%*%C
		    for(i in 1:n){
		      loglik_samples[i,j] <- mvtnorm::dmvnorm(Y[i,], mean=mu[i,], sigma=Sigma, log=TRUE)
		    }
		  }
		}
	# End Gibbs sampler
	}
	
  ###################
  # Discard burn-in #
  ###################
  B_samples <- utils::tail(B_samples, max_steps-burnin)
  if(!is.null(confounders)){
    C_samples <- utils::tail(C_samples, max_steps-burnin)
  }
  Sigma_samples <- utils::tail(Sigma_samples, max_steps-burnin)
  if(model_criteria==TRUE){
    loglik_samples <- loglik_samples[ , (burnin+1):max_steps]
  }
  
  #################################
  # Extract the posterior median, #
  # 2.5th, and 97.5th quantiles   #
  #################################
  B_arr <- array(unlist(B_samples), c(p,q,length(B_samples)))
  B_quantiles <- apply(B_arr, 1:2, function(x) stats::quantile(x, prob=c(.025,.5,.975)))
  # Take posterior median as point estimate for B
  B_est <- B_quantiles[2,,]
  # For marginal credible intervals
  B_CI_lower <- B_quantiles[1,,]
  B_CI_upper <- B_quantiles[3,,]
  
  if(!is.null(confounders)){
    C_arr <- array(unlist(C_samples), c(r,q,length(C_samples)))
    C_quantiles <- apply(C_arr, 1:2, function(x) stats::quantile(x, prob=c(.025,.5,.975)))
    # Take posterior median as point estimate for B
    C_est <- C_quantiles[2,,]
    # For marginal credible intervals
    C_CI_lower <- C_quantiles[1,,]
    C_CI_upper <- C_quantiles[3,,]
  }
  
  Sigma_arr <- array(unlist(Sigma_samples), c(q,q,length(Sigma_samples)))
  Sigma_quantiles <- apply(Sigma_arr, 1:2, function(x) stats::quantile(x, prob=c(.025,.5,.975)))
  # Take posterior median as point estimate for B
  Sigma_est <- Sigma_quantiles[2,,]
  # For marginal credible intervals
  Sigma_CI_lower <- Sigma_quantiles[1,,]
  Sigma_CI_upper <- Sigma_quantiles[3,,]
  
  ################################################
  # Compute DIC and WAIC if model_criteria==TRUE #
  ################################################
  DIC = WAIC = NA
  if(model_criteria==TRUE){
    # Obtain the posterior means
    B_hat <- apply(B_arr, 1:2, mean)
    if(!is.null(confounders)){
      C_hat <- apply(C_arr, 1:2, mean)
    }
    Sigma_hat <- apply(Sigma_arr, 1:2, mean)
  
    # Calculate log-likehood with B_hat, C_hat, Sigma_hat
    ll_theta_hat <- 0
    if(is.null(confounders)){
      mu = X%*%B_hat
      for(i in 1:n){
        ll_theta_hat <- ll_theta_hat + mvtnorm::dmvnorm(Y[i,], mean=mu[i,], sigma=Sigma_hat, log=TRUE)
      }
    } else {
      # If there are confounders
      mu = X%*%B_hat+Z%*%C_hat
      for(i in 1:n){
        ll_theta_hat <- ll_theta_hat + mvtnorm::dmvnorm(Y[i,], mean=mu[i,], sigma=Sigma_hat, log=TRUE)
      }
    }
  
    # Compute DIC
    DIC <- -4*mean(colSums(loglik_samples)) + 2*ll_theta_hat 

    # Compute WAIC
    var_ll <- apply(loglik_samples, 1, stats::var)
    lik_samples <- exp(loglik_samples)
    mean_lik <- apply(lik_samples, 1, mean)
    WAIC <- 2*(sum(var_ll)) - 2*sum(log(mean_lik))
  }
  
  ##############################
  # Perform variable selection #
  ##############################
  
  # Initialize matrix of binary entries: 0 for inactive variable b_ij, 1 for active b_ij
  classification_matrix <- matrix(rep(0,p*q),nrow=p,ncol=q)
  
  # Find the active covariates
  for(t in 1:p){
    for(l in 1:q) {
      if(B_CI_lower[t,l] < 0 && B_CI_upper[t,l] < 0)
        classification_matrix[t,l] <- 1
      else if(B_CI_lower[t,l] > 0 && B_CI_upper[t,l] > 0)
        classification_matrix[t,l] <- 1
    }
  }
  
  # Active predictors according to our estimates
  active_predictors <- which(rowSums(classification_matrix) != 0)
  
  ################################
  # Return:                      #
  # B_est, CI_lower, CI_upper,   #
  # active_predictors, B_samples #   
  ################################
  
  if(save_samples){
    if(is.null(confounders)){
      mbsp_output <- list(B_est = B_est, 
                          B_CI_lower = B_CI_lower,
                          B_CI_upper = B_CI_upper,
                          active_predictors = active_predictors,
                          B_samples = B_samples,
                          Sigma_est = Sigma_est, 
                          Sigma_CI_lower = Sigma_CI_lower,
                          Sigma_CI_upper = Sigma_CI_upper,
                          Sigma_samples = Sigma_samples,
                          DIC = DIC,
                          WAIC = WAIC)
    } else {
      mbsp_output <- list(B_est = B_est,
                          B_CI_lower = B_CI_lower,
                          B_CI_upper = B_CI_upper,
                          active_predictors = active_predictors,
                          B_samples = B_samples,
                          C_est = C_est,
                          C_CI_lower = C_CI_lower,
                          C_CI_upper = C_CI_upper,
                          C_samples = C_samples,
                          Sigma_est = Sigma_est, 
                          Sigma_CI_lower = Sigma_CI_lower,
                          Sigma_CI_upper = Sigma_CI_upper,
                          Sigma_samples = Sigma_samples,
                          DIC = DIC,
                          WAIC = WAIC)
    }
  } else {
    if(is.null(confounders)){
      mbsp_output <- list(B_est = B_est, 
                          B_CI_lower = B_CI_lower,
                          B_CI_upper = B_CI_upper,
                          active_predictors = active_predictors,
                          Sigma_est = Sigma_est, 
                          Sigma_CI_lower = Sigma_CI_lower,
                          Sigma_CI_upper = Sigma_CI_upper,
                          DIC = DIC,
                          WAIC = WAIC)
    } else {
      mbsp_output <- list(B_est = B_est,
                          B_CI_lower = B_CI_lower,
                          B_CI_upper = B_CI_upper,
                          active_predictors = active_predictors,
                          C_est = C_est,
                          C_CI_lower = C_CI_lower,
                          C_CI_upper = C_CI_upper,
                          Sigma_est = Sigma_est, 
                          Sigma_CI_lower = Sigma_CI_lower,
                          Sigma_CI_upper = Sigma_CI_upper,
                          DIC = DIC,
                          WAIC = WAIC)    
    }
  }
  
  # Return list
	return(mbsp_output)
}