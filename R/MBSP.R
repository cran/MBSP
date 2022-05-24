################################################
# Function to implement Gibbs sampler for MBSP #
# with the three parameter beta-normal prior   #          
################################################


######################
# FUNCTION ARGUMENTS #
######################
# X = design matrix
# Y = matrix of responses
# u = first parameter in the TPBN density (default is u=0.5 for horseshoe)
# a = second parameter in the TPBN density (default is a=0.5 for horseshoe) 
# tau = global tuning parameter
#       If not specified by user, defaults to 1/(p*sqrt(n log(n)))
# max_steps = number of total iterations to run MCMC
# burnin = # of samples in burn-in
# save_samples = Boolean variable for whether or not to return all posterior samples

##################
# RETURNS A LIST #
##################
# B_est = posterior median point estimator
# CI_lower = 2.5th percentile of posterior density
# CI_upper = 97.5th percentile of posterior density
# active_predictors = row indices of the active predictors
# B_samples = all posterior samples of B

MBSP = function(X, Y, u=0.5, a=0.5, tau=NA, max_steps = 6000, burnin=1000,
                save_samples=TRUE) {
  
  # Burnin time should not exceed the total number of iterations.
  if (burnin > max_steps){
    stop("ERROR: Burn-in cannot be greater than # of iterations. \n")
  }
  
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
	# if user specified tau, must be between 0 and 1.
	if ( (tau<0) | (tau>1) ){
	  stop("ERROR: tau should be strictly between 0 and 1. \n")
	} 
	
	# Time saving
	XtX <- t(X) %*% X	
	XtY <- t(X) %*% Y
	# List to hold the draws of B and Sigma
	B_samples <- rep(list(matrix(0,p,q)), max_steps)

	#########################
	# Initial guesses for B #
	#########################
	min_sing <- min(svd(XtX)$d)
	delta <- 0.01
	lambda <- min_sing+delta
	
	# Initial guess for B
	if(p <= n){
	  B <- chol2inv(chol(XtX + lambda*diag(p))) %*% XtY
	  
	} else if(p > n) { 
	  # Using Woodbury matrix identity
	  term1 <- (1/lambda)*XtY
	  term2 <- (1/lambda^2)*t(X) %*% chol2inv(chol((1/lambda)*X%*%t(X) + diag(n))) %*% X %*% XtY
	  B <- term1 - term2
	}
	
	###########################
	# Initial guess for Sigma #
	###########################
	resid <- Y - X%*%B
	Sigma <- (n-1)/n * stats::cov(resid) #  Initial guess for Sigma
	
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
		  post_M <- inv_ridge %*% XtY # Conditional posterior mean
		  # Draw from matrix normal density
		  B <- matrix_normal(post_M, inv_ridge, Sigma)
		} else if (p>n){
		  # Use the more efficient sampling algorithm if p>n
		  
		  # Draw U ~ MN(O, D, Sigma)
		  D <- diag(zeta)
		  Zero_1 <- matrix(0, nrow=p, ncol=q)
		  U <- matrix_normal(Zero_1, D, Sigma)
		  # Draw M ~ MN(O, I_n, Sigma)
		  Zero_2 <- matrix(0, nrow=n, ncol=q)
		  M <- matrix_normal(Zero_2, diag(n), Sigma)
		  # Set V = X%*%U + M
		  V <- X%*%U + M
		  # Solve for W
		  W <- chol2inv(chol(X%*%(zeta*t(X))+diag(n)))%*%(Y-V)
		  # Draw B from conditional distribution and return Theta 
		  B <- U + D%*%t(X)%*%W
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
		resid <- Y - X %*% B
		sum1 <- t(resid)%*%resid
		sum2 <- t(B)%*%(B/zeta)
		Sigma <- MCMCpack::riwish(n+d+p, sum1+sum2+k*diag(q))
		
		# Save the most recent estimate of B and Sigma to the list
		B_samples[[j]] <- B
	}
	
  ###################
  # Discard burn-in #
  ###################
  B_samples <- utils::tail(B_samples, max_steps-burnin)
  
  #################################
  # Extract the posterior median, #
  # 2.5th, and 97.5th quantiles   #
  #################################
  arr <- array(unlist(B_samples), c(p,q,length(B_samples)))
  mbsp_quantiles <- apply(arr, 1:2, function(x) stats::quantile(x, prob=c(.025,.5,.975)))
  # Take posterior median as point estimate for B
  B_est <- mbsp_quantiles[2,,]
  # For marginal credible intervals
  CI_lower <- mbsp_quantiles[1,,]
  CI_upper <- mbsp_quantiles[3,,]
  
  ##############################
  # Perform variable selection #
  ##############################
  
  # Initialize matrix of binary entries: 0 for inactive variable b_ij, 1 for active b_ij
  classification_matrix <- matrix(rep(0,p*q),nrow=p,ncol=q)
  
  # Find the active covariates
  for(k in 1:p){
    for(l in 1:q) {
      if(CI_lower[k,l] < 0 && CI_upper[k,l] < 0)
        classification_matrix[k,l] <- 1
      else if(CI_lower[k,l] > 0 && CI_upper[k,l] > 0)
        classification_matrix[k,l] <- 1
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
    mbsp_output <- list(B_est = B_est, 
                        CI_lower = CI_lower,
                        CI_upper = CI_upper,
                        active_predictors = active_predictors,
                        B_samples = B_samples)
  } else {
    mbsp_output <- list(B_est = B_est, 
                        CI_lower = CI_lower,
                        CI_upper = CI_upper,
                        active_predictors = active_predictors)
  }
  
  # Return list
	return(mbsp_output)
}