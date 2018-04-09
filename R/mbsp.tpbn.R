#####################################################
# Function to implement Gibbs sampler for MBSP-TPBN #
#####################################################

######################
# FUNCTION ARGUMENTS #
######################
# X = design matrix
# Y = matrix of responses
# (u,a) = parameters in the TPBN density (default is u=a=0.5 for horseshoe)
# tau = global tuning parameter
#   If not specified by user, defaults to 1/(p*sqrt(n log(n)))
# max.steps = # of times to run MCMC
# burnin = # of samples in burn-in

##################
# RETURNS A LIST #
##################
# B.est = posterior median point estimator
# lower.endpoints = 2.5th percentile of posterior density
# upper.endpoints = 97.5th percentile of posterior density
# active.predictors = row indices of the active predictors

mbsp.tpbn = function(X, Y, u=0.5, a=0.5, tau=NA, max.steps = 15000, burnin=5000) {
  
  # Burnin time should not exceed the total number of iterations.
  if (burnin > max.steps){
    stop("ERROR: Burn-in cannot be greater than # of iterations. \n")
  }
  
  # Extract dimensions n, p, and q
  n <- nrow(X)
	p <- ncol(X)
	q <- ncol(Y)
	
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
	# List to hold the draws of B
	B.samples <- rep(list(matrix(0,p,q)), max.steps)
	
	#########################
	# Initial guesses for B #
	#########################
	min.sing <- min(svd(XtX)$d)
	delta <- 0.01
	# Initial guess for B
	B <- chol2inv(chol(XtX + (delta+min.sing)*diag(p))) %*% XtY
	
	###########################
	# Initial guess for Sigma #
	###########################
	resid <- Y - X%*%B
	Sigma <- (n-1)/n * cov(resid) #  Initial guess for Sigma
	
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
  d <- 3
	resid.vec <- as.vector(resid)
  k <- var(resid.vec)
  
  ###########################
	# Start the Gibbs sampler #
  ###########################
	j <- 0
	while (j < max.steps) {
		j <- j + 1

		if (j %% 1000 == 0) {
			cat("Iteration:", j, "\n")
		}
		
		############
		# Sample B #
		############
		if (p <= n){
		  ridge <- XtX + diag(1/zeta)
		  inv.ridge <- chol2inv(chol(ridge))
		  post.M <- inv.ridge %*% XtY # Conditional posterior mean
		  # Draw from matrix normal density
		  B <- matrix.normal(post.M, inv.ridge, Sigma)
		} else if (p>n){
		  # Use the more efficient sampling algorithm if p>n
		  
		  # Draw U ~ MN(O, D, Sigma)
		  D <- diag(zeta)
		  Zero.1 <- matrix(0, nrow=p, ncol=q)
		  U <- matrix.normal(Zero.1, D, Sigma)
		  # Draw M ~ MN(O, I_n, Sigma)
		  Zero.2 <- matrix(0, nrow=n, ncol=q)
		  M <- matrix.normal(Zero.2, diag(n), Sigma)
		  # Set V = X%*%U + M
		  V <- X%*%U + M
		  # Solve for W
		  W <- chol2inv(chol(X%*%D%*%t(X)+diag(n)))%*%(Y-V)
		  # Draw B from conditional distribution and return Theta 
		  B <- U + D%*%t(X)%*%W
		}
		
		# Calculate square root of Sigma^-1
		inv.Sigma <- solve(Sigma)
		invSigma.eig <- eigen(inv.Sigma)
		invSigma.half <- invSigma.eig$vectors %*% diag(sqrt(invSigma.eig$values)) %*% t(invSigma.eig$vectors) 
		
		##############################
		# Sample zeta_i's and nu_i's #
		##############################
		for (i in 1:p){
		  norm.term <- sum((t(B[i,])%*%invSigma.half)^2) 
		  v <- max(norm.term, .Machine$double.eps) # to prevent chi parameter from collapsing to 0
		  zeta[i] <- rgig(n=1, lambda=u-q/2, chi=v, psi=2*nu[i])
	    nu[i] <- rgamma(n=1, shape=a, scale=1/(tau+zeta[i]))
		}
		
		#####################################		
		# Sample Sigma from inverse-Wishart #
		#####################################
		resid <- Y - X %*% B
		sum1 <- t(resid)%*%resid
		sum2 <- t(B)%*%((1/zeta)*diag(p))%*%B
		Sigma <- riwish(n+d+p, sum1+sum2+k*diag(q))
		
		# Save the most recent estimate of B to the list
		B.samples[[j]] <- B
	}
	
  ###################
  # Discard burn-in #
  ###################
  B.samples <- tail(B.samples,max.steps-burnin)
  
  #######################################
  # Extract the posterior mean, median, #
  # 2.5th, and 97.5th quantiles         #
  #######################################
  arr <- array(unlist(B.samples), c(p,q,length(B.samples)))
  mbsp.quantiles <- apply(arr, 1:2, function(x) quantile(x, prob=c(.025,.5,.975)))
  # Take posterior median as point estimate for B
  B.est <- mbsp.quantiles[2,,]
  # For marginal credible intervals
  lower.endpoints <- mbsp.quantiles[1,,]
  upper.endpoints <- mbsp.quantiles[3,,]
  
  ##############################
  # Perform variable selection #
  ##############################
  
  # Initialize matrix of binary entries: 0 for inactive variable b_ij, 1 for active b_ij
  classification.matrix <- matrix(rep(0,p*q),nrow=p,ncol=q)
  
  # Find the active covariates
  for(k in 1:p){
    for(l in 1:q) {
      if(lower.endpoints[k,l] < 0 && upper.endpoints[k,l] < 0)
        classification.matrix[k,l] <- 1
      else if(lower.endpoints[k,l] > 0 && upper.endpoints[k,l] > 0)
        classification.matrix[k,l] <- 1
    }
  }
  
  # Active predictors according to our estimates
  active.predictors <- which(rowSums(classification.matrix) != 0)

  # B.est = posterior median point estimator
  # lower.endpoints = 2.5th quantile of posterior density
  # upper.endpoints = 97.5th quantile of posterior density
  # active.predictors
  
  ##########################################
  # Return list of B.est, lower.endpoints, #
  # and upper.endpoints, active.predictors #
  ##########################################
  mbsp.output <- list(B.est = B.est, 
                      lower.endpoints = lower.endpoints,
                      upper.endpoints = upper.endpoints,
                      active.predictors = active.predictors)
  # Return list
	return(mbsp.output)
}
