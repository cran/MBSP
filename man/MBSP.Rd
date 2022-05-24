\name{MBSP}
\alias{MBSP}

\title{MBSP Model with Three Parameter Beta Normal (TPBN) Family}

\description{This function provides a fully Bayesian approach for obtaining a (nearly) sparse estimate of the \eqn{p \times q} coefficients matrix \eqn{B} in the multivariate linear regression model, 

\deqn{Y = XB+E,} 

using the three parameter beta normal (TPBN) family. The complete model is described in Bai and Ghosh (2018).
}

\usage{
MBSP(X, Y, u=0.5, a=0.5, tau=NA, 
     max_steps=6000, burnin=1000, save_samples=TRUE)}

\arguments{
  \item{X}{design matrix of \eqn{n} samples and \eqn{p} covariates.}

  \item{Y}{response matrix of \eqn{n} samples and \eqn{q} response variables.}

  \item{u}{The first parameter in the TPBN family. Defaults to \eqn{u=0.5} for the horseshoe prior.}

  \item{a}{The second parameter in the TPBN family. Defaults to \eqn{a=0.5} for the horseshoe prior.}
  
  \item{tau}{The global parameter. If the user does not specify this (tau=NA), the Gibbs sampler will use \eqn{\tau = 1/(p*n*log(n))}. The user may also specify a value for \eqn{\tau} between \eqn{0} and \eqn{1}, otherwise it defaults to \eqn{1/(p*n*log(n))}. }

\item{max_steps}{The total number of iterations to run in the Gibbs sampler. Defaults to \code{6000}.}

  \item{burnin}{The number of burn-in iterations for the Gibbs sampler. Defaults to \code{1000}.}
  
\item{save_samples}{A Boolean variable for whether to save all of the posterior samples of the regression coefficients matrix B. Defaults to \code{"TRUE"}.}

}

\details{
The function performs (nearly) sparse estimation of \eqn{B} and variable selection from the \eqn{p} covariates. The lower and upper endpoints of the 95 percent posterior credible intervals for each of the \eqn{pq} elements of \eqn{B} are also returned so that the user may assess uncertainty quantification.

In the three parameter beta normal (TPBN) family, \eqn{(u,a)=(0.5,0.5)} corresponds to the horseshoe prior, \eqn{(u,a)=(1,0.5)} corresponds to the Strawderman-Berger prior, and \eqn{(u,a)=(1,a), a > 0} corresponds to the normal-exponential-gamma (NEG) prior. This function uses the horseshoe prior as the default shinkrage prior.

}

\value{The function returns a list containing the following components:
\item{B_est}{ the point estimate of the \eqn{p \times q} matrix \eqn{B} (taken as the componentwise posterior median for all \eqn{pq} entries).}
\item{CI_lower}{ The 2.5th percentile of the posterior density (or the lower endpoint of the 95 percent credible interval) for all \eqn{pq} entries of \eqn{B}.}
\item{CI_upper}{ The 97.5th percentile of the posterior density (or the upper endpoint of the 95 percent credible interval) for all \eqn{pq} entries of \eqn{B}. }
\item{active_predictors}{ The row indices of the active (nonzero) covariates chosen by our model from the \eqn{p} total predictors.}
\item{B_samples}{ All \code{max_steps-burnin} samples.}
}

\references{
Armagan, A., Clyde, M., and Dunson, D.B. (2011) Generalized Beta Mixtures of Gaussians. In J. Shawe-Taylor, R. Zemel, P. Bartlett, F. Pereira, and K. Weinberger (Eds.) \emph{Advances in Neural Information Processing Systems 24}, 523-531.

Bai, R. and Ghosh, M. (2018). High-Dimensional Multivariate Posterior Consistency Under Global-Local Shrinkage Priors. \emph{Journal of Multivariate Analysis}, \bold{167}: 157-170.

Berger, J. (1980). A Robust Generalized Bayes Estimator and Confidence Region for a Multivariate Normal Mean. \emph{Annals of Statistics}, \bold{8}(4): 716-761.

Carvalho, C.M., Polson, N.G., and Scott., J.G. (2010). The Horseshoe Estimator for Sparse Signals. \emph{Biometrika}, \bold{97}(2): 465-480. 

Strawderman, W.E. (1971). Proper Bayes Minimax Estimators of the Multivariate Normal Mean. \emph{Annals of Mathematical Statistics}, \bold{42}(1): 385-388.
}

\author{
Ray Bai and Malay Ghosh
}

\examples{

n <- 100 
p <- 40
q <- 3      # number of response variables is 3
p_act <- 5  # number of active (nonzero) predictors is 5

#############################
# Generate design matrix X. #
#############################
set.seed(123)
times <- 1:p
rho <- 0.5
H <- abs(outer(times, times, "-"))
V <- rho^H
mu <- rep(0, p)
# Rows of X are simulated from MVN(0,V)
X <- MASS::mvrnorm(n, mu, V)
# Center X
X <- scale(X, center=TRUE, scale=FALSE)

############################################
# Generate true coefficient matrix B_true. #
############################################
# Entries in nonzero rows are drawn from Unif[(-5,-0.5)U(0.5,5)] 
B_act <- runif(p_act*q,-5,4)
disjoint <- function(x){
         if(x <= -0.5)
            return(x)
         else
           return(x+1)
     }
B_act <- matrix(sapply(B_act, disjoint),p_act,q)

# Set rest of the rows equal to 0
B_true <- rbind(B_act,matrix(0,p-p_act,q))
B_true <- B_true[sample(1:p),] # permute the rows

#########################################
# Generate true error covariance Sigma. #
#########################################
sigma_sq=2
times <- 1:q
H <- abs(outer(times, times, "-"))
Sigma <- sigma_sq * rho^H
  
############################
# Generate noise matrix E. #
############################
mu <- rep(0,q)
E <- MASS::mvrnorm(n, mu, Sigma)

###############################
# Generate response matrix Y. #
###############################
Y <- crossprod(t(X),B_true) + E


#########################################
# Fit the MBSP model on synthetic data. #
#########################################

# Should use default of max_steps=6000, burnin=1000 in practice
mbsp_model = MBSP(X=X, Y=Y, max_steps=1000, burnin=500)

# indices of the true nonzero rows
true_active_predictors <- which(rowSums(B_true)!=0)
true_active_predictors

# variables selected by the MBSP model
mbsp_model$active_predictors

# the true nonzero rows
B_true[true_active_predictors, ]

# the MBSP model's estimates of the nonzero rows
mbsp_model$B_est[true_active_predictors, ]
}