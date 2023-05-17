\name{MBSP}
\alias{MBSP}

\title{MBSP Model with Three Parameter Beta Normal (TPBN) Family}

\description{This function provides a fully Bayesian approach for obtaining a (nearly) sparse estimate of the \eqn{p \times q} regression coefficients matrix \eqn{B} in the multivariate linear regression model, 

\deqn{Y = XB+E,} 

using the three parameter beta normal (TPBN) family. Here \eqn{Y} is the \eqn{n \times q} matrix with \eqn{n} samples of \eqn{q} response variables, \eqn{X} is the \eqn{n \times p} design matrix with \eqn{n} samples of \eqn{p} covariates, and \eqn{E} is the \eqn{n \times q} noise matrix with independent rows. The complete model is described in Bai and Ghosh (2018).

If there are \eqn{r} confounding variables which \emph{must} remain in the model and should \emph{not} be regularized, then these can be included in the model by putting them in a
separate \eqn{n \times r} confounding matrix \eqn{Z}. Then the model that is fit is

\deqn{Y = XB+ZC+E,}

where \eqn{C} is the \eqn{r \times q} regression coefficients matrix corresponding to the confounders. In this case, we put a flat prior on \eqn{C}. By default, confounders are not included.

If the user desires, two information criteria can be computed: the Deviance Information Criterion (DIC) of Spiegelhalter et al. (2002) and the widely applicable information criterion (WAIC) of Watanabe (2010).
}

\usage{
MBSP(Y, X, confounders=NULL, u=0.5, a=0.5, tau=NA, 
     max_steps=6000, burnin=1000, save_samples=TRUE,
     model_criteria=FALSE)}

\arguments{
  \item{Y}{Response matrix of \eqn{n} samples and \eqn{q} response variables.}

  \item{X}{Design matrix of \eqn{n} samples and \eqn{p} covariates. The MBSP model regularizes the regression coefficients \eqn{B} corresponding to \eqn{X}.}

  \item{confounders}{Optional design matrix \eqn{Z} of \eqn{n} samples of \eqn{r} confounding variables. By default, confounders are not included in the model (\code{confounders=NULL}). However, if there are some confounders that \emph{must} remain in the model and should \emph{not} be regularized, then the user can include them here.}

  \item{u}{The first parameter in the TPBN family. Defaults to \eqn{u=0.5} for the horseshoe prior.}

  \item{a}{The second parameter in the TPBN family. Defaults to \eqn{a=0.5} for the horseshoe prior.}
  
  \item{tau}{The global parameter. If the user does not specify this (tau=NA), the Gibbs sampler will use \eqn{\tau = 1/(p*n*log(n))}. The user may also specify any value for \eqn{\tau} strictly greater than 0; otherwise it defaults to \eqn{1/(p*n*log(n))}. }

\item{max_steps}{The total number of iterations to run in the Gibbs sampler. Defaults to \code{6000}.}

  \item{burnin}{The number of burn-in iterations for the Gibbs sampler. Defaults to \code{1000}.}
  
\item{save_samples}{A Boolean variable for whether to save all of the posterior samples of the regression coefficients matrix B and the covariance matrix Sigma. Defaults to \code{"TRUE"}.}

\item{model_criteria}{A Boolean variable for whether to compute the following information criteria: DIC (Deviance Information Criterion) and WAIC (widely applicable information criterion). Can be used to compare models with (for example) different choices of \code{u}, \code{a}, or \code{tau}. Defauls to \code{"FALSE"}.}

}

\details{
The function performs (nearly) sparse estimation of the regression coefficients matrix \eqn{B} and variable selection from the \eqn{p} covariates. The lower and upper endpoints of the 95 percent posterior credible intervals for each of the \eqn{pq} elements of \eqn{B} are also returned so that the user may assess uncertainty quantification.

In the three parameter beta normal (TPBN) family, \eqn{(u,a)=(0.5,0.5)} corresponds to the horseshoe prior, \eqn{(u,a)=(1,0.5)} corresponds to the Strawderman-Berger prior, and \eqn{(u,a)=(1,a), a > 0} corresponds to the normal-exponential-gamma (NEG) prior. This function uses the horseshoe prior as the default shrinkage prior.

The user also has the option of including an \eqn{n \times r} matrix with \eqn{r} confounding variables. These confounders are variables which are included in the model but should \emph{not} be regularized. 

Finally, if the user specifies \code{model_criteria=TRUE}, then the \code{MBSP} function will compute two model selection criteria: the Deviance Information Criterion (DIC) of Spiegelhalter et al. (2002) and the widely applicable information criterion (WAIC) of Watanabe (2010). This permits model comparisons between (for example) different choices of \code{u}, \code{a}, and \code{tau}. The default horseshoe prior and choice of \code{tau} performs well, but the user may wish to experiment with \code{u}, \code{a}, and \code{tau}. In general, models with \emph{lower} DIC or WAIC are preferred.
}

\value{The function returns a list containing the following components:
\item{B_est}{ The point estimate of the \eqn{p \times q} matrix \eqn{B} (taken as the componentwise posterior median for all \eqn{pq} entries).}
\item{B_CI_lower}{ The 2.5th percentile of the posterior density (or the lower endpoint of the 95 percent credible interval) for all \eqn{pq} entries of \eqn{B}.}
\item{B_CI_upper}{ The 97.5th percentile of the posterior density (or the upper endpoint of the 95 percent credible interval) for all \eqn{pq} entries of \eqn{B}.}
\item{active_predictors}{ The row indices of the active (nonzero) covariates chosen by our model from the \eqn{p} total predictors.}
\item{B_samples}{ All \code{max_steps-burnin} samples of \eqn{B}.}
\item{C_est}{ The point estimate of the \eqn{r \times q} matrix \eqn{C} corresponding to the confounders (taken as the componentwise posterior median for all \eqn{rq} entries). This matrix is not returned if there are no confounders (i.e. \code{confounders=NULL}).}
\item{C_CI_lower}{ The 2.5th percentile of the posterior density (or the lower endpoint of the 95 percent credible interval) for all \eqn{rq} entries of \eqn{C}. This is not returned if there are no confounders (i.e. \code{confounders=NULL}).}
\item{C_CI_upper}{ The 97.5th percentile of the posterior density (or the upper endpoint of the 95 percent credible interval) for all \eqn{rq} entries of \eqn{C}. This is not returned if there are no confounders (i.e. \code{confounders=NULL}) }
\item{C_samples}{ All \code{max_steps-burnin} samples of \eqn{C}. This is not returned if there are no confounders (i.e. \code{confounders=NULL})}
\item{Sigma_est}{ The point estimate of the \eqn{q \times q} covariance matrix \eqn{\Sigma} (taken as the componentwise posterior median for all \eqn{q^2} entries).}
\item{Sigma_CI_lower}{ The 2.5th percentile of the posterior density (or the lower endpoint of the 95 percent credible interval) for all \eqn{q^2} entries of \eqn{\Sigma}.}
\item{Sigma_CI_upper}{ The 97.5th percentile of the posterior density (or the upper endpoint of the 95 percent credible interval) for all \eqn{q^2} entries of \eqn{\Sigma}.}
\item{Sigma_samples}{ All \code{max_steps-burnin} samples of \eqn{C}.}
\item{DIC}{ The Deviance Information Criterion (DIC), which can be used for model comparison. Models with smaller DIC are preferred. This only returns a numeric value if \code{model_criteria=TRUE} is specified.}
\item{WAIC}{ The widely applicable information criterion (WAIC), which can be used for model comparison. Models with smaller WAIC are preferred. This only returns a numeric value if \code{model_criteria=TRUE} is specified. The WAIC tends to be more stable than DIC.}
}

\references{
Armagan, A., Clyde, M., and Dunson, D.B. (2011) Generalized beta mixtures of Gaussians. In J. Shawe-Taylor, R. Zemel, P. Bartlett, F. Pereira, and K. Weinberger (Eds.) \emph{Advances in Neural Information Processing Systems 24}, 523-531.

Bai, R. and Ghosh, M. (2018). High-dimensional multivariate posterior consistency under global-local shrinkage priors. \emph{Journal of Multivariate Analysis}, \bold{167}: 157-170.

Berger, J. (1980). A robust generalized Bayes estimator and confidence region for a multivariate normal mean. \emph{Annals of Statistics}, \bold{8}(4): 716-761.

Carvalho, C.M., Polson, N.G., and Scott., J.G. (2010). The horseshoe estimator for sparse signals. \emph{Biometrika}, \bold{97}(2): 465-480. 

Spiegelhalter, D. J., Best, N. G., Carlin, B. P., and Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{64}(4): 583-639.
 
Strawderman, W.E. (1971). Proper Bayes Minimax Estimators of the Multivariate Normal Mean. \emph{Annals of Mathematical Statistics}, \bold{42}(1): 385-388.

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. \emph{Journal of Machine Learning Research}, \bold{11}: 3571-3594.
}

\author{
Ray Bai and Malay Ghosh
}

\examples{

###################################
# Set n, p, q, and sparsity level #
###################################

n <- 100
p <- 40
q <- 3 # number of response variables is 3
p_act <- 5 # number of active (nonzero) predictors is 5

#############################
# Generate design matrix X. #
#############################
set.seed(1234)
times <- 1:p
rho <- 0.5
H <- abs(outer(times, times, "-"))
V <- rho^H
mu <- rep(0, p)
# Rows of X are simulated from MVN(0,V)
X <- mvtnorm::rmvnorm(n, mu, V)
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
E <- mvtnorm::rmvnorm(n, mu, Sigma)

##############################
# Generate response matrix Y #
##############################
Y <- crossprod(t(X),B_true) + E

# Note that there are no confounding variables in this synthetic example

#########################################
# Fit the MBSP model on synthetic data. #
#########################################

# Should use default of max_steps=6000, burnin=1000 in practice.
mbsp_model = MBSP(Y=Y, X=X, max_steps=1000, burnin=500, model_criteria=FALSE)

# Recommended to use the default, i.e. can simply use: mbsp_model = MBSP(Y, X)
# If you want to return the DIC and WAIC, have to set model_criteria=TRUE.


# indices of the true nonzero rows
true_active_predictors <- which(rowSums(B_true)!=0)
true_active_predictors

# variables selected by the MBSP model
mbsp_model$active_predictors

# true regression coeficients in the true nonzero rows
B_true[true_active_predictors, ]

# the MBSP model's estimates of the nonzero rows
mbsp_model$B_est[true_active_predictors, ]
}