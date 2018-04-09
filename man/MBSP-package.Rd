\name{MBSP-package}
\alias{MBSP-package}
\alias{MBSP}
\docType{package}
\title{Multivariate Bayesian Model with Shrinkage Priors}
\description{
This package contains the function, \code{mbsp.tpbn}, for implementing the Multivariate Bayesian Model with Shrinkage Priors (MBSP), with the three parameter beta normal (TPBN) family, as described in Section 4.2 of Bai and Ghosh (2018).

This model gives a row sparse estimate of the \eqn{p \times q} coefficients matrix \eqn{B} in the multivariate regression model,

\deqn{Y = X B + E,}

where \eqn{Y} is the \eqn{n \times q} matrix of responses, \eqn{X} is  the \eqn{n \times p} design matrix with \eqn{p} covariates, and \eqn{E} is the \eqn{n \times q} noise matrix. In \eqn{E}, each row is distributed as multivariate normal \eqn{N(0, \Sigma)}, where \eqn{\Sigma} is estimated from the data using the inverse-Wishart prior.

In addition, the package provides the function, \code{matrix.normal}, for sampling of a random matrix from the matrix-normal distribution.
}

\details{

The DESCRIPTION file:
\packageDESCRIPTION{MBSP}
\packageIndices{MBSP}

This package implements the MBSP model for sparse Bayesian multivariate linear regression using the TPBN family. Our package performs both estimation and variable selection. The \code{mbsp.tpbn} function also returns the endpoints of the credible intervals (i.e. the 2.5th and 97.5th percentiles) for every single coefficient of the matrix \eqn{B} so that uncertainty quantification can be assessed.

}
\author{
Ray Bai and Malay Ghosh

Maintainer: Ray Bai <raybai07@ufl.edu>
}

\references{
Armagan, A., Clyde, M., and Dunson, D.B. (2011) Generalized Beta Mixtures of Gaussians. In J. Shawe-taylor, R. Zemel, P. Bartlett, F. Pereira, and K. Weinberger (Eds.) \emph{Advances in Neural Information Processing Systems 24}, 523-531.

Bai, R. and Ghosh, M. (2018). High-Dimensional Multivariate Posterior Consistency Under Global-Local Shrinkage Priors. Revision at \emph{Journal of Multivariate Analysis}, arXiv:1711.07635.

Carvalho, C.M., Polson, N.G., and Scott., J.G. (2010). The Horseshoe Estimator for Sparse Signals. \emph{Biometrika}, \bold{97}(2): 465-480. 

}