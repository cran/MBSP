\name{mbsp.tpbn}
\alias{mbsp.tpbn}

\title{MBSP Model with Three Parameter Beta Normal (TPBN) Family}

\description{This function provides a fully Bayesian approach for obtaining a sparse estimate of the \eqn{p \times q} coefficients matrix \eqn{B} in the multivariate linear regression model, 

\deqn{Y = XB+E,} 

using the three parameter beta normal (TPBN) family.
}

\usage{
mbsp.tpbn(X, Y, u=0.5, a=0.5, tau=NA, max.steps=15000, burnin=5000)}

\arguments{
  \item{X}{design matrix of \eqn{n} samples and \eqn{p} covariates.}

  \item{Y}{response matrix of \eqn{n} samples and \eqn{q} response variables.}

  \item{u}{The first parameter in the TPBN family. Defaults to \eqn{u=0.5} for the horseshoe prior.}

  \item{a}{The second parameter in the TPBN family. Defaults to \eqn{a=0.5} for the horseshoe prior.}
  
  \item{tau}{The global parameter. If the user does not specify this (tau=NA), the Gibbs sampler will use \eqn{\tau = 1/(p*n*log(n))}. The user may also specify a value for \eqn{\tau} between \eqn{0} and \eqn{1}, otherwise it defaults to \eqn{\tau = 1/(p*n*log(n))}. }

\item{max.steps}{The total number of iterations to run in the Gibbs sampler. Defaults to 15,000.}

  \item{burnin}{The number of burn-in iterations for the Gibbs sampler. Defaults to 5,000.}

}

\details{
The function performs sparse estimation of \eqn{B} and variable selection from the \eqn{p} covariates. The lower and upper endpoints of the 95 percent posterior credible intervals for each of the \eqn{pq} elements of \eqn{B} are also returned so that the user may assess uncertainty quantification.

In the three parameter beta normal (TPBN) family, \eqn{(u,a)=(0.5,0.5)} corresponds to the horseshoe prior, \eqn{(u,a)=(1,0.5)} corresponds to the Strawderman-Berger prior, and \eqn{(u,a)=(1,a), a > 0} corresponds to the normal-exponential-gamma (NEG) prior. This function uses the horseshoe prior as the default shinkrage prior.

}

\value{The function returns a list containing the following components:
\item{B.est}{ the point estimate of the \eqn{p \times q} matrix \eqn{B} (taken as the componentwise posterior median for all \eqn{pq} entries).}
\item{lower.endpoints}{ The 2.5th percentile of the posterior density (or the lower endpoint of the 95 percent credible interval) for all \eqn{pq} entries of \eqn{B}.}
\item{upper.endpoints}{ The 97.5th percentile of the posterior density (or the upper endpoint of the 95 percent credible interval) for all \eqn{pq} entries of \eqn{B}. }
\item{active.predictors}{ The active (nonzero) covariates chosen by our model from the \eqn{p} total predictors.}
}

\references{
Armagan, A., Clyde, M., and Dunson, D.B. (2011) Generalized Beta Mixtures of Gaussians. In J. Shawe-taylor, R. Zemel, P. Bartlett, F. Pereira, and K. Weinberger (Eds.) \emph{Advances in Neural Information Processing Systems 24}, 523-531.

Bai, R. and Ghosh, M. (2018). High-Dimensional Multivariate Posterior Consistency Under Global-Local Shrinkage Priors. Revision at \emph{Journal of Multivariate Analysis}. arXiv:1711.07635

Berger, J. (1980). A Robust Generalized Bayes Estimator and Confidence Region for a Multivariate Normal Mean. \emph{Annals of Statistics}, \bold{8}(4): 716-761.

Carvalho, C.M., Polson, N.G., and Scott., J.G. (2010). The Horseshoe Estimator for Sparse Signals. \emph{Biometrika}, \bold{97}(2): 465-480. 
Strawderman, W.E. (1971). Proper Bayes Minimax Estimators of the Multivariate Normal Mean. \emph{Annals of Mathematical STatistics}, \bold{42}(1): 385-388.
}

\author{
Ray Bai and Malay Ghosh
}

\examples{

#############################################
#############################################
## EXAMPLE ON SYNTHETIC DATA:              ##
## Can change n, p, q, p.act, max.steps,   ##
## and burnin below to reproduce the       ##
## simulations from Section 5.1 of Bai     ##
## and Ghosh                               ##
#############################################
#############################################
n <- 50
p <- 10
q <- 3
# Active number of predictors is 2
p.act <- 2

#############################
# Generate design matrix X. #
#############################
times <- 1:p
rho <- 0.5
H <- abs(outer(times, times, "-"))
V <- rho^H
mu <- rep(0, p)
# Rows of X are simulated from MVN(0,V)
X <- mvrnorm(n, mu, V)
# Center X
X <- scale(X, center=TRUE, scale=FALSE)

#########################################
# Generate true coefficient matrix B_0. #
#########################################
# Entries in nonzero rows are drawn from Unif[(-5,-0.5)U(0.5,5)] 
B.act <- runif(p.act*q,-5,4)
disjoint <- function(x){
         if(x <= -0.5)
            return(x)
         else
           return(x+1)
     }
B.act <- matrix(sapply(B.act, disjoint),p.act,q)

# Set rest of the rows equal to 0
B.true <- rbind(B.act,matrix(0,p-p.act,q))
B.true <- B.true[sample(1:p),] # permute the rows

########################################
# Generate true error covariance Sigma #
########################################
sigma.sq=2
times <- 1:q
H <- abs(outer(times, times, "-"))
Sigma <- sigma.sq * rho^H
  
###########################
# Generate noise matrix E #
###########################
mu <- rep(0,q)
E <- mvrnorm(n, mu, Sigma)

##############################
# Generate response matrix Y #
##############################
Y <- crossprod(t(X),B.true) + E

#########################################
# Run the MBSP model on synthetic data. #
#########################################

# For optimal estimation, change max.steps to 15,000 
# and change burnin to be 5000

mbsp.model = mbsp.tpbn(X=X, Y=Y, max.steps = 2000, burnin=1000)

\dontrun{
################################
# Compute MSE_est and MSE_pred #
################################
MSE.est <- sum((mbsp.model$B.est-B.true)^2)/(p*q)
MSE.est # Note that in the paper it is scaled by factor of 100

MSE.pred <- sum((crossprod(t(X), mbsp.model$B.est - B.true))^2)/(n*q)
MSE.pred # Note that in the paper it is scaled by a factor of 100

################################
# Compute the FDR, FNR, and MP #
################################

# Initialize vector for classification of predictors i = 1,...,p
predictor.classifications <- rep(0, p)
# Initialize vector for true classifications of predictors i = 1,...,p
true.classifications <- rep(0,p)

# True active predictors in B.True
true.active.indices <- which(rowSums(B.true) != 0)
# Rest true signals as 1
true.classifications[true.active.indices] <- 1

# Active predictors according to our estimates
predictor.classifications[mbsp.model$active.predictors] <- 1

# Keep track of false positives and false negatives
false.pos <- length(which(predictor.classifications != 0 & true.classifications == 0))
tot.pos <- length(which(predictor.classifications != 0))
false.neg <- length(which(predictor.classifications == 0 & true.classifications != 0))
tot.neg <- length(which(predictor.classifications == 0))

# Compute FDR, FNR, and MP
fdr.rate <- false.pos/tot.pos
fnr.rate <- false.neg/tot.neg
mis.prob <- (false.pos+false.neg)/p

fdr.rate
fnr.rate
mis.prob
}

\dontshow{
ptm <- proc.time()

##############################
# EXAMPLE ON SYNTHETIC DATA: #
##############################
n <- 30
p <- 20
q <- 4
# Active number of predictors is 3
p.act <- 3

#############################
# Generate design matrix X. #
#############################
times <- 1:p
rho <- 0.5
H <- abs(outer(times, times, "-"))
V <- rho^H
mu <- rep(0, p)
# Rows of X are simulated from MVN(0,V)
X <- mvrnorm(n, mu, V)
# Center X
X <- scale(X, center=TRUE, scale=FALSE)

#########################################
# Generate true coefficient matrix B_0. #
#########################################
# Entries in nonzero rows are drawn from Unif[(-5,-0.5)U(0.5,5)] 
B.act <- runif(p.act*q,-5,4)
disjoint <- function(x){
         if(x <= -0.5)
            return(x)
         else
           return(x+1)
     }
B.act <- matrix(sapply(B.act, disjoint),p.act,q)

# Set rest of the rows equal to 0
B.true <- rbind(B.act,matrix(0,p-p.act,q))
B.true <- B.true[sample(1:p),] # permute the rows

########################################
# Generate true error covariance Sigma #
########################################
sigma.sq=2
times <- 1:q
H <- abs(outer(times, times, "-"))
Sigma <- sigma.sq * rho^H
  
###########################
# Generate noise matrix E #
###########################
mu <- rep(0,q)
E <- mvrnorm(n, mu, Sigma)

##############################
# Generate response matrix Y #
##############################
Y <- crossprod(t(X),B.true) + E

#########################################
# Run the MBSP model on synthetic data. #
#########################################

# For optimal estimation, change max.steps to 15,000 
# and change burnin to be 5000

mbsp.model = mbsp.tpbn(X=X, Y=Y, max.steps = 1000, burnin=500)

mbsp.model$B.est
mbsp.model$active.predictors

sim.time <- (proc.time()-ptm)[3]
sim.time
}


\dontrun{
#
#
#############################################
############################################# 
## MBSP analysis of yeast cell cycle       ##
## data set (Section 5.2 of Bai and Ghosh) ##              
#############################################
#############################################

# Load yeast data set
data(yeast)

# Set seed
set.seed(12345)

# Design matrix X and response matrix Y
X <- yeast$x
X <- scale(X, center=TRUE, scale=FALSE)
Y <- yeast$y
Y <- scale(Y, center=TRUE, scale=FALSE)

# Make sure they are matrices
X <- matrix(X, nrow=nrow(X))
Y <- matrix(Y, nrow=nrow(Y))

###################################
# Run the MBSP model on the yeast #
# cell cycle data set             #
###################################
mbsp.model = mbsp.tpbn(X=X,Y=Y)

# Display the active predictors (correspond to row indices in B)
mbsp.model$active.predictors

# Display names of four of the active TFs
colnames(yeast$x)[2]
colnames(yeast$x)[38]
colnames(yeast$x)[61]
colnames(yeast$x)[94]

# For horizontal axis (time)
time <- seq(from=0, to=119, by=7)

##############################################
# Open pdf to plot 4 of the active TFs.      #
# This reproduces Figure 1 of Bai and  Ghosh #
##############################################
pdf(file=file.path(tempdir(), "significantTFs.pdf"), width=5, height=5)
par(mfrow=c(2,2), mar=c(2,2,2,2))
plot(time, mbsp.model$B.est[2,], type="l", cex.axis = 0.5, xlab="", 
     main="ACE2", ylim=c(-0.6,0.6), ylab="")
abline(h=0)
lines(time, mbsp.model$lower.endpoints[2,], lty=2)
lines(time, mbsp.model$upper.endpoints[2,], lty=2)

plot(time, mbsp.model$B.est[38,], type="l", cex.axis = 0.5, xlab="", 
     main="HIR1", ylim=c(-0.6,0.6), ylab="")
abline(h=0)
lines(time, mbsp.model$lower.endpoints[38,], lty=2)
lines(time, mbsp.model$upper.endpoints[38,], lty=2)

plot(time, mbsp.model$B.est[61,], type="l", cex.axis = 0.5, xlab="",
     main="NDD1", ylim=c(-0.6,0.6), ylab="")
abline(h=0)
lines(time, mbsp.model$lower.endpoints[61,], lty=2)
lines(time, mbsp.model$upper.endpoints[61,], lty=2)

plot(time, mbsp.model$B.est[94,], type="l", cex.axis = 0.5, xlab="", 
     main="SWI6", ylim=c(-0.6,0.6), ylab="")
abline(h=0)
lines(time, mbsp.model$lower.endpoints[94,], lty=2)
lines(time, mbsp.model$upper.endpoints[94,], lty=2)
dev.off()
}

}
