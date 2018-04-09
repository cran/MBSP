\name{matrix.normal}
\alias{matrix.normal}

\title{Matrix-Normal Distribution}

\description{This function provides a way to draw a sample from the matrix-normal distribution, given the mean matrix, the covariance structure of the rows, and the covariance structure of the columns.
}

\usage{
matrix.normal(M, U, V)}

\arguments{
  \item{M}{mean \eqn{a \times b} matrix}

  \item{U}{\eqn{a \times a} covariance matrix (covariance of rows).}

  \item{V}{\eqn{b \times b} covariance matrix (covariance of columns).}

}

\details{
This function provides a way to draw a random \eqn{a \times b} matrix from the matrix-normal distribution,

\deqn{MN(M, U, V),}

where \eqn{M} is the \eqn{a \times b} mean matrix, \eqn{U} is an \eqn{a \times a} covariance matrix, and \eqn{V} is a \eqn{b \times b} covariance matrix.

}

\value{A randomly drawn \eqn{a \times b} matrix from \eqn{MN(M,U,V)}.
}

\author{
Ray Bai and Malay Ghosh
}

\examples{

# Draw a random 50x20 matrix from MN(O,U,V),
# where:
#    O = zero matrix of dimension 50x20
#    U has AR(1) structure,
#    V has sigma^2*I structure

# Specify Mean.mat
p <- 50
q <- 20
Mean.mat <- matrix(0, nrow=p, ncol=q)

# Construct U
rho <- 0.5
times <- 1:p
H <- abs(outer(times, times, "-"))
U <- rho^H

# Construct V
sigma.sq <- 2
V <- sigma.sq*diag(q)

# Draw from MN(Mean.mat, U, V)
mn.draw <- matrix.normal(Mean.mat, U, V)
}
