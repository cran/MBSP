#########################################
# Function to sample from Matrix Normal #
#########################################

# M = mean of matrix
# U = covariance matix of the columns
# V = covariance matrix of the rows

matrix_normal = function(M, U, V){
  a <- dim(M)[1]
  b <- dim(M)[2]
  
  # Draw Z from MN(O, I, I)
  Z <- matrix(stats::rnorm(a*b,0,1), a, b)
  
  # Cholesky decomposition of U and V
  L1 <- chol(U)
  L2 <- chol(V)
  
  # Return draw from MN(M,U,V)
  return(M + crossprod(L1,Z) %*% L2)
}