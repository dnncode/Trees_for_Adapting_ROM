# rSVD
#X = u.combine

rSVD = function(X,r){
  
  n = nrow(X)
  m = ncol(X)
  
  # create the P matrix
  P = diag(m)
  #P = P[, floor(seq(1, m, m/10))]
  P = matrix(rnorm(m*r),nrow=m)
  
  # QR for Z
  Z = X %*% P
  test = qr(Z)
  Q = qr.Q(test , complete = FALSE)
  
  # compute Y
  Y = ginv(Q) %*% X
  
  # svd on Y
  svd.Y = fast.svd(Y)
  U.Y = svd.Y$u
  
  # get U
  U = Q %*% U.Y
  
  return(U)

  
}