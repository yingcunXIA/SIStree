knn.CV.k <- function(X, Y, max.k)
{
  cvi = 1
  p = dim(X)[2]
  n = dim(X)[1]
  Y = as.numeric(Y)

  D = as.matrix(dist(X))
  Iorder = t(apply(D, 2, order))
  Yi = matrix(Y[Iorder], n, n)

  K = seq(3, max.k, 4)
  nK = length(K)
  cv = rep(0, nK)

  for (ik in 1:nK)
  {
      yp = rowMeans(Yi[,2:(K[ik]+1)])
      cv[ik] = mean( (yp>0.5) != Y)
  }
  k= K[which.min(cv)]

  return(list(k=k, cv=min(cv)))
  
}
