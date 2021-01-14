library(class)
source('knnCVk.R')    

knn.cv1 <- function(x, y)
{
  err0 = 1.0e5
  for (k in seq(3, 51, 2))
  {
    err = mean(y != class::knn.cv(x, y, k=k))
    if (err <err0)
      err0 = err
  }
  
  return(err0)
}

knn.cv_ <- function(x, y, fold)
{
  n = length(y)
  y = as.factor(y)
  x = as.matrix(x, nrow=n)
  fold = min(fold, n)
  K = ceiling(n/fold)
  mydata = list(x=x,y=y) 
  err = 0
  
  for (k in 1:fold)
  {
    test = min((k-1)*K +1, n):min(n, k*K)
    train = setdiff(1:n, test)

    A = knn.CV.k(x[train,], y[train], max.k=51)
    k0 = A$k
    kn = knn(x[train,], x[test,], as.factor(y[train]), k=k0) 
    err = err + mean( kn != as.factor(y[test]))/fold
  }
  

  for (k in 1:fold)
  {
    test = sample(1:n,K)
    train = setdiff(1:n, test)
    
    A = knn.CV.k(x[train,], y[train], max.k=51)
    k0 = A$k
    kn = knn(x[train,], x[test,], as.factor(y[train]), k=k0) 
    err = err + mean( kn != as.factor(y[test]))/fold
  }
  
  return(err)  
}

knnCUT <- function(X, y, fold=5, M = -1)
{  
  if (length(unique(y)) > 2)
  {
    stop('THis function only works for two-classes classification')
  }
  
  n = dim(X)[1]
  p = dim(X)[2]

  if (max(M) < 0)
  {
    cut.max = min(p-1,24)
    m = floor(p/2)
    M = 2*1.25^(1:cut.max)
    M = round(M/max(M)*p)
    M = sort(M)
    M0 = c(2:p)[1:cut.max]
    M = (M>=M0)*M + (M0>M)*M0
  }

  cv = c()
  for (jm in M)
  {
    cvj = knn.cv_(X[,1:jm], y, fold=fold)
#    cvj = knn.cv1(X[,1:jm], y)
    cv = c(cv, cvj)
  }
  
  m = M[which(cv==min(cv))]
  m = floor(mean(m))
  
  return(m)
}


