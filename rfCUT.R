library(randomForest)

rfcv_ <- function(x, y, fold)
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
    test = ((k-1)*K +1):min(n, k*K)
    train = setdiff(1:n, test)

    rf0 =  randomForest(x[train,], y[train])
    pred0 = predict(rf0, x[test,])
    err = err + mean(pred0 != y[test])/fold
  }

  return(err)  
}

rfCUT <- function(X, y, fold=5, M = -1)
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
    cvj = rfcv_(X[,1:jm], y, fold=fold)
    cv = c(cv, cvj)
  }
  
  m = M[which(cv==min(cv))]
  m = floor(mean(m))
  
  return(m)
}
