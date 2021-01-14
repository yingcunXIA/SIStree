library(dHSIC)
library(dcov)
library(HHG)


hhg <- function(x, y)
{
  x = as.matrix(x)
  y = as.matrix(y)
  Dy = as.matrix(dist(y))
  Dx = as.matrix(dist(x))
  n = nrow(Dx)
  hhg = hhg.test(Dy, Dx, nr.perm = 0)$sum.chisq /n/(n-2)/(n-3)
  return(hhg)
}  

CvM<-function(X,Y){
  #INPUT
  #X: an n x p matrix
  #Y: an n x 1 vector
  #OUTPUT
  #CvM correlation
  #mbkr = int {corr(I(X <= x), I(Y <= y))}^2 dF(x)dF(y)
  
  n = length(Y)  
  X = matrix(X, nrow=n)
  Y = matrix(Y, ncol=1)
  n=dim(X)[1]
  p=dim(X)[2]
  mbkr = rep(0,p)
  
  m = floor(n^0.8)
  Qpoints = (2:(m-1))/(m+1)
  m = length(Qpoints)
  yknots = quantile(Y, Qpoints)
  
  yknots = Y; m=n
  if (n > 2)
  {
    YY = matrix(Y, n, m) <= matrix(yknots, n, m, byrow=TRUE)
    Ymean = colMeans(YY)
    YY=  YY - matrix(Ymean, n, m, byrow=TRUE)
    YY = YY/matrix(sqrt(colSums(YY^2))+1.0e-10, n, m, byrow=TRUE)
    YY = t(YY)
    
    for (k in 1:p){
      xknots = quantile(X[,k], Qpoints)
      xknots = X[,k]
      Xk =  matrix(X[,k], n, m) <= matrix(xknots, n, m, byrow=TRUE)
      Xmean = colMeans(Xk)
      Xk = Xk- matrix(Xmean, n, m, byrow=TRUE)
      Xk = Xk/matrix(sqrt(colSums(Xk^2))+1.0e-10, n, m, byrow=TRUE)
      mbkr[k] = mean((YY%*%Xk)^2)
    }
  }
  return(mbkr)  
}

dcor3d <- function(x, y)
{
  d = dcor2d(x, y)
  d = max(c(d, dcor2d(abs(x-mean(x)), y)))
  return(d)
}


MI <- function(X, Y)
{
  Y = scale(Y)
  X = as.matrix(scale(X))
  n = length(Y)
  p = ncol(X)
  
  
  
  Y = ((rank(Y)+1)/(n+2))
  h = 2*sd(Y)/n^0.2
  
  dy = as.matrix(dist(Y/h))^2 
  #  fy = exp(-Y^2/2)
  fy = rowMeans(exp(-dy/2))
  
  mi = rep(0, p)
  for (ip in 1:p)
  {
    xi = X[,ip]
    
    xi = ((rank(xi)+1)/(n+2))
    
    dx = as.matrix(dist(xi/h))^2 
    #    fx = exp(-xi^2/2)
    fx = rowMeans(exp(-dx/2))
    
    fxy = rowMeans(exp(-(dx+dy)/2))
    
    mi[ip] = mean(log(fxy/(fx+1.0e-10)/(fy+1.0e-10)))
  }
  
  return(mi)
}


VS0 <- function(X, y, i1=-1, V0=-1, method="dcor")
{
  n = length(y)
  p = ncol(X)
  

  if (i1 < 0)  
    i1 = min(which(V0 == max(V0)))
  
  var1 = X[,i1] - mean(X[,i1])
  
  if (method=="Pearson")
  {
    vcut0 = quantile(var1, 0.5)
    IL = which(var1 <= vcut0)
  }
  else
  {
    var1.sort = sort(unique(var1)) 
    
    m = length(var1.sort)
    if (m > 100)
    {
      cut = quantile(var1.sort, seq(0, 1, length=100)[10:90])
    }else
    {
      if (m>10)
        cut = var1.sort[seq(5, m-5)]
      else
        cut = mean(var1.sort)
    }
    
    vcut0 = 1.0e10
    for (icut in cut)
    {
      I = which(var1<=icut)
      #      print(1)
      #      print(c(length(I), n-length(I)))
      if (method=="dcor")
      {
        d1 = 0
        if (length(I) > 5)
          d1 = dcor3d(var1[I], y[I])*length(I)/n
        
        d2 = 0
        if (n-length(I) > 5)
          d2 = dcor3d(var1[-I], y[-I])*(1-length(I)/n)
        vcut = d1+d2
      }
      
      if (method=="hhg")
      {
        d1 = 0
        if (length(I) > 5)
          d1 = hhg(var1[I], y[I])*length(I)/n
        
        d2 = 0
        if (n-length(I) > 5)
          d2 = hhg(var1[-I], y[-I])*(1-length(I)/n)
        vcut = d1+d2
      }
      
      #      print(2)
      
      if (method=="mi")
      {
        vcut = MI(var1[I], y[I])*length(I)/n
        + MI(var1[-I], y[-I])*(1-length(I)/n)
      }
      
      if (method=="HSIC")
      {
        vcut = dhsic(var1[I], y[I])$dHSIC*length(I)/n
        + dhsic(var1[-I], y[-I])$dHSIC*length(I)/n*(1-length(I)/n)
      }
      
      if (method=="CvM")
      {
        vcut = CvM(var1[I], y[I])*length(I)/n + 
          CvM(var1[-I], y[-I])*(1-length(I)/n)
      }
      
      vcut[!is.finite(vcut)] <- 1.0e5
      if (vcut < vcut0)
      {
        IL = I
        vcut0 = vcut
      }
    }
  }
  IR = setdiff(1:n, IL)
  
  if ((length(IR) < 10)|(length(IL) < 10))
  {
    split = FALSE
    V0=NA
    V1 = NA 
    IL = NA 
    IR = NA
    VL.adj = NA 
    VR.adj = NA
    VL = NA
    VR = NA 
    vsplit=NA
  }else  
  {
    split = TRUE
    
    VL = V0
    VR = V0
    for (i in 1:p)  
    {
      
      if (method=="dcor")
      {
        VL[i] = dcor3d(X[IL,i], y[IL])
        VR[i] = dcor3d(X[IR,i], y[IR])
      }
      
      if (method=="hhg")
      {
        VL[i] = hhg(X[IL,i], y[IL])
        VR[i] = hhg(X[IR,i], y[IR])
      }
      
      if (method=="HSIC")
      {
        VL[i] = dhsic(X[IL,i], y[IL])$dHSIC
        VR[i] = dhsic(X[IR,i], y[IR])$dHSIC
      }
      
    }
    if (method=="CvM")
    {
      VL = CvM(X[IL,], y[IL])
      VR = CvM(X[IR,], y[IR])
    }
    if (method=="mi")
    {
      VL = MI(X[IL,], y[IL])
      VR = MI(X[IR,], y[IR])
    }
    if (method=="Pearson")
    {
      VL = as.numeric(cor(X[IL,], y[IL])^2)
      VL[is.na(VL)] <- 0
      
      VR = VL
      if (length(IR) > 3)
      {
        VR = as.numeric(cor(X[IR,], y[IR])^2)
        VR[is.na(VR)] <- 0
      }
    }
    
    
    VL.adj = VL
    VR.adj = VR
    VL.adj[i1] = max(VL[i1], V0[i1])     # using the maximum value for the splitting variable
    VR.adj[i1] = max(VR[i1], V0[i1])
    
    V1 = (V0 + VL.adj*length(IL)/n + VR.adj*length(IR)/n)/2
  }
  
  
  return(list(split=split, V0=V0, V1 = V1, IL = IL, IR = IR, VL.adj = VL.adj, 
              VR.adj = VR.adj, VL = VL, VR = VR, vsplit=i1))
  
}


SIStree <- function(X, y, depth=5, method="dcor", min.size=20)
{
  # output:   ORDER.importance: px(depth+1) 
  #    each column gives the order of importance of each variable 
  #     x_i, i = 1. ,,,p in each layer k, k = 1,..., depth+1. 
  #     k = 1 corresponding to Dcor 
  #       ORDER.importance.av: px1
  #     averaged rank 
  #  method=c("dcor","CvM", "mi", "HSIC", "Pearson", "hhg")
  
  # Example
  # X = matrix(rnorm(200*2000), 200, 2000)
  # y = X[,1] + X[,2]^2 + 2*X[,4]*X[,5] + rnorm(200)
  # sis = SIStree(X, y, depth=5)
  #
  # # show the first 20 variables ranked by SIS
  # sis$rank.SIS[1:20]
  #
  # # show the first 20 variables ranked by SIStree
  # sis$rank.SIStree[1:20]
  
  if (sum(method == c("dcor","CvM", "mi", "HSIC", "Pearson", "hhg"))< 1)
  {
    warning("method must be one of 'dcor','CvM', 'mi', 'HSIC', 'Pearson'")
    method = "dcor"
  }
  
  n = length(y)
  p = ncol(X)
  
  #  if (depth > 10)
  #  {
  #    print("depth is be smaller than 10")
  #    depth = 10
  #  }
  
  DCOR = list()
  DCOR.adj = list()
  DCOR.adj[[1]] = matrix(0, 1, 2)
  SPLIT = list()
  SIZE = list()
  
  SPLIT[[1]] = matrix(1:n, n, 1)
  SIZE[[1]] = n
  
  dcorALL = matrix(0, p, depth+1)
  split_vars = c()
  
  V0 = rep(0, p)
  if (method=="dcor")
  {
    for (i in 1:p)
      V0[i] = dcor3d(X[,i], y)
  }
  
  if (method=="hhg")
  {
    for (i in 1:p)
      V0[i] = hhg(X[,i], y)
  }
  
  if (method=="mi")
  {
    V0 = MI(X,y)
  }
  if (method=="CvM") V0 = CvM(X, y)
  
  if (method=="HSIC")
  {
    for (i in 1:p)
      V0[i] = dhsic(y, X[,i])$dHSIC
  }
  
  if (method=="Pearson") 
  {
    V0 = as.numeric(cor(y, X)^2)
    V0[is.na(V0)] = 0
  }
  
  for (idepth in 1:depth)
  {
    dcori.adj = c()
    dcori = c()
    spliti = c()
    size = c()
    M = 1; if (idepth > 1) M = ncol(DCOR[[idepth]])
    for (ipart in 1:M)
    {
      Ii = which(SPLIT[[idepth]][,ipart] > 0)
      if (length(Ii) <= min.size) 
      {
        dcori.adj = cbind(dcori.adj, DCOR.adj[[idepth]][,ipart])
        dcori = cbind(dcori, DCOR[[idepth]][,ipart])
        spliti = cbind(spliti, SPLIT[[idepth]][,ipart])
        size = c(size, SIZE[[idepth]][ipart]) 
      }else
      {
        if (idepth >1)
        {
#          V0 = 0
#          for (j in 1:M)
#            V0 = V0 + DCOR[[idepth]][,j]*SIZE[[idepth]][j]/n
          
          V0 = DCOR[[idepth]][,ipart]
        }
        
        if (length(unique(V0))==1)
        {
          dcori.adj = cbind(dcori.adj, DCOR.adj[[idepth]][,ipart])
          dcori = cbind(dcori, DCOR[[idepth]][,ipart])
          spliti = cbind(spliti, SPLIT[[idepth]][,ipart])
          size = c(size, SIZE[[idepth]][ipart]) 
        }else
        {
          
          parti = VS0(X[Ii,], y[Ii], i1=-1, V0, method=method)
          
          if (idepth == 1) 
          {
            DCOR.adj[[1]] = matrix(parti$V0, ncol=1)
            DCOR[[1]] = matrix(parti$V0, ncol=1)
            dcorALL[,1] = parti$V0
          }
          
          
          if (!parti$split)
          {
            dcori.adj = cbind(dcori.adj, DCOR.adj[[idepth]][,ipart])
            dcori = cbind(dcori, DCOR[[idepth]][,ipart])
            spliti = cbind(spliti, SPLIT[[idepth]][,ipart])
            size = c(size, SIZE[[idepth]][ipart]) 
          }else
          {

            split_vars = c(split_vars, parti$vsplit)
            
            dcori = cbind(dcori, parti$VL, parti$VR)
            dcori.adj = cbind(dcori.adj, parti$VL.adj, parti$VR.adj)
            
            IL = rep(0, n)
            IR = rep(0, n)
            
            IL[Ii[parti$IL]] = 1
            IR[Ii[parti$IR]] = 1
            
            spliti = cbind(spliti, IL, IR)
            size = c(size, sum(IL), sum(IR))
          }
        }
      }
    }
    DCOR[[idepth+1]] = dcori
    
    #    print(dcori)
    #    print(split_vars)

    DCOR.adj[[idepth+1]] = dcori.adj
    SPLIT[[idepth+1]] = spliti
    SIZE[[idepth+1]] = size
    

    for (i in 1:length(size))
      dcorALL[,idepth+1] = dcorALL[,idepth+1] + dcori.adj[,i]*size[i]/n 
    
    if (max(size) < min.size)
    {
      depth = idepth
      break
    }
  }
  
  dcorALL.cum = dcorALL
  for (i in 2:(idepth+1))
    dcorALL.cum[,i] = rowMeans(dcorALL.cum[,1:i])
  
  ORDER.importance = matrix(0, p, (idepth+1))
  for (i in 1:(idepth+1))
    ORDER.importance[,i] = order(dcorALL.cum[,i], decreasing = TRUE)
  
  colnames(ORDER.importance) = c('dcor', paste('layer.', 1:depth, sep=""))
  
  DCOR.av = rowMeans(dcorALL.cum)  ## averaged dcor
  ORDER.importance.av = order(DCOR.av, decreasing = TRUE)
  
  split_vars = unique(split_vars)
  
  return(list(rank.SIS =ORDER.importance[,1],
              rank.SIStree=ORDER.importance.av,
              split.vars = split_vars,
              all =ORDER.importance))
}  


