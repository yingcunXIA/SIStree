library(randomForest)
library(class)
library(HHG)
library(dcov)


#Example 1 ###############################
source('SIStree.R')
X = matrix(rnorm(200*1000), 200, 1000)
y = X[,1] + X[,5]^2 + 2*X[,9]*X[,10] + rnorm(200)
ranks = SIStree(X, y, depth=5, method="dcor")

# the first 20 variables ranked by SIS
sis20 =ranks$rank.SIS[1:20]
print(sis20)

# the first 20 variables ranked by SIStree
sistree20 = ranks$rank.SIStree[1:20]
print(sistree20)

#Example 2 ###############################

source('SIStree.R')
source('knnCUT.R')

# please note the knnCUT/rfCUT only applies to binary classification

X = matrix(rnorm(200*100), 200, 100)
y = (X[,1] + X[,5] + rnorm(200) > 0)
ranks = SIStree(X, y, depth=5, method="dcor")
m = knnCUT(X[,ranks$rank.SIStree], y)

# variables remained after screening

print(ranks$rank.SIStree[1:m])


