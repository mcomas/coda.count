library(coda.count)
set.seed(1)
n = 10
MU.orig = c(0,0)
SIGMA.orig = diag(3, 2)
X = rlrnm(MU.orig, SIGMA.orig, 100, n)



# Using a threshold of 30, we have
c(1, 58, 41)
# x1~x2+x3
# x2~x3
# c(?, log(58/41))

c(64,33,3)
# x1+x2~x3
# x1~x2
# c(?, log(64,33))

c(92, 2, 6)
# c(?, ?)


log(58/41)
