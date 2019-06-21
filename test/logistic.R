library(coda.base)
logistic = function(h) sqrt(length(h)+1) * prod(composition(h))

library(cubature)
pcubature(logistic, lowerLimit = c(-15), upperLimit = c(15))
# 1
pcubature(logistic, lowerLimit = c(-15,-15), upperLimit = c(15,15))
# 1/2
pcubature(logistic, lowerLimit = c(-15,-15,-15), upperLimit = c(15,15,15))
# 1/(2*3)
vegas(logistic, lowerLimit = c(-15,-15,-15,-15), upperLimit = c(15,15,15,15))
# 1/(2*3*4)

library(data.table)
df = data.table(a=2,3)
a = df
a
df[,2] = 10
df
a
