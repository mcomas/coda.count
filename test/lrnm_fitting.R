library(dplyr)
library(coda.count)
library(coda.base)
ALPHA = 10+c(100, 10, 5)

MU_ALR = digamma(ALPHA[1:2])-digamma(ALPHA[3])
MU = coordinates(composition(MU_ALR, 'alr'))

SIGMA_ALR = diag(trigamma(ALPHA[1:2])) + trigamma(ALPHA[3])
SIGMA_CLR = alr_basis(3) %*% SIGMA_ALR %*% t(alr_basis(3))
SIGMA = t(ilr_basis(3)) %*% SIGMA_CLR %*% ilr_basis(3)

# set.seed(1)
# X = rlrnm(MU, SIGMA, 100, 1000)
#
# fit = fit_lrnm(X)
# fit$mu
# fit$sigma

X_ = simplex_lattice(30, 3)

df_lrnm = as_tibble(X_, .name_repair = ~paste0('X',1:3)) %>%
  mutate(
    f = apply(., 1, dlrnm, mu = MU, sigma = SIGMA)
  )

library(ggtern)
ggtern(data = df_lrnm) +
  geom_point(aes(x = X2, y = X1, z = X3, col = f), size = 2) +
  scale_color_continuous(low = 'blue', high = 'red') +
  theme_minimal()

#
#
# ############
# coordinates(exp(digamma(ALPHA)))
# coordinates(ALPHA)
#
# inv.digamma = function(y){
#   x = ifelse(y >= -2.22, exp(y) + 0.5, -1/(y + digamma(1)))
#   x <- x - (digamma(x) - y) / trigamma(x)
#   x <- x - (digamma(x) - y) / trigamma(x)
#   x <- x - (digamma(x) - y) / trigamma(x)
#   x <- x - (digamma(x) - y) / trigamma(x)
#   x
# }
# inv.digamma(MU_ALR)
# inv.digamma(digamma(34320000))
# inv.trigamma = function(y){
#   x <- 0.5+1/y
#   tri <- trigamma(x)
#   x <- x+tri*(1-tri/y)/psigamma(x,deriv=2)
#   tri <- trigamma(x)
#   x <- x+tri*(1-tri/y)/psigamma(x,deriv=2)
#   tri <- trigamma(x)
#   x <- x+tri*(1-tri/y)/psigamma(x,deriv=2)
#   tri <- trigamma(x)
#   x <- x+tri*(1-tri/y)/psigamma(x,deriv=2)
#   x
# }
# inv.trigamma(trigamma(c(10,240)))
#
# # In github -> limma/R/ebayes.R
# # inv.trigamma <- function(x) {
# #   #	Solve trigamma(y) = x for y
# #   #	Gordon Smyth
# #   #	8 Sept 2002.  Last revised 12 March 2004.
# #
# #   #	Non-numeric or zero length input
# #   if(!is.numeric(x)) stop("Non-numeric argument to mathematical function")
# #   if(length(x)==0) return(numeric(0))
# #
# #   #	Treat out-of-range values as special cases
# #   omit <- is.na(x)
# #   if(any(omit)) {
# #     y <- x
# #     if(any(!omit)) y[!omit] <- Recall(x[!omit])
# #     return(y)
# #   }
# #   omit <- (x < 0)
# #   if(any(omit)) {
# #     y <- x
# #     y[omit] <- NaN
# #     warning("NaNs produced")
# #     if(any(!omit)) y[!omit] <- Recall(x[!omit])
# #     return(y)
# #   }
# #   omit <- (x > 1e7)
# #   if(any(omit)) {
# #     y <- x
# #     y[omit] <- 1/sqrt(x[omit])
# #     if(any(!omit)) y[!omit] <- Recall(x[!omit])
# #     return(y)
# #   }
# #   omit <- (x < 1e-6)
# #   if(any(omit)) {
# #     y <- x
# #     y[omit] <- 1/x[omit]
# #     if(any(!omit)) y[!omit] <- Recall(x[!omit])
# #     return(y)
# #   }
# #
# #   #	Newton's method
# #   #	1/trigamma(y) is convex, nearly linear and strictly > y-0.5,
# #   #	so iteration to solve 1/x = 1/trigamma is monotonically convergent
# #   y <- 0.5+1/x
# #   iter <- 0
# #   repeat {
# #     iter <- iter+1
# #     tri <- trigamma(y)
# #     dif <- tri*(1-tri/x)/psigamma(y,deriv=2)
# #     y <- y+dif
# #     if(max(-dif/y) < 1e-8) break
# #     if(iter > 50) {
# #       warning("Iteration limit exceeded")
# #       break
# #     }
# #   }
# #   y
# # }
