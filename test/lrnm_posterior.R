library(coda.count)
library(coda.base)
library(dplyr)
library(tidyr)
library(ggtern)

# MU = c(0,0)
# SIGMA = matrix(c(1,-0.8,
#                  -0.8,1), nrow = 2) * 0.25

ALPHA = c(10, 50, 50)

ALR2ILR = MASS::ginv(alr_basis(3)) %*% ilr_basis(3)

MU_ALR = digamma(ALPHA[1:2])-digamma(ALPHA[3])
MU = coordinates(composition(MU_ALR, 'alr')) # equivalent, MU_ALR %*% ALR2ILR


SIGMA_ALR = diag(trigamma(ALPHA[1:2])) + trigamma(ALPHA[3])
SIGMA = t(ALR2ILR) %*% SIGMA_ALR %*% ALR2ILR

x = c(10,5,20)
h = coordinates(x)

ldnormal(matrix(h, nrow = 1), MU, solve(SIGMA))
mvtnorm::dmvnorm(matrix(h, nrow = 1), mean = MU, sigma = SIGMA, log = TRUE)

Pc = simplex_lattice(100, 3)
Pc = Pc[rowSums(Pc == 0) == 0,]
P = Pc/rowSums(Pc)
H = coordinates(P)


df = bind_cols(as_tibble(P, .name_repair = ~paste0('p',1:3)),
               as_tibble(H, .name_repair = ~paste0('h',1:2)))
df$f_prior = scale(exp(ldnormal(H, MU, solve(SIGMA))))[,1]
df$f_posterior = scale(exp(lpnm_join(x, MU, solve(SIGMA), P, H)) / dlrnm(x, MU, SIGMA))[,1]

dfplot =  df %>%
  gather(distribution, density, starts_with('f_')) %>%
  mutate(distribution = factor(distribution, levels = c('f_prior', 'f_posterior')))

ggtern() +
  geom_point(data=dfplot, aes(x=p2,y=p1,z=p3, col=density), size = 2) +
  geom_point(data = data.frame(x = x[2], y = x[1], z = x[3]), aes(x=x,y=y,z=z), col = 'black') +
  scale_color_continuous(low = 'white', high = 'red') +
  facet_wrap(~distribution) + theme(legend.position = 'none')

ggplot() +
  geom_point(data=dfplot, aes(x=h1,y=h2, col=density), size = 1) +
  geom_point(data = data.frame(x = h[1], y = h[2]), aes(x=x,y=y), col = 'black') +
  scale_color_continuous(low = 'white', high = 'red') +
  facet_wrap(~distribution) +
  coord_equal() +
  theme_minimal() + theme(legend.position = 'none')

###########

X = apply(simplex_lattice(10,2), 2, rev)
sum(p <- apply(X, 1, dlrnm, -0.5, 0.1))
names(p) = sprintf("(%d,%d)", X[,1], X[,2])
barplot(p, cex.axis = 0.8, cex.names = 0.8)
