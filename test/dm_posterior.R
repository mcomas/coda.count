library(coda.count)
library(coda.base)
library(dplyr)
library(tidyr)
library(ggtern)

ALPHA = c(10, 50, 50)
x = c(10,5,20)
h = coordinates(x)

Pc = simplex_lattice(100, 3)
Pc = Pc[rowSums(Pc == 0) == 0,]
P = Pc/rowSums(Pc)
H = coordinates(P)

library(gtools)
df = bind_cols(as_tibble(P, .name_repair = ~paste0('p',1:3)),
               as_tibble(H, .name_repair = ~paste0('h',1:2)))
df$f_prior = scale(apply(P, 1, ddirichlet, ALPHA))[,1]
df$f_posterior = scale(apply(P, 1, ddirichlet, ALPHA + x))[,1]

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
  theme_void() + theme(legend.position = 'none')

###########

X = apply(simplex_lattice(10,2), 2, rev)
sum(p <- apply(X, 1, dlrnm, -0.5, 0.1))
names(p) = sprintf("(%d,%d)", X[,1], X[,2])
barplot(p, cex.axis = 0.8, cex.names = 0.8)
