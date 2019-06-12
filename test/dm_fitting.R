library(coda.count)
ALPHA = 10+c(100, 10, 5)

set.seed(1)
X = rdm(ALPHA, 100, 1000)

alpha = fit_dm(X)
alpha

X_ = simplex_lattice(30, 3)

df_dm = as_tibble(X_, .name_repair = ~paste0('X',1:3)) %>%
  mutate(
    f = apply(., 1, ddm, alpha = ALPHA)
  )

library(ggtern)
ggtern(data = df_dm) +
  geom_point(aes(x = X2, y = X1, z = X3, col = f), size = 3) +
  scale_color_continuous(low = 'blue', high = 'red') +
  theme_minimal()

