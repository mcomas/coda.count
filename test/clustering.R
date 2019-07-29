X = parliament2015[,c('jxsi', 'psc', 'pp', 'catsp', 'cs', 'cup')]
labels = parliament2015$name
X$ind = X$jxsi + X$cup
X$esp_e = X$psc
X$esp_d = X$pp + X$cs
X$nsnc = X$catsp
X = as.matrix(df <- X[,c('ind', 'esp_e', 'esp_d', 'nsnc')])

fit = fit_lrnm(X, method = 'hermite', probs = TRUE, B = ilr_basis(4))
fit2 = fit_dm(X)
P = as.data.frame(t(t(X) + fit2[[1]]))
D = count.dist(X, fit$mu, fit$sigma)
hc = hclust(D)
g = cutree(hc, k = 5)

ggtern(data = df) +
  geom_mask() +
  geom_point(aes(x=ind, y=esp_d, z=esp_e, col = factor(g)), alpha = 0.4) +
  theme_minimal()
P = as.data.frame(fit$P)
names(P) = colnames(X)
ggtern(data = P) +
  geom_mask() +
  geom_point(aes(x=ind, y=esp_d, z=esp_e), alpha = 0.4) +
  theme_minimal()

H = coordinates(as.data.frame(fit$P))
ggplot(data = H) +
  geom_point(aes(x=ilr1, y=ilr2, col = factor(g)), alpha = 0.4) +
  theme_minimal()

labels[g==4]
