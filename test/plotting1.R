library(dplyr)
library(coda.base)
source('test/ggbiplot.R')

plot_clrbiplot = function(E, title, subtitle, xlimits, groups=NULL){
  clr = coordinates(E, 'clr', label = names(E))
  names(clr) = names(E)
  pc = prcomp(clr)

  ggbiplot(pc, scale = 1, groups = groups, alpha = 0.2, labels.size = 4#vlength=6,
  ) + xlim(xlimits) +
    #ylim(c(-5,5)) + xlim(c(-10,10)) +
    theme_minimal() +
    theme(legend.position = 'top') +
    ggtitle(title, subtitle = subtitle)
}

library(dplyr)
load('test/fitting1.RData')

plots = lapply(sheets, function(sheet){
  E = as.data.frame(fitting[[sheet]][[3]])
  xlimits = c(-3,3)
  plot_clrbiplot(E, "CLR biplot after LRNM imputation", sheet, xlimits)
})

pdf(file = 'test/fitting1.pdf', width = 6.5, height = 4.5)
for(i in 1:length(plots)){
  print(plots[[i]])
}
dev.off()
