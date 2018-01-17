library(dplyr)
library(WriteXLS)
library(readxl)
load('test/fitting1.RData')

expected = list()
for(i in 1:length(fitting)){
  sheet = sheets[i]
  fit = fitting[[i]]

  data = read_excel('~/donato/conteggi_dic17.xlsx', sheet = sheet, col_names = TRUE)
  vnames = gsub(' ', '_', gsub("'","",names(data), fixed=TRUE), fixed=TRUE)
  nc = c(ncol(data), which(colSums(data) == 0))
  vnames = vnames[-nc]


  expected[[sheet]] = as.data.frame(nth(fit,3))
  names(expected[[sheet]]) = vnames
}

WriteXLS(expected, ExcelFileName = "test/expected1.xlsx")
expected1 = expected[[1]]
expected2 = expected[[2]]
expected3 = expected[[3]]
expected4 = expected[[4]]

save(expected1, expected2, expected3, expected4, file = 'test/expected1.RData')

