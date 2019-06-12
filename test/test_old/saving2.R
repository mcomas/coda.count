library(dplyr)
library(WriteXLS)
library(readxl)
load('test/fitting2.RData')

expected = list()
for(i in 1:length(fitting)){
  sheet = sheets[i]
  fit = fitting[[i]]

  data = read_excel('~/donato/database_e_carote_conteggi_v2_v5_v6.xls', sheet = sheet, col_names = TRUE)
  vnames = LETTERS[1:ncol(data)]
  nc = c(1, ncol(data), which(colSums(data) == 0))
  vnames = vnames[-nc]

  expected[[sheet]] = as.data.frame(nth(fit,3))
  names(expected[[sheet]]) = vnames
}

WriteXLS(expected, ExcelFileName = "test/expected2.xlsx")
expected1 = expected[[1]]
expected2 = expected[[2]]
expected3 = expected[[3]]
expected4 = expected[[4]]
expected5 = expected[[5]]
expected6 = expected[[6]]
expected7 = expected[[7]]
expected8 = expected[[8]]
expected9 = expected[[9]]
expected10 = expected[[10]]
expected11 = expected[[11]]

save(expected1, expected2, expected3, expected4,
     expected5, expected6, expected7, expected8,
     expected9, expected10, expected11, file = 'test/expected2.RData')

