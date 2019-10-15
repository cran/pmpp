library(testthat)
library(pmpp)

test_check("pmpp", filter = "^((?!lint).*)$", perl = TRUE)
