require(tidyverse)
require(abcrf)
forest_list <- readRDS("forest_list.rds")
ref_table <- read_tsv("test.txt")
sum_stats <- ref_table[, 6:76]
param_names <- colnames(ref_table)[1:5]


crown_age <- 40

b1 <- 0.1
c1 <- 1
b2 <- 0.4
mu1 <- mu2 <- 0

true_params <- c(b1, c1, b2, mu1, mu2)
found <- c()
for (r in 1:100) {

  emp_tree <- physim::pbd_sim_rcpp(pars = true_params,
                                   age = crown_age)
  ref_stats <- as.data.frame(t(treestats::calc_all_stats(emp_tree)))
  to_add <- c()
  for (i in 1:5) {
    cat(r, i, "\n")
    param <- ref_table[, i]
    data2 <- data.frame(param, sum_stats)
    colnames(data2)[1] <- "r"
    vv <- predict(object = forest_list[[i]],
                  obs = ref_stats,
                  training = data2, paral = TRUE)
    to_add[i] <- vv$expectation
  }

  cat(r, to_add, "\n")
  found <- rbind(found, to_add)
}
