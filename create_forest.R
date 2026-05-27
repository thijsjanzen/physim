require(tidyverse)
ref_table <- read_tsv("test.txt")
sum_stats <- ref_table[, 6:76]
param_names <- colnames(ref_table)[1:5]

forest_list <- list()

for (i in 1:5) {
  cat(i, "\n")
  param <- ref_table[, i]
  data2 <- data.frame(param, sum_stats)
  colnames(data2)[1] <- "r"
  forest_list[[i]] <- abcrf::regAbcrf(r ~ ., data2, ntree = 1000, paral = TRUE)

}

saveRDS(forest_list, "forest_list.rds")
