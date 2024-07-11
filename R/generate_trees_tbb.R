#' function to perform ABC-SMC
#' @param number_of_trees number of trees to generate
#' @param prior_means means of exponential prior distributions for: c(la_g,
#' la_i, mu_g, mu_i, completion_rate)
#' @param min_tips minimum number of tips (inclusive)
#' @param max_tips maximum number of tips (inclusive)
#' @param crown_age crown age
#' @param file_name file name to write output to
#' @param num_threads number of threads
#' @return a matrix with parameter values and the associated summary statistics
#' @export
generate_trees_tbb <- function(number_of_trees = 1000,
                               prior_means = c(1, 1, 1, 1, 1),
                               min_tips = 50,
                               max_tips = 150,
                               crown_age = NULL,
                               file_name = "test.txt",
                               num_threads = -1) {

  if (is.null(crown_age)) {
    stop("Please either provide a reference tree, or provide the crown age")
  }
  RcppParallel::setThreadOptions(numThreads = num_threads)
  sim_result <- physim::create_ref_table_tbb_par(num_repl = number_of_trees,
                                                 prior_means = prior_means,
                                                 crown_age = crown_age,
                                                 min_lin = min_tips,
                                                 max_lin = max_tips)

  convert_to_phylo <- function(newick_string) {
    phylo_tree <- ape::read.tree(text = newick_string)
    return(phylo_tree)
  }

  phylo_trees <- list()
  for (r in 1:length(sim_result$trees)) {
    phylo_trees[[r]] <- ape::read.tree(text = sim_result$trees[r])
  }


  phylo_trees <- lapply(sim_result$trees, convert_to_phylo)

  cat("simulating trees is done\n")

  # now we calculate stats
  cat("calculating summary statistics for all trees...\n")

 stats <- pbmcapply::pbmclapply(phylo_trees, treestats::calc_all_stats,
                                 mc.cores = num_threads)

  stat_matrix <- matrix(unlist(stats, use.names = FALSE),
                        ncol = 70,
                        byrow = TRUE)

  results <- cbind(sim_result$parameters, stat_matrix)

  test_tree <- ape::rphylo(n = 5, 1, 0)
  test_stats <- treestats::calc_all_stats(test_tree)

  colnames(results) <-
    c("lambda0", "lambda1", "mu0", "mu1", "compl_rate",
      names(test_stats))

  results <- tibble::as_tibble(results)

  readr::write_tsv(results, file = file_name)
  cat(paste("reference table written to:", file_name, "\n"))
  return(results)
}
