% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_trees_tbb.R
\name{generate_trees_tbb}
\alias{generate_trees_tbb}
\title{function to perform ABC-SMC}
\usage{
generate_trees_tbb(
  number_of_trees = 1000,
  prior_means = c(1, 1, 1, 1, 1),
  min_tips = 50,
  max_tips = 150,
  crown_age = NULL,
  file_name = "test.txt",
  num_threads = -1
)
}
\arguments{
\item{number_of_trees}{number of trees to generate}

\item{prior_means}{means of exponential prior distributions for: c(la_g,
la_i, mu_g, mu_i, completion_rate)}

\item{min_tips}{minimum number of tips (inclusive)}

\item{max_tips}{maximum number of tips (inclusive)}

\item{crown_age}{crown age}

\item{file_name}{file name to write output to}

\item{num_threads}{number of threads}
}
\value{
a matrix with parameter values and the associated summary statistics
}
\description{
function to perform ABC-SMC
}
