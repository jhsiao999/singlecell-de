% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulateExpression.R
\name{makeSimCount2groups.filter}
\alias{makeSimCount2groups.filter}
\title{Generate count matrix of all null genes}
\usage{
makeSimCount2groups.filter(counts, Ngenes = NULL,
  sample_method = c("per_gene", "all_genes"))
}
\arguments{
\item{counts}{Gene expression count matrix from a dataset.}

\item{Ngenes}{Number of genes in the simulated dataset.}

\item{Nsample}{Number of samples in each condition.}

\item{pi0}{Proportion of null genes. Default to be 1.}
}
\description{
This is an updated version of makeSimCount2groups which takes in filtered data
              and performs permutation at gene-level and at the sample level.
}
\examples{
library(singleCellRNASeqHumanTungiPSC)
eset <- HumanTungiPSC
counts <- exprs(eset)[,pData(eset)$individual == "NA19101"]

sim_counts <- makeSimCount2groups(counts,
                                  Ngenes = 100,
                                  Nsample = 20,
                                  sample_method = "all_genes")

}
