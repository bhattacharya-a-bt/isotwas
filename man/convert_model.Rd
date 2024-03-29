% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_model.R
\name{convert_model}
\alias{convert_model}
\title{Convert isotwas model list to tibble}
\usage{
convert_model(model, snp_annot = NULL, snp_var = NULL)
}
\arguments{
\item{model, }{isotwas model}

\item{snp_annot, }{annotations for SNPs, typically map obj from bigsnpr}

\item{snp_var, }{column name of SNP marker id}
}
\value{
tibble of the model
}
\description{
The function converts an isotwas model and returns a tibble
}
