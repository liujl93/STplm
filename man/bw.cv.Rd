% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bwselection.R
\name{bw.cv}
\alias{bw.cv}
\title{Bandwidth Selection}
\usage{
bw.cv(h, st, KERN1, KERN2, loctim, P, X, Z, beta.fixed, Sigma.fixed,
  sigma2.fixed, Sigma = NULL)
}
\arguments{
\item{h}{the prespecified bandwidth}

\item{st}{the index of kernel}

\item{KERN1, KERN2}{two candidate kernels}

\item{loctim}{(s,t) the spatiotemporal index}

\item{P}{the dimension of data}

\item{X}{X}

\item{Z}{y}

\item{beta.fixed, Sigma.fixed, sigma2.fixed}{prespecified parameters and covariance matrix}
}
\value{
CV, GCVce criterion
}
\description{
Bandwidth Selection Using Cross-Valication and Generalized Cross-Validation Correlated Criterion
}
