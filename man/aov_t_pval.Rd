% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funciones.R
\name{aov_t_pval}
\alias{aov_t_pval}
\title{P-value of analysis of Variance}
\usage{
aov_t_pval(y, g, pvdigits = 4)
}
\arguments{
\item{y}{Response variable.}

\item{g}{Grouping variable.}

\item{pvdigits}{number of digits for the p-value.}
}
\value{
p-value of t-test when comparing two groups or p-value of the analysis
of variance test aov(y~g) when comparing more than two groups. For comparing two
groups t.test is preferred as it allows for different variances in each group.
}
\description{
P-value of analysis of Variance
}
\examples{
df <- data.frame(g=gl(3,10),y=rnorm(30))
aov_t_pval(df$y,df$g,pvdigits=3)
}
