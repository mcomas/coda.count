% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random.R
\name{rmultinomial}
\alias{rmultinomial}
\title{Multinomial random sample}
\usage{
rmultinomial(n = NULL, size, p)
}
\arguments{
\item{n}{number of random samples. Ignored if p is a matrix}

\item{size}{vector to set the multinomial sampling size. vector is reused to have length equal to the number of rows of probs}

\item{p}{vector or matrix specifiying the probability of each class. If a matrix is provided the number of rows
will be equal to the number of samples to generate. Row number i defines the multinomial probability
to be used when generating sample i. In this last scenario, parameter names (size and p) should be provided when calling the function.}
}
\value{
multinomial random sample
}
\description{
Multinomial random sample
}
\examples{
rmultinomial(10, 4, c(0.3, 0.4, 0.3))
probs = matrix(c(0.2, 0.4, 0.4,
                 0.5, 0.1, 0.4), byrow = TRUE, nrow = 6, ncol = 3)
size = c(10, 500, 10000)
rmultinomial(size = size, p = probs)
}
