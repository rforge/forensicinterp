% Generated by roxygen2 (4.0.2): do not edit by hand
\name{TwoLevelLR_Background}
\alias{TwoLevelLR_Background}
\title{Compute background parameters for the two level likelihood ratio.}
\usage{
TwoLevelLR_Background(data, fac)
}
\arguments{
\item{data}{a \eqn{n x p} numeric matrix, with \eqn{p \geq 2}. The background data containing \eqn{n} measurements on \eqn{p} variables.}

\item{fac}{a factor of length \eqn{p}, indicating the 'population' of each measurement.}
}
\value{
A list containing the overall and group means (\code{all.means} and \code{group.means}), within- and between group variances and covariances matrices (\code{W} and \code{B}).
}
\description{
Compute background parameters (overall and group means, within- and between group variances and covariances matrices) for the two level likelihood ratio.
}
\author{
Silvia Bozza \cr Alexandre Thiery
}
\references{
Bozza S, Taroni F, Marquis R and Schmittbuhl M (2008). "Probabilistic evaluation of handwriting evidence: likelihood ratio for authorship." \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, \strong{57} (3), pp. 329-341.
}
\seealso{
\code{\link{TwoLevelLR}}
}

