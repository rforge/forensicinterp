% Generated by roxygen2 (4.0.2): do not edit by hand
\docType{data}
\name{characterO}
\alias{characterO}
\title{Measurements of Fourier parameters on handwritten characters 'O'.}
\format{This \code{list} contains the following object:
\describe{
  \item{measurements}{A \eqn{554 x 13} \code{matrix}. It contains the extracted fourier parameters (the order 0 harmonic \eqn{A_0} and the first 6ths hamornics \eqn{an_i} and \eqn{bn_i}) for the 554 characters.}
  \item{info}{A \code{data.frame} containing the character's annotation. That is with columns
  \itemize{
  \item \code{writer} the writer's id (from \code{1} to \code{11})
  \item \code{page} the character's page number (from \code{1} to \code{5}).
  \item \code{number} the character's number on the page (from \code{1} to \code{10}).
}
}
}}
\description{
The \code{characterO} is a \code{list} that contains the extracted Fourier parameters and annotation for 650 handwritten characters 'O', written by 11 writers. It is a subset of the data collected by R. Marquis (see \emph{Marquis et al. (2005)}).
}
\author{
Alexandre Thiery \email{alexandre.thiery@unil.ch}
}
\references{
Marquis R, Schmittbuhl M, Mazzela W and Taroni F (2005). "Quantification of the shape of handwritten characters: a step to objective discrimination between writers based on the study of the capital character O." \emph{Forensic Science International}, \bold{150}, pp. 23-32.
}
\keyword{data}
\keyword{document}
\keyword{forensic}
\keyword{handwriting}
\keyword{handwritten}

