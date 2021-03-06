% Generated by roxygen2 (4.0.2): do not edit by hand
\name{ImageGetPoints}
\alias{ImageGetPoints}
\title{Extracts images black pixels' x and y coordinates.}
\usage{
ImageGetPoints(image, black = 0)
}
\arguments{
\item{image}{A numeric matrix, representing the black and white image (a pixel of value \code{0} is from the background).}

\item{black}{A numerical value (\code{0} or \code{1}), indicating the cell value of a black pixel.}
}
\value{
A \code{data.frame} containg the x and y coordinates of the black pixels.
}
\description{
This method allows to extract am image black pixels' x and y coordinates.
}
\author{
Alexandre Thiery
}
\keyword{internal}

