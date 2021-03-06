% Generated by roxygen2 (4.0.2): do not edit by hand
\name{ExtractFourier}
\alias{ExtractFourier}
\title{Extract fourier parameters from close shape}
\usage{
ExtractFourier(files, n.fourier = 8, n.samp = 128, skeletonize = TRUE,
  character_pixel = 0, output, verbose = TRUE)
}
\arguments{
\item{files}{A vector of string, giving the images' filenames to analyse.}

\item{n.fourier}{A integer, giving the number of harmonics (or fourier parameters) to extract.}

\item{n.samp}{A integer, giving the number of points to subsample the closed shape.}

\item{skeletonize}{A logical value, indicating if the character should be skeletonized.}

\item{character_pixel}{A integer value (0 or 1), indicating which pixel value is from the character (see details).}

\item{output}{A string, giving the name of the output file. See details.}

\item{verbose}{A logical value, indicating if progress is to be printed on the console.}
}
\value{
If \code{output} is missing, it returns a list of matrices (one for each \code{files}), containing the extracted An and Bn values. Otherwise, the results are printed in a \emph{csv} file (one record or row for each image).
}
\description{
This method allows to extract fourier parameters (i.e. An and Bn parameters) from closed shapes present in images.
}
\details{
The vector containing the images' filenames, the \code{file} argument, is pre-processed before the analysis. Duplicated file names and files that do not exists that are discared.

\code{character_pixel} is a integer value, either 0 or 1, it indicates which pixel value is from the character. For example, if \code{character_pixel = 1} then pixels that have a value of \code{1} will correspond to the character, and pixels that have a value of \code{0} will correspond to the background.
}
\examples{
\dontrun{
## Example using a sample character O from Marquis et al.
## In this image, black pixels are from the character (value of 0),
## and white pixels are from the bakcground (value of 1). The
## character is not skeletonized, thus we use skeletonize = TRUE.

## We will consider a re-sampling size of the character's skeleton of
## 128 (n.samp = 128). All of the Fourier parameters will be extracted,
## thus n.fourier = n.samp / 2 = 64.
file = system.file("extdata", "fig-O.png", package = "ForensicDocument")
image = ExtractFourier(files = file, n.fourier = 64, n.samp = 128,
         skeletonize = TRUE, character_pixel = 0)
}
}
\author{
Alexandre Thiery
}

