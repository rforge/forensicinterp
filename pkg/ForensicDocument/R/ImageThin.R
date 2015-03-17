#' Extract the skeleton of a black and white image.
#' 
#' Extracts the skeleton of a black and white image (that is, represented by a matrix).
#' 
#' @details The skeletonisation algorithm implemented is based on \emph{Stentiford (1983)}. It is a recursive thinning process that allows to extract a curve with a thickness of 1 pixel.
#' 
#' @details In the original version of this algorithm, endpoints pixels 1 could not be removed (endpoint pixel are defined as a black element that are 8–connected to only one black element). The algorithm used in this package enables the removal of such pixels, in order to extract closed loops only.
#' 
#' @param matrix A numeric matrix, representing the black and white image (a pixel of value \code{0} is from the background).
#' @return A matrix with values of \code{1} for the contour and values of \code{0} for the background. 
#' @references F.W.M. Stentiford, R.G. Mortimer (1983), "Some new heuristics for thinning binary handprinted characters for OCR." \emph{IEEE Transactions on Systems, Man and Cybernetics}, \bold{13}, pp. 81-84.
#' @author Alexandre Thiery
#' @examples 
#' \dontrun{
#' ## We use the 'fig-O.png' example image available in this package
#' file = system.file("extdata", "fig-O.png", package = "ForensicDocument")
#' 
#' ## Read the image using 'ImageRead' method, or (similarly with 'png::readPNG')
#' image1 = ForensicDocument:::ImageRead(file, character_pixel = 0)
#' 
#' ## Create the image skeleton
#' image2 = ForensicDocument:::ImageThin(image1)
#' 
#' par(mfrow = c(1,2))
#' image(t(image1), xlim = c(0,1), ylim = c(1,0), col = c("white", "black"),
#'  main = "Original Image", asp = 1)
#' image(t(image2), xlim = c(0,1), ylim = c(1,0), col = c("white", "black"),
#'  main = "Skeleton Image", asp = 1)
#' )
#' }
#' @keywords internal
#' @useDynLib ForensicDocument
ImageThin = function(matrix) 
{
    # check matrix
    if(!is.matrix(matrix))
        stop(sprintf("%s is not a matrix.", sQuote("matrix")))
    # parameters, width and height
    w = ncol(matrix)
    h = nrow(matrix)
    
    # remove 1 from edges
    matrix[1,] = 0
    matrix[h,] = 0
    matrix[,1] = 0
    matrix[,w] = 0
    # transform into vector    
    mat = do.call("c", lapply(1:h, function(x) matrix[x,]))
    mat = as.vector(as.integer(mat))
    
    # thin
    r = .C("ImageThin", mat, as.vector(h), as.vector(w))
    
    # result
    return(matrix(r[[1]], ncol = w, byrow = TRUE))
}