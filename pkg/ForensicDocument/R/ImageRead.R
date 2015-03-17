#' Read a bitmap image into a binary matrix
#' 
#' Reads an image stored from various file format (see Details) into a matrix. of \code{0}'s and \code{1}'s.
#' 
#' @param file A character: name of the file to read from.
#' @param character_pixel A interger value, indicating which pixel value is from the character.
#' @details
#' This method allows to read images of format type:
#' \enumerate{
#'   \item PNG, using the \pkg{png} package.
#' }
#' 
#' @details The image is read as a matrix. of the dimensions height x width, with real values between 0 and 1. If there is more than one channel in the image, the mean values across the channels is taken.
#' @import png
#' @seealso the method \code{readPNG} from the package \code{png}.
#' @author Alexandre Thiery
#' @keywords internal
ImageRead <- function(file, character_pixel = 1) 
{    
#     file = "inst/extdata/fourier.png"
    m = regexpr("\\.[^\\.]+$", file)
    ext = regmatches(file, m)
    
    result = NULL
    if(ext == ".png") {
        result = png::readPNG(file)
    } 
    # else if(ext %in% c(".jpeg", ".jpg")) {
    #     result = jpeg::readJPEG(file)
    # }
    
    if(is.null(result))
    {
        stop(sprintf("File '%s' not in correct format", file))
    }
    if(length(dim(result)) == 3)
        result = (result[,,1]+result[,,2]+result[,,3])/3
    
    result = round(result,0)
    if(character_pixel == 0)
        result = abs(result - 1)
    
    return(result)
}