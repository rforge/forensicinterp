#' Extracts images black pixels' x and y coordinates.
#' 
#' This method allows to extract am image black pixels' x and y coordinates.
#' 
#' @param image A numeric matrix, representing the black and white image (a pixel of value \code{0} is from the background).
#' @param black A numerical value (\code{0} or \code{1}), indicating the cell value of a black pixel.
#' @return A \code{data.frame} containg the x and y coordinates of the black pixels. 
#' @author Alexandre Thiery
#' @keywords internal
ImageGetPoints <- function(image, black = 0)
{
    ncol = ncol(image)
    nrow = nrow(image)
    
    xy = lapply(1:nrow, function(y) {
        x = which(image[y,] == black)
        if(length(x) == 0) {
            return(data.frame(x = NULL, y = NULL))
        } 
        return(data.frame(x = x, y = rep(y, length(x))))
    })
    xy = do.call("rbind", xy)
    xy
}

