#' Extracts the radius of a closed shape, relative to its barycenter.
#' 
#' This method allows to extract  the radius of a closed shape, contained in a image, relative to its barycenter.
#' 
#' @param image A numeric matrix, representing the black and white image (a pixel of value \code{0} is from the background).
#' @param seq A numeric vector, giving the angles to subsample the closed shape.
#' @return A numeric vector containing the radius values.
#' @author Alexandre Thiery
#' @keywords internal
ImageGetClosedShapeRadius = function(image, seq)
{
    ncol = ncol(image)
    nrow = nrow(image)
    
    xy = ImageGetPoints(image = image, black = 1)
    
    bary = colMeans(xy)
    xy = data.frame(x = xy$x - bary['x'], y = xy$y - bary['y'])
    #xy = data.frame(t(t(xy) - (bary)))
    
    
    angles = atan2(xy$y, xy$x)
    angles[angles < 0] = angles[angles < 0] + 2*pi 
    
    angle.sort = sapply(seq, function(x) {
        temp1 = abs(angles - x)
        temp2 = abs(angles - x + 2*pi)
        m1 = min(temp1)
        m2 = min(temp2)
        if(m1 < m2) {
            return(c(m1, which(temp1 == m1)[1]))
        } else {
            return(c(m2, which(temp2 == m2)[1]))
        }
    })
    radius = sqrt(rowSums(xy[angle.sort[2,],]^2))
    return(as.numeric(radius))
}