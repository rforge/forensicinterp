#' Extract fourier parameters from close shape
#' 
#' This method allows to extract fourier parameters (i.e. An and Bn parameters) from closed shapes present in images.
#' 
#' @param files A vector of string, giving the images' filenames to analyse.
#' @param n.fourier A integer, giving the number of harmonics (or fourier parameters) to extract.
#' @param n.samp A integer, giving the number of points to subsample the closed shape.
#' @param skeletonize A logical value, indicating if the character should be skeletonized.
#' @param character_pixel A integer value (0 or 1), indicating which pixel value is from the character (see details).
#' @param output A string, giving the name of the output file. See details.
#' @param verbose A logical value, indicating if progress is to be printed on the console.
#' @return If \code{output} is missing, it returns a list of matrices (one for each \code{files}), containing the extracted An and Bn values. Otherwise, the results are printed in a \emph{csv} file (one record or row for each image).
#' @details The vector containing the images' filenames, the \code{file} argument, is pre-processed before the analysis. Duplicated file names and files that do not exists that are discared. 
#' 
#' \code{character_pixel} is a integer value, either 0 or 1, it indicates which pixel value is from the character. For example, if \code{character_pixel = 1} then pixels that have a value of \code{1} will correspond to the character, and pixels that have a value of \code{0} will correspond to the background. 
#' @examples 
#' \dontrun{
#' ## Example using a sample character O from Marquis et al. 
#' ## In this image, black pixels are from the character (value of 0),
#' ## and white pixels are from the bakcground (value of 1). The 
#' ## character is not skeletonized, thus we use skeletonize = TRUE.
#' 
#' ## We will consider a re-sampling size of the character's skeleton of 
#' ## 128 (n.samp = 128). All of the Fourier parameters will be extracted,
#' ## thus n.fourier = n.samp / 2 = 64.

#' file = system.file("extdata", "fig-O.png", package = "ForensicDocument")
#' image = ExtractFourier(files = file, n.fourier = 64, n.samp = 128, 
#'          skeletonize = TRUE, character_pixel = 0)
#' }
#' @author Alexandre Thiery
#' @export
ExtractFourier = function(files, n.fourier = 8, n.samp = 128, skeletonize = TRUE, character_pixel = 0, output, verbose = TRUE)
{
    if(missing(files))
        stop("No files given!")
    # get unique files
    files = unique(files)
    # length of files
    len.f0 = length(files) 
    # check if exists
    if(any(!file.exists(files)))
        stop("Some files don't exist!")
    # normalise
    files = normalizePath(files)
    # length of files
    len.f1 = length(files) 
    
    # check n.samp
    if(!is.numeric(n.samp) | n.samp <= 0)
        stop("Incorrect subsample number.")
    # check n.fourier
    if(!is.numeric(n.fourier) | (n.fourier < 0 | n.fourier > n.samp/2))
        stop("Incorrect number of fourier coefficiants to extract.")
    
    # check output
    writefile = TRUE
    if(missing(output) || is.null(output))
    {
        writefile = FALSE
    } else {
        if(!is.character(output))
            stop("'output' parameter incorrect format.")
    }
    
    
    
    seq = seq(0, to = 2*pi, length.out = n.samp + 1)[-(n.samp + 1)]
    result = lapply(1:len.f1, function(i) {
        
        if(verbose) message(sprintf(" - %s (%d/%d)", files[i],  i, len.f1))
        A = ImageRead(files[i], character_pixel)
        
        if(skeletonize)
            A = ImageThin(matrix = A)
        
        if(sum(A) <= 100)
        {
            return(NA)
        }
        radius = ImageGetClosedShapeRadius(A, seq)
        
        # get the fourier coefficiants
        result1 = stats::fft(radius)[1:(n.fourier)]
        result1 = result1/pi
        result1[1] = result1[1]/2
        result1 = cbind(Re(result1),Im(result1))
        
        #result1 = do.call("rbind", lapply(stats::fft(radius)[1:(n.fourier+1)], function(x) c(Re(x), Im(x))))
        #result1 = result1/pi
        #result1[1,1] = result1[1,1]/2
        gc()
        colnames(result1) = c("an", "bn")
        rownames(result1) = 0:(n.fourier-1)
        result1
    })
    
    if(writefile)
    {
        result = do.call("rbind", lapply(1:length(result), function(i) {
            if(class(result[[i]]) != "matrix" && is.na(result[[i]]))
                return(rep(NA, n.fourier * 2 + 1))
            return(unlist(result[[i]])[-(n.fourier+2)])
        }))
        
        result = as.data.frame(result)
        colnames(result) = c(paste("a",0:n.fourier,sep = ""),paste("b",1:n.fourier, sep = ""))
        write.csv(x = result, file = output, row.names = files)
        
        return(invisible())
    } else {
        names(result) = files
        return(result)        
    }
}
