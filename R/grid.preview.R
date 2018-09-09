#' Preview of Uniform points
#'
#' Function generates a pairs plot to show the hypothetical distribution of 
#' uniform points that might be generated for \code{\link{achc}}.
#' 
#' @param Y A matrix or data frame of data
#' @param pts The number of points (in all dimensions)
#' @param pt.scale The relative distance between points, as a fraction of the expanse 
#' of representation of data (by variable).  the default, 0.05, means a 5 percent
#' increment along axes.
#' @param ... other arguments for graphical parameters
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @examples 
#' library(RRPP)
#' data("Pupfish")
#' P <- prcomp(Pupfish$coords)$x[,1:3] # first 3 PCs
#' 
#' grid.preview(P, pts = 100, pt.scale = 0.1, pch = 19, cex = 0.3, asp = 1)
#' grid.preview(P, pts = 500, pt.scale = 0.05, pch = 19, cex = 0.3, asp = 1)
grid.preview <- function(Y, pts = 500, pt.scale = 0.05, ...){
  Y <- as.matrix(Y)
  Z <- uniform.grid.sample(Y, pts, pt.scale)
  pairs(Z, panel = "points", ...)
}