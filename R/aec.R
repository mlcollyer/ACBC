#' Analysis of Ellipsoid Coverage
#'
#' Function obtains relative frequencies of ellipsoid coverage in the represented data space,
#' and uses a permutation procedure to generate confidence limits for random assignmnet of observations to groups.
#'
#' A description is needed here

#' @param dat Either a data frame or matrix of data to be analyzed.  Because of computational
#' limitations, the number of variables is currently limited to 5.
#' @param std A logical value that if TRUE finds standard deviates of the data 
#' (data are both centered and scaled by variable standard deviations).
#' @param group A factor or vector coerible to factor for defining groups.
#' @param grid.points The desired number of points, sampled from a uniform distribution of points
#' in the data space, within a convex hull for all observed points.  This number might be
#' less than the maximum possible number of points (which could be huge).  
#' @param grid.space The approximate spacing of uniform points along each axis.  For example, 0.05 means points will
#' be placed at increments that are 5 percent of the expanse of data, per axis.
#' @param confidence The confidence level for ellipsoids, based on the covariance matrix of 
#' data by groups.  Multivariate normality is assumed in estimation.  If NULL, then ellipsoids
#' merely reflect the span of eigenvalues for the covariance matrix of the data.  Otherwise, the value should be between
#' 0.01 and 1.
#' @param ellipse.density A numeric value to indicate how many discrete points (in a circle)
#' are used to approximate the continuous ellipse function.  More points mean a more precise curve, 
#' but increase computation time.  The default, 120 points, is the same as 3 degrees (pi/60 radians) increments.
#' @param iter The number of iterations (permutations) to run for the test.  Because the observed case counts as one
#' iteration, this should be the number desired, minus one.
#' @param seed Change the random seed, if desired.  If NULL, the seed will equal the number of permutations.
#' @param print.progress A logical value to indicate if permutation progress should be printed to the screen.
#' This is useful for analyses that will run a long time.

#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class \code{achc} is a list containing the following
#' \item{grid}{The grid points obtained.}
#' \item{analysis}{A matrix of 0s and 1s for convex hull presenece at grid points (rows)
#' by groups (columns) for each permutation.}
#' \item{group}{The factor of group levels used; useful for downstream functions.}
#' \item{perms}{The number of permutations.}
#' \item{perm.schedule}{The sampling frames in each permutation.}
#' \item{std}{Whether data were standardized.}

#' @examples 
#' 
#' # 3 dimensions of data
#' library(RRPP)
#' data("Pupfish")

#' group <- interaction(Pupfish$Sex, Pupfish$Pop)
#' P <- prcomp(Pupfish$coords)$x[,1:3] # first 3 PCs
#' 
#' pupCHC <- aec(P, std = FALSE, group = group, iter = 99, confidence = 0.95,
#'  ellipse.density = 60, grid.points = 100, grid.space = 0.1)
#' pupCHC
#' summary(pupCHC, confidence = 0.95)
#' plot(pupCHC, lwd = 2)
#' plot(pupCHC, lwd = 2, confidence = 0.99)

aec <- function(dat, std = FALSE, group, grid.points = 500, grid.space = 0.05,
                 confidence = NULL, ellipse.density = 120,
                 iter = 99, seed = NULL, print.progress = TRUE){

  if(!inherits(dat, c("data.frame", "matrix")))
    stop("\nData must be a data frame or matrix.")
  dat <- as.data.frame(dat)
  if(std) data <- scale(dat)
  Y <- as.matrix(dat)
  if(!is.numeric(Y))
    stop("\nNot all data are numeric.")
  pca <- prcomp(Y)
  d <- pca$sdev^2
  k <- which(zapsmall(d) > 0)
  Y <- pca$x[, k]
  p <- length(k)
  group <- as.factor(group)
  
  # Make grid
  if(print.progress)
    cat("\nSampling grid points This might take a moment.\n")
  grid <- uniform.grid.sample(Y, pts = grid.points, 
                              pt.scale = grid.space)
  
  # Random outcomes
  perms <- iter + 1
  if(print.progress) {
    cat("\nRandom hull calculations and assignment of groups:\n")
    pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
  }
  
  ind <- perm.schedule(nrow(Y), iter = iter, seed = seed)
  aec.args <- list(Y = Y, group = group, grid = grid, confidence = confidence,
                    ellipse.density = ellipse.density)
  analysis <- lapply(1:perms, function(j){
    step <- j
    if(print.progress) setTxtProgressBar(pb,step)
    aec.args$Y <- Y[ind[[j]], ]
    do.call(groups.at.points, aec.args)
  })
  
  step <- perms + 1
  if(print.progress) {
    setTxtProgressBar(pb,step)
    close(pb)
  }
  
  names(analysis) <- c("obs", paste("iter", 1:iter, sep = "."))
  
  out <- list(grid = grid, analysis = analysis, group = group, perms = perms,
              perm.schedule = ind, std = std)
  class(out) <- "aec"
  out
}
