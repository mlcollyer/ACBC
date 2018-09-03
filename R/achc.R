#' Analysis of Convex Hull Coverage
#'
#' Function obtains relative frequencies of convex hull coverage in the represented data space,
#' and uses a permutation procedure to generate confidence limits for random assignmnet of observations to groups.
#'
#' A description is needed here

#' @param dat Either a data frame or matrix of data to be analyzed.  Because of computational
#' limitations, the number of variables is currently limited to 5.
#' @param group A factor or vector coerible to factor for defining groups.
#' @param grid.points The targeted number of uniform points to use to defined the expanse of the data space.  Emphasis 
#' is placed on maintaining uniformity, so the optimal solution will likely not match exactly the desired number of points.
#' Optimization involves fitting a spectrum of point desnities and choosing the solution that most closely matches the 
#' number of grid points chosen.  Disparity between targeted and obtained numbers are likely due to edge effects
#' (as small changes in density can add or remove many landmarks at the edges of data space occupation).
#' @param iter The number of iterations (permutations) to run for the test.  Because the observed case counts as one
#' iteration, this should be the number desired, minus one.
#' @param seed AChange the random seed, if desired.  If NULL, the seed will equal the number of permutations.
#' @param print.progress A logical value to indicate if permutation progress should be printed to the screen.
#' This is useful for analyses that will run a long time.

#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class \code{achc} is a list containing the following
#' \item{grid}{grid points obtained.}
#' \item{analysis}{A matrix of 0s and 1s for convex hull presenece at grid points (rows)
#' by groups (columns) for each permutation.}
#' \item{group}{The factor of group levels used; useful for downstream functions.}

#' @examples 
#' 
#' library(RRPP)
#' data("Pupfish")

#' group <- interaction(Pupfish$Sex, Pupfish$Pop)
#' P <- prcomp(Pupfish$coords)$x[,1:3] # first 3 PCs
#'
#' pupCHC <- achc(P, group, iter = 99, grid.points = 1000)
#' pupCHC
#' summary(pupCHC, confidence = 0.95)
#' plot(pupCHC, lwd = 2)
#'
#' # The grid used
#' library(rgl)
#' plot3d(pupCHC$grid)
#' aspect3d("iso")
#' 
achc <- function(dat, group, grid.points = 500,
                 iter = 99, seed = NULL, print.progress = TRUE){
  
  require(geometry)
  
  if(!inherits(dat, c("data.frame", "matrix")))
    stop("\nData must be a data frame or matrix.")
  Y <- as.matrix(dat)
  dat <- as.data.frame(dat)
  p <- ncol(Y)
  if(p > 5) stop("\nToo many variables; not computationally feasible.")
  if(!is.numeric(Y))
    stop("\nNot all data are numeric.")
  group <- as.factor(group)
  
  # Make grid
  if(print.progress)
    cat("\nOptimizing grid.  This might take a moment.\n")
  grid <- uniform.grid(Y, pts = grid.points)
  upts <- nrow(grid)
  if(print.progress)
    cat("\n", upts, "points produced. (This is the optimal solution compared to the desired number of points.\n")
  
  # Observed Hulls
  glev <- levels(group)
  ng <- nlevels(group)
  hulls.obs <- lapply(1:ng, function(j){
    y <- Y[group == glev[j],]
    convhulln(y, options = "FA")
  })
  
  # Random outcomes
  perms <- iter + 1
  if(print.progress) {
    cat("\nRandom hull calculations and assignment of groups:\n")
    pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
  }
  
  ind <- RRPP:::perm.index(nrow(Y), iter = iter, seed = seed)
  achc.args <- list(Y = Y, group = group, grid = grid)
  analysis <- lapply(1:perms, function(j){
    step <- j
    if(print.progress) setTxtProgressBar(pb,step)
    achc.args$Y <- Y[ind[[j]], ]
    do.call(groups.at.points, achc.args)
  })
  
  step <- perms + 1
  if(print.progress) {
    setTxtProgressBar(pb,step)
    close(pb)
  }
  
  names(analysis) <- c("obs", paste("iter", 1:iter, sep = "."))
  
  out <- list(grid = grid, analysis = analysis, group = group)
  class(out) <- "achc"
  out
}
