#' @name ACHC-package
#' @docType package
#' @aliases ACHC
#' @title Analysis of convex hull coverage
#' @author Michael Collyer
#' @return Key functions for this package:
#' \item{\code{\link{achc}}}{Permutation procedure to generate many random hulls and measure
#' the amount of space covered by individual hulls.}
#' \item{\code{\link{plot.achc}}}{Generates plots associated with ACHC.}

#' The name, "ACHC", is an acronym for, "Analysis of Convex Hull Coverage."  Through
#' the various functions in this package, one can determine if the distribution of convex hulls in a data
#' space, represented by groups to be compared, is consistent with or different than expectation from
#' random assembly of observations.  ACHC fits a uniform grid of points to the total expanse of multidimensional space
#' covered by observations.  Hulls for groups are described, and the frequency of 0, 1, ..., g hulls for 
#' g groups are measured at every point.  A relative frequency distribution is generated from these frequencies.  Next,
#' the procedure is repeated many times with random assigment of observations to groups.  Confidence bands
#' are generated for the 0, 1, 2, ..., g groups.  By comapring the observed distribution to the confidence limits
#' generated, one can ascertain if the data space is meaningfully partitioned in some way.
#' 

#' @import stats
#' @import graphics
#' @import utils
#' @import geometry
#' @import RRPP
#' @import rgl

NULL


uniform.cube <- function(Y, pts){
  p <- NCOL(Y)
  tar <- ceiling(pts^(1/p))
  scope <- sapply(1:NCOL(Y), function(j) {
    y <- Y[,j]
    max(y) - min(y)
  })
  coords <- lapply(1:p, function(j){
    y <- Y[, j]
    seq(min(y), max(y), (max(y) - min(y))/(tar -1))
  })
  cube <- expand.grid(coords)
  as.matrix(cube)
}


uniform.grid <- function(Y, pts) {
  p <- ncol(Y)
  pt.scale <- seq(1,(p+2),0.05)
  hull.obs <- convhulln(Y, options="Fa")
  uni.cubes <- lapply(1:length(pt.scale), function(j){
    uniform.cube(P, pts = 1000*pt.scale[j])
  })
  
  uclen <- sapply(uni.cubes, nrow)
  
  uni.cubes.red <- list()
  uni.cubes.red[[1]] <- uni.cubes[[1]]
  for(i in 2:length(uclen)) {
    if(uclen[i] != uclen[i-1]) uni.cubes.red <- c(uni.cubes.red, list(uni.cubes[[i]]))
  }
  
  Result <- lapply(1:length(uni.cubes.red), function(j){
    uc <- uni.cubes.red[[j]]
    result <- list()
    for(i in 1:nrow(uc)) {
      ar <- convhulln(rbind(uc[i,], Y), options="Fa")$vol
      if(ar == hull.obs$vol) result$keep = c(result$keep, i) 
    }
    result
  })
  
  res.check <- sapply(Result, function(x) length(x$keep))
  best <- which.min(abs(pts - res.check))
  
  uc.best <- uni.cubes.red[[best]]
  uc.best[Result[[best]]$keep,]
}

groups.at.points <- function(Y = Y, group = group, grid = grid) {
  require(geometry)
  glev <- levels(group)
  ng <- nlevels(group)
  g.a.p <- sapply(1:nrow(grid), function(j){
    g.point <- grid[j,]
    hd <- sapply(1:ng, function(jj){
      y <- Y[group == glev[jj],]
      yy <- rbind(g.point, y)
      v1 <- convhulln(y, options = "FA")$vol
      v2 <- convhulln(yy, options = "FA")$vol
      res <- ifelse(v1==v2, 1, 0)
      res
    })
    
    hd
  })
  
  rownames(g.a.p) <- glev
  t(g.a.p)
}
