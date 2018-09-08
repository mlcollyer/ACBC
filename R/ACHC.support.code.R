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


uniform.cube <- function(Y, pt.space){
  p <- NCOL(Y)
  scope <- sapply(1:NCOL(Y), function(j) {
    y <- Y[,j]
    max(y) - min(y)
  })
  coords <- lapply(1:p, function(j){
    y <- Y[, j]
    seq(min(y), max(y), (max(y) - min(y))* pt.space)
  })
  cube <- expand.grid(coords)
  as.matrix(cube)
}


uniform.grid <- function(Y, pt.space) {
  p <- ncol(Y)
  Y.hull <- hull.pts.by.planes(Y)
  dim.index <- combn(p, 2)
  coords <- as.list(array(NA, p))
  for(i in 1: ncol(dim.index)) {
    a.match <- dim.index[1, i]
    b.match <- dim.index[2, i]
    yh <- Y.hull[[i]]
    x <- seq(min(yh[,1]), max(yh[,1]), (max(yh[,1]) - min(yh[,1]))*pt.scale)
    y <- seq(min(yh[,2]), max(yh[,2]), (max(yh[,2]) - min(yh[,2]))*pt.scale)
    z <- expand.grid(x, y)
    ch <- sort(chull(yh))
    check <- sapply(1:nrow(z), function(j){
      z.i <- z[j,]
      names(z.i) <- colnames(yh)
      chp <- sort(chull(rbind(yh, z.i)))
      identical(ch, chp)
    })
    res <- z[check,]
    coords[[a.match]] <- c(coords[[a.match]], z[,1])
    coords[[b.match]] <- c(coords[[b.match]], z[,2])
  }
  coords <- lapply(coords, unique)
  coords <- lapply(coords, na.omit)
  coords
}

uniform.grid.sample <- function(Y, pts, pt.scale){
  U <- uniform.grid(Y, pt.scale)
  p <- ncol(Y)
  tpts <- pts
  pts <- 0
  Y.hull<- hull.pts.by.planes(Y)
  coords <- list()
  while(pts < tpts) {
    u <- sapply(1:p, function(j) sample(U[[j]], size = 1))
    ar <- in.hull(Y.hull, u, p)
    if(ar) {
      coords <- c(coords, list(u))
      pts <- pts + 1
      }
  }
  
  res <- t(simplify2array(coords))
  colnames(res) <- colnames(Y)
  res
}


hull.by.planes <- function(Y){
  Y <- as.matrix(Y)
  p <- ncol(Y)
  n <- nrow(Y)
  dim.list <- combn(p, 2)
  lapply(1:ncol(dim.list), function(j){
    dims <- dim.list[, j]
    sort(chull(Y[, dims]))
  })
}

hull.pts.by.planes <- function(Y){
  Y <- as.matrix(Y)
  p <- ncol(Y)
  n <- nrow(Y)
  dim.list <- combn(p, 2)
  hbp <-   lapply(1:ncol(dim.list), function(j){
    dims <- dim.list[, j]
    sort(chull(Y[, dims]))
  })
  lapply(1:length(hbp), function(j){
    dims <- dim.list[, j]
    y <- Y[, dims]
    y[hbp[[j]],]
  })
}

in.hull <- function(Y.hull, pt, p) { # assumes hull and Y.hull correspond
  if(length(pt) != p) stop("unequal dimensions between point and matrix")
  dim.list <- combn(p, 2)
  pt.i <- lapply(1:ncol(dim.list), function(j) pt[dim.list[,j]])
  for(i in 1:ncol(dim.list)){
    dims <- dim.list[, i]
    ch <- sort(chull(Y.hull[[i]]))
    chp <- sort(chull(rbind(Y.hull[[i]], pt.i[[i]])))
    res <- identical(ch, chp)
    if(!res) break
  }
  res
}

groups.at.points <- function(Y = Y, group = group, grid = grid) {
  pca <- prcomp(Y)
  d <- pca$sdev^2
  d <- d[which(zapsmall(d) > 0)]
  p <- length(d)
  P <- pca$x[,1:p]
  n <- nrow(P)
  glev <- levels(group)
  ng <- nlevels(group)
  gp.hull.pts.by.plane <- lapply(1:ng, function(j){
    gp.j <- which(group == glev[[j]])
    y <- P[gp.j,]
    hull.pts.by.planes(y)
  })
  g.a.p <- sapply(1:nrow(grid), function(j){
    g.point <- grid[j,]
    hd <- sapply(1:ng, function(j){
      y.hull <- gp.hull.pts.by.plane[[j]]
      in.hull(y.hull, g.point, p)
    })
    as.numeric(hd)
  })
  
  rownames(g.a.p) <- glev
  t(g.a.p)
}

grid.preview <- function(Y, pts = 500, pt.scale = 0.05){
  Z <- uniform.grid.sample(Y, pts, pt.scale)
  pairs(Z, pch = 19, cex = 0.3, asp = 1)
}
