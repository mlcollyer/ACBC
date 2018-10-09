#' Principal Component Plot with Ellipsoids
#'
#' Function generates a principal component plot to best view ellipsoids in low dimensions.
#' 
#' Principal components are obtained from eigen analysis of an among-group
#' covariance or correlation matrix.  (Correlation matrix is used if std = TRUE.)
#' Typical plot parameters should be available but parameters for ellipsoids are controlled
#' with specific function arguments.  One can add points, lines, arrows, or a legend 
#' to the plot canvas, among other plot tricks, if desired.  (See \code{\link{points}},
#' \code{\link{lines}}, \code{\link{arrows}}, and \code{\link{legend}}).
#' 
#' Note that the PC plot is 2-dimensional, so the ellipsoids are ellipses.  However, the
#' ellipses are 2-dimensional projections of ellipsoids in the full data space.  If original
#' data are 2-dimensional, the plot simply rotates the data and ellispes to best align
#' with variation among grup centroids (means).
#' 
#' @param dat A matrix or data frame of data
#' @param std A logical value that if TRUE finds standard deviates of the data 
#' (data are both centered and scaled by variable standard deviations).
#' @param PC A vector of length 2 to indicate which PCs to view.  The default is c(1, 2).
#' This can be changed to view alternative dimensions.  Illogical requests will defualt to
#' c(1, 2) with a warning
#' @param confidence The confidence level for ellipsoids, based on the covariance matrix of 
#' data by groups.  Multivariate normality is assumed in estimation.  If NULL, then ellipsoids
#' merely reflect the span of eigenvalues for the covariance matrix of the data.  Otherwise, the value should be between
#' 0.01 and 1.
#' @param ellipse.density A numeric value to indicate how many discrete points (in a circle)
#' are used to approximate the continuous ellipse function.  More points mean a more precise curve, 
#' but increase computation time.  The default, 120 points, is the same as 3 degrees (pi/60 radians) increments.
#' @param group A factor of vector coercible to factor.  If null, the a single convex
#' hull with be returned.
#' @param group.cols A optional vector with length equal to the number of group levels
#' to describe the colors of the hulls.  If NULL, R standard colors will be used.
#' @param group.lwd A optional vector with length equal to the number of group levels
#' to describe the line weight (magnification) of the hulls.  If NULL, 
#' group.lwd = 1 will be used.
#' @param group.lty A optional vector with length equal to the number of group levels
#' to describe the line type of the hulls.  If NULL, 
#' group.lty = 1 (solid line) will be used.
#' @param ... other arguments for graphical parameters
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @examples 
#' library(RRPP)
#' data("Pupfish")
#' group <- interaction(Pupfish$Sex, Pupfish$Pop)
#' pc.ellipse.plot(Pupfish$coords, group = group)
pc.ellipse.plot <- function(dat, std = FALSE, PC = c(1,2), confidence = NULL,
                            ellipse.density = 120,
                            group = NULL, group.cols = NULL, 
                         group.lwd = NULL, group.lty = NULL, ...){
  if(NCOL(dat) < 2) stop("Cannot generate hulls in fewer than 2 dimensions")
  n <- NROW(dat)
  Y <- as.matrix(dat)
  if(std) Y <- scale(Y) 
  if(is.null(group)) group <- rep(1, n)
  group <- as.factor(group)
  ng <- nlevels(group)
  glev <- levels(group)
  if(is.null(group.cols)) group.cols <- 1:ng
  if(is.null(group.lwd)) group.lwd <- rep(1, ng)
  if(is.null(group.lty)) group.lty <- rep(1, ng)
  if(length(group.cols) != ng) stop("Number of requested group colors does not match the number of groups")
  if(length(group.lwd) != ng) stop("Number of requested group widths does not match the number of groups")
  if(length(group.lty) != ng) stop("Number of requested group line types does not match the number of groups")
  
  if(ng == 1) formula <- Y ~ 1 else 
    formula <- Y ~ group
  dat <- data.frame(Y=Y, group = group)
  X <- model.matrix(formula, data = dat)
  fitted <- lm.fit(X, Y)$fitted.values
  pca <- prcomp(fitted)
  d <- pca$sdev^2
  k <- which(zapsmall(d) > 0)
  Rot <- pca$rotation[, k]
  Y <- Y %*% Rot
  pc <- match(PC, k)
  if(any(is.na(pc))) {
    cat("\nIllogical PC dimensions chosen; thus the first two are shown\n")
    cat("These are the possible dimensions you could chose:", k, "\n")
    pc <- c(1, 2)
  }
  P <- Y[, pc]
  plot(P, asp = 1,
       xlab = paste("PC1", round(d[1]/sum(d)*100, 1), "%"),
       ylab = paste("PC1", round(d[2]/sum(d)*100, 1), "%"),
       ...)
  
  ellipse.pts <- lapply(1:ng, function(j){
    g <- which(group == glev[j])
    Yp <- P[g,] 
    ellipsoid.pts.by.planes(Yp, confidence, ellipse.density)
  })
  
  for(i in 1:ng){
    y <- ellipse.pts[[i]][[1]]
    points(y[,1], y[,2], type = "l", lty = group.lty[i],
           lwd = group.lwd[i], col = group.cols[i])
  }
}