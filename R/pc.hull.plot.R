#' Principal Component Plot of Convex Hulls
#'
#' Function generates a principal component plot to best view convex hulls in low dimensions.
#' 
#' Principal components are obtained from eigen analysis of an among-group
#' covariance or correlation matrix.  (Correlation matrix is used if std = TRUE.)
#' Typical plot parameters should be available but parameters for hulls are controlled
#' with specific function arguments.  One can add points, lines, arrows, or a legend 
#' to the plot canvas, among other plot tricks, if desired.  (See \code{\link{points}},
#' \code{\link{lines}}, \code{\link{arrows}}, and \code{\link{legend}}).
#' 
#' @param dat A matrix or data frame of data
#' @param std A logical value that if TRUE finds standard deviates of the data 
#' (data are both centered and scaled by variable standard deviations).
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
#' pc.hull.plot(Pupfish$coords, group = group)
pc.hull.plot <- function(dat, std = FALSE, group = NULL, group.cols = NULL, 
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
  plot(Y, asp = 1,
       xlab = paste("PC1", round(d[1]/sum(d)*100, 1), "%"),
       ylab = paste("PC1", round(d[2]/sum(d)*100, 1), "%"),
       ...)
  hull.pts <- lapply(1:ng, function(j){
    g <- which(group == glev[j])
    Z <- Y[g,]
    hp <- chull(Z)
    hp <- c(hp, hp[1])
    Z[hp, ]
  })
  
  for(i in 1:ng){
    y <- hull.pts[[i]]
    points(y[,1], y[,2], type = "l", lty = group.lty[i],
           lwd = group.lwd[i], col = group.cols[i])
  }
}