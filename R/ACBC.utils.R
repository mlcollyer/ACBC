## achc

#' Print/Summary Function for ACHC
#'
#' @param x print/summary object (from \code{\link{achc}})
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.achc <- function(x, ...){
  grid <- x$grid
  n <- nrow(grid)
  group <- x$group
  ng <- nlevels(group)
  ref <- seq(0, ng, 1)
  ref.mat <- rm <- matrix(0, n, length(ref), byrow = TRUE)
  analysis <- x$analysis
  obs <- analysis[[1]]
  obs.sum <- rowSums(obs)
  for(i in 1:n) {
    loc <- match(obs.sum[i], ref)
    rm[i, loc] <- 1
  }
  nrh <- colSums(rm)/n
  names(nrh) <- paste(0:(length(nrh) - 1), "hulls", sep=".")
  cat("\nRelative frequencies of hull coverage:\n\n")
  print(nrh)
  cat("\n")
}

#' Print/Summary Function for ACHC
#'
#' @param object print/summary object (from \code{\link{achc}})
#' @param confidence Confidence level for ACHC test
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.achc <- function(object, confidence = 0.95, ...){
  x <- object
  if(confidence <= 0) stop("Confidence level must be between 0 and 1")
  if(confidence > 1) stop("Confidence level must be between 0 and 1")
  alpha <- 1 - confidence
  grid <- x$grid
  n <- nrow(grid)
  group <- x$group
  ng <- nlevels(group)
  ref <- seq(0, ng, 1)
  ref.mat <- rm <- matrix(0, n, length(ref), byrow = TRUE)
  analysis <- x$analysis
  obs <- analysis[[1]]
  obs.sum <- rowSums(obs)
  for(i in 1:n) {
    loc <- match(obs.sum[i], ref)
    rm[i, loc] <- 1
  }
  nrh <- colSums(rm)/n
  names(nrh) <- paste(0:(length(nrh) - 1), "hulls", sep=".")
  perms <- length(analysis)
  rf <- sapply(1:perms, function(j){
    rs <- rowSums(analysis[[j]])
    rm <- ref.mat
    for(i in 1:n) {
      loc <- match(rs[i], ref)
      rm[i, loc] <- 1
    }
    res <- colSums(rm)/n
    names(res) <- ref
    res
  })
  
  lcl <- apply(rf, 1, function(x) quantile(x, alpha/2))
  ucl <- apply(rf, 1, function(x) quantile(x, (1 - alpha/2)))
  
  
  df <- data.frame(obs = nrh, lcl = lcl, ucl = ucl)
  colnames(df)[2:3] <- paste(c(alpha/2 *100, (1-alpha/2)*100),
                             c("%lcl", "%ucl"), sep = "")
  cat("\nObserved distribution of relative frequencies\n",
      "and", confidence*100, "% confidence intervals\n\n")
  print(df)
  invisible(df)
}

#' Plot Function for ACHC
#' 
#' @param x plot object (from \code{\link{achc}})
#' @param confidence Confidence level for ACHC test
#' @param ... other arguments passed to plot parameters (helpful to employ
#' different colors or symbols for different groups).  See
#' \code{\link{plot.default}} and \code{\link{par}}.  Some limitations
#' occur because of polygon use.
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.achc <- function(x, confidence = 0.95, ...){
  df <- summary.achc(x, confidence = confidence)
  rfmax <- apply(df, 1, max)
  xvar <- seq(0, (nrow(df) - 1), 1)
  plot(xvar, rfmax, type = "n", ylim = c(0, max(rfmax)),
       xlab = "Number of hulls", ylab = "Relative frequency", ...)
  
  pxvar <- c(xvar, xvar[order(xvar, decreasing = TRUE)])
  y3 <- df[,3]
  y2 <- df[,2]
  pyvar <- c(y3, y2[order(xvar, decreasing = TRUE)])
  polygon(pxvar, pyvar, col = "gray", border = "gray", ...)
  points(xvar, df[,1], type = "l", ...)
}

## aec

#' Print/Summary Function for AEC
#'
#' @param x print/summary object (from \code{\link{achc}})
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
print.aec <- function(x, ...){
  grid <- x$grid
  n <- nrow(grid)
  group <- x$group
  ng <- nlevels(group)
  ref <- seq(0, ng, 1)
  ref.mat <- rm <- matrix(0, n, length(ref), byrow = TRUE)
  analysis <- x$analysis
  obs <- analysis[[1]]
  obs.sum <- rowSums(obs)
  for(i in 1:n) {
    loc <- match(obs.sum[i], ref)
    rm[i, loc] <- 1
  }
  nrh <- colSums(rm)/n
  names(nrh) <- paste(0:(length(nrh) - 1), "hulls", sep=".")
  cat("\nRelative frequencies of hull coverage:\n\n")
  print(nrh)
  cat("\n")
}

#' Print/Summary Function for AEC
#'
#' @param object print/summary object (from \code{\link{achc}})
#' @param confidence Confidence level for AEC test
#' @param ... other arguments passed to print/summary
#' @export
#' @author Michael Collyer
#' @keywords utilities
summary.aec <- function(object, confidence = 0.95, ...){
  x <- object
  if(confidence <= 0) stop("Confidence level must be between 0 and 1")
  if(confidence > 1) stop("Confidence level must be between 0 and 1")
  alpha <- 1 - confidence
  grid <- x$grid
  n <- nrow(grid)
  group <- x$group
  ng <- nlevels(group)
  ref <- seq(0, ng, 1)
  ref.mat <- rm <- matrix(0, n, length(ref), byrow = TRUE)
  analysis <- x$analysis
  obs <- analysis[[1]]
  obs.sum <- rowSums(obs)
  for(i in 1:n) {
    loc <- match(obs.sum[i], ref)
    rm[i, loc] <- 1
  }
  nrh <- colSums(rm)/n
  names(nrh) <- paste(0:(length(nrh) - 1), "hulls", sep=".")
  perms <- length(analysis)
  rf <- sapply(1:perms, function(j){
    rs <- rowSums(analysis[[j]])
    rm <- ref.mat
    for(i in 1:n) {
      loc <- match(rs[i], ref)
      rm[i, loc] <- 1
    }
    res <- colSums(rm)/n
    names(res) <- ref
    res
  })
  
  lcl <- apply(rf, 1, function(x) quantile(x, alpha/2))
  ucl <- apply(rf, 1, function(x) quantile(x, (1 - alpha/2)))
  
  
  df <- data.frame(obs = nrh, lcl = lcl, ucl = ucl)
  colnames(df)[2:3] <- paste(c(alpha/2 *100, (1-alpha/2)*100),
                             c("%lcl", "%ucl"), sep = "")
  cat("\nObserved distribution of relative frequencies\n",
      "and", confidence*100, "% confidence intervals\n\n")
  print(df)
  invisible(df)
}

#' Plot Function for AEC
#' 
#' @param x plot object (from \code{\link{achc}})
#' @param confidence Confidence level for AEC test
#' @param ... other arguments passed to plot parameters (helpful to employ
#' different colors or symbols for different groups).  See
#' \code{\link{plot.default}} and \code{\link{par}}.  Some limitations
#' occur because of polygon use.
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @keywords visualization
plot.aec <- function(x, confidence = 0.95, ...){
  df <- summary.achc(x, confidence = confidence)
  rfmax <- apply(df, 1, max)
  xvar <- seq(0, (nrow(df) - 1), 1)
  plot(xvar, rfmax, type = "n", ylim = c(0, max(rfmax)),
       xlab = "Number of hulls", ylab = "Relative frequency", ...)
  
  pxvar <- c(xvar, xvar[order(xvar, decreasing = TRUE)])
  y3 <- df[,3]
  y2 <- df[,2]
  pyvar <- c(y3, y2[order(xvar, decreasing = TRUE)])
  polygon(pxvar, pyvar, col = "gray", border = "gray", ...)
  points(xvar, df[,1], type = "l", ...)
}



