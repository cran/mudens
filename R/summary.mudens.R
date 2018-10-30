#'
#' Display the most important input parameters used in calling the 'mudens' function
#'
#' Display the most important input parameters used in calling the 'mudens' function
#'
#' @description It also displays some of the output data. Common to all three methods:
#'              (1) number of observations, (2) number of censored observations, (3) bandwidth method
#'              used (global, local or nearest neighbor), (4) boundary correction type (none, left only,
#'              both left and right), (5) kernel type (rectangle, Epanechnikov, biquadradic, triquadratic),
#'              (6) minimum time, (7) maximum time, (8) number of points in MSE minimization grid,
#'              (9) number of points in estimation grid, (10) pilot bandwidth,
#'              (11) estimated IMSE for optimal bandwidth.
#'
#'  This function will also report the following two results
#'  for different selections of methods: (12) Smoothing Bandwidth and
#'  (13) Optimal Number of Nearest Neighbor.
#'  If \code{bw.method="global"}, this function will display the optimal global bandwidth.
#'  If \code{bw.method="knn"}, it will show optimal number of nearest neighbors.
#'  If \code{bw.method="local"} and \code{bw.method="knn"}, the
#'  summary result will include the smoothing bandwidth used to smooth the optimal local bandwidths.
#'
#' @usage
#' \method{summary}{mudens}(object, ...)
#'
#' @param object object of class \code{mudens} (output from calling \code{mudens(.)} function)
#' @param ... Additional arguments to be passed along.
#'
#' @return The summary result for mudens estimation
#'
#' @references Hess, K.R. and Zhong, M. Density Function Estimation for Possibly Right-Censored Data Using Kernel Functions. Submitted.
#'
#' @author Kenneth R. Hess
#'
#' @seealso \link{mudens}
#' @import survival
#'
#' @method summary mudens
#' @export
summary.mudens <- function (object, ...) {
if (object$pin$dens == 0) {
  cat("\n\t\tHazard Function Estimates")
}
else {
  cat("\n\t\tDensity Estimates")
}
cat("\nNumber of Observations ..........", object$pin$nobs)
cat("\nCensored Observations ...........", object$pin$nobs - sum(object$pin$delta))
cat("\nMethod used ..................... ")
y <- object$pin$method
if (y == 1) {
  cat("Global")
}
else if (y == 2) {
  cat("Local")
}
else if (y == 3) {
  cat("Nearest Neighbor")
}
else {
  cat("Unknown method")
}
cat("\nBoundary Correction Type ........ ")
y <- object$pin$b.cor
if (y == 0) {
  cat("None")
}
else if (y == 1) {
  cat("Left Only")
}
else if (y == 2) {
  cat("Left and Right")
}
else {
  cat("Unknown")
}
cat("\nkernel .......................... ")
alpha <- object$pin$alpha
beta <- object$pin$beta
if (alpha == 0 && beta == 0) {
  cat("rectangle")
}
else if (alpha == 0 && beta == 1) {
  cat("Epanechnikov")
}
else if (alpha == 1 && beta == 2) {
  cat("biquadratic")
}
else if (alpha == 2 && beta == 3) {
  cat("triquadratic")
}
cat("\nMinimum Time ....................", round(object$pin$min.time,
                                                 2))
cat("\nMaximum Time ....................", round(object$pin$max.time,
                                                 2))
cat("\nNumber of minimization points ...", object$pin$n.min.grid)
cat("\nNumber of estimation points .....", object$pin$n.est.grid)
cat("\nPilot Bandwidth .................", round(object$pin$bw.pilot,
                                                 2))
cat("\nSmoothing Bandwidth .............", round(object$pin$bw.smooth,
                                                 2))
y <- object$pin$method
if (y == 1) {
  cat("\nOptimal Global Bandwidth ........", round(object$bw.glob,
                                                   2))
}
else if (y == 3) {
  cat("\nOptimal Nearest Neighbor ........", object$kopt)
}
cat("\nMinimum IMSE ....................", round(object$imse.opt,
                                                 2))
cat("\n")
return(invisible())
}
#'
#' @examples
#' data(ovarian,package="survival")
#' attach(ovarian)
#' fit <- mudens(futime, fustat)
#' summary.mudens(fit)
#'
