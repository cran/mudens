#'
#' Plots estimated density function from an object of class 'mudens'.
#'
#' Plots estimated density function from an object of class 'mudens'.
#'
#' @description Plot estimated density function from an object of class 'mudens'.
#'
#' @usage
#' \method{plot}{mudens}(x,...)
#'
#' @param x Object of class mudens (output from calling mudens function)
#' @param ... Additional arguments of basic \code{plot(.)} function to be passed along.
#'
#' @details Default time limits are those provided to mudens,
#' which default to zero and the time corresponding to
#' the maximum value of the times. Default y-axis limits are 0 and
#' the maximum estimated density.
#'
#' @references Hess, K.R. and Zhong, M. Density Function Estimation for Possibly Right-Censored Data Using Kernel Functions. Submitted.
#'
#' @author Kenneth R. Hess
#'
#' @seealso \link{mudens}
#' @import graphics
#' @method plot mudens
#' @export
#'
plot.mudens <- function (x, ...) {
    y <- x$haz.est
    if (x$pin$dens == 0) {
        plot(x$est.grid, y, type = "l", ylim = c(0, max(y)),
            xlab = "Follow-up Time", ylab = "Hazard Rate", ...)
    }
    else {
        plot(x$est.grid, y, type = "l", ylim = c(0, max(y)),
            xlab = "Follow-up Time", ylab = "Density", ...)
    }
    return(invisible())
}
