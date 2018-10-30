#'
#' Plots estimated density function from an object of class 'mudens'.
#'
#' Plots estimated density function from an object of class 'mudens'.
#'
#' @description Add connected segments to a plot of estimated density function
#'              from an object of class 'mudens'.
#'
#' @usage
#' \method{lines}{mudens}(x,...)
#'
#' @param x Object of class mudens (output from calling mudens function)
#' @param ... Additional arguments for basic \code{lines(.)} function to be passed along.
#'
#' @details Plots estimated density function from an object of class 'mudens'.
#' Default time limits are those provided to mudens,
#'           which default to zero and the max observed time.
#'           Default y-axis limits are 0 and the maximum estimated density.
#'
#' @author Kenneth R. Hess
#'
#' @references Hess, K.R. and Zhong, M. Density Function Estimation for Possibly Right-Censored Data Using Kernel Functions. Submitted.
#'
#' @seealso \link{mudens}
#' @import graphics
#' @method lines mudens
#' @export
#'
#'
lines.mudens <- function (x, ...) {
    y <- x$haz.est
    if (x$pin$dens == 0) {
        lines(x$est.grid, y, type = "l", ylim = c(0, max(y)),
            xlab = "Follow-up Time", ylab = "Hazard Rate", ...)
    }
    else {
        lines(x$est.grid, y, type = "l", ylim = c(0, max(y)),
            xlab = "Follow-up Time", ylab = "Density", ...)
    }
    return(invisible())
}
