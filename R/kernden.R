#'
#' Calculate global bandwidth kernel estimates of density function
#'
#' Calculate global bandwidth kernel estimates of density function for survival times
#'
#' @description Estimates global bandwidth kernel from right-censored data using the Epanechnikov kernel described in Silverman BW (1986).
#'
#' @usage kernden(times, status, estgrid, bandwidth)
#'
#' @param times A vector of survival times. It does not need to be sorted.
#' @param status A vector indicating censoring: 0 - censored (alive), 1 - uncensored (dead).
#'               If status is missing, all the observations are assumed uncensored.
#' @param estgrid A vector of time points at which the estimation will be made.
#' @param bandwidth Bandwidth used to determine the degree of smoothing.
#'                  Larger values of bandwidth will result in smoother of mean estimates.
#'                  It is suggested to start with a value of approximately 20\% of the range
#'                  of the survival times.
#'
#' @return Returns an object containing the time points of estimations (estgrid) and corresponding density estimates (denest).
#'
#' @references Hess, K.R. and Zhong, M. Density Function Estimation for Possibly Right-Censored Data Using Kernel Functions. Submitted.
#'
#' @import survival
#' @author R. Herrick and Dan Serachitopol
#'
#'
kernden = function (times, status, estgrid, bandwidth) {
    e <- length(estgrid)
    denest <- rep(NA, e)
    sfit <- survfit(Surv(times, status))
    s <- 0 - diff(c(1, sfit$surv))
    nt <- length(s)
    udt <- sfit$time
    for (i in 1:e) {
        sum <- 0
        x0 <- estgrid[i]
        for (j in 1:nt) {
            arg <- (udt[j] - x0)/bandwidth
            arg <- arg^2
            if (arg < 1)
                sum <- sum + (1 - arg) * s[j]
        }
        denest[i] <- sum * 0.75/bandwidth
    }
    return(list(estgrid, denest))
}
