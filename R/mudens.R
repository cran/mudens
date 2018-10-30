#'
#' Estimate density function from right-censored survival data
#'
#' Estimate density function from a vector of right-censored survival times.
#'
#' @description Estimate density function from a vector of right-censored survival times using kernel functions.
#' Options include three types of bandwidth functions, three types of boundary correction, and four shapes
#' for the kernel function. Uses the global and local bandwidth selection algorithms and the boundary kernel
#' formulations described in Mueller and Wang (1994). The nearest neighbor bandwidth formulation is based
#' on that described in Gefeller and Dette (1992). The statistical properties of many of these estimators
#' are reported and compared in Hess et al. The \code{mudens(.)} function is an R wrapper around C code
#' and returns an object of class 'mudens' based on the density estimation in the HADES program developed
#' by H.G. Mueller.
#'
#' @usage mudens(times, delta, subset, min.time, max.time, bw.grid, bw.pilot,
#' bw.smooth, bw.method="local", b.cor="both", n.min.grid=51,
#'n.est.grid = 101, kern="epanechnikov")
#'
#' @param times A vector of survival times. It does not need to be sorted.
#' @param delta A vector indicating censoring: \code{0} - censored (alive), \code{1} - uncensored (dead).
#'              If delta is missing, all the observations are assumed uncensored.
#' @param subset A logical vector indicating the observations used in analysis.
#'               \code{TRUE} - observation is used, \code{FALSE} - observation is not used.
#'               If missing, all the observations will be used.
#' @param min.time Left bound of the time domain used in analysis. If missing, \code{min.time} is set to \code{0}.
#' @param max.time Right bound of the time domain used in analysis.
#'                 If missing, \code{max.time} is the maximum value of times.
#' @param bw.grid Bandwidth grid used in the MSE minimization.
#'                If \code{bw.method="global"} and \code{bw.grid} has one component only, no MSE minimization is performed.
#'                The hazard estimates are computed for the value of \code{bw.grid}.
#'                If \code{bw.grid} is missing, then a bandwidth grid of 21 components is built, having as bounds:
#'                \code{[0.2*bw.pilot, 20*bw.pilot]}
#' @param bw.pilot Pilot bandwidth used in the MSE minimization.
#'                 If missing, the default value is the one recommended by Mueller and Wang (1994):
#'
#'                 \code{bw.pilot = (max.time-min.time)/(8*nz^0.2)},
#'
#'                  where \code{nz} is the number of uncensored observations.
#' @param bw.smooth Bandwidth used in smoothing the local bandwidths. Not used if
#'
#'                  \code{bw.method="global"}.
#'                  If missing:   \code{bw.smooth=5*bw.pilot}.
#' @param bw.method Algorithm to be used. Possible values are: \code{"global"} - same bandwidth for all grid points.
#'                  In this case, the optimal bandwidth is obtained by minimizing the IMSE.
#'                  \code{"local"} - different bandwidths at each grid point, and the optimal bandwidth at a grid point
#'                  is obtained by minimizing the local MSE. \code{"knn"} - k nearest neighbors distance bandwidth,
#'                  and the optimal number of neighbors is obtained by minimizing the IMSE.
#'                  Note: The default value is \code{"local"}. Only the first letter needs to be given (e.g. "g", instead of \code{"global"}).
#' @param b.cor Boundary correction type. Possible values are:
#'              \code{"none"} - no boundary correction, \code{"left"} - left only correction,
#'              \code{"both"} - left and right corrections.
#'              The default value is set to \code{"both"}. Only the first letter needs to be given (e.g. b.cor="n").
#' @param n.min.grid Number of points in the minimization grid. This value greatly influences the computing time. Default value is \code{51}.
#' @param n.est.grid Number of points in the estimation grid, where hazard estimates are computed. Default value is \code{101}.
#' @param kern Boundary kernel function to be used. Possible values are:
#'            \code{"rectangle"}, \code{"epanechnikov"}, \code{"biquadratic"}, \code{"triquadratic"}.
#'            The default value is \code{"epanechnikov"}. Only the first letter needs to be given (e.g. kern="b").
#'
#' @details The mudens object contains a list of the input data and parameter values as well as a variety of output data.
#' The density function estimate is contained in the \code{haz.est} element and the corresponding time points are in \code{est.grid}.
#' The unsmoothed and smoothed local bandwidths are in \code{bw.loc} and \code{bw.loc.sm}, respectively.
#'
#' When setting \code{bw.method='local'} or \code{'knn'}, to check the shape of the bandwidth function used in the estimation,
#' use \code{plot(fit$pin$min.grid, fit$bw.loc)} to plot the unsmoothed bandwidths and
#' use \code{lines(fit$est.grid, fit$bw.loc.sm)} to superimpose the smoothed bandwidth function.
#' We can also use \code{bw.smooth} to change the amount of smoothing used on the bandwidth function.
#'
#' For \code{bw.method='global'}, use \code{plot(fit$bw.grid, fit$globlmse)} to check the minimization process, and
#' plot the estimated IMSE values over the bandwidth search grid; while for
#' \code{bw.method='k'}, use \code{plot(fit$k.grid, fit$k.imse)}.
#'
#' You may want to repeat the search using a finer grid over a shorter interval to fine-tune the optimization or if the observed
#' minimum is at the extreme of the grid you should specify a different grid.
#'
#' @return Returns an object of class \code{'mudens'}, containing input and output values.
#'         Methods working on such an object are:
#'         \code{plot}, \code{lines}, \code{summary}. For a detailed description of its components,
#'         see \code{object.mudens} in the mudens package.
#'
#' @references Hess, K.R. and Zhong, M. Density Function Estimation for Possibly Right-Censored Data Using Kernel Functions. Submitted.
#'
#' H.G. Mueller and J.L. Wang. Hazard Rates Estimation Under Random Censoring with Varying Kernels and Bandwidths. Biometrics 50:61-76, March, 1994.
#'
#' O.Gefeller and H. Dette. Nearest Neighbor Kernel Estimation of the Hazard Function From Censored Data. J. Statist. Comput. Simul., Vol.43:93-101, 1992.
#'
#' @author Kenneth R. Hess
#'
#' @importFrom  Rcpp sourceCpp
#' @useDynLib mudens, .registration = TRUE
#' @export mudens
#'
#' @examples
#' time <- rexp(1000)
#' stat <- sample(c(0,1), 1000, 0.5)
#' fit <- mudens(time, stat)
#' summary(fit)
#'
#'
mudens = function (times, delta, subset, min.time, max.time, bw.grid, bw.pilot, bw.smooth, bw.method = "local",
                   b.cor = "both", n.min.grid = 51, n.est.grid = 101, kern = "epanechnikov") {
    kernidx <- pmatch(kern, c("rectangle", "epanechnikov", "biquadratic",
        "triquadratic"))
    if (is.na(kernidx)) {
        stop("\nkern MUST be one of: 'rectangle', 'epanechnikov',\n               'biquadratic', or 'triquadratic'\n")
    }
    else if (kernidx == 1) {
        alpha <- 0
        beta <- 0
        nu <- 0
    }
    else if (kernidx == 2) {
        alpha <- 0
        beta <- 1
        nu <- 0
    }
    else if (kernidx == 3) {
        alpha <- 1
        beta <- 2
        nu <- 0
    }
    else if (kernidx == 4) {
        alpha <- 2
        beta <- 3
        nu <- 0
    }
    dens <- 1
    method <- pmatch(bw.method, c("global", "local", "knn"))
    if (is.na(method))
        stop("\nbw.method MUST be one of: 'global', 'local', or 'knn'\n")
    b.cor <- pmatch(b.cor, c("none", "left", "both"))
    if (is.na(b.cor))
        stop("\nb.cor MUST be one of: 'none', 'left', or 'both'\n")
    b.cor <- b.cor - 1
    if (missing(times))
        stop("\nParameter times is missing\n")
    nobs <- length(times)
    if (missing(delta)) {
        delta <- rep(1, nobs)
    }
    else {
        if (length(delta) != nobs) {
            stop("\ntimes and delta MUST have the same length\n")
        }
    }
    if (missing(subset)) {
        subset <- rep(T, nobs)
    }
    else {
        if (is.logical(subset)) {
            if (length(subset) != nobs) {
                stop("\ntimes and subset MUST have the same length\n")
            }
        }
        else {
            stop("\nsubset MUST contain ONLY logical values (T or F)\n")
        }
    }
    times <- times[subset]
    delta <- delta[subset]
    ix <- order(times)
    times <- times[ix]
    delta <- delta[ix]
    nobs <- length(times)
    if (missing(min.time)) {
        startz <- 0
    }
    else {
        if (min.time > times[1]) {
            warning("minimum time > minimum Survival Time\n\n")
        }
        startz <- min.time
    }
    if (missing(max.time)) {
        endz <- times[nobs]
    }
    else {
        if (max.time > times[nobs]) {
            warning("maximum time > maximum Survival Time\n\n")
        }
        endz <- max.time
    }
    if (startz > endz)
        stop("\nmin.time MUST be < max.time\n")
    nz <- sum(delta)
    if (missing(bw.pilot)) {
        bw.pilot <- endz/8/(nz^0.2)
    }
    if (missing(bw.smooth)) {
        bw.smooth <- 5 * bw.pilot
    }
    z <- seq(startz, endz, len = n.min.grid)
    zz <- seq(startz, endz, len = n.est.grid)
    endl <- startz
    endr <- endz
    pin.common <- list(times = times, delta = delta, nobs = nobs,
        min.time = startz, max.time = endz, n.min.grid = n.min.grid,
        min.grid = z, n.est.grid = n.est.grid, bw.pilot = bw.pilot,
        bw.smooth = bw.smooth, method = method, b.cor = b.cor,
        alpha = alpha, beta = beta, nu = nu, dens = dens)
    if (method < 3) {
        if (missing(bw.grid)) {
            bw.grid <- seq(0.2 * bw.pilot, 20 * bw.pilot, len = 25)
        }
        gridb <- length(bw.grid)
        m1 <- method - 1
        ans <- .C("newhad", as.integer(nobs), as.double(times),
            as.double(delta), as.integer(m1), as.double(z), as.integer(n.min.grid),
            as.integer(nu), as.integer(alpha), as.integer(beta),
            as.integer(dens), as.double(zz), as.integer(n.est.grid),
            as.double(bw.pilot), as.double(bw.grid), as.integer(gridb),
            as.double(endl), as.double(endr), as.double(bw.smooth),
            as.integer(b.cor), fzz = double(n.est.grid), bopt = double(n.min.grid),
            bopt1 = double(n.est.grid), msemin = double(n.min.grid),
            biasmin = double(n.min.grid), varmin = double(n.min.grid),
            imsemin = double(1), globlb = double(1), globlmse = double(gridb))
        if (method == 1) {
            ans <- list(pin = pin.common, est.grid = zz, haz.est = ans$fzz,
                imse.opt = ans$imsemin, bw.glob = ans$globlb,
                glob.imse = ans$globlmse, bw.grid = bw.grid)
        }
        else {
            ans <- list(pin = pin.common, est.grid = zz, haz.est = ans$fzz,
                bw.loc = ans$bopt, bw.loc.sm = ans$bopt1, msemin = ans$msemin,
                bias.min = ans$biasmin, var.min = ans$varmin,
                imse.opt = ans$imsemin)
        }
    }
    else if (method == 3) {
        kmin <- 2
        kmax <- as.integer(nz/2)
        m1 <- 2
        ans <- .C("knnhad", as.integer(nobs), as.double(times),
            as.double(delta), as.integer(nu), as.integer(alpha),
            as.integer(beta), as.integer(dens), as.integer(m1),
            as.integer(n.min.grid), as.double(z), as.integer(n.est.grid),
            as.double(zz), as.double(bw.pilot), as.double(endl),
            as.double(endr), as.double(bw.smooth), as.integer(b.cor),
            fzz = double(n.est.grid), kopt = as.integer(kmin),
            as.integer(kmax), bopt = double(n.min.grid), bopt1 = double(n.est.grid),
            kimse = double(kmax - kmin + 1))
        ans <- list(pin = pin.common, est.grid = zz, haz.est = ans$fzz,
            bw.loc = ans$bopt, bw.loc.sm = ans$bopt1, k.grid = kmin:kmax,
            k.imse = ans$kimse, imse.opt = ans$kimse[ans$kopt -
                kmin + 1], kopt = ans$kopt)
    }
    else {
        stop("\nmethod MUST be 1, 2, or 3\n")
    }
    class(ans) <- "mudens"
    ans
}
