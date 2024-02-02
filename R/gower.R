
is.wholenumber = function(x, tol = .Machine$double.eps^0.5)
{
    abs(x - round(x)) < tol
}

#' Compute a highest posterior density (HPD) interval.
#'
#' @details This function uses a given posterior sample to compute an HPD interval at a given significance level.
#'
#' @param x the posterior sample.
#' @param alpha the desired significance level.
#'
#' @return A 2-vector containing the lower endpoint and the upper endpoint, respectively.
#'
#' @export

hpd = function(x, alpha = 0.05)
{
    n = length(x)
    m = round(n * alpha)
    x = sort(x)
    y = x[(n - m + 1):n] - x[1:m]
    z = min(y)
    k = which(y == z)[1]
    c(x[k], x[n - m + k])
}

#' Apply the discrete metric to two scores.
#'
#' @details This function applies the discrete metric to two scores. This is an appropriate distance function for nominal data.
#'
#' @param x a score.
#' @param y a score.
#'
#' @return 0 if \code{x} is equal to \code{y}, 1 if \code{x} is not equal to \code{y}, \code{NA} if either score is \code{NA}.
#'
#' @seealso \code{\link{L1.dist}}
#'
#' @export

discrete.dist = function(x, y)
{
    d = as.numeric(x != y)
    d
}

#' Apply the L1 distance function to two scores.
#'
#' @details This function applies the L1 distance function to two scores. This is an appropriate distance function for ordinal data.
#'
#' @param x a score.
#' @param y a score.
#' @param range the range of the scores, i.e., the difference between the maximum score and the minimum score.
#'
#' @return \eqn{|x - y| / r}, where \eqn{r} is the range of the scores; or \code{NA} if either score is \code{NA}.
#'
#' @seealso \code{\link{discrete.dist}}
#'
#' @export

L1.dist = function(x, y, range)
{
    d = as.numeric(abs(x - y) / range)
    d
}

# The following function computes all pairwise distances for a given row, i.e., unit, in the dataset.
# The distances are returned in a vector. NAs are removed.

row.dist = function(x, dist, ...)
{
    d = NULL
    n = length(x)
    for (i in 1:(n - 1))
        for (j in (i + 1):n)
            d = c(d, dist(x[i], x[j], ...))
    d = d[! is.na(d)]
    d
}

# The following function computes all of the row statistics for a given dataset, using the mean distance for
# each row. The row statistics are returned in a vector. NAs are removed. This function is used when the data
# are nominal, or when the data are ordinal and dist.type is equal to "mean".

row.stat = function(data, dist, ...)
{
    n = nrow(data)
    g = numeric(n)
    for (i in 1:n)
        g[i] = 1 - mean(row.dist(data[i, ], dist, ...))
    g = g[! is.na(g)]
    g
}

# The following function computes all of the row statistics for a given dataset, using the maximum distance for
# each row. The row statistics are returned in a vector. NAs are removed. This function is used when the data
# are ordinal and dist.type is equal to "max".

row.stat.max = function(data, dist, ...)
{
    n = nrow(data)
    g = numeric(n)
    for (i in 1:n)
        g[i] = 1 - max(row.dist(data[i, ], dist, ...))
    g = g[! is.na(g)]
    g
}
    
#' Apply the Bayesian Gower agreement methodology to nominal or ordinal data.
#'
#' @details This is the package's flagship function. It applies the Bayesian Gower methodology to nominal or ordinal data, and provides an estimate of the posterior mean along with a credible interval.
#'
#' 
#'
#' @param data a matrix of scores. Each row corresponds to a unit, each column to a coder.
#' @param data.type the type of scores to be analyzed, either \code{"nominal"} or \code{"ordinal"}.
#' @param dist.type for ordinal data, whether the row statistics are computed using the mean of the pairwise distances or the maximum pairwise distance.
#' @param design the sampling design, either \code{"one-way"} or \code{"two-way"}. For the former, the units are random and the coders are fixed. For the latter, both the units and the coders are random.
#' @param iter the desired size of the posterior sample. The default value is 10,000.
#' @param \dots additional arguments for the distance function. These are ignored for nominal data. For ordinal data the range of the scores must be provided via argument \code{range}.  
#'
#' @return Function \code{gower.agree} returns an object of class \code{"gower"}, which is a list comprising the following elements.
#'         \item{mu.hat}{the estimate of the posterior mean.}
#'         \item{mu.sample}{the posterior sample.}
#'         \item{call}{the matched call.}
#'         \item{units}{the number of units.}
#'         \item{coders}{the number of coders.}
#'         \item{data}{the data matrix, sans rows that were removed due to missigness.}
#'         \item{data.type}{the type of scores, nominal or ordinal.}
#'         \item{dist.type}{for ordinal data, the manner in which the row statistics were computed.}
#'         \item{design}{the sampling design, one-way or two-way.}
#'         \item{row.stats}{the vector of row statistics.}
#'         \item{del}{the number of rows that were deleted due to missingness.}
#'
#' @export
#'
#' @examples
#'
#' # Fit the liver data, using the mean distance for each row of the data matrix.
#' # The range (which is equal to 4) must be passed to \code{\link{gower.agree}}
#' # since these data are ordinal and the L1 distance function is used. We assume
#' # a one-way sampling design for these data, i.e., units are random and coders
#' # are fixed.
#'
#' data(liver)
#' liver = as.matrix(liver)
#' fit = gower.agree(liver, data.type = "ordinal", range = 4)
#' summary(fit)

gower.agree = function(data, data.type = c("nominal", "ordinal"), dist.type = c("mean", "max"),
                       design = c("one-way", "two-way"), iter = 10000, ...)
{
    call = match.call()
    if (missing(data) || ! is.matrix(data) || ! is.numeric(data))
        stop("You must supply a numeric data matrix.")
    n = nrow(data)
    m = ncol(data)
    del = NULL
    for (i in 1:n)
    {
        present = sum(! is.na(data[i, ]))
        if (present < 2)
            del = c(del, i)
    }
    if (length(del) > 0)
        data = data[-del, ]
    n = nrow(data)
    data.type = match.arg(data.type)
    if (data.type == "nominal")
        dist = discrete.dist
    else
        dist = L1.dist
    dist.type = match.arg(dist.type)
    design = match.arg(design)
    if (! is.wholenumber(iter) || iter < 1)
    {
        warning("Argument 'iter' must be a positive integer. Using the default value of 10,000.")
        iter = 10000
    }
    if (design == "one-way")
    {
        if (dist.type == "mean")
            g = row.stat(data, dist, ...)
        else
            g = row.stat.max(data, dist, ...)
        mu.sample = numeric(iter)
        for (j in 1:iter)
        {
            u = runif(n - 1)
            w = diff(c(0, sort(u), 1))
            mu.sample[j] = w %*% g
        }
    }
    else
    {
        mu.sample = numeric(iter)
        for (j in 1:iter)
        {
            rows = sample(1:n, n)
            dat = data[rows, ]
            cols = sample(1:m, m)
            dat = dat[, cols]
            if (dist.type == "mean")
                g = row.stat(data, dist, ...)
            else
                g = row.stat.max(data, dist, ...)
            u = runif(n - 1)
            w = diff(c(0, sort(u), 1))
            mu.sample[j] = w %*% g
        }
        if (dist.type == "mean")
            g = row.stat(data, dist, ...)
        else
            g = row.stat.max(data, dist, ...)
    }
    mu.hat = mean(mu.sample)
    names(mu.hat) = "mu"
    object = list(call = call, units = n, coders = m, data = data, data.type = data.type, row.stats = g,
                  dist.type = dist.type, design = design, mu.hat = mu.hat, mu.sample = mu.sample, del = length(del))
    class(object) = "gower"
    object
}

#' Compute DFBETAs for units and/or coders.
#'
#' @details This function computes DFBETAs for one or more units and/or one or more coders.
#'
#' @param model a fitted model object, the result of a call to \code{\link{gower.agree}}.
#' @param units a vector of integers. A DFBETA will be computed for each of the corresponding units.
#' @param coders a vector of integers. A DFBETA will be computed for each of the corresponding coders.
#' @param ... additional arguments. These are ignored.
#
#' @return A list comprising at most four elements.
#'         \item{dfbeta.units}{a vector containing DFBETAs for the units specified via argument \code{units}.}
#'         \item{dfbeta.coders}{a vector containing DFBETAs for the coders specified via argument \code{coders}.}
#'         \item{fits.units}{a list containing fit objects for the omitted units specified via argument \code{units}.}
#'         \item{fits.coders}{a list containing fit objects for the omitted coders specified via argument \code{coders}.}
#'
#' @method influence gower
#'
#' @export
#'
#' @examples
#'
#' # Analyze nominal data previously considered by Krippendorff.
#' # Assume a one-way design. Compute a DFBETA for unit 6, which
#' # should be rather influential.
#'
#' kripp = matrix(c(1,2,3,3,2,1,4,1,2,NA,NA,NA,
#'                  1,2,3,3,2,2,4,1,2,5,NA,3,
#'                  NA,3,3,3,2,3,4,2,2,5,1,NA,
#'                  1,2,3,3,2,4,4,1,2,5,1,NA), 12, 4)
#' kripp
#'
#' set.seed(12)
#' fit = gower.agree(kripp)
#' summary(fit)
#' inf = influence(fit, units = 6)
#'
#' # Report the DFBETA for unit 6 and the estimate of mu when unit 6
#' # is ommitted.
#'
#' inf$dfbeta.units
#' fit$mu.hat - inf$dfbeta.units

influence.gower = function(model, units, coders, ...)
{
    mod.call = model$call
    dfbeta.units = NULL
    fits.units = NULL
    if (! missing(units))
    {
        dfbeta.units = numeric(length(units))
        k = 1
        for (j in units)
        {
            mod.call$data = model$data[-j, ]
            fit = eval(mod.call)
            fits.units = c(fits.units, fit)
            dfbeta.units[k] = model$mu.hat - fit$mu.hat
            k = k + 1
        }
        names(dfbeta.units) = units
    }
    dfbeta.coders = NULL
    fits.coders = NULL
    if (! missing(coders))
    {
        dfbeta.coders = numeric(length(coders))
        k = 1
        for (j in coders)
        {
            mod.call$data = model$data[, -j]
            fit = eval(mod.call)
            fits.coders = c(fits.coders, fit)
            dfbeta.coders[k] = model$mu.hat - fit$mu.hat
            k = k + 1
        }
        names(dfbeta.coders) = coders
    }
    result = list()
    if (! is.null(dfbeta.units))
    {
        result$dfbeta.units = dfbeta.units
        result$fits.units = fits.units
    }
    if (! is.null(dfbeta.coders))
    {
        result$dfbeta.coders = dfbeta.coders
        result$fits.coders = fits.coders
    }
    result        
}

#' Print a summary of a Bayesian Gower fit.
#'
#' @details This function prints a summary of the fit.
#'
#' @param object an object of class \code{"gower"}, the result of a call to \code{\link{gower.agree}}.
#' @param conf.level the confidence level for the credible interval. The default is 0.95.
#' @param digits the number of significant digits to display. The default is 4.
#' @param \dots additional arguments. These are currently ignored.
#'
#' @return A list containing various elements of the summary is returned invisibly, but this function should be called for its side effects.
#'
#' @seealso \code{\link{gower.agree}}
#'
#' @method summary gower
#'
#' @export
#'
#' @examples
#'
#' # Fit the liver data, using the mean distance for each row of the data matrix.
#' # The range (which is equal to 4) must be passed to \code{\link{gower.agree}}
#' # since these data are ordinal and the L1 distance function is used. We assume
#' # a one-way sampling design for these data, i.e., units are random and coders
#' # are fixed.
#'
#' data(liver)
#' liver = as.matrix(liver)
#' fit = gower.agree(liver, data.type = "ordinal", range = 4)
#' summary(fit)

summary.gower = function(object, conf.level = 0.95, digits = 4, ...)
{
    cat("\nBayesian Gower agreement\n========================\n")
    cat("\nData:", object$units, "units x", object$coders, "coders\n")
    ans = NULL
    ans$units = object$units
    ans$coders = object$coders
    if (object$del > 0)
        cat("\n", object$del, " rows removed due to missingness\n", sep = "")
    cat("\nData type:", object$data.type, "\n")
    ans$data.type = object$data.type
    if (object$data.type == "ordinal")
    {
        cat("\nDistance type:", object$dist.type, "\n")
        ans$dist.type = object$dist.type
    }
    cat("\nSampling design:", object$design, "\n")
    ans$design = object$design
    cat("\nPosterior sample size:", length(object$mu.sample), "\n")
    cat("\nCall:\n\n")
    print(object$call)
    ans$call = object$call
    coef.table = cbind(object$mu.hat, t(confint(object, level = conf.level)))
    colnames(coef.table) = c("Estimate", "Lower", "Upper")
    rownames(coef.table) = names(object$mu.hat)
    ans$coef.table = coef.table
    cat("\nResults:\n\n")
    print(signif(coef.table, digits = digits))
    cat("\n")
    invisible(ans)
}

#' Compute a credible interval for a Bayesian Gower fit.
#'
#' @details This function computes a credible interval for mu, the agreement parameter.
#'
#' This function uses the bootstrap sample to compute a credible interval. If \code{type = "HPD"} is not an element of \dots, the lower and upper limits of the interval are appropriate quantiles of the bootstrap sample. Otherwise the limits are for an HPD interval.
#'
#' @param object an object of class \code{"gower"}, the result of a call to \code{\link{gower.agree}}.
#' @param parm always ignored since there is only one parameter, mu.
#' @param level the desired confidence level for the interval. The default is 0.95.
#' @param \dots additional arguments. These are passed to \code{\link{quantile}}.
#'
#' @return A vector with entries giving the lower and upper limits of the interval. These will be labelled as (1 - level) / 2 and 1 - (1 - level) / 2 if the quantile method was used. They will be labeled "Lower" and "Upper" if the interval is an HPD interval.
#'
#' @seealso \code{\link{gower.agree}}
#'
#' @method confint gower
#'
#' @export
#'
#' @examples
#'
#' # Fit the liver data, using the mean distance for each row of the data matrix.
#' # The range (which is equal to 4) must be passed to \code{\link{gower.agree}}
#' # since these data are ordinal and the L1 distance function is used. We assume
#' # a one-way sampling design for these data, i.e., units are random and coders
#' # are fixed. Also compute a 95\% credible interval using the quantile method
#' # and then the HPD method.
#'
#' data(liver)
#' liver = as.matrix(liver)
#' fit = gower.agree(liver, data.type = "ordinal", range = 4)
#' confint(fit)
#' confint(fit, type = "HPD")

confint.gower = function(object, parm = "mu", level = 0.95, ...)
{
    n = object$units
    args = list(...)
    if ("type" %in% names(args) && args$type == "HPD")
    {
        ci = hpd(object$mu.sample, 1 - level)
        names(ci) = c("Lower", "Upper")
    }
    else
    {
        alpha = pnorm(sqrt(n / (n - 1)) * qt((1 - level) / 2, df = n - 1))
        ci = quantile(object$mu.sample, c(alpha, 1 - alpha))
        names(ci) = c((1 - level) / 2, 1 - (1 - level) / 2)
    }
    ci
}

