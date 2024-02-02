
#' @importFrom utils packageDescription
#' @importFrom stats quantile
#' @importFrom stats influence
#' @importFrom stats confint
#' @importFrom stats density
#' @importFrom stats runif
#' @importFrom stats pnorm
#' @importFrom stats qt
#' @importFrom stats var
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @importFrom graphics hist
#' @importFrom graphics lines

.onAttach = function(libname, pkgname)
{
    temp = packageDescription("goweragreement")
    msg = paste("\n", temp$Package, ": ", temp$Title, "\n", "Version ", temp$Version,
                " created on ", temp$Date, ".\n", sep = "")
    msg = paste(msg, "copyright (c) 2024 John Hughes\n",
                sep = "")
    msg = paste(msg, 'For citation information, type citation("goweragreement").\n', sep = "")
    msg = paste(msg, 'Type help(package = goweragreement) to get started.\n', sep = "")
    packageStartupMessage(msg)
}

