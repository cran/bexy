% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bexy.R
\name{plot.bexy}
\alias{plot.bexy}
\title{Plotting a bexy object. To customize the plots, please call the respective plotting functions directly (e.g. plotTernary, plotCounts, plotBars...).}
\usage{
\method{plot}{bexy}(x, ...)
}
\arguments{
\item{x}{A bexy object.}

\item{...}{Other parameters used for plotting.}
}
\value{
No return value, called for side effects.
}
\description{
Plotting a bexy object. To customize the plots, please call the respective plotting functions directly (e.g. plotTernary, plotCounts, plotBars...).
}
\examples{
bex <- bexy()
plot(bex)
}
\seealso{
\code{\link{bexy}}
}
