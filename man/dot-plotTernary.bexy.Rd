% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bexy.R
\name{.plotTernary.bexy}
\alias{.plotTernary.bexy}
\title{Function to plot a ternary plot of the probabilities for being aneuploid, XX or XY}
\usage{
.plotTernary.bexy(object, colors, sampleNamesToHighlight, label, ...)
}
\arguments{
\item{object}{A bexy object.}

\item{colors}{A vector of length three with the colors for aneuploid, XX and XY samples, respectively.}

\item{sampleNamesToHighlight}{A vector of sample names that should be highlighted. If empty, samples are automatically highlighted based on the posterior probabilities.}

\item{label}{Boolean indicating whether samples should be labeled.}

\item{...}{More parameters passed to TernaryPlot()}
}
\value{
The upper boundary of the plot.
}
\description{
Function to plot a ternary plot of the probabilities for being aneuploid, XX or XY
}
\keyword{internal}
