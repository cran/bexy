% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bexy.R
\name{.plotXVsY.bexy}
\alias{.plotXVsY.bexy}
\title{Function to plot X vs Y counts and color them according to BeXY posterior mode}
\usage{
.plotXVsY.bexy(object, ix_X, ix_Y, colors, sampleNamesToHighlight, label, ...)
}
\arguments{
\item{object}{A bexy object.}

\item{ix_X}{An integer representing the index of the scaffold to be plotted on the x-axis.}

\item{ix_Y}{An integer representing the index of the scaffold to be plotted on the y-axis.}

\item{colors}{A vector of length three with the colors for aneuploid, XX and XY samples, respectively.}

\item{sampleNamesToHighlight}{A vector of sample names that should be highlighted. If empty, samples are automatically highlighted based on the posterior probabilities.}

\item{label}{Boolean indicating whether samples should be labeled.}

\item{...}{Other parameters used for plotting}
}
\value{
No return value, called for side effects.
}
\description{
Function to plot X vs Y counts and color them according to BeXY posterior mode
}
\keyword{internal}
