% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bexy.R
\name{plotBars}
\alias{plotBars}
\title{Plotting a barplot with the posterior probabilities of each sex karyotype per sample}
\usage{
plotBars(
  x,
  colors = c("lightblue", "lightpink", "khaki", "maroon1", "dodgerblue3", "red4", "navy"),
  maxNumSamplesPerPlot = 100,
  sortByKaryotype = TRUE,
  sampleNames = c(),
  ...
)
}
\arguments{
\item{x}{A bexy object.}

\item{colors}{A vector of length 7 with the colors for XY, XX, X0, XXY, XYY, XXX and XXYY, respectively.}

\item{maxNumSamplesPerPlot}{How many samples to place next to each other into one barplot. If there are more samples, multiple barplots will be generated. This is to ensure that the labels are readable.}

\item{sortByKaryotype}{If TRUE, the samples (i.e. the bars) are ordered by the posterior mode of the sex karyotype.}

\item{sampleNames}{A vector of sample names that represents the order in which the bars are plotted. Overrides any argument given to sortByKaryotype.}

\item{...}{Other parameters used for plotting.}
}
\value{
No return value, called for side effects.
}
\description{
Plotting a barplot with the posterior probabilities of each sex karyotype per sample
}
\examples{
bex <- bexy()
plotBars(bex)
}
\seealso{
\code{\link{bexy}}
}
