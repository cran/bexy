% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bexy.R
\name{writePosteriorModeAutosomalTrisomies}
\alias{writePosteriorModeAutosomalTrisomies}
\title{Write a file with the posterior mode for each autosomal trisomy}
\usage{
writePosteriorModeAutosomalTrisomies(
  object,
  scaffoldName,
  file,
  threshold_certainty = 0.9
)
}
\arguments{
\item{object}{A bexy object.}

\item{scaffoldName}{The scaffold name, must be an autosome.}

\item{file}{The name of the output file.}

\item{threshold_certainty}{The threshold for certainty on the posterior probabilities. All samples that have a posterior probability less than this threshold are classified as 'uncertain'.}
}
\value{
No return value, called for side effects.
}
\description{
Write a file with the posterior mode for each autosomal trisomy
}
\examples{
bex <- bexy()
tf <- tempfile(fileext = ".txt")
writePosteriorModeAutosomalTrisomies(bex, "21", tf)
}
\seealso{
\code{\link{bexy}}
}
