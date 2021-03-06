% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/joinSpectraData.R
\name{joinSpectraData}
\alias{joinSpectraData}
\title{Join Spectra Data}
\usage{
joinSpectraData(x, y, by.x = "spectrumId", by.y, suffix.y = ".y")
}
\arguments{
\item{x}{A \code{\link[=Spectra]{Spectra()}} object.}

\item{y}{A \code{DataFrame} with spectra data to be merged with
\code{spectraData(x)}.}

\item{by.x}{A \code{character(1)} specifying of the spectra variable
used for merging. Default is \code{"spectrumId"}.}

\item{by.y}{A \code{character(1)} specifying of the column used for
merging. Set to \code{by.x} if missing.}

\item{suffix.y}{A \code{character(1)} specifying the suffix to be used
for making the names of columns in the merged spectra variables
unique. This suffix will be used to amend \code{names(y)}, while
\code{spectraVariables(x)} will remain unchanged.}
}
\description{
Individual spectra data variable can be directly added with the
\verb{$<-} or \verb{[[<-} syntax. The \code{joinSpectraData()} function allows to
merge a \code{DataFrame} to the existing spectra data.

This function diverges from the \code{merge()} method in two main ways:
\itemize{
\item The \code{by.x} and \code{by.y} column names must be of length 1.
\item If variable names are shared in \code{x} and \code{y}, the spectra
variables of \code{x} are not modified. It's only the \code{y} variables
that are appended the suffix defined in \code{suffix.y}. This is to
avoid modifying any core spectra variables that would lead to an
invalid object.
}
}
\examples{

library("Spectra")
library("msdata")
library("magrittr")
library("PSM")

## Creat a Spectra object
ms <- Spectra(msdata::proteomics(pattern = "2014", full.names = TRUE)) \%>\%
    filterMsLevel(2L) \%>\%
    dropNaSpectraVariables()
spectraVariables(ms)

## Additional spectra variables
id <- readPSMs(msdata::ident(full.names = TRUE)) \%>\%
    filterPSMs()
id$modLocation <- NULL
id <- unique(id)
names(id)

## Add the new spectra variables
ms <- joinSpectraData(ms, id,
                      by.x = "spectrumId",
                      by.y = "spectrumID")
spectraVariables(ms)
}
\author{
Laurent Gatto.
}
