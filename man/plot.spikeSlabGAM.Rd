% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.spikeSlabGAM}
\alias{plot.spikeSlabGAM}
\title{Generates graphical summaries of a fitted model}
\usage{
\method{plot}{spikeSlabGAM}(
  x,
  labels = NULL,
  cumulative = TRUE,
  commonEtaScale = FALSE,
  aggregate = mean,
  quantiles = c(0.1, 0.9),
  gridlength = 20,
  base_size = 12,
  ggElems = list(),
  nrow = min(ceiling(sqrt(length(plotList))), 3),
  ncol = min(ceiling(length(plotList)/nrow), 3),
  interactive.dev = c("x11", "quartz", "windows"),
  ...
)
}
\arguments{
\item{x}{a fitted \code{spikeSlabGAM} model}

\item{labels}{a character vector of names of model terms to be plotted}

\item{cumulative}{Defaults to TRUE, in which case all lower order terms that
are involved in an interaction are cumulated and then plotted (e.g, if a
model contains 2 smooth effects and their interaction, ONLY the sum of the
marginal smooth and linear terms and all their interactions are plotted.) If
FALSE, a separate plot for every term in the model is produced.}

\item{commonEtaScale}{use the same limits for all vertical axes of the
different panels? Defaults to FALSE. Can be useful to compare effect sizes
more easily between panels, but tends to mess up the scales.}

\item{aggregate}{(function) which summary statistic to use for the posterior
of the model terms. Defaults to the mean.}

\item{quantiles}{which quantiles to use for the borders of credible regions.
Defaults to 10\% and 90\% percentiles. Set to NULL to omit plotting credible
regions.}

\item{gridlength}{length of the (univariate) grids for the covariates on which
to evaluate the posterior. Defaults to 20.}

\item{base_size}{default base font size for plot (see e.g.
\code{\link[ggplot2]{theme_gray}})}

\item{ggElems}{a list of plot elements to give to \code{ggplot}. Use this to
supply custom themes or colors, fortify the plot(s) with partial residuals
etc. Defaults to an empty list. Unless specified differently here, the
default ggplot-theme (\code{\link[ggplot2]{theme_gray}}) is changed to a
white background with major gridlines in gray (\code{'grey95'}), no minor
grid lines, and smaller text for the legends.}

\item{nrow}{number of rows per page, defaults to min(sqrt(no. of plots), 3).
See \code{\link[gridExtra]{marrangeGrob}}.}

\item{ncol}{number of columns per page, defaults to min((no. of plots)/nrow,
3). See \code{\link[gridExtra]{marrangeGrob}}.}

\item{interactive.dev}{alternative device to use in interactive mode in
RStudio if output needs to be spread on multiple pages, since the RStudio
graphical device does not support opening multiple displays.}

\item{...}{arguments passed to \code{\link[gridExtra]{marrangeGrob}}.}
}
\value{
a list of \code{\link[ggplot2]{ggplot}}-objects (invisible)
}
\description{
This function plots the estimated linear predictors of the terms in a model on
a grid of values. By default displays all 3-way, 2-way interactions and main
effects present in the model. Starting with ggplot-0.9.2 these are no longer
aligned by their axes due to internal changes in grid and ggplot2. Uses
\pkg{gridExtra}'s  \code{\link[gridExtra]{marrangeGrob}} to arrange the plots
for the terms, also over multiple pages if necessary. This means the graphics
device type is temporarily set to the value of \code{interactive.dev} in
interactive use in RStudio if necessary since the \code{RStudioGD} graphical
device does not support opening multiple pages.
}
\note{
Note that \code{cumulative = TRUE} will only find all relevant terms to
 accumulate if, for all numeric covariates that have a smooth term, the
 smooth term is specified \emph{after} the linear term in the formula.
}
\examples{
#see ?spikeSlabGAM
}
\seealso{
\code{\link{plotTerm}} for more details on the specific plots
}
\author{
Fabian Scheipl
}
