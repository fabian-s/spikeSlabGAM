.onAttach <- function(libname, pkgname) {
  vers <- packageDescription("spikeSlabGAM")[["Version"]]
  packageStartupMessage("## ----- This is spikeSlabGAM ", vers, " ----- ##\n",
    "Please note that a recent update to gridExtra has made it necessary to change the interface for plot.spikeSlabGAM in version 1.1-9.\n",
    " Instead of arguments 'rows', 'cols', 'widths', 'heights', 'maxPlotsPerPage', it now accepts only 'nrow' and 'ncol'.\n",
    " Arguments 'widths' & 'heights' can still be defined and are handed over to gridExtra:::marrangeGrob.\n",
    " Sorry for the inconvenience.")
}
