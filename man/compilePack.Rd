\name{compilePack}
\alias{compilePack}
\title{Compile a mbPackWorkspace}
\usage{
  compilePack(mb, pack, table, name = "packed")
}
\description{
  Compiles the contents of a mbPackWorkspace such that all
  spectra selected in the \code{\link{table}} are converted
  into MassBank record data. The compiled
  \code{mbPackWorkspace} then can be used like an usual
  \code{mbWorkspace} to finish the generation of MassBank
  records with \code{\link{mbWorkflow}} steps 5:8.
}

