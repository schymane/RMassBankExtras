\name{collectPeaks}
\alias{collectAnnotationInfo}
\alias{collectCompoundInfo}
\alias{collectPeaks}
\alias{collectSpectraInfo}
\alias{combineTables}
\alias{functions:}
\title{txWorkflow}
\usage{
  collectPeaks(mb, pack, table)

  collectCompoundInfo(mb, pack, table)

  collectSpectraInfo(mb, pack, table)

  collectAnnotationInfo(mb, peakTable)

  combineTables(peakTable, compoundInfo, spectraInfo,
    annotationInfo)
}
\description{
  Functions to export packed MassBank datasets in table
  format.
}
\examples{
\dontrun{
# Preface: load a pack, see the packing functions for info
loadList("E:/RMassBank/Compoundlist.csv")
arcfiles <- c("111028_noNA-pH_wsettings.RData", "111116_pH_RF_wsettings.RData" )
pack <- initializePack()
pack <- addToPack(pack, arcfiles)
# load pack table, edited before
pack.table <- readPackTable("pack_include_pH.csv")
# make the mbPackWorkspace
mb <- new("mbPackWorkspace")
mb <- loadInfolists(mb, "E:/RMassBank/Namelists/")
# then add additional peaks if required, e.g.:
# mb <- addPeaks(mb, "111028_pH_Failpeaks_wOKs.csv")
mb <- compilePack(mb, pack, pack.table)

### here the table export workflow starts:
# Export to table:
# Collect all the peaks in table form
peakTable <- collectPeaks(mb, p2, pack.table)
# collect all basic compound info in table form
compoundInfo <- collectCompoundInfo(mb, p2, pack.table)
# collect spectra-specific info in table form
spectraInfo <- collectSpectraInfo(mb, p2, pack.table)
# collect compound annotation info in table form
annotationInfo <- collectAnnotationInfo(mb, peakTable)
# merge the tables together
finalTable <- combineTables(peakTable, compoundInfo, spectraInfo, annotationInfo)
# export a CSV
write.csv(finalTable, file="tx.csv")
}
}
\author{
  Michael Stravs, Eawag
}

