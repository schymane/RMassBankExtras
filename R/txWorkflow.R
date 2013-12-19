# txWorkflow: export spectra to table

#' txWorkflow
#' 
#' Functions to export packed MassBank datasets in table format.
#' 
#' @aliases collectPeaks collectCompoundInfo collectSpectraInfo collectAnnotationInfo combineTables
#' 
#' functions:
#' @usage
#'  collectPeaks(mb, pack, table)
#'  
#'  collectCompoundInfo(mb, pack, table)
#'  
#'  collectSpectraInfo(mb, pack, table)
#'  
#'  collectAnnotationInfo(mb, peakTable)
#'  
#'  combineTables(peakTable, compoundInfo, spectraInfo, annotationInfo)
#'  
#' @examples
#' 
#' \dontrun{
#' # Preface: load a pack, see the packing functions for info
#' loadList("E:/RMassBank/Compoundlist.csv")
#' arcfiles <- c("111028_noNA-pH_wsettings.RData", "111116_pH_RF_wsettings.RData" )
#' pack <- initializePack()
#' pack <- addToPack(pack, arcfiles)
#' # load pack table, edited before
#' pack.table <- readPackTable("pack_include_pH.csv")
#' # make the mbPackWorkspace
#' mb <- new("mbPackWorkspace")
#' mb <- loadInfolists(mb, "E:/RMassBank/Namelists/")
#' # then add additional peaks if required, e.g.:
#' # mb <- addPeaks(mb, "111028_pH_Failpeaks_wOKs.csv")
#' mb <- compilePack(mb, pack, pack.table)
#' 
#' ### here the table export workflow starts:
#' # Export to table:
#' # Collect all the peaks in table form
#' peakTable <- collectPeaks(mb, p2, pack.table)
#' # collect all basic compound info in table form
#' compoundInfo <- collectCompoundInfo(mb, p2, pack.table)
#' # collect spectra-specific info in table form
#' spectraInfo <- collectSpectraInfo(mb, p2, pack.table)
#' # collect compound annotation info in table form
#' annotationInfo <- collectAnnotationInfo(mb, peakTable)
#' # merge the tables together
#' finalTable <- combineTables(peakTable, compoundInfo, spectraInfo, annotationInfo)
#' # export a CSV
#' write.csv(finalTable, file="tx.csv")
#' }
#' @author Michael Stravs, Eawag
#' @export
collectPeaks <- function(mb, pack, table)
{
  table_ok <- table[table$OK != "",]
  if(any(duplicated(table_ok$cpd)))
    stop("Duplicate entries are present! Choose only one entry per compound.")
  
  # cut out all compounds which do not really occur
  pack.croppedpeaks <- lapply(pack, function(archive)
  {


    # copy additional peaks and mbdata_archive to the pack environment
    additionalPeaks <- mb@additionalPeaks
    refiltered <- archive@refilteredRcSpecs
    # cut the archive to pieces so only the selected files stay in there
    cpdList <- as.character(table_ok[table_ok$archive==archive@name,"cpd"])

    peaks <- refiltered$peaksOK
    reanPeaks <- refiltered$peaksReanOK
    
    # if cpdIDs are numeric, force everything to numeric 
    if(is.numeric(additionalPeaks$cpdID))
    {
      peaks$cpdID <- as.numeric(as.character(peaks$cpdID))
      reanPeaks$cpdID <- as.numeric(as.character(reanPeaks$cpdID))
    }
    
    # get all the peaks for the given cpdIDs
    
    peaks <- peaks[(peaks$cpdID %in% cpdList),,drop=FALSE]
    peaks$intrel <- numeric(nrow(peaks))
    reanPeaks <- reanPeaks[(reanPeaks$cpdID %in% cpdList),,drop=FALSE]
    reanPeaks$intrel <- numeric(nrow(reanPeaks))
    
    addPeaks <- additionalPeaks[ 
      (additionalPeaks$OK == 1) & 
        (additionalPeaks$cpdID %in% cpdList),
      ,drop=FALSE]
    addPeaks$intrel <- numeric(nrow(addPeaks))
    #addPeaks$parentScan <- integer(nrow(addPeaks))
    addPeaks$formula <- character(nrow(addPeaks))
    addPeaks$dppm <- numeric(nrow(addPeaks))
    addPeaks$mzCalc <- numeric(nrow(addPeaks))
    addPeaks$formulaCount <- integer(nrow(addPeaks))
    # convert cpdIDs to character

#     
#     # recompute the parentScan for addPeaks
#     addPeaks <- merge(addPeaks,
#                       unique(peaks[,c("cpdID", "scan", "parentScan")]),
#                       by=c("cpdID", "scan"),
#                       all.x = TRUE,
#                       all.y = FALSE,
#                       suffixes=c("", ".peaks")
#                       )
    
    # peak tables:
    # Fragment Peak m/z Fragment Peak intensity Fragment 
    # Peak rel intensity Fragment Peak annotated formula Fragment
    # Peak annotated formula - mass Fragment Peak annotated formula - mass error 
    # Fragment Peak annotated formula - number
    peaks.select <- peaks[,
                          c("cpdID",
                            "scan",
                            "mzFound",
                            "int",
                            "intrel",
                            "formula",
                            "mzCalc",
                            "dppm",
                            "formulaCount"), drop=FALSE]
    reanPeaks.select <- reanPeaks[,
                          c("cpdID",
                            "scan",
                            "mzFound",
                            "int",
                            "intrel",
                            "reanalyzed.formula",
                            "reanalyzed.mzCalc",
                            "reanalyzed.dppm",
                            "reanalyzed.formulaCount"), drop=FALSE]
    addPeaks.select <- addPeaks[,
                                c("cpdID",
                                  "scan",
                                  "mzFound",
                                  "int",
                                  "intrel",
                                  "formula",
                                  "mzCalc",
                                  "dppm",
                                  "formulaCount"), drop=FALSE]
    
    
    colnames(reanPeaks.select) <- colnames(peaks.select)
    colnames(addPeaks.select) <- colnames(peaks.select)
    return(rbind(peaks.select, reanPeaks.select, addPeaks.select))
  })
  peaks <- do.call(rbind,pack.croppedpeaks)
  # recompute relative intensities
  # find spectrum maximum for each peak, and merge into table
  maxInt <- as.data.frame(aggregate(peaks$int, 
                                    by=list(as.numeric(as.character(peaks$cpdID)),
                                            as.numeric(as.character(peaks$scan))), max)
  )
  colnames(maxInt) <- c("cpdID", "scan", "maxint")
  maxInt$cpdID <- as.numeric(as.character(maxInt$cpdID))
  maxInt$scan <- as.numeric(as.character(maxInt$scan))
  peaks <- merge(peaks, maxInt, by=c("cpdID", "scan"))
  peaks$intrel <- peaks$int / peaks$maxint
  peaks$maxint <- NULL
  
  
  
  colnames(peaks) <- c("cpdID",
                       "scan",
                       "FragmentMZ",
                       "FragmentIntensity",
                       "FragmentRelativeIntensity",
                       "FragmentAnnotatedFormula",
                       "FragmentFormulaMZ",
                       "FragmentFormulaErrorPpm",
                       "FragmentFormulaCount")
  
  
  
  return(peaks)
}


gatherCompoundInfo <- function(spec)
{
  id <- spec$id
  imode <- spec$mode
  
  # define positive or negative, based on processing mode.
  ion_modes <- list(
    "pH" = "POSITIVE",
    "pNa" = "POSITIVE",
    "mH" = "NEGATIVE",
    "mFA" = "NEGATIVE",
    "pM" = "POSITIVE",
    "mM" = "NEGATIVE")
  mode <- ion_modes[[imode]]
  
  # for format 2.01
  ac_ms <- list();
  ac_ms[['MS_TYPE']] <- getOption("RMassBank")$annotations$ms_type
  ac_ms[['IONIZATION']] <- getOption("RMassBank")$annotations$ionization
  ac_ms[['ION_MODE']] <- mode
  
  # This list could be made customizable.
  ac_lc <- list();
  rt  <- spec$parentHeader[1,"retentionTime"] / 60
  ac_lc[['COLUMN_NAME']] <- getOption("RMassBank")$annotations$lc_column
  ac_lc[['FLOW_GRADIENT']] <- getOption("RMassBank")$annotations$lc_gradient
  ac_lc[['FLOW_RATE']] <- getOption("RMassBank")$annotations$lc_flow
  ac_lc[['rttable']] <- sprintf("%.1f", rt)  
  ac_lc[['SOLVENT A']] <- getOption("RMassBank")$annotations$lc_solvent_a
  ac_lc[['SOLVENT B']] <- getOption("RMassBank")$annotations$lc_solvent_b
  
  mbdata <- list();
  mbdata[['AC$INSTRUMENT']] <- getOption("RMassBank")$annotations$instrument
  mbdata[['AC$INSTRUMENT_TYPE']] <- getOption("RMassBank")$annotations$instrument_type
  
  
  return(list(
    id = id,
    mode = mode,
    rt = rt,
    ac_ms = ac_ms,
    ac_lc = ac_lc,
    mbdata = mbdata
  ))
}

#' @export
collectCompoundInfo <- function(mb, pack, table)
{
  table_ok <- table[table$OK != "",]
  if(any(duplicated(table_ok$cpd)))
    stop("Duplicate entries are present! Choose only one entry per compound.")
  
  # cut out all compounds which do not really occur
  pack.compoundinfo <- lapply(pack, function(archive)
  {
    archive@mbdata_archive <- mb@mbdata_archive
    # load the right settings for the archive
    loadRmbSettings(archive@settings)
    # cut the archive to pieces so only the selected files stay in there
    cpdList <- as.character(table_ok[table_ok$archive==archive@name,"cpd"])
    specFound <- archive@aggregatedRcSpecs$specFound[
      unlist(lapply(
        archive@aggregatedRcSpecs$specFound, function(s)
          # this makes "0192" == "192" == 192
          as.character(as.numeric(as.character(s$id))) %in% cpdList
      ))]
    info.all <- lapply(specFound, gatherCompoundInfo)
    # flatten the info to a table structure
    # Instrument Intrument type MS-Type (MS2) Ion_Mode (pos/neg) #
    # Ionization Type (ESI) Fragmentation_Mode Collision Energy Resolution Retention Time
    info.rows <- lapply(info.all, function(info) c(
      cpdID = info$id,
      Instrument = info$mbdata[["AC$INSTRUMENT"]],
      InstrumentType = info$mbdata[["AC$INSTRUMENT_TYPE"]],
      MSType = info$ac_ms[["MS_TYPE"]],
      IonMode = info$ac_ms[["ION_MODE"]],
      IonizationType = info$ac_ms[["IONIZATION"]],
      RT = info$ac_lc[["rttable"]]
      ))
    info.table <- do.call(rbind, info.rows)
    
    return(info.table)
  })
  
  pack.ci.table <- do.call(rbind, pack.compoundinfo)
  return(pack.ci.table)
}

# get info for individual spectra, ex settings and spectra scan etc
#' @export
collectSpectraInfo <- function(mb, pack, table)
{
  table_ok <- table[table$OK != "",]
  if(any(duplicated(table_ok$cpd)))
    stop("Duplicate entries are present! Choose only one entry per compound.")
  
  # cut out all compounds which do not really occur
  pack.spectrainfo <- lapply(pack, function(archive)
  {
    archive@mbdata_archive <- mb@mbdata_archive
    # load the right settings for the archive
    loadRmbSettings(archive@settings)
    # cut the archive to pieces so only the selected files stay in there
    cpdList <- as.character(table_ok[table_ok$archive==archive@name,"cpd"])
    specFound <- archive@aggregatedRcSpecs$specFound[
      unlist(lapply(
        archive@aggregatedRcSpecs$specFound, function(s)
          as.character(as.numeric(as.character(s$id))) %in% cpdList
      ))]
    
    spectraInfo <- lapply(specFound, function(spec)
      {
      info <- lapply(spec$msmsdata, function(m) {
        i <- gatherSpectrumInfo(spec, m)
        if(is.null(i)) return(NULL)
        i.row <- c(
          cpdID = i$id,
          scan = i$scan,
          FragmentationMode = i$ac_ms[["FRAGMENTATION_MODE"]],
          CollisionEnergy = i$ac_ms[["cetable"]],
          Resolution = i$ac_ms[["RESOLUTION"]],
          PrecursorType = i$ms_fi[["PRECURSOR_TYPE"]],
          PrecursorMZ = i$ms_fi[["PRECURSOR_M/Z"]]
          )
        return(i.row)
      })
      return(do.call(rbind, info))                    
    })
    return(do.call(rbind, spectraInfo))
  })
  return(do.call(rbind, pack.spectrainfo))
}

# Fragmentation_Mode Collision Energy Resolution 
# Precursor Type Precursor m/z 
gatherSpectrumInfo <- function(spec, msmsdata)
{
  # If the spectrum is not filled, return right now. All "NA" spectra will
  # not be treated further.
  if(msmsdata$specOK == FALSE)
    return(NULL)
  # get data
  scan <- msmsdata$scan
  id <- spec$id
  # Further fill the ac_ms datasets, and add the ms$focused_ion with spectrum-specific data:
  precursor_types <- list(
    "pH" = "[M+H]+",
    "pNa" = "[M+Na]+",
    "mH" = "[M-H]-",
    "mFA" = "[M+HCOO-]-",
    "pM" = "[M]+",
    "mM" = "[M]-")
  ac_ms <- list()
  ac_ms[['FRAGMENTATION_MODE']] <- msmsdata$info$mode
  #ac_ms['PRECURSOR_TYPE'] <- precursor_types[spec$mode]
  ac_ms[['cetable']] <- as.numeric(sub("%","",msmsdata$info$ces))
  ac_ms[['RESOLUTION']] <- msmsdata$info$res
  
  # Calculate exact precursor mass with Rcdk, and find the base peak from the parent
  # spectrum. (Yes, that's what belongs here, I think.)
  precursorMz <- findMz(spec$id, spec$mode)
  ms_fi <- list()
  ms_fi[['BASE_PEAK']] <- round(spec$parentMs[which.max(spec$parentMs[,"int"]),"mz"],4)
  ms_fi[['PRECURSOR_M/Z']] <- round(precursorMz$mzCenter,4)
  ms_fi[['PRECURSOR_TYPE']] <- precursor_types[[spec$mode]]
  
  return(list(
    id = id,
    scan = scan,
    ac_ms = ac_ms,
    ms_fi = ms_fi
    ))
}

#' @export
collectAnnotationInfo <- function(mb, peakTable)
{
  mbdata_archive <- mb@mbdata_archive
  
  # Substanz_ID (kommt aus der Tabelle Alles unserer Access Datenbank und steht hoffentlich im Massbank Record, da ich diese Id zum Verknüpfen nehmen will) 
  # Name Formula CAS INCHI Key
  
  # if the id in the peaktable is numeric, convert the archive id to numeric
  if(is.numeric(peakTable$cpdID))
  {
    mbdata_archive$id <- as.numeric(as.character(mbdata_archive$id))
  }
  cpdList <- unique(peakTable$cpdID)
  
  archive <- mbdata_archive[
      mbdata_archive$id %in% cpdList,
      c("id","CH.NAME1", "CH.FORMULA", "CH.LINK.CAS", "CH.LINK.INCHIKEY")
      ]
  colnames(archive) <- c("cpdID", "Name", "Formula", "CAS", "InChIKey")
  return(archive)
}

# Complete the peak table with compound info
#' @export
combineTables <- function(peakTable, compoundInfo, spectraInfo, annotationInfo)
{
  compoundInfo <- as.data.frame(compoundInfo)
  compoundInfo$cpdID <- as.numeric(as.character(compoundInfo$cpdID))
  spectraInfo <- as.data.frame(spectraInfo)
  spectraInfo$cpdID <- as.numeric(as.character(spectraInfo$cpdID))
  spectraInfo$scan <- as.numeric(as.character(spectraInfo$scan))
  peakTable$cpdID <- as.numeric(as.character(peakTable$cpdID))
  annotationInfo$cpdID <- as.numeric(as.character(annotationInfo$cpdID))
  t1 <- merge(peakTable, compoundInfo, by="cpdID")
  t2 <- merge(t1, annotationInfo, by="cpdID")
  t3 <- merge(t2, spectraInfo, by=c("cpdID", "scan"))
  # cut the scan column out and rename cpdID
  colnames(t3)[[1]] <- "CompoundID"
  t3 <- t3[,c("CompoundID", "FragmentMZ", "FragmentIntensity",  "FragmentRelativeIntensity", 
              "FragmentAnnotatedFormula", "FragmentFormulaMZ",  "FragmentFormulaErrorPpm", 
              "FragmentFormulaCount", "Instrument",  "InstrumentType", "MSType", "IonMode", 
              "IonizationType", "RT",  "Name", "Formula", "CAS", "InChIKey", "FragmentationMode",
              "CollisionEnergy",  "Resolution", "PrecursorType", "PrecursorMZ")]
  
  return(t3)
}