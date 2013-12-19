# TODO: Add comment
# 
# Author: stravsmi
###############################################################################


# Run the entire script before starting with the packing workflow
# as explained below.

# This workflow allows to collect many spectra from different archives
# and choose which ones to export as MassBank records.

# Limitations of this workflow:
# * All settings files must have the same annotations$internal_id_fieldname
#	  setting, since the compound ID for the molfile filename is extracted
#   from there in the molfile creation step.
# * The records themselves and the files will have the
#   correct prefix as specified in the individual RMassBank settings.

# First, cumbersome step: ###for old data (pre version 0.99.0)
# go into every folder with an archive you want to use, 
# load the _RF archive,
# and re-save the archive including the RMassBank options:
# (RMassBank >= 0.99.0 will save the options into the archive by default.)
.step1 <- function()
{
	setwd("myArchiveFolder")
	w <- loadMsmsWorkspace("myArchive_RF.RData")
	loadRmbSettings("mySettingsFile.R")
	archiveResults(w, "myArchive_wsettings.RData")
}
# Second step: copy all archives with options into one folder.
# If required, also copy the hand-corrected additional-peaks lists into that folder.

# Third step: go into the folder with all archives.
# Load the compound list.
# Load the files, pack them, make a table for one specific 
# ion type (e.g. all [M+H]+ spectra) and export it.
.step3 <- function()
{
	setwd("allArchivesFolder")
	arcfiles <- list.files(,".RData")
	pack <- initializePack()
	pack <- addToPack(pack, arcfiles)
	pack.summary <- summarizePack(pack, "pH")
	pack.table <- makePackTable(pack.summary)
	exportPackTable(pack.table, "pack.csv")
}

# Fourth step: Open the exported table in Excel and go through all the entries.
# Set the OK column to a value, for example "X", (anything not empty) if you want 
# to use the spectra from that line.
# Save the file.

# Fifth step: Reload the table in R.
# Load your infolists.
# Load all additional peak lists if required
# Compile the pack.
# This executes the mbWorkflow until step 4.
# Finally, complete the task by executing steps 5:8.
.step4 <- function()
{
	pack.table <- readPackTable("pack_edited.csv")
	#pack.table <- readPackTable("pack.csv")
	mb <- new("mbPackWorkspace")
	mb <- loadInfolists(mb, "myInfolistPath")
	#mb <- loadInfolists(mb, "H:/Namelists/")
	mb <- addPeaks(mb, "myadditionalpeaks1.csv")
	mb <- addPeaks(mb, "myadditionalpeaks2.csv")
	#...
	mb <- compilePack(mb, pack, pack.table)
    mb <- mbWorkflow(mb, steps=c(5:8))
}


#' @exportClass mbPackWorkspace
#' @export
setClass("mbPackWorkspace", representation(
				"mbWorkspace", name="character",
				settings="list"
				))


#' Initialize a pack of spectra
#' 
#' Creates a new empty container to add archives into
#' 
#' @return An empty list of class \code{RMBPack}. 
#' 
#' @author Michael Stravs
#' @export
initializePack <- function()
{
	pack <- list()
	class(pack) <- "RMBPack"
	return(pack)
}

#' Add archives to a spectra pack
#' 
#' @author Michele Stravs, Eawag
#' @export
addToPack <- function(pack, archiveName)
{
	for(name in archiveName)
	{
		w <- loadMsmsWorkspace(name, TRUE)
		mb <- new("mbPackWorkspace",
					aggregatedRcSpecs = w@aggregatedRcSpecs,
					refilteredRcSpecs = w@refilteredRcSpecs,
					name = name,
					settings = w@settings)
		pack[[name]] <- mb
	}
	return(pack)
}

#' Summarize compound spectra in a pack
#' 
#' @author Michele Stravs, Eawag
#' @export
summarizePack <- function(pack, mode = c())
{
	cpdIdList <- list()
	for(mb in pack)
	{
		for(spec in mb@aggregatedRcSpecs$specFound)
		{
			# If the spectrum is in the mode we are interested in:
			# add an entry into the table
			if((length(mode) == 0) || (spec$mode %in% mode))
			{
				id <- as.character(as.numeric(spec$id))
				if(is.null(cpdIdList[[id]]))
					cpdIdList[[id]] <- list()
				# No. spectra found
				countSpectra <- length(which(unlist(lapply(
												spec$msmsdata, function(m) m$specOK
										))))
				
				# Base peak intensity range
				intensity <- unlist(lapply(spec$msmsdata, function(m)
						{
							if(!m$specOK)
								return(NA)
							else
								return(max(m$childFilt$int))
						}))
				intmax <- max(intensity, na.rm=T)
				intmin <- min(intensity, na.rm=T)
				intmed <- median(intensity, na.rm=T)
				
				# Precursor intensity
				precursor <- spec$msmsdata[[1]]$header$precursorIntensity
				
				cpdIdList[[id]][[mb@name]] <- list(
						cpd = id,
						nspecs = countSpectra,
						precursor = precursor,
						max = intmax,
						med = intmed,
						min = intmin
				)
			}
		}
	}
	return(cpdIdList)
}

#' Create and export pack table
#' 
#' \code{makePackTable} creates a table with a list of spectra, their summary properties, and 
#' which pack they are in. \code{exportPackTable} then exports the table to CSV.
#' The user can then select the appropriate spectra
#' to compile by editing the table and adding OK markers to the spectra they
#' want to select.
#' 
#' @aliases makePackTable exportPackTable
#' @author Michele Stravs, Eawag
#' @seealso readPackTable
#' @export
makePackTable <- function(sumtable)
{
	datablocks <- lapply(sumtable, function(cpd)
			{
				block <- as.data.frame(do.call(rbind, cpd))
				block$archive <- names(cpd)
				block$OK <- ""
				block$name <- ""
				block[which.max(block$precursor),"OK"] <- "X"
				block <- block[,c(
								"cpd", "name", "OK", "archive", "nspecs", "precursor", "max", "med", "min")]
				
			})
	block <- do.call(rbind, datablocks)
	for(colname in colnames(block))
		block[[colname]] <- unlist(block[[colname]])
	block$name <- unlist(lapply(block$cpd, findName))
	return(block)
}

#' @export
exportPackTable <- function(table, fileName)
{
	table$precursor <- sprintf("%.1e", table$precursor)
	table$min <- sprintf("%.1e", table$min)
	table$max <- sprintf("%.1e", table$max)
	table$med <- sprintf("%.1e", table$med)
	write.csv(table, fileName, row.names=F, na="")
}

#' Read pack table
#' 
#' Reads an edited pack table from CSV (as exported by \code{\link{exportPackTable}}).
#' 
#' @export
readPackTable <- function(fileName)
{
	return(read.csv(fileName))
}

#' Compile a mbPackWorkspace
#' 
#' Compiles the contents of a mbPackWorkspace such that all spectra selected
#' in the \code{\link{table}} are converted into MassBank record data. The compiled
#' \code{mbPackWorkspace} then can be used like an usual \code{mbWorkspace} to finish
#' the generation of MassBank records with \code{\link{mbWorkflow}} steps 5:8.
#' 
#' @export
compilePack <- function(mb, pack, table, name="packed")
{
	table_ok <- table[table$OK != "",]
	if(any(duplicated(table_ok$cpd)))
		stop("Duplicate entries are present! Choose only one entry per compound.")
	
	mb@name <- name
	#env$files <- as.character(table_ok$cpd)
	
	compiled_all <- lapply(pack, function(archive)
	{
		# copy additional peaks and mbdata_archive to the pack environment
		archive@additionalPeaks <- mb@additionalPeaks
		archive@mbdata_archive <- mb@mbdata_archive
		# load the right settings for the archive
		loadRmbSettings(archive@settings)
		# cut the archive to pieces so only the selected files stay in there
		cpdList <- as.character(table_ok[table_ok$archive==archive@name,"cpd"])
		archive@aggregatedRcSpecs$specFound <- archive@aggregatedRcSpecs$specFound[
				unlist(lapply(
								archive@aggregatedRcSpecs$specFound, function(s)
									as.character(as.numeric(s$id)) %in% cpdList
								))] 
		# process steps 1:4
		archive <- mbWorkflow(archive, steps=c(1:4))
		return(archive@compiled_ok)

	})
	# Set the folder name (i.e. the prefix at this point) to the "name" parameter
	RmbSettings <- getOption("RMassBank")
	RmbSettings$annotations$entry_prefix <- name
	loadRmbSettings(RmbSettings)
	
	mb@compiled_ok <- do.call(c, compiled_all)
	return(mb)
}

