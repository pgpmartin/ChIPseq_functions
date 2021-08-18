#!/usr/bin/env Rscript

# The function expects the following arguments:
	# covFile: path to an rds file containing a coverage (RleList)
	# sampleName: a character string for the name of the sample 
	# outdir: an output directory for the results
#It returns:
	# a bigwig file with the coverage 
	

#Libraries
require(GenomicAlignments)
require(GenomeInfoDb)
require(rtracklayer)
require(R.utils)
require(BiocParallel)


###-----------------
### Timing
###-----------------
cat("Processing started:", date(),"\n")
ptm<-proc.time()

###-----------------
### Arguments
###-----------------

args=R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE,
                          defaults=list(covFile  = NULL,
                                        sampleName = "MySample",
                                        outdir = getwd(),
                                        Ncores = 1))

#Verification
if (!file.exists(args$covFile)) {
    stop(args$covFile, " is absent")
}

sampleName <- args$sampleName
covFile <- args$covFile
outdir <- args$outdir

#Verbose
cat("Analysis of sample:", sampleName, "\n")
cat("List of arguments for the call:\n")
args

###-----------------
### Set number of cores
###-----------------
register(MulticoreParam(args$Ncores))

###-----------------
### Import coverage
###-----------------
covr <- readRDS(covFile)

###-----------------
### Export as bigwig
###-----------------
#file name
outn <- file.path(outdir, paste0(sampleName, '.bigWig'))
if (file.exists(outn)) {
warning(outn, " already exists, it will be overwritten")
}

#export
rtracklayer::export(covr, outn)


#Verbose
cat("\n")
cat("Sample name:", sampleName,"\n")
cat("Exported bigWig file:", outn,"\n")
cat("\n")


###-----------------
### Timing
###-----------------
cat("Processing Ended:", date(),"\n")
cat("Timing:\n")
proc.time()-ptm

