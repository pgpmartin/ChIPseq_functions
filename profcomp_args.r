#!/usr/bin/env Rscript

# The function expects the following arguments:
	# inputcov: path to a rds file containing a coverage (RleList) or a list of coverages (list of RleLists)
	# windowset: path to a rds file containing the windows (GRanges) on which profiles should be extracted
	# outputdir: Output directory
	# outname: Sample Name (will be the name of the output rds file)
	# Ncores (may not be provided): Number of cores to use

#It returns:
	# an RleList of profiles if inputcov is a coverage
	# a list of RleLists of profiles if inputcov is a list of coverages

#---------------
#Session Setup
#---------------

#Libraries
require(R.utils)
require(GenomicRanges)
require(BiocParallel)

#Functions
profcomp_tx = function(covr, gr)
#Compute profiles from a genome-wide coverage and a set of windows
#covr is a genome-wide coverage (typically obtained with the coverage function)
#gr is a GRanges object containing the windows over which to compute the profiles
#Note that names of covr and seqnames of gr must match
	{
	require(GenomicRanges)
	prof <- covr[gr]
	prof[strand(gr)=='-'] <- lapply(prof[strand(gr)=='-'], rev)
	if (!is.null(names(gr)))
		(names(prof)=names(gr))
	return(prof)
	}


#---------------
#Arguments
#---------------
###-----------------
### Timing
###-----------------
cat("Processing started:", date(),"\n")
ptm <- proc.time()

###-----------------
### Arguments
###-----------------
args=R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE,
                            defaults=list(inputcov=NULL,
                                            windowset=NULL,
                                            outputdir=getwd(),
                                            outname=NULL,
                                            Ncores=1))


#Test arguments
isnullarg <- sapply(args[c('inputcov','windowset','outname')], is.null)
if (any(isnullarg)) {
    stop("missing argument(s): ",
         paste(c('inputcov', 'windowset', 'outname')[isnullarg],
               collapse=', '))}

#Prepare Parallel environment
register(MulticoreParam(args$Ncores))

#Get other arguments
covr = readRDS(args$inputcov) #path to normalized coverage (.rds) or a list of normalized coverages
windowset = readRDS(args$windowset) #path to Windows on which coverage should be extracted (GRanges saved as rds file)
outdir = args$outputdir #output directory
SampleName = args$outname #SampleName


#Output file
outfile = file.path(outdir,
                    paste0(gsub('.rds$', '', SampleName), ".rds"))

#Test if class of arguments are OK
if (!is(windowset, "GRanges"))
	{stop(paste(args$windowset, "must be a GRanges"))}
if (!is(covr, "list") & !is(covr, "RleList"))
	{stop(paste(args$inputcov, "must be an RleList or a list of RleLists"))}
if (!dir.exists(outdir))
	{stop(paste(args$outputdir, "is not a valid directory"))}
if (file.exists(outfile))
	{stop(paste(outfile, "already exists"))}


#---------------
#Compute profiles
#---------------

# Case 1: covr is a coverage
if (is(covr,"RleList"))
	{
	profs = profcomp_tx(covr,windowset)
	}

# Case 2: covr is a list of coverages
if (is(covr,"list"))
	{
	if (!all(sapply(covr, is, "RleList")))
		(stop(paste(args$inputcov, "must be a list of RleLists")))	
	profs = bplapply(covr, profcomp_tx, windowset)
	}

# Save the results as rds
saveRDS(profs, outfile)

###-----------------
### Timing
###-----------------
cat("Processing Ended:", date(),"\n")
cat("Timing:\n")
proc.time() - ptm

#Done


