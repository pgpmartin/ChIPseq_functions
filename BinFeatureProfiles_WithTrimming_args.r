#!/usr/bin/env Rscript

# The function expects the following arguments:
	# prof: path to a rds file containing an RleList with profiles on features (e.g. genes) or a list of RleLists
	# nbins: (default: 100) number of bins to define on each feature
		#OR binwidth (default: NULL): size of bins to define on each feature
	# aggregFUN: (default: 'mean') one of 'min', 'max', 'mean' or 'median' or the full path to an R function used for aggregation
	# returnRemovedFeatures: (default: FALSE) logical indicating if the names of genes removed by the binning procedure should be returned
	# asMatrix: (default: FALSE) logical indicating if the result should be returned as a matrix instead of an RleList
	# outputdir: (default: getwd()) Output directory
	# outFileName: (optional) Sample Name (will be the name of the output rds file)
	# Ncores: (optional) Number of cores to use

#It returns:
	# an Rle with the binned profiles (aggregation done with the aggregFUN function) or a matrix of binned profiles when asMatrix=TRUE

#Libraries
suppressMessages(require(IRanges))
suppressMessages(require(R.utils))
suppressMessages(require(BiocParallel))
suppressMessages(require(BiocGenerics))
suppressMessages(require(parallel))
suppressMessages(require(S4Vectors))
suppressMessages(require(stats4))


###------------------------------------------------------------------------------------------------------------------------------
### BinFeatureProfiles FUNCTION
###------------------------------------------------------------------------------------------------------------------------------

BinFeatureProfiles <- function(FeatureProfiles, 
                                nbins = 100, 
                                binwidth = NULL, 
                                aggregFUN = 'mean', 
                                returnRemovedFeatures = FALSE, 
                                asMatrix = FALSE)
{
#------------------------------------------------------------------
# Function that takes a matrix (Features in row) or an RleList as input, performs a binning on the object and returns an aggregated value for each bin
# The binning is defined either by the number of bins (nbins, DEFAULT) or by the width of the bins (binwidth)
# Features of length < nbins or <binwidth are removed
# The aggregation function aggregFUN can be one of 'min','max','mean' or 'median' or a user-defined function (some arguments may wot work though)
# The function returns an RleList, except if asMatrix=T where it returns, if possible, a matrix
# If returnRemovedFeatures=TRUE, the function returns a list with 1) the RleList/Matrix result and 2) names of removed features
#------------------------------------------------------------------
suppressMessages(require(IRanges))

#Some controls
if (!missingArg(nbins) && !missingArg(binwidth))
    warning("Both nbins and binwidth were provided. Using nbins only")

if (missingArg(nbins) && missingArg(binwidth))
    cat("Working with 100 bins (default). Note that different features may have bins of different size")

if (!missingArg(nbins) && length(nbins)!=0 && !is.na(nbins) && (!is.numeric(nbins) || nbins!=as.integer(nbins)))
    stop("nbins must be an integer")

if (!missingArg(binwidth) && length(binwidth)!=0 && !is.na(binwidth) && (!is.numeric(binwidth) || binwidth!=as.integer(binwidth)))
    stop("binwidth must be an integer")

if (!is.function(aggregFUN) && !(is.character(aggregFUN) && (tolower(aggregFUN) %in% c('min','max','mean','median'))))
    stop("aggregFUN must be a function or one of 'min', 'max', 'mean' or 'median'")

if(!is.logical(returnRemovedFeatures))
    stop("returnRemovedFeatures should be a logical")

if (!is.matrix(FeatureProfiles) && !is(FeatureProfiles,"RleList"))
    stop("FeatureProfiles should be a matrix or an RleList")

if (is.function(aggregFUN) && (length(aggregFUN(1:10))!=1 || !is.numeric(aggregFUN(1:10))))
    stop("aggregFUN should return a single numeric value")

#If FeatureProfile is a matrix, convert to an RleList
if (is.matrix(FeatureProfiles)){
    FeatureProfiles <- as(apply(FeatureProfiles,1,Rle),"RleList")
}

#type of binning: FixedNumBin=TRUE (nbins is given) or FALSE (binwith is given)
FixedNumBin <- TRUE

if (!missingArg(binwidth) && missingArg(nbins)){
	FixedNumBin <- FALSE
}

#Remove Features not compatible with nbins or binsize
FeatureLengths = elementNROWS(FeatureProfiles)

if (FixedNumBin){
	FeatureIsOK <- FeatureLengths >= nbins
	cat("Removing", sum(!FeatureIsOK), 
        "features of size lower than", nbins, "bp\n")
} else {
	FeatureIsOK <- FeatureLengths >= binwidth
	cat("Removing", sum(!FeatureIsOK), 
        "features of size lower than", binwidth, "bp\n")
}

#if returnRemovedFeatures=TRUE, store the names of genes that were removed
if (returnRemovedFeatures){
	RemovedFeatures <- names(FeatureProfiles)[!FeatureIsOK]
}

FeatureProfiles <- FeatureProfiles[FeatureIsOK]
FeatureLengths <- FeatureLengths[FeatureIsOK]


#Create the tiling
if (FixedNumBin){

    tileFeatures = Views(FeatureProfiles,
                         tile(IRanges(start=1,
                                      end=FeatureLengths),
                              n=nbins))
    cat("Bin size is", 
        round(mean(mean(width(tileFeatures))),2),
        "+/-",
        round(sd(mean(width(tileFeatures))),2),
        "bp (mean +/- sd)\n")
} else {
    tileFeatures=Views(FeatureProfiles,
                       tile(IRanges(start=1,
                                    end=FeatureLengths),
                            width=binwidth))
    cat("Average bin size is", 
        round(mean(mean(width(tileFeatures))),2), "bp\n")
}

#If aggregFUN is 'min', 'max', 'mean' or 'median'
if (is.character(aggregFUN) && aggregFUN=='min') {
    res <- as(viewMins(tileFeatures), "CompressedRleList")
}

if (is.character(aggregFUN) && aggregFUN=='max') {
    res <- as(viewMaxs(tileFeatures), "CompressedRleList")
}

if (is.character(aggregFUN) && aggregFUN=='mean') {
    res <- as(viewMeans(tileFeatures), "CompressedRleList")
}

if (is.character(aggregFUN) && aggregFUN=='median') {
    res <- as(viewApply(tileFeatures,median), "CompressedRleList")
}

#If aggregFUN is something else
#Note that some arguments (such as trim in the mean function) may not work

if (!is.character(aggregFUN)) {
    res <- as(viewApply(tileFeatures, aggregFUN), "CompressedRleList")
}

#Convert to a matrix if asMatrix=TRUE and if all binned features have the same length
if (asMatrix) {
    if(any(elementNROWS(res)!=length(res[[1]]))) {
        cat("Binned features have different lengths. Cannot return a matrix. Returning an RleList")
    } else {
        res <- matrix(as.numeric(unlist(res,use.names=F)),
                      nrow=length(res),
                      byrow=TRUE,
                      dimnames=list(names(res), NULL))
    }
}

#Create a list with the result and a vector of removed features
if (returnRemovedFeatures){
    res <- list(binnedFeatures = res, 
                removedFeatures = RemovedFeatures)
}

return(res)
}
###------------------------------------------------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------------------------------------------------


###------------------------------------------------------------------------------------------------------------------------------
### TrimRleList FUNCTION
###------------------------------------------------------------------------------------------------------------------------------

TrimRleList <- function(rlelist, startprof=1, endprof=100) {
require(GenomicRanges)

if (!is(rlelist,"RleList")) {
    stop("must provide an RleList for trimming")
}

if (is.null(names(rlelist))) {
    warning('Names of RleList are NULL')
}

ng <- length(rlelist)

if (!is.integer(startprof)) {startprof  = as.integer(startprof)}
if (!is.integer(endprof)) {endprof = as.integer(endprof)}

if (length(startprof)!=1) {
  if (length(startprof)!=ng) {
    stop('length of startprof is not equal to length of rlelist')} } else {
  startprof <- rep(startprof, ng)
  }

if (length(endprof)!=1) {
  if (length(endprof)!=ng) {
    stop('length of endprof is not equal to length of rlelist')} } else {
  endprof <- rep(endprof, ng)
  }

gr <- GRanges(seqnames = names(rlelist),
              ranges = IRanges(start=startprof, end=endprof))

return(rlelist[gr])
}



###-----------------
### Timing
###-----------------
cat("Processing started:", date(),"\n")
ptm<-proc.time()

###-----------------
### Arguments
###-----------------

args=R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE,
                            defaults=list(nbins = 100,
                                          binwidth = as.integer(NULL),
                                          aggregFUN = 'mean',
                                          returnRemovedFeatures = FALSE, 
                                          asMatrix = FALSE,
                                          trimCol = FALSE,
                                          trimStart = 1,
                                          trimEnd = 100,
                                          outputdir = getwd(),
                                          outFileName = paste0("BinnedProfile_",
                                                              gsub(" ", "_", Sys.time()),
                                                              ".rds"),
                                          Ncores=1))

#Import aggregFUN if it is not 'min','max','mean' or 'median'
if (!(tolower(args$aggregFUN) %in% c('min','max','mean','median'))){
    stopifnot(file.exists(args$aggregFUN))
    aggregFUN <- source(args$aggregFUN)$value
## Be carefull not to import a function that has a name used in this script !
    if (!exists(aggregFUN) || !is.function(aggregFUN))
        stop("Cannot import",args$aggregFUN)
        stopifnot(length(aggregFUN(1:5))==1,
                  is.numeric(aggregFUN(1:5)))
} else {
    aggregFUN = args$aggregFUN
}

#Check nbins and binwidth argument
if (length(args$binwidth)==0 || is.na(args$binwidth)) args$binwidth=NULL
if (length(args$nbins)==0 || is.na(args$nbins)) args$nbins=NULL

#Verbose
cat("List of arguments for the call:\n")
args

cat("aggregFUN used to aggregate the data within bins:\n")
aggregFUN

###-----------------
### Set number of cores
###-----------------
register(MulticoreParam(args$Ncores))

###-----------------
### read-in input file
###-----------------

inputProfiles = readRDS(args$prof)

if (!is(inputProfiles,"RleList") && !(is(inputProfiles,"list") && all(sapply(inputProfiles,is,"RleList")))){
    stop("prof must be an RleList or a list of RleLists")
}

###-----------------
### Trim the input file if necessary
###-----------------

if (args$trimCol) {
require(GenomicRanges)

  if (is(inputProfiles,"RleList")) {
  inputProfiles <- TrimRleList(inputProfiles,
                               startprof  = args$trimStart,
                               endprof = args$trimEnd)
  }

  if (is(inputProfiles,"list") && all(sapply(inputProfiles,is,"RleList"))) {
  inputProfiles <- lapply(inputProfiles,
                          TrimRleList,
                          startprof  = args$trimStart,
                          endprof =  args$trimEnd)
   }

cat("Before binning, the profiles are trimmed based on the provided trimStart and trimEnd values\n")

}

###-----------------
### Run the binning
###-----------------

if (is(inputProfiles,"RleList")){
    res <- BinFeatureProfiles(inputProfiles,
                                nbins = args$nbins, 
                                binwidth = args$binwidth, 
                                aggregFUN = aggregFUN, 
                                returnRemovedFeatures = args$returnRemovedFeatures, 
                                asMatrix = args$asMatrix)
}

if (is(inputProfiles,"list") && all(sapply(inputProfiles,is,"RleList"))){
    if (args$Ncores>=8 & args$Ncores>=length(inputProfiles)){
        res <- bplapply(inputProfiles,
                        BinFeatureProfiles, 
                        nbins = args$nbins, 
                        binwidth = args$binwidth, 
                        aggregFUN = aggregFUN, 
                        returnRemovedFeatures = args$returnRemovedFeatures, 
                        asMatrix = args$asMatrix)
        } else {
            res <- lapply(inputProfiles,
                          BinFeatureProfiles, 
                          nbins = args$nbins, 
                          binwidth = args$binwidth, 
                          aggregFUN = aggregFUN, 
                          returnRemovedFeatures = args$returnRemovedFeatures, 
                          asMatrix = args$asMatrix)
        }
    names(res) = names(inputProfiles)
}

###-----------------
### Save Extracted profiles
###-----------------

saveRDS(res,
        file.path(args$outputdir,
                  paste0(gsub(".rds$", "", args$outFileName),
                        ".rds")))

###-----------------
### Timing
###-----------------
cat("Processing Ended:", date(),"\n")
cat("Timing:\n")
proc.time()-ptm

#Done
