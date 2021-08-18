#!/usr/bin/env Rscript

# The function (FOR PAIRED-END DATA ONLY) expects the following arguments:
    # bampath: path to a BAM file containing paired-end reads on which to compute the coverage (note that the path to an rds file containing a GAlignment object also works)
    # sampleName: a character string for the name of the sample under study
    # outdir: an output directory for the results
    # resizeLength: the length at which the reads must be resized to compute the coverage on resized fragments
    # minFragSize: an integer indicating the minimum size of a fragment. Fragments of length < minFragSize are dropped
    # maxFragSize: an integer indicating the maximum size of a fragment. Fragments of length > maxFragSize are dropped
    # chromSizes: path to an rds file containing a 2 column table with chromosome names and lengths (Only these chromosomes will be used for coverage computation and RPM normalization)
    # MultiplicationFactor: The function calculates "Reads Per Million" by default but MultiplicationFactor allows you to adjust this. For example, if you want Reads per 10 millions then use MultiplicationFactor=10
    # Ncores: number of core to use
    # RemoveDuplicateReads: logical indicating if duplicate reads should be removed (TRUE) or not (FALSE, default)
#The function generates the following files:
    # an rds file (bamGAp_*.rds) with the result of the readGAlignmentPairs function (if a bam file and not an rds file is provided for bampath)
    # an rds file (FragGR_*.rds) containing a GRanges oject with the coordinates of the fragments (possibly filtered based on minFragSize/maxFragSize and resized to resizeLength).
        # This file is filtered to keep only fragments that are in the chromSizes file
    # an rds file (covrRaw_*.rds) containing an RleList with the raw coverage (without normalization for library size) from the fragments (allows to apply a different normalization if necessary)
    # an rds file (normcovr_*.rds) containing an RleList with the normalized coverage (reads per million -RPM- * MultiplicationFactor)

#Libraries
suppressPackageStartupMessages(require(GenomicAlignments))
suppressPackageStartupMessages(require(GenomeInfoDb))
suppressPackageStartupMessages(require(R.utils))
suppressPackageStartupMessages(require(BiocParallel))


###-----------------
### Timing
###-----------------
cat("Processing started:", date(),"\n")
ptm <- proc.time()

###-----------------
### Arguments
###-----------------

args=R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE,
                          defaults = list(bampath = NULL,
                                          sampleName = "MySample",
                                          outdir = getwd(),
                                          resizeLength = NULL,
                                          minFragSize = NULL,
                                          maxFragSize = NULL,
                                          chromSizes = "/<path_to>/STDchromSizes.rds",
                                          MultiplicationFactor = 1,
                                          Ncores=1,
                                          RemoveDuplicateReads=FALSE))

#Verification
if (!file.exists(args$chromSizes)) {
    stop(args$chromSizes, " is absent")
}

if (!file.exists(args$bampath)) {
    stop(args$bampath, " is absent")
}

sampleName <- args$sampleName
bampath <- args$bampath
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
### Parameters for importing BAM files
###-----------------
param=if (args$RemoveDuplicateReads) {
        ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE,
                                        hasUnmappedMate = FALSE,
                                        isNotPassingQualityControls = FALSE,
                                        isProperPair = TRUE,
                                        isDuplicate = FALSE))} else {
        ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE,
                                        hasUnmappedMate = FALSE,
                                        isNotPassingQualityControls = FALSE,
                                        isProperPair = TRUE))
        }

###-----------------
#Import BAM file
# if it is a .rds file it is imported directly
# If it is a bam file it is imported using param (see above) and saved as rds
###-----------------
if (length(grep(".rds$", basename(bampath)))==1) {
    bam <- readRDS(bampath)
    cat("BAM file imported from:", bampath, "\n")
} else {
    bam=readGAlignmentPairs(bampath, param=param)
    saveRDS(bam, file=file.path(outdir, paste0("bamGAp_", sampleName, ".rds")))
    cat("BAM saved as rds file\n")
}


###-----------------
#Convert to fragments (resize fragments if necessary) and save as rds
###-----------------
#Fragments as a GRanges:
	grbam <- granges(bam)
#output file name:
	outn <- file.path(outdir, paste0("FragGR_", sampleName, ".rds"))

#Filter fragments based on min/max size if necessary:
if (any(c(!is.null(args$minFragSize), 
          !is.null(args$maxFragSize)))) {
  if (!is.null(args$minFragSize)) {
    if (!is.numeric(args$minFragSize))
      {warning("minFragSize is not numeric. No filtering done based on minFragSize")} else {
  minFragSize <- as.integer(args$minFragSize)
  okMinFragSize <- width(grbam)>=minFragSize
  grbam <- grbam[okMinFragSize]
  outn <- paste0(gsub(".rds$","",outn), "_min", minFragSize, ".rds")
  cat("minFragSize Filter:\n")
  cat("Fragments analyzed:", 
        length(okMinFragSize),"\n")
  cat("Fragments removed (size<minFragSize):", 
        sum(!okMinFragSize), 
        paste0("(", round(mean(!okMinFragSize),1), "%)"),
        "\n")
  cat("Fragments kept (size>=minFragSize):", 
        sum(okMinFragSize), 
        paste0("(", round(mean(okMinFragSize),1), "%)"),
        "\n")
      }
  }
  if (!is.null(args$maxFragSize)) {
    if (!is.numeric(args$maxFragSize))
      {warning("maxFragSize is not numeric. No filtering done based on maxFragSize")} else {
  maxFragSize <- as.integer(args$maxFragSize)
  okMaxFragSize <- width(grbam)<=maxFragSize
  grbam <- grbam[okMaxFragSize]
  outn <- paste0(gsub(".rds$","",outn), "_max", maxFragSize, ".rds")
  cat("maxFragSize Filter:\n")
  cat("Fragments analyzed:", 
        length(okMaxFragSize), "\n")
  cat("Fragments removed (size>maxFragSize):", 
        sum(!okMaxFragSize), 
        paste0("(", round(mean(!okMaxFragSize),1), "%)"),
        "\n")
  cat("Fragments kept (size<=maxFragSize):", 
        sum(okMaxFragSize), 
        paste0("(", round(mean(okMaxFragSize),1), "%)"),
        "\n")
      }
  }
}

#resize if necessary
if (!is.null(args$resizeLength)) {
    if (!is.numeric(args$resizeLength))
      {warning("resizeLength is not numeric. No fragment resizing done")} else {
      resizeLength <- as.integer(args$resizeLength)
      grbam <- resize(grbam, resizeLength, fix="center")
      outn <- paste0(gsub(".rds$","",outn), "_resized", resizeLength, ".rds")
  cat("Fragments resized to:", paste0(resizeLength, "bp"), "\n")
      }
}

###-----------------
# Filter grbam based on chromSizes and set the seqlengths of grbam
###-----------------
chromSizes <- readRDS(args$chromSizes)
NonSTDchrom <- setdiff(seqlevels(grbam), chromSizes[,1])

grbam <- keepSeqlevels(grbam, chromSizes[,1], pruning.mode="coarse")
seqlengths(grbam) <- chromSizes[match(names(seqlengths(grbam)), chromSizes[,1]),2]


###-----------------
# save Fragment file (only if the bam file used is not and rds file)
###-----------------

if (length(grep(".rds$",basename(bampath)))==0) {
    saveRDS(grbam,file=outn) 
    if (!is.null(args$resizeLength)) {
      cat("Fragment (RESIZED!!) file saved\n")
    } else {cat("Fragment file saved\n")}
}

###-----------------
#Get normalization factor (reads per MILLION)
###-----------------

Numfrag=length(grbam)
Knorm=1e6/Numfrag

cat("\n\n")
cat("Sample name:", sampleName,"\n")
cat("Using the following standard chromosomes only:", chromSizes[,1], "\n")
cat("The following chromosomes are not considered for normalization:", NonSTDchrom, "\n")
cat("Number of fragments (standard chromosomes only):", Numfrag, "\n")
cat("Normalization factor (reads per million):", Knorm, "\n")
cat("\n\n")


###-----------------
#Get raw coverage from fragments
###-----------------
#Get coverage
	rawcov <- coverage(grbam)
#Output file name
	outn <- file.path(outdir, paste0("covrRaw_", sampleName, ".rds"))

if (!is.null(args$minFragSize) & is.numeric(args$minFragSize)) {
  outn <- paste0(gsub(".rds$","",outn), "_min", as.integer(args$minFragSize), ".rds")
}
if (!is.null(args$maxFragSize) & is.numeric(args$maxFragSize)) {
  outn <- paste0(gsub(".rds$","",outn), "_max", as.integer(args$maxFragSize),".rds")
}
if (!is.null(args$resizeLength) & is.numeric(args$resizeLength)) {
  outn <- paste0(gsub(".rds$", "", outn), "_resized", as.integer(args$resizeLength), ".rds")
}
#Save raw coverage
	saveRDS(rawcov, file=outn) #coverage with raw fragments
cat("Raw coverage saved\n")


###-----------------
#Get normalized coverage (default is reads per MILLION) from fragments
###-----------------
#Get normalized coverage
    MultiplicationFactor <- as.numeric(args$MultiplicationFactor)
    normcov <- rawcov * Knorm * MultiplicationFactor
#output file name
	outn <- file.path(outdir, paste0("normcovr_", sampleName, ".rds"))

if (!is.null(args$minFragSize) & is.numeric(args$minFragSize)) {
  outn <- paste0(gsub(".rds$","",outn), "_min", as.integer(args$minFragSize), ".rds")
}
if (!is.null(args$maxFragSize) & is.numeric(args$maxFragSize)) {
  outn <- paste0(gsub(".rds$","",outn), "_max", as.integer(args$maxFragSize), ".rds")
}
if (!is.null(args$resizeLength) & is.numeric(args$resizeLength)) {
  outn <- paste0(gsub(".rds$","",outn), "_resized", as.integer(args$resizeLength), ".rds")
}
#Save Normalized coverage
	saveRDS(normcov, file=outn)
cat("Normalized coverage saved\n")


###-----------------
### Timing
###-----------------
cat("Processing Ended:", date(), "\n")
cat("Timing:\n")
proc.time() - ptm

#Done
