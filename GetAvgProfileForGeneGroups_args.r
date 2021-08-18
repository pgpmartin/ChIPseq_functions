#!/usr/bin/env Rscript

# The function expects the following arguments:
    # prof: path to an rds file containing an RleList or a matrix of profiles
    # featuregroups: path to an rds file containing a named list of featuregroups for which the average profiles should be computed
    # pos: (optional) an rds file containing a vector of numeric positions for the columns/elements of each profile / or 2 integers for the start and end positions
	# outputdir: (default: getwd()) Output directory
	# outFileName: (optional) Sample Name (will be the name of the output rds file)
    # Ncores: number of cores to use

###-----------------
### Libraries
###-----------------
suppressMessages(require(R.utils))
suppressMessages(require(IRanges))
suppressMessages(require(BiocParallel))

###-----------------
### Functions
###-----------------

## Convert RleList to matrix
RleList2matrix <- function(rlelist)
{
matrix(as.numeric(unlist(rlelist, use.names = FALSE)),
        nrow = length(rlelist),
        byrow = TRUE,
        dimnames = list(names(rlelist), NULL))
}

## Perform winsorization on the matrix
FilterMatrixForQuantiles <- function(mat, 
                                     low = 1, 
                                     high = 99, 
                                     verbose = FALSE,
                                     warningThreshold = 20)
    {
#Function to filter the uppest and lowest values of a numeric matrix (winsorization)
#All the values that are above quantile(mat,high/100) are set to quantile(mat,high/100)
#All the values that are below quantile(mat,low/100) are set to quantile(mat,low/100)
#If verbose=TRUE, the function will print: 
    # the total number of values replaced by the corresponding quantiles,
    # a summary of the percentage of values replaced in each column
    # a warning if the % of values replaced by a quantile in one column exceeds warningThreshold

#Some controls
    if (length(low)!=1 | length(high)!=1) {
        stop("high and low must be of length 1")
    }

    if (is.numeric(low) & is.numeric(high) & low>high) {
        stop("low must be lower than high")
    }

#Initialize
    res = as.matrix(mat)
    quantlow = quanthigh=NA

#Get lower quantile
    if (is.numeric(low)){
        if (low<0 | low>100){
            stop("lower threshold must be in [0:100]")
        } else {
            quantlow <- quantile(res,probs=low/100,na.rm=T)
        }
    }

#Get upper quantile
    if (is.numeric(high)){
        if (high<0 | high>100){
            stop("higher threshold must be in [0:100]")
        } else {
        quanthigh <- quantile(res,probs=high/100,na.rm=T)
        }
    }

#Replace extreme values by lower/upper quantile
    if (is.numeric(quantlow)) {
        res[res<quantlow & !is.na(res)] <- quantlow
    }

    if (is.numeric(quanthigh)) {
        res[res>quanthigh & !is.na(res)] <- quanthigh
    }

#Verbose
    if (verbose){
    cat("Number of values replaced by the ", 
        high, "th percentile: ",
        length(mat[mat>quanthigh & !is.na(mat)]),
        "\n", sep="")
    cat("Number of values replaced by the ", 
        low, "th percentile: ",
        length(mat[mat<quantlow & !is.na(mat)]),
        "\n", sep="")
    diffmats <- abs(mat-res) > 1e-15
    percValuesChanged <- colMeans(diffmats,na.rm=T)*100
    cat("Summary for the percentage of replaced values by columns:\n")
    print(summary(percValuesChanged))
    AboveWarningThresh <- percValuesChanged>=warningThreshold
    if (any(AboveWarningThresh))
        cat("Columns with more than ", warningThreshold, 
            "% of the values replaced by quantiles: ", 
            paste0(which(AboveWarningThresh), collapse=", "),
            "\n")
	}
	return(res)
	}


GetAverageProfileWithCI <- function(FeatureProfiles,
                                    selFeatures = NA,
                                    xlim = NA,
                                    conf = 0.95,
                                    conftype = c("boot","basic"),
                                    resample = 500,
                                    ncores = NA,
                                    normtype = c("none","freq","avg"),
                                    filterOutliers = TRUE,
                                    LowerQuantile = 0,
                                    UpperQuantile = 99.99)
{
#FeatureProfiles can be provided as a matrix or an RleList
#conftype: bootstrap or basic
#conf is the confidence interval (default to 95%)
#resample gives the number of samples in the bootstrap procedure
#The function returns how many values are filtered by FilterOUtliers
#normtype freq = divide colsums by total sum
#normtype avg = divide colsums by total sum and multiply by ncol
#selFeatures allows to select a specific group of genes


#libraries
	suppressMessages(require(IRanges))
	suppressMessages(require(boot))
	suppressMessages(require(parallel))
	suppressMessages(require(BiocParallel))

#Set the number of cores
    if (is.na(ncores)) {
        availcores <- detectCores()
        ncores <- if (!is.na(availcores) && availcores>1) availcores-1 else 1
    }

    register(MulticoreParam(ncores))

#Tests on arguments
    if (!is(FeatureProfiles,"RleList") && !is.matrix(FeatureProfiles)) {
        stop("FeatureProfiles should be an RleList or a matrix")
    }

    if (length(normtype)!=1) {
        stop ("normtype should be ONE of 'none','freq' or 'avg'")
    }

    npos = if (is(FeatureProfiles,"RleList")) length(FeatureProfiles[[1]]) else ncol(FeatureProfiles)
    if (any(is.na(xlim))) {
        xlim <- c(1,npos)
    }
    
    if ((length(xlim)==2 && (xlim[2]-xlim[1]+1)!=npos) || (length(xlim)>2 && length(xlim)!=npos)) {
        stop("xlim range is not equal to the number of columns in FeatureProfiles")
    }

# Data selection and formatting
    ## If selFeatures is NA or NULL, select all features
    if (is.na(selFeatures) || is.null(selFeatures)) {
        selFeatures <- if (is(FeatureProfiles,"RleList")) 1:length(FeatureProfiles) else 1:nrow(FeatureProfiles)
    }

    ## Select the Features if FeatureProfiles is a matrix
    if (is.matrix(FeatureProfiles)) {
        FeatureProfiles <- FeatureProfiles[selFeatures,]
    }

    ## Convert RleList to a matrix and select the features
    if (is(FeatureProfiles,"RleList")){
        nc <- elementNROWS(FeatureProfiles[selFeatures])
        if (sd(nc)!=0) {
            stop("Features profiles in RleList have different lengths")
        }
        FeatureProfiles <- RleList2matrix(FeatureProfiles[selFeatures])
	}

# (Optional) Filtering for extreme values
    #if filterOUtliers=T replace the outliers by the corresponding quantiles
    #We could also choose to remove the corresponding values and introduce NA but this would complicate the avgFUN below
    if (filterOutliers) {
        cat("Filtering data matrix to remove extreme values\n")
        FeatureProfiles <- FilterMatrixForQuantiles(FeatureProfiles,
                                                    low = LowerQuantile,
                                                    high = UpperQuantile,
                                                    verbose = TRUE,
                                                    warningThreshold = 20)
    }


# Define the function applied to the columns
    #mean by columns (normtype="none")
    if (normtype=="none"){
        avgFUN <- function(data,idx){
        d <- data[idx,]
        colMeans(d, na.rm = TRUE) #faster than colSums(d)/nrow(d)
        }
    }

    #Frequence (Proportion of signal in each column)
    #identical to ChIPseeker getSgn function (see https://github.com/GuangchuangYu/ChIPseeker/blob/master/R/utilities.R)
    #note that if we change the number of columns (i.e. we change the window width) the values will change only because the grand sum changes
    #it is OK to use this to compare profiles of identical width
    if (normtype=="freq"){
        avgFUN <- function(data,idx){
        d <- data[idx,]
        ss <- colSums(d, na.rm = TRUE)
        ss <- ss / sum(ss, na.rm = TRUE)
        }
    }

    #Average signal by column divided by grand mean of the signal (analog to "freq" but independent from the size of the window)
    #use this rather than "freq" to compare profiles of different width
    if (normtype=="avg"){
        avgFUN <- function(data,idx){
        d <- data[idx,]
        ss <- colSums(d, na.rm = TRUE)
        ss <- ss * ncol(d) / sum(ss, na.rm = TRUE)
        }
    }


# Apply the function on the whole dataset
    avgprof <- avgFUN(FeatureProfiles,
                      1:nrow(FeatureProfiles))

# Get confidence interval
    ## Basic confidence interval based on a normal assumption 
    if (conftype=="basic") {
        avgprof_sem=apply(FeatureProfiles, 2, sd, na.rm = TRUE) / sqrt(nrow(FeatureProfiles))

            if (normtype=="freq"){
                avgprof_sem=avgprof_sem*nrow(FeatureProfiles) / sum(FeatureProfiles, na.rm = TRUE)
            }

            if (normtype=="avg"){
                avgprof_sem=avgprof_sem / mean(FeatureProfiles, na.rm=TRUE)
            }

        alpha <- 1-conf
        CI=cbind.data.frame(Upper = avgprof + qnorm(1-alpha/2) * avgprof_sem,
                            Lower = avgprof + qnorm(alpha/2) * avgprof_sem)
    }

    ## Confidence intervals based on bootstrap (same approach as ChIPseeker)
    if (conftype=="boot"){
        bootsamples <- boot(FeatureProfiles, statistic = avgFUN,
                            R = resample,
                            parallel = "multicore", ncpus = ncores)
        CI <- t(sapply(1:ncol(FeatureProfiles),
                       function(i){
                           rev(boot.ci(bootsamples,
                                       conf=conf,type="perc",
                                       index=i)$percent)[1:2]
                        }))
        colnames(CI)=c("Upper","Lower")
    }

# Format The result as a data frame

    res <- data.frame(Position = if (length(xlim)==2) (xlim[1]:xlim[2]) else xlim,
                      AverageProfile = avgprof,
                      CI)

return(res)
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
                          defaults=list(prof = NULL,
                                        featuregroups = NULL,
                                        posi = NULL,
                                        outputdir = getwd(),
                                        outFileName = paste0("GeneGroupProfiles_",
                                                            gsub(" ", "_", Sys.time()),
                                                            ".rds"),
                                        Ncores = 1))


#Evaluate prof argument
if (is.null(args$prof)) {
  stop('No Profiles provided')
} else if (!file.exists(args$prof)) {
  stop(args$prof," not found")
}

#Verbose
cat("List of arguments for the call:\n")
args

###-----------------
### Set number of cores
###-----------------
register(MulticoreParam(args$Ncores))


###-----------------
### Import profiles
###-----------------
cat("Extracting average profiles from: ",args$prof,"\n")

prof <- readRDS(args$prof) #import prof

if (!is(prof,"RleList") && !is.matrix(prof)) {
    stop("prof should be an RleList or a matrix")
}

#Convert prof to a matrix
if (is(prof,"RleList")) {
    prof <- RleList2matrix(prof)
}

###-----------------
### Import featuregroups
###-----------------

if (is.null(args$featuregroups)) {
    gg <- list(All=rownames(prof))
}

if (!is.null(args$featuregroups)) {
  if (!file.exists(args$featuregroups)) {
    stop(args$featuregroups, " not found")
  } else {
    gg <- readRDS(args$featuregroups)
}
}

if (!is(gg,'list') || is.null(names(gg))) {
  stop("featuregroups should be a named list (possibly of length 1)")
}

if (length(unique(names(gg))) < length(gg)) {
  warning("names of groups in featuregroups are not unique")
}

#Verbose
cat("Computing Average profiles for the lists in", 
    args$featuregroups, "\n") 


###-----------------
### Filter featuregroups for features with profiles
###-----------------

#Number of features provided
cat("Initial number of features in each group (as provided):\n", 
    sapply(gg,length), "\n\n") 

#Remove genes that have no profile
ggf <- lapply(gg, intersect, rownames(prof))

#Number of features in each group
ggfn <- sapply(ggf, length)

if (any(ggfn==0)) {
    warning("After filtering, the following groups have 0 features:\n",
            names(ggf)[ggfn==0], "\n\n")
    ggf <- ggf[ggfn>0] #Remove groups with 0 features
    ggfn <- ggfn[ggfn>0]
}

#Verbose
cat("Final number of features in each group (after filtering):\n", ggfn,"\n\n") 

###-----------------
### Define xlims
###-----------------

if (is.null(args$posi)) {
  xlims <- NA
} else {
    if (!file.exists(args$posi)) {
      cat(args$posi, " not found. Working with default position values")
      xlims <- NA
      } else {
          xlims <- readRDS(args$posi)
        }
  }

cat('ONE: length of xlims is: ',length(xlims),'\n')
cat('ONE: head of xlims is: ', head(xlims),'\n')

if (any(is.na(xlims))) {
  cat(args$posi," contains missing values. Working with default position values")
  xlims <- NA
}

if (!any(is.na(xlims)) & !is.numeric(xlims)) {
  cat(args$posi," is not a numeric vector. Working with default position values")
  xlims <- NA
}

if (!any(is.na(xlims))) {
  npos=ncol(prof)
  if ((length(xlims)==2 && (xlims[2]-xlims[1]+1)!=npos) || (length(xlims)>2 && length(xlims)!=npos)) {
    cat("Range of ", args$posi, " is not equal to the number of positions in ",args$prof,
            "\n  Working with default position values")
    xlims <- NA
  }
}

cat('length of xlims is: ',length(xlims),'\n')
cat('head of xlims is: ', head(xlims),'\n')

###-----------------
### Get average profiles for each group
###-----------------

avgprofs <- list()
avgprofs <- lapply(ggf,function(x){GetAverageProfileWithCI(prof[x,],
                                                            xlim = xlims,
                                                            conf = 0.95,
                                                            conftype = 'basic',
                                                            normtype = 'none',
                                                            ncores = args$Ncores)})
names(avgprofs) <- names(ggf)


###-----------------
### Save average profiles
###-----------------

saveRDS(avgprofs,
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

