# Scope

This repository contains functions used to analyze ChIP-seq data. They were used in the following articles:  

  - Yu X\*, Martin PGP\*, Michaels SD. (\* co-first authors). BORDER proteins protect expression of neighboring genes by promoting 3' Pol II pausing in plants. [Nat Commun](https://rdcu.be/cux4Z). 2019 Sep 25;10(1):4359. doi: 10.1038/s41467-019-12328-w. 
  - Yu X\*, Martin PGP\*, Zhang Y, Trinidad JC, Xu F, Huang J, Thum KE, Li K, Shao S, Gu Y, Wang X, Michaels SD. (\* co-first authors). *Manuscript in preparation*  

We usually prepare our sequencing libraries with [NEBNext Ultra II](https://international.neb.com/products/e7645-nebnext-ultra-ii-dna-library-prep-kit-for-illumina) DNA library prep kits (NEB E7645) and sequence them in paired end mode in order to have information on the fragment sizes.

# Workflow

## Initial QC
QC on raw reads are obtained with [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).  
They can be nicely combined across multiple samples using [multiQC](https://multiqc.info/)

## Trimming
To multiplex samples we often use [NEBNext single index oligos](https://international.neb.com/tools-and-resources/selection-charts/nebnext-multiplex-oligos-selection-chart).  
For each sample, I build a file containing the sequence of potential adapters that may contaminate the reads.  
For single index, the adapter file is built from the [`oligoBase_SingleIndex.fa`](oligoBase_SingleIndex.fa) in which we replace the first oligo with the appropriate index oligo obtained from [`NEBNext_Multiplex_Oligos.fasta`.](NEBNext_Multiplex_Oligos.fasta)  
<br/>
Residual adapter sequences are removed using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)  
For a paired-end run the typical command is:  

    java -jar /<path_to>/trimmomatic-0.36.jar PE \
        -threads ${threadNumber} \
        -phred33 \
        ${fastqfile1} \
        ${fastqfile2} \
        -baseout ${output_directory}/${SampleName}_trimmo.fq.gz \
        ILLUMINACLIP:${adapterFile}:2:25:7:1:true \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:20 \
        ;

It's a good idea to generate QC reports also after trimming.

## Alignment
For ChIP-seq I generally use [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).  
The typical command is:  

    bowtie2 \
        --very-sensitive \
        -p ${threadNumber} \
        -X ${MaxInsertSize} \
        --dovetail \
        -x ${BowtieIndexPrefix} \
        -1 ${trimmed_fastqfile1} \
        -2 ${trimmed_fastqfile2} \
        -S ${output_directory}/${output_Prefix}.sam \
        ;

I generally set `MaxInsertSize=1000`.  
The `--dovetail` argument is especially important with trimmed reads from libraries of small insert size.

## Filtering
I generally filter the aligned reads using [Picard Markduplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) and [samtools](http://www.htslib.org/).  
Typical commands include:

  1. Sort and convert to BAM:

    samtools sort \
        -@ ${threadNumber} \
        -m 3G \
        -o ${output_bamfile} \
        ${input_samfile} \
        ;


  2. Mark duplicates :

    java -jar \
        -Xmx4g \
        /<path_to>/picard.jar MarkDuplicates \
        INPUT=${input_bamfile} \
        OUTPUT=${output_bamfile} \
        METRICS_FILE=${log_directory}/${SampleName}_MarkedDupMetrics.log \
        ASSUME_SORT_ORDER=coordinate \
        ;

  3. Filter reads :

    samtools view \
        -bh \
        -f 3 \
        -F 1804 \
        -q ${MinMapQ} \
        ${input_bamfile_WithMarkedDuplicates} | \
    samtools sort \
        -@ ${threadNumber} \
        -o ${output_filtered_bamfile} -


Read pairs with MapQ < `${MinMapQ}` are dropped (`-q ${MinMapQ}`). I generally set `MinMapQ=2` or `MinMapQ=10` which, in my experience, only minimally affects the filtering and efficiently removes reads that map to multiple locations in the genome.  
Only paired reads and reads mapped in proper pairs are kept retained (`-f 3`)  
Reads with the following characteristics are removed (`-F 1804`):  

  - fails platform/vendor quality checks are removed
  - not primary alignment
  - PCR or optical duplicate
  - read is unmapped
  - mate is unmapped

To adapt these filters to your needs, see [Explain SAM flags](https://broadinstitute.github.io/picard/explain-flags.html)
<br/>
Various QC can be generated for example with:

  - [multiQC](https://multiqc.info/)
  - `samtools idxstats`
  - `samtools flagstat`
  - `samtools stats` followed by the `plot-bamstats` function

The number of reads before (or after) filtering can be easily obtained with:

    ReadsBefore=$(samtools view -c ${bamfile_before_filtering})


## Calculate coverage and generate genome browser tracks

The script [`GetCoveragesFromBAMPE_args.r`](GetCoveragesFromBAMPE_args.r) calculates the coverage of ChIP fragments along the chromosomes.  
Furthermore, it allows:

  - to filter the fragments based on their size (`minFragSize` and `maxFragSize` arguments)  
  - to resize the fragments to a fixed length (`resizeLength`)
  - to select specific chromosomes (e.g. excluding mitochondrial or plastid chromosomes) via the `chromSizes` argument. [Ath_STDchromSizes.rds](Ath_STDchromSizes.rds) is an example of such rds file for Arabidopsis thaliana (use `readRDS` function in R to import this file).
  - to obtain both raw coverage and normalized coverage ("Read-per-million" normalization possibly adjusted by a multiplication factor)

The obtained coverage can be easily exported as a [bigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html) track file using the [rtracklayer](http://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html) package.  
With multiple files, I use the script [`ExportCoverageAsBigWig_args.r`](ExportCoverageAsBigWig_args.r) to automate the task.


## Extract signal at specific windows

The [GenomicRanges](http://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html) package from [Biocconductor](http://www.bioconductor.org/) makes it straightforward to extract the signal from a genome-wide coverage at specific windows (e.g. around the TSS, along the gene bodies, etc.).  
Examples of such "windows of interest" (WOI) for Arabidopsis thaliana are provided in the [Ath_TAIR10_woi](Ath_TAIR10_woi) folder.  
The script [`profcomp_args.r`](profcomp_args.r) allows to automate the process for one or several signal files and takes care to reverse the signal based on the strand indicated in the WOI file.  


## Signal binning

We the signal is extracted from windows of different sizes (e.g. gene bodies), it may be necessary to calculate the average signal over bins of variable sizes (e.g. 100 bins along the gene body).  
Even for windows of fixed sizes (e.g. TSS +/- 2kb), it can be useful to summarize the signal over bins of fixed sizes (e.g. 10bp bins) when bp-level signal is not necessary (e.g. to plot heatmaps).  
The script [BinFeatureProfiles_args.r](BinFeatureProfiles_args.r) performs such "signal binning".  
<br/>
In addition the script [BinFeatureProfiles_WithTrimming_args.r](BinFeatureProfiles_WithTrimming_args.r) is a slighty modified version that performs the binning only on a specific part of the windows.  
For example if we extracted the signal at TSS +/- 2kb, we have 4001 values for each gene. We may want to bin only the first 2000 values corresponding to the signal before the TSS.  
<br/>
After binning, it is possible to reassemble profiles in order to plot e.g. metagene profiles or heatmaps. For example, we can assemble :

  1. The signal (without binning) for 2kb upstream of the TSS up to the TSS
  2. The binned signal along the gene body
  3. The signal (without binning) for the TES (transcript end site sometimes refered to as PAS for polyadenylation site) up to 2kb after the TES

This generates data that can be stored in a matrix (genes in row and positions in columns) or in an RleList.


## Calculate average profiles (and confidence intervals) for specific groups of genes

The script [`GetAvgProfileForGeneGroups_args.r`](GetAvgProfileForGeneGroups_args.r) calculates average profiles and the associated confidence intervals for different groups of genes (provided as a named list of gene names).  
The scripts accepts a limited number of arguments for simplicity but the functions within the script accept more options.  
In particular:

  - By default the data matrix is winsorized to reduce the effect of extreme values above the 99.99th percentile. This can be adjusted as needed to limit (or not) the impact of extremely high or low values.  
  - Confidence intervals can be calculated either based on a normal assumption or using bootstrap as in the [ChIPseeker](http://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html) package. Practically, they're expected to give very simiar results when the number of genes is large enough but the bootstrap CI are much longer to calculate.  
  - By default, the profiles are not normalized within the defined window but two normalization options are available which may be necessary to compare profiles from different ChIPs (`"freq"`) or for windows of different width (`"avg"`)


## Plotting

The data is generally plotted as:

  - "metagene profiles" representing the average signal for groups of genes along with confidence intervals. I use [ggplot2](https://ggplot2.tidyverse.org/) to do this
  - "heatmaps". I use [EnrichedHeatmap](https://bioconductor.org/packages/release/bioc/html/EnrichedHeatmap.html) / [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html) to do this

