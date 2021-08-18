# Scope

This repository contains functions used to analyze ChIP-seq data. They were used in the following articles:  

  - Yu X\*, Martin PGP\*, Michaels SD. (\* co-first authors). BORDER proteins protect expression of neighboring genes by promoting 3' Pol II pausing in plants. [Nat Commun](https://rdcu.be/cux4Z). 2019 Sep 25;10(1):4359. doi: 10.1038/s41467-019-12328-w. 
  - Yu X\*, Martin PGP\*, Zhang Y, Trinidad JC, Xu F, Huang J, Thum KE, Li K, Shao S, Gu Y, Wang X, Michaels SD. (\* co-first authors). *Manuscript in preparation*  

We usually prepare our sequencing libraries with [NEBNext Ultra II](https://international.neb.com/products/e7645-nebnext-ultra-ii-dna-library-prep-kit-for-illumina) DNA library prep kits (NEB E7645) and sequence them in paired end mode in order to have information on the fragment sizes.

# General workflow

## Initial QC
QC on raw reads are obtained with [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

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

