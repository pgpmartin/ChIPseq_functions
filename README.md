# Scope

This repository contains functions used to analyze ChIP-seq data. They were used in the following articles:  

  - Yu X**\***, Martin PGP**\***, Michaels SD. (**\*** co-first authors). BORDER proteins protect expression of neighboring genes by promoting 3' Pol II pausing in plants. [Nat Commun](https://rdcu.be/cux4Z). 2019 Sep 25;10(1):4359. doi: 10.1038/s41467-019-12328-w. 
  - Yu X\*, Martin PGP\*, Zhang Y, Trinidad JC, Xu F, Huang J, Thum KE, Li K, Shao S, Gu Y, Wang X, Michaels SD. (\* co-first authors). *Manuscript in preparation*  

We usually prepare our sequencing libraries with [NEBNext Ultra II](https://international.neb.com/products/e7645-nebnext-ultra-ii-dna-library-prep-kit-for-illumina) DNA library prep kits (NEB E7645) and sequence them in paired end mode in order to have information on the fragment sizes.

# General workflow

## Initial QC
QC on raw reads are obtained with [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

## Trimming
To multiplex samples we often use [NEBNext single index oligos](https://international.neb.com/tools-and-resources/selection-charts/nebnext-multiplex-oligos-selection-chart).  
For each sample, we build a file containing the sequence of potential adapters that may contaminate the reads.  
For single index, the adapter file is built from the `oligoBase_SingleIndex.fa` in which we replace the first oligo with the appropriate index oligo obtained from `NEBNext_Multiplex_Oligos.fasta`.  
<br/>
Residual adapter sequences are removed using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)  
For a paired-end run the typical commandd is:  

    java -jar /<path_to>/trimmomatic-0.36.jar PE \
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




 


