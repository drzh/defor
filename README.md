# DEFOR - Depth and Frequency Based Copy Number Alternation Detector


## Table of content
* [Introduction](#introduction)
* [Download](#download)
* [Installing](#installing)
* [Tutorial](#tutorial)
* [Input files](#input-files)
* [Programs and scripts](#programs-and-scripts)
* [Contact](#contact)

## Introduction

DEFOR is a set of programs to identify copy number alternations from tumor/normal pair samples.


## Installing

* Prerequisites

    [samtools](http://www.htslib.org/) was required if the bam files were used as input.

* Compile the source code

    ````
    cd defor-xxx/
    make
    ```

* (Optional) Copy the programs under 'bin' to any folder you like

    ```
    cp bin/* ${INSATLL_DIRECOTRY}
    ```

## Tutorial

1. Change to the testing directory

    ````
    cd test
    ````

2. Download genome reference file and bam files

    ````
    wget hs37d5.fa
    wget test_normal.bam
    wget test_tumor.bam
    ````
    
3. Estimate depth ratio for tumor/normal pair

    ````
    ../bin/calc_deprat <(samtools mpileup -q 10 -f hs37d5.fa test_normal.bam) <(samtools mpileup -q 10 -f hs37d5.fa test_tumor.bam) > test_normal_tumor.dep
    ````


4. Calculate the allele frequency for both tumor and normal samples

    ````
    samtools mpileup -q 10 -f hs37d5.fa test_normal.bam | ../bin/calc_freq.pl -d 30 -f 0.01 -F 0.99 > test_normal.freq
    samtools mpileup -q 10 -f hs37d5.fa test_tumor.bam | ../bin/calc_freq.pl -d 30 -f 0.01 -F 0.99 > test_tumor.freq
    ````

5. Estimate allele frequency clusters

    ````
    ../bin/calc_nclust test_normal.freq test_tumor.freq > test_normal_tumor.nclust
    cat test_normal_tumor.nclust | ../script/extract_nclust.pl > test_normal_tumor.nclust.seg
    ````
    
6. Estimate the copy number alterations

    ````
    ../bin/calc_cna.pl -r test_normal_tumor.dep -c test_normal_tumor.nclust.seg > test_normal_tumor.cna
    ````

## Input files

1. mpileup format

    DEFOR can take the input from mpileup files directly

2. bam or sam files

    Bam or sam files can also be used as input files. samtools was required for converting the bam or sam files

## Programs and scripts

#### calc_deprat  --  Calculate the depth ratio between tumor and normal

    Usage: calc_deprat [options] mpileup_file_normal mpileup_file_tumor > normal_tumor.dep

         Options:
           -w : window size [1000000]
           -s : step size [1000]
           -d : maximum allele frequency [0.95]
           -n : maximum iteration [10]

#### calc_freq.pl  --  Calculate the allele frequency

    Usage: cat mpileup_file | calc_freq.pl [options] > [normal/tumor].freq

         Options:
           -d : minimum depth to report [30]
           -f : minimum allele frequency [0.01]
           -F : maximum allele frequency [0.99]

#### calc_nclust  --  Identify the allele frequency clusters

    Usage: calc_nclust [options] normal.freq tumor.freq > normal_tumor.nclust

         Options:
           -w : window size [1000000]
           -f : minimum allele frequency [0.05]
           -F : maximum allele frequency [0.95]
           -I : maximum iteration [10]

#### extract_nclust.pl  -- Extract the allele frequency clusters

    Usage: extract_nclust.pl [options] -i normal_tumor.nclust > test_normal_tumor.nclust.seg
    
         Options:
           -w : window size [1000000]
           -d : minimum allele frequency difference to distinguish two clusters [0.1]
           -f : low allele frequency cutoff [0.1]
           -F : high allele frequency cutoff [0.9]
           -m : lower cutoff of mean frequency for the middle cluster [0.45]
           -M : higher cutoff of mean frequency for the middle cluster [0.6]
           -v : verbose level [0]

## Authors

* He Zhang (he.zhang@utsouthwestern.edu)

## Acknowledgments

* Xiaowei Zhan
* Yang Xie
* University of Texas Southwestern Medical Center

## Contact

Any questions or comments can be sent to He Zhang (he.zhang@utsouthwestern.edu) 

