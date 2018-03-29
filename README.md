# DEFOR - Depth and Frequency Based Copy Number Alternation Detector


## Table of content
* [Introduction](#introduction)
* [Installing](#installing)
* [Tutorial](#tutorial)
* [Input files](#input-files)
* [Programs and scripts](#programs-and-scripts)
* [Contact](#contact)

## Introduction

DEFOR is a set of programs to identify copy number alternations from tumor/normal pair samples. In this method, both read depth and allele frequency were considered, and it doesn’t rely on the assumption that there are no large-scale CNAs in tumor cells. In the evaluation using a published dataset, DEFOR has better accuracy than the other methods especially in the situation where there are large-scale CNAs in the tumor genome.


## Installing

* Download
    ```
    git clone https://github.com/drzh/defor.git
    ```

* Compile the source code

    ```
    cd defor
    make
    ```

* (Optional) Copy the programs under 'bin' to any folder you like

    ```
    cp bin/* ${INSATLL_DIRECOTRY}
    ```

## Tutorial

1. Download the test dataset and decompress it

    ```
    wget https://cloud.biohpc.swmed.edu/index.php/s/NguCFR0bPP6T198/download
    tar zxf test.tar.gz
    cd test
    ```

3. Estimate depth ratio for tumor/normal pair

    ```
    ../bin/calc_deprat test_normal.mpileup test_tumor.mpileup > test_normal_tumor.dep
    ```

4. Calculate the allele frequency for both tumor and normal samples

    ```
    cat test_normal.mpileup | ~/script/package/cna/script/calc_freq.pl -d 30 -f 0.01 -F 0.99 > test_normal.freq
    cat test_tumor.mpileup | ~/script/package/cna/script/calc_freq.pl -d 30 -f 0.01 -F 0.99 > test_tumor.freq 
    ```

5. Estimate allele frequency clusters

    ```
    ../bin/calc_nclust test_normal.freq test_tumor.freq > test_normal_tumor.nclust
    cat test_normal_tumor.nclust | ../script/extract_nclust.pl > test_normal_tumor.nclust.seg
    ```
    
6. Estimate the copy number alterations

    ```
    ../bin/calc_cna.pl -r test_normal_tumor.dep -c test_normal_tumor.nclust.seg > test_normal_tumor.cna
    ```

## Input files

1. mpileup format

    DEFOR can take the mpileup files as input directly

2. bam or sam files
    To use bam or sam files as input, you have two options:
    1. Convert the bam or sam files to mpileup files using samtools, and then use the mpileup file as input. Here is an example:
           
        ```
        samtools mpileup -q 10 -d 200 -f hs37d5.fa test_normal.bam > test_normal.mpileup
        ```
            
    2. Pipe the output from samtools to DEFOR programs. Here is an example:
 
        ```
        ../bin/calc_deprat <(samtools mpileup -q 10 -f hs37d5.fa test_normal.bam) <(samtools mpileup -q 10 -f hs37d5.fa test_tumor.bam) > test_normal_tumor.dep
        ```
            
    
## Output file

The final output file from the whole pipeline is composed of six columns. Here is an example:

    2       212989532       243175930       loh       1       -0.44
    3       6521842         94773340        loh       1       -0.45
    3       94773341        96718082        loss      0       -0.64
    3       96718083        197947755       normal    2       -0.00
    4       143643221       171526951       loh-amp   2       -0.01
    4       176649734       191024530       normal    2       0.04
    5       5038810         180903064       amp       3       0.38
    6       304501          58779245        normal    2       -0.00
    
The meaning of each column:
* The 1st column: chromosome
* The 2nd column: start position (1-based)
* The 3rd column: end position (1-based)
* The 4th column: copy number event occured in this region
    * normal - normal status, copy number is two, no CNA
    * loss - loss of both alleles, copy number is zero
    * loh - loss of heterozygosity, loss of one allele, copy number is one
    * amp - amplification of one or both alleles, copy nubmber is greater than two
    * loh-amp - loss of heterozygosity followed by amplification, copy number is two
* The 5th column: the expected copy number if the purity of tumor is 100%
* The 6th column: the observed copy number, the value was normalized as log(copy number / 2)

## Programs and scripts

#### calc_deprat  --  Calculate the depth ratio between tumor and normal

    Usage: calc_deprat [options] mpileup_file_normal mpileup_file_tumor > normal_tumor.dep

         Options:
           -w : window size [1000000]
           -s : step size [1000]
           -d : minimal depth [10]
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

#### calc_cna.pl  -- Estimate the copy numbers 

    Usage: calc_cna.pl [options] -r normal_tumor.dep -c normal_tumor.nclust.seg > normal_tumor.cna
    
         Options:
           -D : minimum allele frequency difference to distinguish two clusters [0.15]
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

