# DEFOR - Depth and Frequency Based Copy Number Detector


## Table of content
* [Introduction](#introduction)
* [Download](#download)
* [Installing](#installing)
* [Tutorial](#tutorial)
* [Input files](#input-files)
* [Programs](#programs)

## Introduction

DEFOR is a set of programs to identify copy number alternations from tumor/normal pair samples.


## Installing

### Prerequisites
[samtools](http://www.htslib.org/) was required to process the bam files

### Compile the source code

````
cd defor-xxx/
make
````

### (Optional) Move the programs under 'bin' to any folder you like

````
mv bin/* ${INSATLL_DIRECOTRY}
````

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
    calc_depratio -d 10 -w 1000000 -s 1000 <(samtools mpileup -q 10 -f hs37d5.fa XP108N.bam) <(samtools mpileup -q 10 -f ~/genome/human.stripe4/hs37d5.fa /home2/s167968/project/kidney/natgen2012/wes/recal/XP108M.bam) | gzip -c > XP108N.XP108M.gz
    ````

4. Estract allele frequency clusters

5. Calculate the allele frequency for both tumor and normal samples

    ````
    Give an example
    ````




## Programs

Add additional notes about how to deploy this on a live system

## Authors

* **He Zhang**

## Acknowledgments

* Xiaowei Zhan
* Yang Xie
* University of Texas Southwestern Medical Center


