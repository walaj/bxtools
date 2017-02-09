[![Build Status](https://travis-ci.org/walaj/bxtools.svg?branch=master)](https://travis-ci.org/walaj/bxtools)

## *bxtools* - Tools for analyzing 10X genomics data

**License:** [MIT][license]

## Note: *bxtools* is an emerging project. If you find an operation that you need that may be in the scope of *bxtools*, please submit an issue report or pull request with the suggested functionality. We are looking for community suggestions for what we might include.

Table of contents
=================

  * [Installation](#installation)
  * [Description](#description)
  * [Components](#components)
    * [Split](#split)
    * [Stats](#stats)
    * [Tile](#tile)
    * [Relabel](#relabel)
    * [Mol](#mol)
    * [Convert](#convert)
  * [Example Recipes](#examples-recipes)
  * [Attributions](#attributions)

Installation
------------

```
git clone --recursive https://github.com/walaj/bxtools
cd bxtools
./configure
make 
make install
```

Description
-----------
*bxtools* is a set of light-weight command line tools for analyzing 10X genomics data. It is built to 
take care of low-level type operations in a 10X-specific way by accounting for the BX tag in 10X data.

Components
----------

#### Split

Split a BAM file by the BX tag.

```
## split a BAM into individual BAMs (called test.<bx>.bam). Don't output tags with < 10 reads
bxtools split $bam -a test -m 10 > counts.tsv

## split a portion of a BAM 
samtools view -h $bam 1:1,000,000-2,000,000 | bxtools split - -a test > counts.tsv

## just get the BX counts and sort by prevalence
bxtools split $bam -x | sort -n -k 2,2 > counts.tsv
```

#### Stats

Collect BX-level statistics from a 10X BAM

```
bxtools stats $bam > stats.tsv
## output is BX  count  median_isize  median_mapq
```

#### Tile

Collect BX-level read counts on a tiled genome
```
## default is 1kb tiles, across entire genome
bxtools tile $bam > counts.bed

## input bed to check (e.g. chr1 only)
samtools view -h $bam 1:1-250,000,000 | bxtools tile - -b chr1.tiles.bed > chr1.tiles.counts.bed
```

#### Relabel
Move the BX barcodes from the ``BX`` tag (e.g. ``BX:ACTTACCGA``) to the read name (e.g. ``qname_ACTTACCGA``)
```
VERBOSE=-v ## print progress
bxtools relabel $bam $VERBOSE > relabeled.bam
```

#### Mol
Get the minimum molecular footprint on the genome as BED file for each MI tag. The 
minimal footprint is defined from the minimum start position to the maximum end position of 
all reads sharing an MI tag. Throws an error message if detects the same MI tag on multiple chromosomes.

The output BED format is chr, start, end, MI, BX, read_count
```
bxtools mol $bam > mol_footprint.bed
```

#### Convert
Switch the alignment chromosome with the BX tag. This is a hack to allow a 10X BAM to be sorted and indexed by BX tag, rather than coordinate. 
Useful for rapid lookup of all BX reads from a particular BX. Note that this switches "-" for "_" to make query possible with ``samtools view``.
This also requires a two-pass solution. The first loop is to get all of the unique BX tags to build the new BAM header. The second makes the switches.
This means that streaming from ``stdin`` is not available.
``
bxtools convert $bam | samtools sort - -o bx_sorted.bam
samtools index bx_sorted.bam
samtools view AGTCCAAGTCGGAAGT_1
``

Example recipes
---------------
#### Get BX level coverage in 2kb bins across genome, ignore low-frequency tags

```
## make a list of bad tags (freq < 100)
samtools view -h $bam 1:1-10,000,000 | bxtools split - -x | awk '$2 < 100' | cut -f1 > excluded_list.txt

## get the coverage, while excluding bad tags (grep: -F literal, -f file, -v exclude)
samtools view -h $bam 1:1-10,000,000 | grep -v -F -f excluded_list.txt | bxtools tile - -w 2000 > bxcov.bed
```

Attributions
------------

This project is developed and maintained by Jeremiah Wala (jwala@broadinstitute.org).

Analysis suggestions and 10X support
* Gavin Ha - Postdoctoral Fellow, Broad Institute
* Srinivas Viswanathan - Oncology Fellow, Dana Farber Cancer Institute
* Chris Whelan - Computational Biologist, Broad Institute
* Tushar Kamath - MD-PhD Student, Harvard Medical School
* Cheng-Zhong Zhang - Assistant Professor, Dana Farber Cancer Institute
* Marcin Imielinski - Assistant Professor, Weill Cornell Medical College
* Rameen Beroukhim - Assistant Professor, Dana Farber Cancer Institute
* Matthew Meyerson - Professor, Dana Farber Cancer Institute

[license]: https://github.com/walaj/bxtools/blob/master/LICENSE
