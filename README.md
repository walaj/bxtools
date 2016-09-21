[![Build Status](https://travis-ci.org/walaj/bxtools.svg?branch=master)](https://travis-ci.org/walaj/bxtools)

## *bxtools* - Tools for analyzing 10X genomics data

Installation
------------

```
git clone --recursive https://github.com/walaj/bxtools
cd bxtools
./configure
make 
make install
```

Split
-----

Split a BAM file by the BX tag.

```
## split a BAM into individual BAMs test.<bx>.bam. Don't output tags with < 10 reads
bxtools split -b $bam -a test -m 10 > counts.tsv

## split a portion of a BAM 
samtools view -h $bam 1:1,000,000-2,000,000 | bxtools split -b - -a test > counts.tsv

## just get the BX counts and sort by prevalence
bxtools split -b $bam - -x | sort -n -k 2,2 > counts.tsv
```

Stats
-----

Collect BX-level statistics from a 10X BAM

```
bxtools stats $bam > stats.tsv
## output is BX  count  median_isize  median_mapq
```
