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
bxtools split $bam - -x | sort -n -k 2,2 > counts.tsv
```

Stats
-----

Collect BX-level statistics from a 10X BAM

```
bxtools stats $bam > stats.tsv
## output is BX  count  median_isize  median_mapq
```

Tile
----

Collect BX-level read counts on a tiled genome
```
bxtools tile $bam > counts.bed

## input bed to check (e.g. chr1 only)
samtools view -h $bam 1:1-250,000,000 | bxtools tile - -b chr1.tiles.bed > chr1.tiles.counts.bed
```

### Example recipes
##### Get BX level coverage in 1kb bins across genome, ignore low-frequency tags
```
## make a list of bad tags (freq < 100)
samtools view -h $bam 1:1-10,000,000 | bxtools split - -x | awk '$2 < 100' | cut -f1 > excluded_list.txt

## get the coverage, while excluding bad tags (grep: -F literal, -f file)
samtools view -h $bam 1:1-10,000,000 | grep -F -f excluded_list.txt | bxtools tile - -w 1000 > bxcov.bed
```

```
