# REMOCON

REMOve CONtaminant reads before variant calling


## Method


## Requirements

1. Perl - https://www.perl.org
2. BWA - http://bio-bwa.sourceforge.net
3. SAMtools 1.x - http://www.htslib.org
4. Common linux commands: bash, gzip, ...


## Install

If you already have Git (https://git-scm.com) installed, you can get the latest development version using Git.
```
git clone https://github.com/jiwoongbio/REMOCON.git
```


## Usages

1. Prepare BWA index files
  ```
  bwa index <genome.fasta>
  bwa index <genome.contaminant.fasta>
  ```

2. Remove contaminant reads
  * Use **remocon.sh**
  ```
  ./remocon.sh <output.prefix> <genome.fasta> <genome.contaminant.fasta> <threads> <input.1.fastq> [input.2.fastq]
  ```
  * Use **remocon.pl**
  ```
  bwa mem <genome.fasta>             <input.1.fastq> [input.2.fastq] | gzip > output.sam.gz
  bwa mem <genome.contaminant.fasta> <input.1.fastq> [input.2.fastq] | gzip > output.contaminant.sam.gz
  perl remocon.pl output.sam.gz output.contaminant.sam.gz | gzip > output.contaminant_removed.sam.gz
  ```
  * Use **remocon.sort.sh** - take SAM files as input instead of FASTQ files
  ```
  ./remocon.sort.sh <output.prefix> <input.sam> <input.contaminant.sam>
  ```
