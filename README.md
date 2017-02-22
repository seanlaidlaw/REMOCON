# REMOCON

REMOve CONtaminant reads before variant calling


## Method


## Requirements

1. Perl - https://www.perl.org
2. BWA - http://bio-bwa.sourceforge.net
3. SAMtools - http://samtools.sourceforge.net, http://www.htslib.org
4. Common linux commands: bash, gzip, ...


## Install

If you already have Git (https://git-scm.com) installed, you can get the latest development version using Git.
```
git clone https://github.com/jiwoongbio/REMOCON.git
```


Usages
------

1. Prepare BWA index files
```
bwa index <genome.fasta>
bwa index <genome.contaminant.fasta>
```

2. Remove contaminant reads
```
./remocon.sh <output.prefix> <genome.fasta> <genome.contaminant.fasta> <threads> <input.1.fastq> [input.2.fastq]
```
