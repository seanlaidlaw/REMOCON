# Author: Jiwoong Kim (jiwoongbio@gmail.com)
#!/bin/bash

directory=`dirname $0`
outputPrefix=$1
index=$2
indexContaminant=$3
threads=$4
fastqFiles=${@:5}

if [ -z "$fastqFiles" ]; then
	echo 'Usage: ./remocon.sh <output.prefix> <bwa.index> <bwa.index.contaminant> <threads> <input.1.fastq> [input.2.fastq]' 1>&2
	echo 1>&2
	exit 1
fi

bwa mem -t $threads $index            $fastqFiles | gzip > $outputPrefix.sam.gz
bwa mem -t $threads $indexContaminant $fastqFiles | gzip > $outputPrefix.contaminant.sam.gz
perl $directory/remocon.pl $outputPrefix.sam.gz $outputPrefix.contaminant.sam.gz 2> $outputPrefix.contaminant_removed.log | gzip > $outputPrefix.contaminant_removed.sam.gz
