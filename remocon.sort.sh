# Author: Jiwoong Kim (jiwoongbio@gmail.com)
#!/bin/bash

directory=`dirname $0`
outputPrefix=$1
sam=$2
samContaminant=$3

if [ -z "$samContaminant" ]; then
	echo 'Usage: ./remocon.sort.sh <output.prefix> <sam> <sam.contaminant>' 1>&2
	echo 1>&2
	exit 1
fi

samtools view -H $sam > $outputPrefix.header.sam
samtools view $sam | LC_ALL=C sort --field-separator=$'\t' -k1,1 | cat $outputPrefix.header.sam - | gzip > $outputPrefix.namesorted.sam.gz

samtools view -H $samContaminant > $outputPrefix.contaminant.header.sam
samtools view $samContaminant | LC_ALL=C sort --field-separator=$'\t' -k1,1 | cat $outputPrefix.contaminant.header.sam - | gzip > $outputPrefix.contaminant.namesorted.sam.gz

perl $directory/remocon.pl -n $outputPrefix.namesorted.sam.gz $outputPrefix.contaminant.namesorted.sam.gz 2> $outputPrefix.contaminant_removed.log | gzip > $outputPrefix.contaminant_removed.sam.gz
