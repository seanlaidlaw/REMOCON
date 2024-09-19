# REMOCON

REMOve CONtaminant reads before variant calling

This is a fork of [REMOCON](https://github.com/jiwoongbio/REMOCON.git).

I have modified the perl code from the original tool to produce output similar
to that from [Xenome](https://github.com/data61/gossamer/blob/master/docs/xenome.md),
so it will produce a separate sam file for each potential source of reads from:
- human.sam
- mouse.sam
- ambiguous.sam

For complete readme, see the original repo, but with my changes the workflow is
as follows:

## Usage

1. Before use of REMOCON, align your sample's reads to both human and mouse reference genomes. You should also have a VCF of the SNPs present in the sample.
```
# assuming you have these already present in directory
sample.SNPs.vcf
sample.humanaligned.sam
sample.mousealigned.sam
```

2. calculate an alignment score to measure the quality of the alignment of each read to its reference, in human and in mouse. This is the basis on which the read classification will take place
```{bash}
perl remocon_alignment_score.pl -b sample.humanaligned.sam sample.SNPs.vcf > sample.humanaligned.alignmentedscore.sam
perl remocon_alignment_score.pl -b sample.mousealigned.sam sample.SNPs.vcf > sample.mousealigned.alignmentedscore.sam
```

3. Classify the reads into mouse, human, and ambiguous.
```{bash}
perl remocon.pl "sample.human.alignmentedscore.sam" "sample.mouse.alignmentedscore.sam"
```

This will produce 3 output files, each containing the reads assigned to that origin:
```
human.sam
mouse.sam
ambiguous.sam
```
