#!/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

my %scoreHash = ();
GetOptions(
	'h' => \(my $help = ''),
	's' => \(my $sort = ''),
	'n' => \(my $namesort = ''),
	'R=f' => \($scoreHash{'match'} = 1),
	'P=f' => \($scoreHash{'mismatch'} = 4),
	'G=s' => \($scoreHash{'gap open'} = '6,6'),
	'E=s' => \($scoreHash{'gap extension'} = '1,1'),
	'b' => \(my $bismark = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl remocon_alignment_score.pl input.sam > output.alignment_score_added.sam
         perl remocon_alignment_score.pl input.sorted.sam input.sorted.vcf > output.alignment_score_added.sam

Options: -h       display this help message
         -s       sort input SAM file by coordinate
         -n       output SAM file sorted lexicographically by query name
         -R FLOAT match reward [$scoreHash{'match'}]
         -P FLOAT mismatch penalty [$scoreHash{'mismatch'}]
         -G STR   comma-separated gap open penalties for deletions and insertions [$scoreHash{'gap open'}]
         -E STR   comma-separated gap extension penalties for deletions and insertions [$scoreHash{'gap extension'}]

EOF
}
@scoreHash{'deletion open', 'insertion open'} = split(/,/, $scoreHash{'gap open'});
@scoreHash{'deletion extension', 'insertion extension'} = split(/,/, $scoreHash{'gap extension'});

my ($samFile, $vcfFile) = @ARGV;
my %chromosomeIndexHash = ();
{
	my $reader;
	my @columnList = ();
	my %tokenHash = ();
	sub setTokenHash {
		while(my $line = <$reader>) {
			chomp($line);
			@tokenHash{@columnList} = split(/\t/, $line);
			return if($tokenHash{'FILTER'} eq '.' || $tokenHash{'FILTER'} eq 'PASS');
		}
		%tokenHash = ();
	}
	if(defined($vcfFile)) {
		open($reader, $vcfFile);
		while(my $line = <$reader>) {
			chomp($line);
			last if($line !~ /^##/ && $line =~ s/^#// && (@columnList = split(/\t/, $line)));
		}
		setTokenHash;
	}
	my $chromosomeIndex = 0;
	my $firstPosition = 0;
	my @variantHashList = ();
	sub setChromosomePosition {
		my ($chromosome, $position) = @_;
		if($chromosomeIndexHash{$chromosome} != $chromosomeIndex) {
			$chromosomeIndex = $chromosomeIndexHash{$chromosome};
			$firstPosition = 0;
			@variantHashList = ();
		}
		if($chromosomeIndexHash{$chromosome} == $chromosomeIndex && $position > $firstPosition) {
			@variantHashList = splice(@variantHashList, $position - $firstPosition);
			$firstPosition = $position;
		}
		while(%tokenHash && $chromosomeIndexHash{$tokenHash{'CHROM'}} < $chromosomeIndex) {
			setTokenHash;
		}
	}
	sub setVariantPosition {
		my ($variantPosition) = @_;
		while(%tokenHash && $chromosomeIndexHash{$tokenHash{'CHROM'}} == $chromosomeIndex && $tokenHash{'POS'} <= $variantPosition) {
			foreach my $altBase (split(/,/, $tokenHash{'ALT'})) {
				my ($position, $refBase) = @tokenHash{'POS', 'REF'};
				if($refBase ne $altBase) {
					while(substr($refBase, -1, 1) eq substr($altBase, -1, 1)) {
						substr($refBase, -1, 1, '');
						substr($altBase, -1, 1, '');
					}
					while(substr($refBase, 0, 1) eq substr($altBase, 0, 1)) {
						$position = $position + 1;
						substr($refBase, 0, 1, '');
						substr($altBase, 0, 1, '');
					}
					my $index = $position - $firstPosition;
					$variantHashList[$index]->{join("\t", $position, $refBase, $altBase)} = 1 if($index >= 0);
				}
			}
			setTokenHash;
		}
	}
	sub hasVariant {
		my ($position, $refBase, $altBase) = @_;
		setVariantPosition($position);
		my $index = $position - $firstPosition;
		return $index >= 0 && defined($variantHashList[$index]->{join("\t", $position, $refBase, $altBase)});
	}
	sub closeVcfFileReader {
		close($reader) if(defined($vcfFile));
	}
}
open(my $reader, $sort ? "samtools sort $samFile | samtools view -h - |" : "samtools view -h $samFile |");
{
	my $line = <$reader>;
	while(defined($line) && $line =~ /^\@/) {
		chomp($line);
		my @tokenList = split(/\t/, $line);
		if($tokenList[0] eq '@HD') {
			@tokenList = map {$_ =~ /^SO:/ ? 'SO:queryname' : $_} @tokenList if($namesort);
		}
		if($tokenList[0] eq '@SQ') {
			my %tokenHash = map {$_->[0] => $_->[1]} map {[split(/:/, $_, 2)]} @tokenList[1 .. $#tokenList];
			$chromosomeIndexHash{$tokenHash{'SN'}} = scalar(keys %chromosomeIndexHash);
		}
		print join("\t", @tokenList), "\n";
		$line = <$reader>;
	}
	open(my $writer, $namesort ? "| LC_ALL=C sort --field-separator='\t' -k1,1" : '>-');
	while(defined($line)) {
		chomp($line);
		my @tokenList = split(/\t/, $line);
		my %tokenHash = ();
		(@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = @tokenList;
		$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
		if(($tokenHash{'flag'} & 4) == 0) {
			my $alignmentScore = getAlignmentScore(%tokenHash);
			push(@tokenList, "Za:f:$alignmentScore");
		}
		print $writer join("\t", @tokenList), "\n";
		$line = <$reader>;
	}
	close($writer);
}
close($reader);
closeVcfFileReader;

sub getAlignmentScore {
	my (%tokenHash) = @_;
	setChromosomePosition(@tokenHash{'rname', 'pos'}) if(defined($vcfFile));
	my $alignmentScore = 0;
	my ($cigar, $md, $position, $index) = (@tokenHash{'cigar', 'MD:Z', 'pos'}, 0);
	while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//) {
		my ($length, $operation) = ($1, $2);
		if($operation eq 'M') {
			if($md =~ s/^([0-9]+)//) {
				if($length > $1) {
					$cigar = join('', $length - $1, 'M', $cigar);
					$length = $1;
				} elsif($1 > $length) {
					$md = join('', $1 - $length, $md);
				}
				$alignmentScore += $scoreHash{'match'} * $length;
				$index += $length;
				$position += $length;
			} elsif($md =~ s/^([A-Z]+)0*//) {
				my $refBase = $1;
				if($length > length($refBase)) {
					$cigar = join('', $length - length($refBase), 'M', $cigar);
					$length = length($refBase);
				} elsif(length($refBase) > $length) {
					$md = join('', substr($refBase, $length), $md);
					$refBase = substr($refBase, 0, $length);
				}
				my $altBase = substr($tokenHash{'seq'}, $index, $length);
				if($bismark) {
					foreach my $call (split(//, substr($tokenHash{'XM:Z'}, $index, $length))) {
						if($call ne '.') {
							$alignmentScore += $scoreHash{'match'};
						} else {
							$alignmentScore -= $scoreHash{'mismatch'};
						}
					}
				} elsif(defined($vcfFile) && hasVariant($position, $refBase, $altBase)) {
					$alignmentScore += $scoreHash{'match'} * $length;
				} else {
					$alignmentScore -= $scoreHash{'mismatch'} * $length;
				}
				$index += $length;
				$position += $length;
			} else {
				die "Failed to parse CIGAR/MD";
			}
		} elsif($operation eq 'I') {
				my ($refBase, $altBase) = ('', substr($tokenHash{'seq'}, $index, $length));
				if(defined($vcfFile) && hasVariant($position, $refBase, $altBase)) {
					$alignmentScore += $scoreHash{'match'} * $length;
				} else {
					$alignmentScore -= $scoreHash{'insertion open'} + $scoreHash{'insertion extension'} * $length;
				}
				$index += $length;
		} elsif($operation eq 'D') {
			if($md =~ s/^\^([A-Z]+)0*// && length($1) == $length) {
				my ($refBase, $altBase) = ($1, '');
				if(defined($vcfFile) && hasVariant($position, $refBase, $altBase)) {
#					$alignmentScore += $scoreHash{'match'} * $length;
				} else {
					$alignmentScore -= $scoreHash{'deletion open'} + $scoreHash{'deletion extension'} * $length;
				}
				$position += $length;
			} else {
				die "Failed to parse CIGAR/MD";
			}
		} elsif($operation eq 'N') {
				$position += $length;
		} elsif($operation eq 'S') {
				$index += $length;
		} elsif($operation eq 'H') {
		} else {
				die "Failed to parse CIGAR/MD";
		}
	}
	die "Failed to parse CIGAR/MD" unless($cigar eq '' && $md eq '' && $index == length($tokenHash{'seq'}));
	return $alignmentScore;
}
