# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum);
use Getopt::Long qw(:config no_ignore_case);

GetOptions('h' => \(my $help = ''),
	'n' => \(my $namesorted = ''),
	's=s' => \(my $scores = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl remocon.pl target.sam contaminant.sam > target.contaminant_removed.sam

Options: -h       display this help message
         -n       input SAM files sorted lexicographically by read name
         -s STR   comma-separated following values of
                  1. match score
                  2. mismatch penalty
                  3. colon-separated gap open penalties for deletion and insertion
                  4. colon-separated gap extension penalties for deletion and insertion
                  5. clipping penalty
                  6. unpaired read pair penalty
                  e.g. 1,4,6:6,1:1,5,17

EOF
}
my %scoreHash = ();
if($scores ne '') {
	@scoreHash{'match', 'mismatch', 'gap open', 'gap extension', 'clipping', 'unpaired'} = split(/,/, $scores);
	@scoreHash{'deletion open', 'insertion open'} = split(/:/, $scoreHash{'gap open'});
	@scoreHash{'deletion extension', 'insertion extension'} = split(/:/, $scoreHash{'gap extension'});
}

my ($samFile, $samFileContaminant) = @ARGV;
system("samtools view -H $samFile");
open(my $reader,            "samtools view -F 2304 $samFile |");
open(my $readerContaminant, "samtools view -F 2304 $samFileContaminant |");
my %numberHash = ();
my ($readName, @tokenListList) = ('');
my $lineContaminant = <$readerContaminant>;
while(my $line = <$reader>) {
	chomp($line);
	my @tokenList = split(/\t/, $line);
	if($tokenList[0] ne $readName) {
		printLines() if(scalar(@tokenListList) > 0);
		($readName, @tokenListList) = ($tokenList[0]);
	}
	push(@tokenListList, \@tokenList);
}
printLines() if(scalar(@tokenListList) > 0);
close($reader);
close($readerContaminant);
print STDERR "$_: $numberHash{$_}\n" foreach('contaminant', 'non-contaminant', 'ambiguous');

sub printLines {
	my $sum            = sum(0, map {getAlignmentScore(@$_)} grep {($_->[1] & 4) == 0} @tokenListList);
	my $sumContaminant = sum(0, map {getAlignmentScore(@$_)} grep {($_->[1] & 4) == 0} getTokenListListContaminant($readName));
	if($sumContaminant > $sum) {
		$numberHash{'contaminant'} += 1;
	} else {
		print join("\t", @$_), "\n" foreach(@tokenListList);
		if($sumContaminant > 0) {
			if($sum > $sumContaminant) {
				$numberHash{'non-contaminant'} += 1;
			} else {
				$numberHash{'ambiguous'} += 1;
			}
		}
	}
}

sub getTokenListListContaminant {
	my ($readName) = @_;
	my @tokenListListContaminant = ();
	for(; $lineContaminant; $lineContaminant = <$readerContaminant>) {
		chomp($lineContaminant);
		my @tokenList = split(/\t/, $lineContaminant);
		next if($tokenList[0] lt $readName and $namesorted);
		last if($tokenList[0] ne $readName);
		push(@tokenListListContaminant, \@tokenList);
	}
	return @tokenListListContaminant;
}

sub getAlignmentScore {
	my %tokenHash = ();
	(@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = @_;
	my $alignmentScore = 0;
	if($scores ne '') {
		my $matchLength = 0;
		{
			my $cigar = $tokenHash{'cigar'};
			while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//) {
				my ($length, $operation) = ($1, $2);
				if($operation eq 'M') {
					$alignmentScore += $length * $scoreHash{'match'};
					$matchLength += $length;
				} elsif($operation eq 'I') {
					$alignmentScore -= $length * $scoreHash{'insertion extension'} + $scoreHash{'insertion open'};
				} elsif($operation eq 'D') {
					$alignmentScore -= $length * $scoreHash{'deletion extension'} + $scoreHash{'deletion open'};
				} elsif($operation eq 'S') {
					$alignmentScore -= $scoreHash{'clipping'};
				} elsif($operation eq 'H') {
					$alignmentScore -= $scoreHash{'clipping'};
				}
			}
		}
		if((my $md) = map {$_->[2]} grep {$_->[0] eq 'MD' && $_->[1] eq 'Z'} map {[split(/:/)]} @tagTypeValueList) {
			$alignmentScore -= ($matchLength - sum(0, $md =~ /([0-9]+)/g)) * ($scoreHash{'match'} + $scoreHash{'mismatch'});
		}
		if($tokenHash{'flag'} & 1 and not $tokenHash{'flag'} & 2) {
			$alignmentScore -= $scoreHash{'unpaired'};
		}
	} else {
		if((my $as) = map {$_->[2]} grep {$_->[0] eq 'AS' && $_->[1] eq 'i'} map {[split(/:/)]} @tagTypeValueList) {
			$alignmentScore += $as;
		}
	}
	return $alignmentScore;
}
