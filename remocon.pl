# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
use List::Util qw(sum);
use Getopt::Long;

GetOptions('fastq=s' => \(my $fastqFile = ''));
{
	my $fastqFileReader;
	sub openFastqFileReader {
		open($fastqFileReader, ($fastqFile =~ /\.gz$/) ? "gzip -dc $fastqFile |" : $fastqFile);
	}
	sub getReadNameList {
		my ($readName) = @_;
		my @readNameList = ();
		while(my $line = <$fastqFileReader>) {
			<$fastqFileReader>;
			<$fastqFileReader>;
			<$fastqFileReader>;
			chomp($line);
			$line =~ s/^\@//;
			$line =~ s/\s.*$//;
			$line =~ s/\/[12]$//;
			if($line eq $readName) {
				@readNameList = @readNameList[0, grep {$readNameList[$_ - 1] ne $readNameList[$_]} 1 .. $#readNameList] if(scalar(@readNameList) > 0);
				return @readNameList;
			} else {
				push(@readNameList, $line);
			}
		}
	}
	sub closeFastqFileReader {
		close($fastqFileReader);
	}
}

my ($samFile, $samFileContaminant) = @ARGV;
system("samtools view -S -H $samFile");
open(my $reader,            "samtools view -S -F 2304 $samFile |");
open(my $readerContaminant, "samtools view -S -F 2304 $samFileContaminant |");
openFastqFileReader() if($fastqFile ne '');
my ($numberContaminant, $numberNotContaminant, $numberAmbiguous) = (0, 0, 0);
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
closeFastqFileReader() if($fastqFile ne '');
print STDERR "$_\n" foreach($numberContaminant, $numberNotContaminant, $numberAmbiguous);

sub printLines {
	if($fastqFile ne '') {
		foreach my $readName (getReadNameList($readName)) {
			my @tokenListListContaminant = getTokenListListContaminant($readName);
			my $sumContaminant = sum(0, map {getAlignmentScore(@$_)} grep {($_->[1] & 4) == 0} getTokenListListContaminant($readName));
			$numberContaminant += 1 if($sumContaminant > 0);
		}
	}
	my $sum            = sum(0, map {getAlignmentScore(@$_)} grep {($_->[1] & 4) == 0} @tokenListList);
	my $sumContaminant = sum(0, map {getAlignmentScore(@$_)} grep {($_->[1] & 4) == 0} getTokenListListContaminant($readName));
	if($sumContaminant > $sum) {
		$numberContaminant += 1;
	} else {
		print join("\t", @$_), "\n" foreach(@tokenListList);
		if($sumContaminant > 0) {
			if($sum > $sumContaminant) {
				$numberNotContaminant += 1;
			} else {
				$numberAmbiguous += 1;
			}
		}
	}
}

sub getTokenListListContaminant {
	my ($readName) = @_;
	my @tokenListListContaminant = ();
	while($lineContaminant) {
		chomp($lineContaminant);
		my @tokenList = split(/\t/, $lineContaminant);
		last if($tokenList[0] ne $readName);
		push(@tokenListListContaminant, \@tokenList);
		$lineContaminant = <$readerContaminant>;
	}
	return @tokenListListContaminant;
}

sub getAlignmentScore {
	my %tokenHash = ();
	(@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = @_;
	my @alignmentScoreList = map {$_->[2]} grep {$_->[0] eq 'AS' && $_->[1] eq 'i'} map {[split(/:/)]} @tagTypeValueList;
	return @alignmentScoreList;
}
