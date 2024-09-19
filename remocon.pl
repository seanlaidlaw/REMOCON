#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
# Modified by: Sean Laidlaw
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(sum max);
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
    'h' => \(my $help = ''),
    'n' => \(my $namesort = ''),
    's' => \(my $strict = ''),
    'a' => \(my $useAS = ''),
    'm=s' => \(my $moleculeNameDelimiter = ''),
);

if($help || scalar(@ARGV) == 0) {
    die <<EOF;

Usage:   perl remocon.pl input.sam input.contaminant.sam [...] > output.contaminant_removed.sam

Options: -h       display this help message
         -n       input SAM files sorted lexicographically by query name
         -s       remove ambiguous reads
         -a       use AS, alignment score generated by aligner
         -m STR   molecule name delimiter

EOF
}

# Open output files for human, mouse, both, and ambiguous
open(my $human_fh, '>', 'human.sam') or die "Cannot open human.sam: $!";
open(my $mouse_fh, '>', 'mouse.sam') or die "Cannot open mouse.sam: $!";
open(my $ambiguous_fh, '>', 'ambiguous.sam') or die "Cannot open ambiguous.sam: $!";

my ($samFile, @contaminantSamFileList) = @ARGV;
open(my $reader, "samtools view -h $samFile |");
my @contaminantReaderList = ();
foreach my $contaminantSamFile (@contaminantSamFileList) {
    open(my $contaminantReader, "samtools view $contaminantSamFile |");
    push(@contaminantReaderList, $contaminantReader);
}
my %numberHash = ();
$numberHash{$_} = 0 foreach('contaminant', 'non-contaminant', 'ambiguous');
my ($currentMoleculeName, @tokenListList) = ('');
my @contaminantTokenListList = ();
foreach my $index (0 .. $#contaminantReaderList) {
    $contaminantTokenListList[$index] = [getTokenList($contaminantReaderList[$index])];
}

while(my @tokenList = getTokenList($reader)) {
    if($tokenList[0] =~ /^@..$/) {
        # Write header lines to all files
        print $human_fh join("\t", @tokenList), "\n";
        print $mouse_fh join("\t", @tokenList), "\n";
        print $ambiguous_fh join("\t", @tokenList), "\n";
        next;
    }
    my $moleculeName = getMoleculeName($tokenList[0]);
    if($moleculeName ne $currentMoleculeName) {
        printLines() if(scalar(@tokenListList) > 0);
        ($currentMoleculeName, @tokenListList) = ($moleculeName);
    }
    push(@tokenListList, \@tokenList);
}
printLines() if(scalar(@tokenListList) > 0);
close($reader);
close($_) foreach(@contaminantReaderList);

# Print the final statistics to stderr
print STDERR "$_: $numberHash{$_}\n" foreach('contaminant', 'non-contaminant', 'ambiguous');

# Close all filehandles
close($human_fh);
close($mouse_fh);
close($ambiguous_fh);

sub printLines {
    my %readNameAlignmentScoreHash = getReadNameAlignmentScoreHash(@tokenListList);
    my @contaminantReadNameAlignmentScoreHashList = ();

    # Gather contaminant read alignment scores
    foreach my $index (0 .. $#contaminantReaderList) {
        my @tokenListList = ();
        for(; @{$contaminantTokenListList[$index]}; $contaminantTokenListList[$index] = [getTokenList($contaminantReaderList[$index])]) {
            my $moleculeName = getMoleculeName($contaminantTokenListList[$index]->[0]);
            next if($moleculeName lt $currentMoleculeName and $namesort);
            last if($moleculeName ne $currentMoleculeName);
            push(@tokenListList, $contaminantTokenListList[$index]);
        }
        my %readNameAlignmentScoreHash = getReadNameAlignmentScoreHash(@tokenListList);
        push(@contaminantReadNameAlignmentScoreHashList, \%readNameAlignmentScoreHash);
    }

    # Determine maximum alignment score for both human and contaminant reads
    my $alignmentScore = max(values %readNameAlignmentScoreHash) // 0;
    my $contaminantAlignmentScore = @contaminantReadNameAlignmentScoreHashList
        ? max(map {values %$_} @contaminantReadNameAlignmentScoreHashList) // 0
        : 0;

    # Write to the appropriate file based on alignment scores
    if($contaminantAlignmentScore > $alignmentScore) {
        # Write to the mouse file if the contaminant score is higher
        print $mouse_fh map {join("\t", @$_), "\n"} @tokenListList;
        $numberHash{'contaminant'} += 1;
    } elsif($alignmentScore > $contaminantAlignmentScore) {
        # Write to the human file if the human score is higher
        print $human_fh map {join("\t", @$_), "\n"} @tokenListList;
        $numberHash{'non-contaminant'} += 1;
    } elsif($alignmentScore == $contaminantAlignmentScore) {
        # Write to the ambiguous file if the scores are equal
        print $ambiguous_fh map {join("\t", @$_), "\n"} @tokenListList;
        $numberHash{'ambiguous'} += 1;
    }
}

sub getMoleculeName {
    my ($readName) = @_;
    my $moleculeName = $readName;
    $moleculeName =~ s/$moleculeNameDelimiter.*$// if($moleculeNameDelimiter ne '');
    return $moleculeName;
}

sub getReadNameAlignmentScoreHash {
    my @tokenListList = @_;
    my %readNameNumberAlignmentScoreListHash = ();
    push(@{$readNameNumberAlignmentScoreListHash{$_->[0]}->{($_->[1] & 192) / 64}}, getAlignmentScore(@$_)) foreach(@tokenListList);
    my %readNameAlignmentScoreHash = ();
    foreach my $readName (keys %readNameNumberAlignmentScoreListHash) {
        $readNameAlignmentScoreHash{$readName} = sum(map {max(@$_)} values %{$readNameNumberAlignmentScoreListHash{$readName}});
    }
    return %readNameAlignmentScoreHash;
}

sub getAlignmentScore {
    my @tokenList = @_;
    my %tokenHash = ();
    (@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = @tokenList;
    $tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
    return 0 if($tokenHash{'flag'} & 4); # Return 0 if the read is unmapped
    if($useAS) {
        return $tokenHash{'AS:i'} if(defined($tokenHash{'AS:i'}));
        return 0;
    }
    return $tokenHash{'Za:f'} if(defined($tokenHash{'Za:f'}));
    return $tokenHash{'AS:i'} if(defined($tokenHash{'AS:i'}));
    return 0;
}

sub getTokenList {
    my ($reader) = @_;
    while(my $line = <$reader>) {
        chomp($line);
        next if($line =~ /^#/); # Skip comment lines
        return split(/\t/, $line, -1);
    }
    return ();
}
