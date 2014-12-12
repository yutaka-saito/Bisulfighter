#!/usr/bin/perl

# Bisulfighter http://epigenome.cbrc.jp/bisulfighter
# by National Institute of Advanced Industrial Science and Technology (AIST)
# is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# http://creativecommons.org/licenses/by-nc-sa/3.0/

use strict;
use warnings;

# format ComMet's output to Bisulfighter's output
# (dmrs are placed to both strands)
#
# usage:
# ./this_program ComMetInput > diff-dmr-1-2.txt

# input format 
=pod
chrX  1446851   1447062   UP    8.62835
chrX  193077    193276    UP    3.69225
chrX  25031411  25031919  DOWN  219.966
chrX  9982882   9983394   DOWN  154.95
=cut

=pod
chrX  1446851   1447062   UP    8.62835
chrX  193077    193276    UP    3.69225
chrX  25031411  25031919  DOWN  219.966
chrX  9982882   9983394   DOWN  154.95
=cut

# output format
=pod
chrX  193077    193276    DMR0000001  +  1  0  0  0  NA  0
chrX  193077    193276    DMR0000002  -  1  0  0  0  NA  0
chrX  1446851   1447062   DMR0000003  +  1  0  0  0  NA  0
chrX  1446851   1447062   DMR0000004  -  1  0  0  0  NA  0
chrX  9982882   9983394   DMR0000005  +  0  0  1  0  NA  0
chrX  9982882   9983394   DMR0000006  -  0  0  1  0  NA  0
chrX  25031411  25031919  DMR0000007  +  0  0  1  0  NA  0
chrX  25031411  25031919  DMR0000008  -  0  0  1  0  NA  0
=cut

my $id = 1;

my $infile = $ARGV[0];
open(IN, "sort -k1,1 -k2,2n < $infile | ") or die "couldn't open input file $infile\n";
while (my $line=<IN>) {
    chomp $line;
    my @cells = split(/\s+/, $line);
    my $nm = $cells[0];
    my $start = $cells[1];
    my $stop = $cells[2];
    my $dir = $cells[3];
    my $score = $cells[4];

    my $sid;

    if ($dir eq "UP") {
	$sid = sprintf("DMR%07d", $id++);
	print join("\t", ($nm,$start,$stop,$sid,"+","1","0","0","0","NA","0")), "\n";
	$sid = sprintf("DMR%07d", $id++);
	print join("\t", ($nm,$start,$stop,$sid,"-","1","0","0","0","NA","0")), "\n";
    }
    elsif ($dir eq "DOWN") {
	$sid = sprintf("DMR%07d", $id++);
	print join("\t", ($nm,$start,$stop,$sid,"+","0","0","1","0","NA","0")), "\n";
	$sid = sprintf("DMR%07d", $id++);
	print join("\t", ($nm,$start,$stop,$sid,"-","0","0","1","0","NA","0")), "\n";
    }
    else {
	die "wrong direction\n";
    }
}
