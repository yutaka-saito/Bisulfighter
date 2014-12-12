#!/usr/bin/perl
# 
# Simulate DMRs so that neighbor CpGs are coordinately changed to 
# the same direction (hyper- or hypo- methylation)   
#
# usage:
# ./this_program chrX.met.fa 
# 
# Bisulfighter (http://epigenome.cbrc.jp/bisulfighter)
# by National Institute of Advanced Industrial Science and Technology (AIST)
# is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# http://creativecommons.org/licenses/by-nc-sa/3.0/

use strict;
use warnings;

our @NDMR = (
  100,
  200,
  500,
);
our @LEN = (
  50,
  500,
  5000,
  50000,
);

our $NAME = "chrX";

srand(time^$$);
for (my $i=0; $i<@NDMR; $i++) {
    for (my $j=0; $j<@LEN; $j++) {
	&region_simulate($NDMR[$i], $LEN[$j], $ARGV[0]);
    }
}


sub region_simulate {
    my ($ndmr, $len, $infile) = @_;
    my $outfile = join(".", ($NAME, "met", "dmr", $ndmr, $len, "fa"));  
    my $dmbfile = join(".", ($ndmr, $len, "dmb"));  
    my $dmrfile = join(".", ($ndmr, $len, "dmr"));  
    my $relax = 2.0;

    my $name;
    my $seq;
    open(IN, $infile) or die "couldn't open input file $infile\n";
    $name = <IN>;
    $seq = join("", <IN>);
    $seq =~ s/\s//g;
    close(IN);

    my @pos_cpg;
    my %dmr_cpg; # from start to stop

    while ($seq =~ /[Ccdvt][Gghba]/g) {
	push(@pos_cpg, pos($seq)-2);
#	print substr($seq, pos($seq)-2, 2), " "; # should be CpG
    }
    for (my $i=0; $i<2*$ndmr; $i++) {
	my ($start, $stop) = (0, 0);
	my $j = int(rand @pos_cpg);
	$start = $pos_cpg[$j];
	$stop = $start;
	while (++$j < @pos_cpg) {
	    $stop = $pos_cpg[$j];
	    $stop - $start + 1 > $len and last;
	}
	if ($j==@pos_cpg) {
	    print STDERR "reached to the end\n";
	    redo;
	}
	if (! ($stop - $start + 1 < $len * $relax)) { 
	    print STDERR "could not find a moderate-size DMR\n";
	    redo;
	}
	my $overlap = 0;
	foreach my $st (keys %dmr_cpg) { 
	    my $sp = $dmr_cpg{$st};
	    ($start<=$sp and $stop>=$st) and $overlap=1 and last;
	}
	if ($overlap) {
#	    print STDERR "overlapped with previous DMRs\n";
	    print STDERR ".";
	    redo;
	}
	$dmr_cpg{$start} = $stop;
	print STDERR "DMR registered: $start ", $stop+1, " (", $stop-$start+1, ")\n";
    }

    open(DMB, "> $dmbfile");
    open(DMR, "> $dmrfile");
    foreach my $start (sort {$a <=> $b} keys(%dmr_cpg)) {
	my $stop = $dmr_cpg{$start};
	my $s = substr($seq, $start, ($stop-$start+1)+1);
#	print STDERR $s, "\n";
	if (rand(1000) % 2) { # UP
	    while ($s =~ /[Ccdvt][Gghba]/g) {
		my $p = $start + pos($s)-2;
		print DMB join("\t", ($NAME, $p, $p+1, "UP", "CpG")), "\n";
	    }
	    print DMR join("\t", ($NAME, $start, $stop+1, "UP", "CpG")), "\n";
	    $s =~ s/[Ccdvt][Gghba]/ta/g;
	    substr($seq, $start, ($stop-$start+1)+1) = $s;
	}
	else { # DOWN
	    while ($s =~ /[Ccdvt][Gghba]/g) {
		my $p = $start + pos($s)-2;
		print DMB join("\t", ($NAME, $p, $p+1, "DOWN", "CpG")), "\n";
	    }
	    print DMR join("\t", ($NAME, $start, $stop+1, "DOWN", "CpG")), "\n";
	    $s =~ s/[Ccdvt][Gghba]/CG/g;
	    substr($seq, $start, ($stop-$start+1)+1) = $s;
	}
    }
    close(DMB);
    close(DMR);

    open(OUT, "> $outfile");
    $seq =~ s/(.{50})/$1\n/g;
    print OUT $name, $seq, "\n";
    close(OUT);
}
