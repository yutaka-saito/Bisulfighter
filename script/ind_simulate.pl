#!/usr/bin/perl
# 
# Simulate DMRs so that neighbor CpGs are independently changed to 
# random direction (hyper- or hypo- methylation)   
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

our @RATE = (
  10,
  20,
  40,
);

our $NAME = "chrX";

srand(time^$$);
for (my $i=0; $i<@RATE; $i++) {
    &base_simulate($RATE[$i], $ARGV[0]);
}


sub base_simulate {
    my ($rate, $infile) = @_;
    my $r = sprintf("%02d", $rate);
    my $outfile = join(".", ($NAME, "met", "dmr", $r, "fa"));  
    my $dmbfile = join(".", ($r, "dmb"));  
    my $dmrfile = join(".", ($r, "dmr"));  
    my $ndmb = 0;

    my $name;
    my $seq;
    open(IN, $infile) or die "couldn't open input file $infile\n";
    $name = <IN>;
    $seq = join("", <IN>);
    $seq =~ s/\s//g;
    close(IN);

    my @pos_cpg;
    my %dmb_cpg; # from start to UP (1) / DOWN (0)

    while ($seq =~ /[Ccdvt][Gghba]/g) {
	push(@pos_cpg, pos($seq)-2);
#	print substr($seq, pos($seq)-2, 2), " "; # should be CpG
    }

    $ndmb = int(@pos_cpg * $rate / 100);
    for (my $i=0; $i<$ndmb; $i++) {
	my $pos = $pos_cpg[ int(rand @pos_cpg) ];
	if (exists($dmb_cpg{$pos})) {
#	    print STDERR "overlapped with previous DMBs\n";
	    print STDERR ".";
	    redo;
	}

	my $ctx = substr($seq, $pos, 2);
	$ctx =~ /^[Ccdvt][Gghba]$/ or die;
	if (rand(1000) % 2) { # UP
	    $ctx eq "ta" and redo;
	    $ctx = "ta";
	    $dmb_cpg{$pos} = 1;
	}
	else { # DOWN
	    $ctx eq "CG" and redo;
	    $ctx = "CG";
	    $dmb_cpg{$pos} = 0;
	}
	substr($seq, $pos, 2) = $ctx;
	print STDERR "DMB:$pos ";
    }
    print STDERR "\n\n";

    open(DMB, "> $dmbfile");
    foreach my $pos (sort {$a <=> $b} keys %dmb_cpg) {
	if ($dmb_cpg{$pos}) {
	    print DMB join("\t", ($NAME, $pos, $pos+1, "UP", "CpG")), "\n";
	}
	else {
	    print DMB join("\t", ($NAME, $pos, $pos+1, "DOWN", "CpG")), "\n";
	}
    }	
    close(DMB);

    open(DMR, "> $dmrfile");
    my ($st_cpg, $sp_cpg) = (0, 0);
    my ($up_cpg, $dn_cpg) = (0, 0);
    foreach my $pos (sort {$a <=> $b} keys %dmb_cpg) {
	my $s = substr($seq, $sp_cpg, $pos-$sp_cpg);
	if ($s =~ /[Ccdvt][Gghba]/) {
	    if ($sp_cpg-$st_cpg >= 2) {
		$up_cpg and print DMR join("\t", ($NAME, $st_cpg, $sp_cpg, "UP", "CpG")), "\n";
		$dn_cpg and print DMR join("\t", ($NAME, $st_cpg, $sp_cpg, "DOWN", "CpG")), "\n";
		print STDERR "DMR: $st_cpg $sp_cpg (", $sp_cpg-$st_cpg, ")\n";
	    }
	    $st_cpg = $pos;
	    $up_cpg = 0;
	    $dn_cpg = 0;
	}

	if ($dmb_cpg{$pos}) {
	    if ($dn_cpg) {
		if ($sp_cpg-$st_cpg >= 2) {
		    print DMR join("\t", ($NAME, $st_cpg, $sp_cpg, "UP", "CpG")), "\n";
		    print STDERR "DMR: $st_cpg $sp_cpg (", $sp_cpg-$st_cpg, ")\n";
		}
		$st_cpg = $pos;
	    }
	    $up_cpg = 1;
	    $dn_cpg = 0;
	}
	else {
	    if ($up_cpg) {
		if ($sp_cpg-$st_cpg >= 2) {
		    print DMR join("\t", ($NAME, $st_cpg, $sp_cpg, "DOWN", "CpG")), "\n";
		    print STDERR "DMR: $st_cpg $sp_cpg (", $sp_cpg-$st_cpg, ")\n";
		}
		$st_cpg = $pos;
	    }
	    $up_cpg = 0;
	    $dn_cpg = 1;
	}
	$sp_cpg = $pos + 1;
    }

    open(OUT, "> $outfile");
    $seq =~ s/(.{50})/$1\n/g;
    print OUT $name, $seq, "\n";
    close(OUT);
}
