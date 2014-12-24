#!/usr/bin/perl

# Bisulfighter http://epigenome.cbrc.jp/bisulfighter
# by National Institute of Advanced Industrial Science and Technology (AIST)
# is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# http://creativecommons.org/licenses/by-nc-sa/3.0/

use strict;
use warnings;

# format Bisulfighter's input to ComMet's input
# (data from both strands are integrated)
#
# usage:
# ./this_program sample1.tsv sample2.tsv > ComMetInput

# input format 
=pod
# sample1
chr1  123      +  CG   0.8    10
chr1  128      +  CHG  0.7    8
chr1  1000000  +  CG   0.234  9
chr2  987      +  CG   0.654  12
chr2  1000     +  CHG  0.788  10
chr2  1057     +  CG   0.826  50
=cut

=pod
# sample2
chr1  123      +  CG   0.3    30
chr1  128      +  CHG  0.24   20
chr1  1000000  +  CG   0.36   8
chr2  987      +  CG   0.64   12
chr2  1000     +  CHG  0.820  9
chr2  1057     +  CG   0.0    32
=cut

# output format 
=pod
chr1  123      10*(0.8)    10*(1-0.8)    30*(0.3)   30*(1-0.3)
chr1  1000000  9*(0.234)   9*(1-0.234)   8*(0.36)   8*(1-0.36)
chr2  987      12*(0.654)  12*(1-0.654)  12*(0.64)  12*(1-0.64)
chr2  1057     50*(0.826)  50*(1-0.826)  32*(0.0)   32*(1-0.0)
=cut

my $command = "";
my $tmpout1 = "tmp.$$.1";
my $tmpout2 = "tmp.$$.2";

&uniq_add($ARGV[0], $tmpout1);
&uniq_add($ARGV[1], $tmpout2);

my ($nm1, $pos1, $m1, $u1) = ("", -1, -1, -1);
my ($nm2, $pos2, $m2, $u2) = ("", -1, -1, -1);
open(IN1, $tmpout1) or die "couldn't open input file $tmpout1\n";
open(IN2, $tmpout2) or die "couldn't open input file $tmpout2\n";
while (1) {
    my ($line1, $line2);

    if (&chreq($nm1,$pos1,$nm2,$pos2)) {
	(defined($line1=<IN1>) and defined($line2=<IN2>)) or last;
	chomp $line1;
	($nm1,$pos1,$m1,$u1) = split(/\s+/, $line1);
	chomp $line2;
	($nm2,$pos2,$m2,$u2) = split(/\s+/, $line2);
    }
    elsif (&chrlt($nm1,$pos1,$nm2,$pos2)) {
	defined($line1=<IN1>) or last;
	chomp $line1;
	($nm1,$pos1,$m1,$u1) = split(/\s+/, $line1);
    }
    elsif (&chrgt($nm1,$pos1,$nm2,$pos2)) {
	defined($line2=<IN2>) or last;
	chomp $line2;
	($nm2,$pos2,$m2,$u2) = split(/\s+/, $line2);
    }
    else {
	die "wrong position\n";
    }

    if (&chreq($nm1,$pos1,$nm2,$pos2)) {
	print STDOUT join("\t", ($nm1, $pos1, $m1, $u1, $m2, $u2)), "\n";
    }
}
close(IN1);
close(IN2);

&invoke_command("\\rm -f $tmpout1");
&invoke_command("\\rm -f $tmpout2");

sub invoke_command {
    my ($command) = @_;
    print STDERR $command, "\n";
    system $command;
}

sub chrgt {
    my ($nm1, $pos1, $nm2, $pos2) = @_;

    if ($nm1 gt $nm2 || ($nm1 eq $nm2 && $pos1 > $pos2)) {
	return 1;
    }
    else {
	return 0;
    }
}
sub chrlt {
    my ($nm1, $pos1, $nm2, $pos2) = @_;

    if ($nm1 lt $nm2 || ($nm1 eq $nm2 && $pos1 < $pos2)) {
	return 1;
    }
    else {
	return 0;
    }
}
sub chreq {
    my ($nm1, $pos1, $nm2, $pos2) = @_;

    if ($nm1 eq $nm2 && $pos1 == $pos2) {
	return 1;
    }
    else {
	return 0;
    }
}

sub uniq_add {
    my ($infile, $outfile) = @_;
    my ($nm, $str, $pos, $m, $u) = ("", "", -1, -1, -1);
    my ($pre_nm, $pre_str, $pre_pos, $pre_m, $pre_u) = ("", "", -1, -1, -1);
    my $mind = 2; # minimum distance between two CpGs
    my $thsh = 0;

    open(IN, $infile) or die "couldn't open input file $infile\n";
    open(OUT, "> $outfile") or die "couldn't open output file $outfile\n";
    while (my $line=<IN>) {
	chomp $line;
	$line =~ /^\#/ and next;
	$line =~ /CG/ or next;
	my @cells = split(/\s+/, $line);
	$nm = $cells[0];
	$str = $cells[2];
	$pos = $cells[1];
	$m = $cells[5] * $cells[4];
	$u = $cells[5] * (1.0-$cells[4]);
	
	if ($str eq "+") {
	}
	elsif ($str eq "-") {
	    $str = "+";
	    $pos -= $mind - 1; 
	}
	else {
	    print STDERR "warning: abnormal strand $str was removed\n";
	    next;
	}


	if ($nm eq $pre_nm && $str eq $pre_str && $pos==$pre_pos) {
	    $pre_m += $m;
	    $pre_u += $u;
	}
	else {
	    if ($nm eq $pre_nm && $pos-$pre_pos < $mind) {
		print STDERR "warning: abnormal position difference $pre_pos $pos, ";
		if ($m + $u > $pre_m + $pre_u) {
		    print STDERR "$pre_pos was removed.\n";
		    $pre_nm = $nm;
		    $pre_str = $str;
		    $pre_pos = $pos;
		    $pre_m = $m;
		    $pre_u = $u;
		}
		else {
		    print STDERR "$pos was removed.\n"
		    }
		next;
	    }
	    if ($pre_m + $pre_u > $thsh) {
		print OUT join("\t", ($pre_nm, $pre_pos, $pre_m, $pre_u)), "\n";
	    }
	    $pre_nm = $nm;
	    $pre_str = $str;
	    $pre_pos = $pos;
	    $pre_m = $m;
	    $pre_u = $u;
	}
    }
    if ($pre_m + $pre_u > $thsh) {
	print OUT join("\t", ($pre_nm, $pre_pos, $pre_m, $pre_u)), "\n";
    }
    close(IN);
    close(OUT);
}

