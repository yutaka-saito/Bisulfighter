#!/usr/bin/perl

# Bisulfighter http://epigenome.cbrc.jp/bisulfighter
# by National Institute of Advanced Industrial Science and Technology (AIST)
# is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# http://creativecommons.org/licenses/by-nc-sa/3.0/

use strict;
use warnings;
use Getopt::Long qw(:config posix_default no_ignore_case gnu_compat);

# convert bsf-call output to ComMet input
# (data from the reverse strand are integrated)
#
# usage:
# ./this_program --cpg --sample1 BsfOutput1.tsv --sample2 BsfOutput2.tsv > ComMetInput

# input format
=pod
# sample1
chr1  123      +  CG   0.8    10
chr1  128      +  CHG  0.7    8
chr1  1000000  -  CG   0.234  9
chr2  987      -  CG   0.654  12
chr2  1000     +  CHG  0.788  10
chr2  1057     +  CG   0.826  50
=cut

=pod
# sample2
chr1  123      +  CG   0.3    30
chr1  128      +  CHG  0.24   20
chr1  1000000  -  CG   0.36   8
chr2  987      -  CG   0.64   12
chr2  1000     +  CHG  0.820  9
chr2  1057     +  CG   0.0    32
=cut

# output format
=pod
chr1  1000000  9*(0.234)   9*(1-0.234)   8*(0.36)   8*(1-0.36)
chr2  987      12*(0.654)  12*(1-0.654)  12*(0.64)  12*(1-0.64)
=cut


my $context = "";
my ($cpg, $chg, $chh) = (0 ,0, 0);
my ($sample1, $sample2) = ("", "");
my ($tmpout1, $tmpout2) = ("tmp.$$.1", "tmp.$$.2");

GetOptions(
    "--cpg" => \$cpg,
    "--chg" => \$chg,
    "--chh" => \$chh,
    "--sample1=s" => \$sample1,
    "--sample2=s" => \$sample2,
);

if ($cpg) {
    $context = "CG";
} 
elsif ($chg) {
    $context = "CHG";
}
elsif ($chh) {
    $context = "CHH";
}
else {
    die "specify the context with --cpg, --chg, or --chh\n";
}
$sample1 eq "" and die "specify the file name with --sample1\n";
$sample2 eq "" and die "specify the file name with --sample2\n";

&bsf_unstrand($sample1, $tmpout1, $context);
&bsf_unstrand($sample2, $tmpout2, $context);

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

sub bsf_unstrand {
    my ($infile, $outfile, $context) = @_;
    my $thsh = 0;
    my $Table = [];
    my @Print;

    open(IN, $infile) or die "couldn't open input file $infile\n";
    while (my $line=<IN>) {
	chomp $line;
	my ($nm, $pos, $str, $cxt, $r, $d) = split(/\s+/, $line);;
	my ($m, $u) = ($r * $d, (1.0-$r) * $d);
	push(@$Table, [$nm, $pos, $str, $cxt, $m, $u]);
	push(@Print, 1);
    }
    close(IN);

    open(OUT, "> $outfile") or die "couldn't open output file $outfile\n";
    for (my $i=0; $i<@$Table; $i++) {
	$Print[$i]==1 or next;
	my ($nm, $pos, $str, $cxt, $m, $u) = @{$Table->[$i]};
	$str eq "-" or next;
	$cxt eq $context or next;
	$m+$u > $thsh or next;
	print OUT join("\t", ($nm, $pos, $m, $u)), "\n";
    }
    close(OUT);
}

