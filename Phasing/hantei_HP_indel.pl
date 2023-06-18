#!/usr/bin/perl

use strict ;
use warnings ;

#ARGV[0] --- ID

#input: output file of extract_HP_indel.pl
open (IN, "${ARGV[0]}_T_HP_indel.txt") or die ;

open (OUT, ">${ARGV[0]}_T_HP_hantei_indel.txt") or die ;

my $line ;
my @LINE ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    if ($line =~ /^Chr\t/) {
	print OUT "$line\tHP\n" ;
	next ;
    }
    @LINE = () ;
    @LINE = split (/\t/, $line) ;
    if ($LINE[10] <= 1 and $LINE[11] >= 3) {
	print OUT "$line\tHP2\n" ;
    } elsif ($LINE[11] <= 1 and $LINE[10] >= 3) {
	print OUT "$line\tHP1\n" ;
    } else {
	print OUT "$line\tHPunk\n" ;
    }
}

close (IN) ;
close (OUT) ;
