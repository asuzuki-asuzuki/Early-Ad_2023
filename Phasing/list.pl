#!/usr/bin/perl

use strict ;
use warnings ;

#ARGV[0] --- ID

#input: the filtered VCF file from GATK Mutect-2
open (IN, "${ARGV[0]}_T.vcf") or die ; 

#output: the list of positions of SNVs, indels or MNVs
open (OUT, ">${ARGV[0]}_T.SNV.list") or die ;
open (OUT2, ">${ARGV[0]}_T.indel.list") or die ;
open (OUT3, ">${ARGV[0]}_T.MNV.list") or die ;

my $line ;
my @LINE ;
my $p ;
my $i ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    if ($line =~ /^\#/) {
	next ;
    }
    @LINE = () ;
    @LINE = split (/\t/, $line) ;
    if (length ($LINE[3]) == 1 and length ($LINE[4]) == 1) { # SNV
	print OUT "$LINE[0]\t$LINE[1]\n" ;
    } elsif (length ($LINE[3]) == length ($LINE[4])) { # MNV
	for ($i = 1 ; $i <= length($LINE[3]) ; $i++) {
	    $p = $LINE[1] + $i - 1 ;
	    print OUT3 "$LINE[0]\t$p\n" ;
	}
    } else { # indel
	print OUT2 "$LINE[0]\t$LINE[1]\n" ;
    }
}

close (OUT) ;
close (OUT2) ;
close (OUT3) ;
close (IN) ;
