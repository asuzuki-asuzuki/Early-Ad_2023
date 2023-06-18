#!/usr/bin/perl

use strict ;
use warnings ;

#ARGV[0] --- ID

#input: block regions from WhatsHap
open (IN, "$ARGV[0]\_T.block-list.tsv") or die ;

my $line ;
my @LINE ;
my %c ;
my %s ;
my %e ;
my %phased ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    if ($line =~ /^\#/) {
	next ;
    }
    @LINE = () ;
    @LINE = split (/\t/, $line) ;
    if ($LINE[4] - $LINE[3] + 1 == 1) {
	next ; # block len == 1
    }
    $c{$LINE[1]}++ ;
    $s{$LINE[1]}{$c{$LINE[1]}} = $LINE[3] ;
    $e{$LINE[1]}{$c{$LINE[1]}} = $LINE[4] ;
    $phased{$LINE[1]}{$c{$LINE[1]}} = $LINE[2] ;
}

close (IN) ;

#input: the filtered VCF file from GATK Mutect-2
open (IN, "${ARGV[0]}_T.vcf") or die ;

#output
open (OUT, ">${ARGV[0]}_T.phased.list") or die ;

my $p ;
my $i ;
my $k ;
my %PS ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    if ($line =~ /^\#/) {
	next ;
    }
    @LINE = () ;
    @LINE = split (/\t/, $line) ;
    if (length ($LINE[3]) == 1 and length ($LINE[4]) == 1) { # SNV
	print OUT "$LINE[0]\t$LINE[1]\t$LINE[3]\t$LINE[4]\tSNV" ;
	%PS = () ;
	unless (exists ($c{$LINE[0]})) {
	    print OUT "\tun\n" ;
	    next ;
	}
	for ($i = 1 ; $i <= $c{$LINE[0]} ; $i++) {
	    if ($s{$LINE[0]}{$i} <= $LINE[1] and $LINE[1] <= $e{$LINE[0]}{$i}) {
		$PS{$phased{$LINE[0]}{$i}}++ ;
	    }
	}
	if (keys %PS == 1) { # in the PS
	    foreach my $item (keys %PS) {
		print OUT "\t$item" ;
	    }
	    print OUT "\n" ;
	} elsif (keys %PS == 0) { # not in PS
	    print OUT "\tun\n" ;
	} else {
	    die ;
	}
    } elsif (length ($LINE[3]) == length ($LINE[4])) { # MNV
	print OUT "$LINE[0]\t$LINE[1]\t$LINE[3]\t$LINE[4]\tMNV" ;
	unless (exists ($c{$LINE[0]})) {
            print OUT "\tun\n" ;
            next ;
        }
	%PS = () ;
	for ($k = 1 ; $k <= length($LINE[3]) ; $k++) { # check all bases for MNV
	    $p = $LINE[1] + $k - 1 ;
	    for ($i = 1 ; $i <= $c{$LINE[0]} ; $i++) {
		if ($s{$LINE[0]}{$i} <= $p and $p <= $e{$LINE[0]}{$i}) {
		    $PS{$phased{$LINE[0]}{$i}}++ ;
		}
	    }
	}
	if (keys %PS == 1) {
	    foreach my $item (keys %PS) {
		if ($PS{$item} == length ($LINE[3])) { # in the PS
		    print OUT "\t$item" ; # all MNV in the same PS
		} else {
		    print OUT "\tun" ; # not all MNV bases in the same PS
		}
	    }
	    print OUT "\n" ;
        } elsif (keys %PS <= length ($LINE[3])) { # not in PS or in multiple PS
	    print OUT "\tun\n" ;
	} else {
	    die ;
	}
    } else { # INDEL
	print OUT "$LINE[0]\t$LINE[1]\t$LINE[3]\t$LINE[4]\tINDEL" ;
	unless (exists ($c{$LINE[0]})) {
            print OUT "\tun\n" ;
            next ;
        }
	%PS = () ;
        for ($i = 1 ; $i <= $c{$LINE[0]} ; $i++) {
            if ($s{$LINE[0]}{$i} <= $LINE[1] and $LINE[1] <= $e{$LINE[0]}{$i}) {
		$PS{$phased{$LINE[0]}{$i}}++ ;
	    }
        }
        if (keys %PS == 1) { # in the PS
            foreach my $item (keys %PS) {
                print OUT "\t$item" ;
            }
            print OUT "\n" ;
        } elsif (keys %PS == 0) { # not in PS
            print OUT "\tun\n" ;
        } else {
            die ;
        }
    }
}

close (OUT) ;
close (IN) ;

