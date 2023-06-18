#!/usr/bin/perl

use strict ;
use warnings ;

#ARGV[0] --- ID

#input: output file of list_phased.pl
open (IN, "${ARGV[0]}.phased.list") or die ; 

my $line ;
my @LINE ;
my %ref ;
my %alt ;
my %phased ;
my @P ;


while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    if ($line =~ /^\#/) {
	next ;
    }
    @LINE = () ;
    @LINE = split (/\t/, $line) ;
    if (length ($LINE[2]) == 1 and length ($LINE[3]) == 1) { #SNV
	unless ($LINE[4] eq "SNV") {
	    die ;
	}
 	$ref{$LINE[0]}{$LINE[1]} = $LINE[2] ;
	$alt{$LINE[0]}{$LINE[1]} = $LINE[3] ;
	$phased{$LINE[0]}{$LINE[1]} = $LINE[5] ;
	push (@P, "$LINE[0]\t$LINE[1]\t$LINE[2]\t$LINE[3]") ;
    }
}

close (IN) ;

my @HP = (1, 2, "*", "otherPS") ; # HP1, HP2, unphased, other phased block

#input: output file of samtoools mpileup
open (IN, "${ARGV[0]}.SNV.pileup") or die ;

my @hp ;
my @ps ;
my $k ;
my $pileup ;
my $base ;
my %REF ;
my %ALT ;
my %OTHERS ;

my %ok ;

while (<IN>) {
    $line = $_ ;
    chomp ($line) ;
    @LINE = () ;
    @LINE = split (/\t/, $line) ;
    #HP tag
    unless ($LINE[7] =~ /^[12,\*]+$/) {
	print "$line\n" ;
	die ;
    }
    @hp = () ;
    @hp = split (/,/, $LINE[7]) ;
    #PS tag
    unless ($LINE[8] =~ /^[\d,\*]+$/) {
	print "$line\n" ;
	die ;
    }
    @ps = () ;
    @ps = split (/,/, $LINE[8]) ;
    #pileup
    $k = -1 ;
    $pileup = $LINE[4] ;
    %REF = () ;
    %ALT = () ;
    %OTHERS = () ;
    while ($pileup =~ /^./g) {
	$base = $& ;
	$pileup = "$'" ;
	if ($base =~ /^[\.,]$/) { # skip (N)
	    print "$line\n" ;
	    die ;
	} elsif ($base =~ /^[ATCGNatcgn]$/) {
	    if ($base =~ /^$ref{$LINE[0]}{$LINE[1]}$/i) { # REF
		$k++ ;
		if ($ps[$k] eq $phased{$LINE[0]}{$LINE[1]} or $ps[$k] eq "*") { # PS check
		    $REF{$hp[$k]}++ ;
		} else {
		    $REF{"otherPS"}++ ;
		}
            } elsif ($base =~ /^$alt{$LINE[0]}{$LINE[1]}$/i) { # ALT
		$k++ ;
		if ($ps[$k] eq $phased{$LINE[0]}{$LINE[1]} or $ps[$k] eq "*") { # PS check
		    $ALT{$hp[$k]}++ ;
		} else {
		    $ALT{"otherPS"}++ ;
		}
	    } else {
		$k++ ;
		if ($ps[$k] eq $phased{$LINE[0]}{$LINE[1]} or $ps[$k] eq "*") { #PS check
		    $OTHERS{$hp[$k]}++ ;
		} else {
		    $OTHERS{"otherPS"}++ ;
		}
	    }
        } elsif ($base =~ /^[\*]$/) {
	    $k++ ;
	    if ($ps[$k] eq $phased{$LINE[0]}{$LINE[1]} or $ps[$k] eq "*") { # PS check
		$OTHERS{$hp[$k]}++ ;
	    } else {
		$OTHERS{"otherPS"}++ ;
	    }
	} elsif ($base eq "+") { # information insertion after a base
	    if ($pileup =~ /^(\d+)/) {
		$pileup = substr ("$'", $&) ;
	    } else {
		print "$line\n" ;
		die ;
	    }
	} elsif ($base eq "-") { # information insertion after a base
            if ($pileup =~ /^(\d+)/) {
                $pileup = substr ("$'", $&) ;
	    } else {
                print "$line\n" ;
                die ;
            }
	} elsif ($base =~ /^[\^]$/) { # read head
	    $pileup = substr ($pileup, 1) ;
	} elsif ($base =~ /^[\$]$/) { # read tail 
	    next ;
	} else {
	    print "$line\n" ;
	    die ;
	}
    }
    unless (@hp == $k + 1) {
	print "$line\n" ;
	die ;
    }
    $ok{$LINE[0]}{$LINE[1]}{$ref{$LINE[0]}{$LINE[1]}}{$alt{$LINE[0]}{$LINE[1]}} = $LINE[3]."\t".$phased{$LINE[0]}{$LINE[1]} ;
    foreach my $item (@HP) {
	$REF{$item} += 0 ;
	$ok{$LINE[0]}{$LINE[1]}{$ref{$LINE[0]}{$LINE[1]}}{$alt{$LINE[0]}{$LINE[1]}} .= "\t$REF{$item}" ;
    }
    foreach my $item (@HP) {
        $ALT{$item} += 0 ;
	$ok{$LINE[0]}{$LINE[1]}{$ref{$LINE[0]}{$LINE[1]}}{$alt{$LINE[0]}{$LINE[1]}} .= "\t$ALT{$item}" ;
    }
    foreach my $item (@HP) {
        $OTHERS{$item} += 0 ;
	$ok{$LINE[0]}{$LINE[1]}{$ref{$LINE[0]}{$LINE[1]}}{$alt{$LINE[0]}{$LINE[1]}} .= "\t$OTHERS{$item}" ;
    }
}

close (IN) ;

#output: count reference/variant tags in each haplotype 
open (OUT, ">${ARGV[0]}_HP_SNV.txt") or die ;

print OUT "Chr\tPosition\tRef\tAlt\tTotal\tPS\t" ;
print OUT "REF_1\tREF_2\tREF_un\tREF_otherPS\t" ;
print OUT "ALT_1\tALT_2\tALT_un\tALT_otherPS\t" ;
print OUT "OTHERS_1\tOTHERS_2\tOTHERS_un\tOTHERS_otherPS\n" ;

foreach my $item (@P) {
    @LINE = split (/\t/, $item) ;
    if (exists ($ok{$LINE[0]}{$LINE[1]}{$LINE[2]}{$LINE[3]})) {
	print OUT "$LINE[0]\t$LINE[1]\t$LINE[2]\t$LINE[3]\t" ;
	print OUT "$ok{$LINE[0]}{$LINE[1]}{$LINE[2]}{$LINE[3]}\n" ;
    } else {
	print OUT "$LINE[0]\t$LINE[1]\t$LINE[2]\t$LINE[3]\t0\t$phased{$LINE[0]}{$LINE[1]}\t" ;
	print OUT "0\t0\t0\t0\t" ;
        print OUT "0\t0\t0\t0\t" ;
	print OUT "0\t0\t0\t0\n" ;
    }
}

close (OUT) ;
